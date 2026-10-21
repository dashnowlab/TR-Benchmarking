import sys
import pysam
import timeit as ti
import argparse as ap
from collections import Counter


def parse_args():
    parser = ap.ArgumentParser(prog="aam")
    parser._action_groups.pop()

    print("\nAAM - Assembly Allele Measurer\n")

    required = parser.add_argument_group('Required arguments')
    required.add_argument('-a', '--aln', required=True, metavar='<FILE>', help='Input BAM file')
    required.add_argument('-l', '--loci', required=True, metavar='<FILE>', help='Input loci file in tabix format')
    optional = parser.add_argument_group('Optional arguments')
    optional.add_argument('-o', '--out', type=str, metavar='<STR>', default='assembly_allele.bed', help='name of the output file, output is in bed format.')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    return parser.parse_args()
      

def match_jump(rpos, repeat_index, loci_coords, match_len, tracked, locus_qpos_range, qpos):
    """
    Calculate the jump in the read position due to a match.

    :param rpos: Current reference position
    :param repeat_index: Current index of the repeat in loci_coords
    :param loci_coords: List of tuples containing the start and end coordinates of the repeats
    :param match_len: Length of the match
    :param tracked: List indicating whether each repeat has been tracked
    :param locus_qpos_range: List of lists containing the query position range for each repeat
    :param qpos: Current query position
    :return: The number of repeats jumped over due to the match
    """

    r = 0 
    for r,coord in enumerate(loci_coords[repeat_index:]):
        ridx = r + repeat_index

        coord_start, coord_end = coord    
        if rpos < coord_start: break
        if rpos - match_len > coord_end: continue

        if not tracked[ridx]:
            if coord_start <= rpos:
                locus_qpos_range[ridx][0] = qpos - (rpos - coord_start)
            if coord_end < rpos:
                locus_qpos_range[ridx][1] = qpos - (rpos - coord_end)

            tracked[ridx] = True

        elif coord_end <= rpos:
            locus_qpos_range[ridx][1] = qpos - (rpos - coord_end)

    jump = 0
    if loci_coords[repeat_index + r - 1][1] < rpos:
        for coord in loci_coords[repeat_index:]:
            if coord[1] < rpos: jump += 1
            else: break

    return jump


def deletion_jump(deletion_length, rpos, repeat_index, loci_keys, tracked, loci_coords, loci_variations,
                  locus_qpos_range, qpos):
    """
    Calculate the jump in the read position due to a deletion.
    
    :param deletion_length: Length of the deletion
    :param rpos: Current reference position
    :param repeat_index: Current index of the repeat in loci_coords
    :param loci_keys: List of keys for the loci
    :param tracked: List indicating whether each repeat has been tracked
    :param loci_coords: List of tuples containing the start and end coordinates of the repeats
    :param loci_variations: Dictionary of variations for each locus
    :param locus_qpos_range: List of lists containing the query position range for each repeat
    :param qpos: Current query position
    :return: The number of repeats jumped over due to the deletion
    """

    for r, coord in enumerate(loci_coords[repeat_index:]):
        ridx = r + repeat_index
        if rpos < coord[0]: break

        del_pos = rpos - deletion_length
        if del_pos > coord[1]: continue

        locus_key = loci_keys[ridx]

        if not tracked[ridx]:
            # if the locus is not tracked
            # deletion is encountered beyond
            if coord[0] <= rpos:
                locus_qpos_range[ridx][0] = qpos        
                tracked[ridx] = True    # set tracked as true
            if coord[1] < rpos:
                locus_qpos_range[ridx][1] = qpos
        
        elif coord[1] < rpos:
            locus_qpos_range[ridx][1] = qpos

        del_len = min(coord[1], rpos) - max(coord[0], del_pos)
        loci_variations[locus_key]['alen'] -= del_len

    jump = 0
    if loci_coords[repeat_index + r - 1][1] < rpos:
        for f in loci_coords[repeat_index:]:
            if f[1] < rpos: jump += 1
            else: break
    
    return jump


def insertion_jump(insertion_length, read_seq, qpos, rpos, repeat_index, loci_keys, tracked, loci_coords,
                   loci_variations, locus_qpos_range):
    """
    Calculate the jump in the read position due to an insertion.

    :param insertion_length: Length of the insertion
    :param read_seq: The read sequence
    :param qpos: Current query position
    :param rpos: Current reference position
    :param repeat_index: Current index of the repeat in loci_coords
    :param loci_keys: List of keys for the loci
    :param tracked: List indicating whether each repeat has been tracked
    :param loci_coords: List of tuples containing the start and end coordinates of the repeats
    :param loci_variations: Dictionary of variations for each locus
    :param locus_qpos_range: List of lists containing the query position range for each repeat
    :return: The number of repeats jumped over due to the insertion
    """

    r=0
    for r, coord in enumerate(loci_coords[repeat_index:]):
        ridx = r + repeat_index
        if rpos < coord[0] - 1: break

        if rpos > coord[1]: continue

        locus_key = loci_keys[ridx]

        if not tracked[ridx]:
            # if the locus is not tracked
            # deletion is encountered beyond
            if coord[0] <= rpos:
                locus_qpos_range[ridx][0] = qpos-insertion_length
                tracked[ridx] = True    # set tracked as true
            if coord[1] == rpos:
                locus_qpos_range[ridx][1] = qpos

        elif coord[1] == rpos:
            locus_qpos_range[ridx][1] = qpos

        insert = read_seq[qpos-insertion_length:qpos]
        loci_variations[locus_key]['alen'] += insertion_length

    jump = 0
    if loci_coords[repeat_index + r - 1][1] < rpos:
        for f in loci_coords[repeat_index:]:
            if f[1] < rpos: jump += 1
            else: break
    
    return jump



def parse_cigar(cigar_tuples, read_start, loci_keys, loci_coords, loci_variations, read_seq):
    """
    Parse the CIGAR string of a read and update the loci variations accordingly.

    :param cigar_tuples: List of tuples representing the CIGAR string
    :param read_start: Start position of the read
    :param loci_keys: List of keys for the loci
    :param loci_coords: List of tuples containing the start and end coordinates of the repeats
    :param loci_variations: Dictionary of variations for each locus
    :param read_seq: The read sequence
    :return: The number of repeats jumped over due to the insertion
    """

    rpos = read_start   # NOTE: The coordinates are 1 based
    qpos = 0            # starts from 0 to sub string the read sequence in python

    repeat_index = 0
    tracked = [False]*len(loci_coords)
    locus_qpos_range = []
    for _ in loci_coords:
        locus_qpos_range.append([0,0])

    for c, cigar in enumerate(cigar_tuples):

        if cigar[0] == 4:       # soft clipped
            # soft clipping or hard clipping do not consume the reference position
            # This should be appropriately handled before
            qpos += cigar[1]

        elif cigar[0] == 2:     # deletion
            rpos += cigar[1]
            repeat_index += deletion_jump(cigar[1], rpos, repeat_index, loci_keys, tracked, loci_coords,
                                          loci_variations, locus_qpos_range, qpos)

        elif cigar[0] == 1:     # insertion
            qpos += cigar[1]
            repeat_index += insertion_jump(cigar[1], read_seq, qpos, rpos, repeat_index, loci_keys, tracked, loci_coords,
                                           loci_variations, locus_qpos_range)

        elif cigar[0] == 0 or cigar[0] == 7 or cigar[0] == 8:     # match
            qpos += cigar[1]; rpos += cigar[1]; match_len = cigar[1]
            repeat_index += match_jump(rpos, repeat_index, loci_coords, match_len, tracked, locus_qpos_range, qpos)

        elif cigar[0] == 5:
            pass

        if rpos > loci_coords[-1][1]:
            # if the reference position moves beyond the end position of the last repeat
            break

    for idx, each_key in enumerate(loci_keys):
        loci_variations[each_key]['seq'] = read_seq[locus_qpos_range[idx][0]:locus_qpos_range[idx][1]]


def cooper(bamfile, tbx, mapq_threshold = 5, flank_len = 0):
    """
    Process the assembly alignment file to extract allele from the assembly of the locus

    :param bamfile: Input BAM file
    :param tbx: Tabix file containing loci information
    :param mapq_threshold: Minimum mapping quality threshold for reads to be considered
    :param flank_len: Length of the flanking region to consider for each locus
    :return: Dictionary containing allele information for each locus
    """

    aln_alleles = {}
    # IMPORTANT NOTE: The repeat coordinate should be 1 based
    contigs     = sorted(tbx.contigs)

    # iterating through the reads in a sorted bamfile
    for contig in contigs:
        loci_vars       = {}
        loci_info       = {}
        last_coordinate = 0
        for row in tbx.fetch(contig):
            last_coordinate = int(row.split('\t')[2])

        for read in bamfile.fetch(contig):
            read_loci_vars = {}
            # skip read with low mapping quality
            if read.mapping_quality < mapq_threshold:
                continue

            read_chrom = read.reference_name
            read_start = read.reference_start
            read_end   = read.reference_end

            if read_start > last_coordinate: break
            # repeat loci covered by the read
            loci_coords = []
            loci_keys   = []
            for row in tbx.fetch(read_chrom, read_start, read_end):
                # adjust read start and end based on soft and hard clippings
                # soft and hard clippings do not consume the reference bases
                row = row.split('\t')
                locus_start = int(row[1])
                locus_end   = int(row[2])

                # if only the read completely covers the repeat
                if ( (locus_start-flank_len) >= read_start ) & ( (locus_end+flank_len) <= read_end ):
                    loci_coords.append((locus_start, locus_end))
                    locus_key = f'{read_chrom}-{locus_start}-{locus_end}'
                    loci_keys.append(locus_key)

                    if locus_key not in read_loci_vars:
                        read_loci_vars[locus_key] = {'alen': locus_end-locus_start, 'rlen': locus_end-locus_start, 'seq':''}
                    if locus_key not in loci_vars:
                        loci_vars[locus_key] = []
                    if locus_key not in loci_info:
                        loci_info[locus_key] = row

            # if no repeats are covered by the read
            if len(loci_coords) == 0: continue
            
            cigar_tuples = read.cigartuples
            read_seq = read.query_sequence

            parse_cigar(cigar_tuples, read_start, loci_keys, loci_coords, read_loci_vars, read_seq)

            for locus_key in read_loci_vars:
                loci_vars[locus_key].append(read_loci_vars[locus_key])

        
        for loc in loci_vars:
            alen_counter = Counter()
            seq_counter = Counter()
            for a in loci_vars[loc]:
                if a['alen'] in alen_counter:
                    alen_counter[a['alen']] += 1
                else: alen_counter[a['alen']] = 1

                if a['seq'] in seq_counter:
                    seq_counter[a['seq']] += 1
                else: seq_counter[a['seq']] = 1
            
            if len(alen_counter) == 1:
                alen       = alen_counter.most_common()[0][0]
                alen_count = alen_counter.most_common()[0][1]
                seq = seq_counter.most_common()[0][0]

                locus_info = loci_info[loc]
                aln_alleles[loc] = (locus_info, alen, seq, alen_count)

    return aln_alleles


if __name__ == "__main__":
    
    t1 = ti.default_timer()
    args = parse_args()
    aln_file =  args.aln
    loci_file = args.loci
    out_file = args.out

    bamfile=pysam.AlignmentFile(aln_file, 'rb')
    tbx = pysam.Tabixfile(loci_file)

    cooper(bamfile, tbx, 5, 0)
    tbx.close()
    bamfile.close()
    t2 = ti.default_timer()
    sys.stderr.write('CPU time: {} seconds\n'.format(t2 - t1))
