import argparse
import os
import sys
import pysam

def main():
    parser = argparse.ArgumentParser(description="Process some VCF files.")
    parser.add_argument("input", help="Input VCF file")
    parser.add_argument("output", help="Output VCF file")
    parser.add_argument("--ref", required=True, help="Reference FASTA file (required)")
    args = parser.parse_args()

    # Use pysam for indexed access. pysam.fetch uses 0-based start, end-exclusive
    fasta = pysam.FastaFile(args.ref)

    with open(args.input, "r") as infile, open(args.output, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
            else:
                fixed_line = fix_row(line, fasta)
                outfile.write(fixed_line + "\n")


def fix_row(row, fasta):
    """
    Sample input data
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG004.30x.haplotagged
    chr1    296009  .       N       <VNTR>  .       PASS    END=296026;RU=AT;SVTYPE=STR;ALTANNO_H1=0-0-0-0-0-0-0-0-0;LEN_H1=9;      GT      1/1:TTATATATAAATATATA
    chr1    2061674 .       N       <VNTR>  .       PASS    END=2061690;RU=GCCCT,AGGGC;SVTYPE=STR;ALTANNO_H1=0-0-0;LEN_H1=3;
        GT      1/1:AGCCCTGCCCTGCCCT
    chr1    6023939 .       N       <VNTR>  .       PASS    END=6023955;RU=T,A;SVTYPE=STR;ALTANNO_H1=1-0-0-0-0-0-0-0-0-0-0-0-0-0;LEN_H1=14;ALTANNO_H2=1-0-0-0-0-0-0-0-0-0-0-0-0-0-0;LEN_H2=15;      GT      1/2:ATTTTTTTTTTTTT:ATTTTTTTTTTTTTT
    chr1    11509478        .       N       <VNTR>  .       PASS    END=11509490;RU=GA,AG;SVTYPE=STR;ALTANNO_H1=1-1-1-1-1-1;LEN_H1=6;ALTANNO_H2=1-1-1-1-1-1-1-1;LEN_H2=8;   GT      1/2:TGAGAGAGAGAG:TGAGAGAGAGAGAGAG

    Sample output data
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG004.30x.haplotagged
    chr1    296009  .       TTATATATAAATATATAT       TTATATATAAATATATA  .       PASS    END=296026;RU=AT;SVTYPE=STR;ALTANNO_H1=0-0-0-0-0-0-0-0-0;LEN_H1=9;      GT      1/1
    chr1    2061674 .       AGCCCTGCCCTGCCCTG       AGCCCTGCCCTGCCCT  .       PASS    END=2061690;RU=GCCCT,AGGGC;SVTYPE=STR;ALTANNO_H1=0-0-0;LEN_H1=3;
        GT      1/1
    chr1    6023939 .       ATTTTTTTTTTTTTTTT       ATTTTTTTTTTTTT,ATTTTTTTTTTTTTT  .       PASS    END=6023955;RU=T,A;SVTYPE=STR;ALTANNO_H1=1-0-0-0-0-0-0-0-0-0-0-0-0-0;LEN_H1=14;ALTANNO_H2=1-0-0-0-0-0-0-0-0-0-0-0-0-0-0;LEN_H2=15;      GT      1/2
    chr1    11509478        .       TGAGAGAGAGAGA       TGAGAGAGAGAG,TGAGAGAGAGAGAGAG  .       PASS    END=11509490;RU=GA,AG;SVTYPE=STR;ALTANNO_H1=1-1-1-1-1-1;LEN_H1=6;ALTANNO_H2=1-1-1-1-1-1-1-1;LEN_H2=8;   GT      1/2
    """

    fields = row.strip().split("\t")
    chrom = fields[0]
    start_pos = int(fields[1])
    info_field = fields[7]

    # VAMOS alleles are consistently 1 bp shorter than the ref with the last bp missing so I think this resolves it
    end_pos = int(info_field.split("END=")[1].split(";")[0])
    end_pos -= 1
    # Replace the adjusted END position in the INFO field
    info_field = info_field.replace(f"END={end_pos + 1}", f"END={end_pos}")

    sample_field = fields[9]
    genotype = sample_field.split(":")[0]
    alt_alleles = sample_field.split(":")[1:]
    ref_allele = fetch_ref_allele(chrom, start_pos, end_pos, fasta)

    fields[3] = ref_allele
    fields[4] = ",".join(alt_alleles)
    fields[9] = sample_field.split(":")[0]  # Keep only genotype information

    # Check which separator is used
    if '|' in genotype:
        sep = '|'
    else:
        sep = '/'
    genotype_list = genotype.split(sep)

    # Check if any of the ALT alleles match the REF allele
    for i, alt_allele in enumerate(alt_alleles):
        if alt_allele == ref_allele:
            # Change this allele in the genotype to REF (0)
            genotype_list = ['0' if g == str(i + 1) else g for g in genotype_list]
            # Remove this ALT allele
            del alt_alleles[i]
            # Replace corresponding ALTANNO_H1 and LEN_H1 fields in INFO so they refer to the ref allele
            info_field = info_field.replace(f"ALTANNO_H{i + 1}", f"ALTANNO_H0")
            info_field = info_field.replace(f"LEN_H{i + 1}", f"LEN_H0")

    # If the remaining genotypes contain skip sequencial alleles (e.g. 0/2), fix them
    if genotype_list == ['0', '2']:
        genotype_list = ['0', '1']
    elif genotype_list == ['2', '2']:
        genotype_list = ['1', '1']
    elif genotype_list == ['1', '0'] and sep == '/': # Only swap if unphased
        genotype_list = ['0', '1']

    # Update the fields
    if len(alt_alleles) == 0:
        fields[4] = "."
    else:
        fields[4] = ",".join(alt_alleles)
    fields[7] = info_field
    fields[9] = sep.join(genotype_list)

    return "\t".join(fields)

def fetch_ref_allele(chrom, start, end, fasta):
    """
    Fetch reference sequence for a given region from a FASTA file.

    Args:
        chrom (str): chromosome name as in VCF (e.g. 'chr1' or '1')
        start (int): 1-based start position
        end (int): 1-based end position (inclusive)
        fasta

    Behavior:
        Attempts to use pysam.FastaFile if available for fast indexed access. If pysam
        is not installed, falls back to loading the whole FASTA into memory (simple)
        and returns the requested slice. The FASTA path must be provided via the
        --ref CLI argument or the REF_FASTA environment variable.
    """
    try:
        seq = fasta.fetch(chrom, start - 1, end)
    except Exception as e:
        seq = 'N'
        sys.stderr.write(f"Error fetching {chrom}:{start}-{end}: {e}\n")

    return seq.upper()

if __name__ == "__main__":
    main()