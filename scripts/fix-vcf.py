#!/usr/bin/env python3
import argparse
import sys
import pysam


def main():
    parser = argparse.ArgumentParser(description="Fix STR/VNTR VCF rows using sample alleles + reference FASTA.")
    parser.add_argument("input", help="Input VCF file")
    parser.add_argument("output", help="Output VCF file")
    parser.add_argument("--ref", required=True, help="Reference FASTA file (indexed .fai required)")
    args = parser.parse_args()

    fasta = pysam.FastaFile(args.ref)

    with open(args.input, "r") as infile, open(args.output, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
                continue

            fixed = fix_row(line, fasta)
            outfile.write(fixed + "\n")


def fix_row(row: str, fasta: pysam.FastaFile) -> str:
    """
    Expects VCF with at least 10 columns:
      CHROM POS ID REF ALT QUAL FILTER INFO FORMAT SAMPLE

    Sample field assumed to be like:
      GT:ALLELE1:ALLELE2...
    where ALLELE1/2 are the actual sequences to use as ALT alleles.
    Sample input data
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG004.30x.haplotagged
    chr1	2795086	.	N	<VNTR>	.	PASS	END=2795101;RU=A;SVTYPE=STR;ALTANNO_H1=0-0-0-0-0-0-0-0-0-0-0-0-0-0-0;LEN_H1=15;ALTANNO_H2=0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0;LEN_H2=16;	GT:RS	1/2:AAAAAAAAAAAAAAA,AAAAAAAAAAAAAAAA

    Sample output data
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG004.30x.haplotagged
    chr1	2795086	.	AAAAAAAAAAAAAAA	AAAAAAAAAAAAAAAA	.	PASS	END=2795101;RU=A;SVTYPE=STR;ALTANNO_H1=0-0-0-0-0-0-0-0-0-0-0-0-0-0-0;LEN_H1=15;ALTANNO_H2=0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0;LEN_H2=16;	GT:RS	0/1
    """

    fields = row.rstrip("\n").split("\t")
    if len(fields) < 10:
        return row.rstrip("\n")  # not a standard VCF data line; leave unchanged

    chrom = fields[0]
    start_pos = int(fields[1])  # Convert from 0-based to 1-based
    info_field = fields[7]

    # VAMOS alleles are consistently 1 bp shorter than the ref with the last bp missing so I think this resolves it
    end_pos = int(info_field.split("END=")[1].split(";")[0])
    # This adjustment no longer needed as of vamos v3.0.6, which outputs correct END positions
    #end_pos -= 1
    # Replace the adjusted END position in the INFO field
    #info_field = info_field.replace(f"END={end_pos + 1}", f"END={end_pos}")

    sample_field = fields[9]
    genotype = sample_field.split(":")[0]
    alt_alleles = sample_field.split(":")[1]
    if "," in alt_alleles:
        alt_alleles = alt_alleles.split(",")
    else:
        alt_alleles = [alt_alleles]

    # Fetch reference allele (inclusive end; pysam fetch is 0-based start, end-exclusive)
    ref_allele = fetch_ref_allele(chrom, start_pos + 1, end_pos, fasta)

    # Set REF and ALT from sequences
    fields[3] = ref_allele
    fields[4] = "." if len(alt_alleles) == 0 else ",".join(alt_alleles)

    # Determine genotype separator
    sep = "|" if "|" in genotype else "/"
    genotype_list = genotype.split(sep)

    # If any ALT allele equals REF, convert that allele index to 0 and remove that ALT
    # NOTE: removing elements while iterating => use while loop
    i = 0
    while i < len(alt_alleles):
        if alt_alleles[i] == ref_allele:
            # remap genotype: allele i+1 becomes 0
            old_idx = str(i + 1)
            genotype_list = ["0" if g == old_idx else g for g in genotype_list]

            # remove the ALT allele
            del alt_alleles[i]

            # update INFO fields ALTANNO_Hx / LEN_Hx to ALTANNO_H0 / LEN_H0
            # (kept as your original intent; simple string replace)
            info_field = info_field.replace(f"ALTANNO_H{i + 1}", "ALTANNO_H0")
            info_field = info_field.replace(f"LEN_H{i + 1}", "LEN_H0")

            # do NOT increment i because list shifted
            continue
        i += 1

    # Fix skipped sequential alleles after removals (your specific cases)
    if genotype_list == ["0", "2"]:
        genotype_list = ["0", "1"]
    elif genotype_list == ["2", "2"]:
        genotype_list = ["1", "1"]
    elif genotype_list == ["1", "0"] and sep == "/":  # only swap if unphased
        genotype_list = ["0", "1"]

    # Update ALT field after potential removals
    fields[4] = "." if len(alt_alleles) == 0 else ",".join(alt_alleles)

    # Write updated INFO and SAMPLE (GT only)
    fields[7] = info_field
    fields[9] = sep.join(genotype_list)

    return "\t".join(fields)


def fetch_ref_allele(chrom: str, start: int, end: int, fasta: pysam.FastaFile) -> str:
    """
    Fetch reference sequence from FASTA.
      start: 1-based inclusive
      end:   1-based inclusive

    pysam.fetch(chrom, start0, end0) uses:
      start0: 0-based inclusive
      end0:   0-based exclusive
    so we pass (start-1, end) to include end base.
    """
    try:
        seq = fasta.fetch(chrom, start - 1, end)
        return seq.upper()
    except Exception as e:
        sys.stderr.write(f"Error fetching {chrom}:{start}-{end}: {e}\n")
        return "N"


if __name__ == "__main__":
    main()
