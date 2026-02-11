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
    """

    fields = row.rstrip("\n").split("\t")
    if len(fields) < 10:
        return row.rstrip("\n")  # not a standard VCF data line; leave unchanged

    chrom = fields[0]
    start_pos = int(fields[1])  # 1-based POS
    info_field = fields[7]

    # Parse END from INFO (required for your logic)
    # Adjust END by -1 (your VAMOS-specific fix)
    # and update INFO END accordingly
    try:
        end_str = info_field.split("END=")[1].split(";")[0]
        end_pos = int(end_str)
    except Exception:
        # If END isn't present/parsable, we cannot fetch reference properly; leave unchanged
        return row.rstrip("\n")

    # Apply your correction: END -= 1
    end_pos_adj = end_pos - 1
    # Replace only the first END occurrence safely
    info_field = info_field.replace(f"END={end_pos}", f"END={end_pos_adj}", 1)

    sample_field = fields[9]
    parts = sample_field.split(":")
    genotype = parts[0]
    alt_alleles = parts[1:]  # sequences from sample columns (may be empty)

    # Fetch reference allele (inclusive end; pysam fetch is 0-based start, end-exclusive)
    ref_allele = fetch_ref_allele(chrom, start_pos, end_pos_adj, fasta)

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
