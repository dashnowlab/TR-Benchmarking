from cyvcf2 import VCF, Writer

def main(vcf: str, output: str):
    """
    Annotate TR VCF with pathogenicity using STRchive loci JSON.

    :param vcf: Input VCF file
    :param output: Output file in extended BED format
    """
    vcf_reader = VCF(vcf)

    with open(output, 'w') as outfile:
        # Write header
        outfile.write("#chrom\tstart\tend\tref_len\tsample\tallele1_seq\tallele2_seq\tallele1_len\tallele2_len\n")

        for record in vcf_reader:
            # Get VCF record coordinates
            chrom = record.CHROM
            pos = record.POS  # 1-based position

            # Calculate END position from REF length
            end_pos = pos + len(record.REF) - 1
            ref_len = len(record.REF)

            # Get all alleles sequences as a list e.g. [REF, ALT1, ALT2]
            # These can then be indexed by genotype indices e.g. 0 for REF, 1 for ALT1, etc.
            alleles = [record.REF, *record.ALT]

            # Iterate over samples in case of multi-sample VCF
            for sample_idx in range(len(vcf_reader.samples)):
                # gt = record.gt_bases[sample_idx] # alternative way to get genotype bases
                gt_idx = record.genotypes[sample_idx]
                gt = [alleles[gt_idx[0]], alleles[gt_idx[1]]] # get allele sequences as a list (typically two per sample for autosomal diploid)
                gt_len = [len(x) for x in gt]
                outfile.write(f"{chrom}\t{pos-1}\t{end_pos}\t{ref_len}\t{vcf_reader.samples[sample_idx]}\t{'\t'.join(gt)}\t{'\t'.join(map(str, gt_len))}\n")

if __name__ == '__main__':
    import defopt
    defopt.run(main)