import argparse
import sys
import glob
import gzip

def main():
    parser = argparse.ArgumentParser(description="Extract pathogenic loci from VCF files.")
    parser.add_argument("vcf", nargs='+', help="Path to the input VCF file(s).")
    parser.add_argument("bed", help="Path to the input BED file of pathogenic loci with locus ID in the 4th column.")
    parser.add_argument("--metadata", help="Path to the output metadata file (optional).")
    args = parser.parse_args()

    #args.vcf = ["NA13515.atarva.vcf"] # for testing
    # all vcfs in the directory /Users/dashnowh/alpine/dashnowlab/projects/TR-benchmarking/benchmark-catalog-v2-Deveson/atarva
    # args.vcf = glob.glob("/Users/dashnowh/alpine/dashnowlab/projects/TR-benchmarking/benchmark-catalog-v2-Deveson/atarva/*.vcf")
    # Files ending in .vcf or .vcf.gz
    #args.vcf = glob.glob("/Users/dashnowh/alpine/dashnowlab/projects/TR-benchmarking/benchmark-catalog-v2-Deveson/*/*.vcf")
    args.vcf = glob.glob("/Users/dashnowh/alpine/dashnowlab/projects/TR-benchmarking/benchmark-catalog-v2-Deveson/*/*.vcf.gz")

    args.bed = "STRchive-disease-loci-v2.15.0.hg38.general.bed" # for testing
    args.metadata = "metadata_pathogenic_sciadv.abm5386.tsv" # for testing

    pathogenic_loci = parse_bed_file(args.bed)
    metadata = parse_metadata_file(args.metadata) if args.metadata else None

    with open("pathogenic_results.tsv", "w") as out:
        out.write("\t".join(["tool", "sample", "gene", "molec_type", "molec_prec_allele1", "molec_prec_allele2", "molec_allele1_len", "molec_allele2_len", "allele1_len", "allele2_len", "allele1_seq", "allele2_seq"]) + "\n")
        for vcf in args.vcf:
            # Get sample ID from the start of the filename and tool from the directory name
            sample_id = vcf.split("/")[-1].split(".")[0]  # Assuming sample ID is the first part of the filename
            tool = vcf.split("/")[-2]  # Assuming tool is the parent directory name
            # if tool != 'atarva' or sample_id != 'NA06890':
            #     continue
            print(f"Processing {vcf} for sample {sample_id}, tool {tool}...")
            # check if sample_id is in metadata  
            if metadata and sample_id in metadata:
                try:
                    loci = filter_vcf_by_bed(vcf, pathogenic_loci)
                except Exception as e:
                    if tool == 'medaka': # treat VCF as a directory and look inside for a file medaka_to_ref.TR.vcf
                        loci = filter_vcf_by_bed(vcf + "/medaka_to_ref.TR.vcf", pathogenic_loci)
                    else:
                        print(f"Error processing VCF {vcf} for sample {sample_id}, skipping: {e}")
                        continue
                # Filter to loci annotated with the gene in metadata
                for locus in loci:
                    if "LOCID=" in locus:
                        loc_id = locus.split("LOCID=")[1].split()[0]
                        #print(f"Found locus ID -{loc_id}- in VCF for sample {sample_id}")
                        #print(f"Known expansion -{metadata[sample_id].get('Known expansion')}-")
                        if loc_id == metadata[sample_id].get("Known expansion"):
                            # Check which kind of molecular testing was done
                            if metadata[sample_id].get("Southern") == "Y":
                                molec_type = "Southern"
                                if metadata[sample_id].get("RP-PCR") == "Y":
                                    molec_type += " and RP-PCR"
                            elif metadata[sample_id].get("RP-PCR") == "Y":
                                molec_type = "RP-PCR"
                            else:
                                molec_type = "Unknown"
                            # Get alleles using GT
                            if tool == "straglr":
                                gt_allele_lens = parse_straglr_vcf(locus)
                                gt_alleles = ['NA', 'NA']  
                            else:
                                locus_fields = locus.split("\t")
                                
                                gt_field_index = locus_fields[8].split(":").index("GT")
                                sample_field = locus_fields[9]
                                # split on / or | and convert to numeric genotype
                                gt_field = sample_field.split(":")[gt_field_index]
                                if gt_field == "." or gt_field == "./." or gt_field == ".|.":
                                    gt_alleles = ['NA', 'NA']
                                else:
                                    gt_indices_raw = gt_field.replace("|", "/").split("/")
                                    gt_indices = [int(x) for x in gt_indices_raw if x.isdigit()]
                                    alleles = [locus_fields[3]] + locus_fields[4].split(",")  # REF + ALT alleles
                                    gt_alleles = [alleles[i] for i in gt_indices if i is not None]
                                    gt_allele_lens = [len(allele) for allele in gt_alleles if allele is not None]
                                    if len(gt_indices) < 2:
                                        gt_alleles.append('NA')  # If one allele is missing, add NA to keep the length consistent
                                        gt_allele_lens.append('NA')
                            # Sort alleles and allele lengths by length (shortest first)
                            molec_alleles = [metadata[sample_id].get("Molecular_allele1_bp"), metadata[sample_id].get("Molecular_allele2_bp")]
                            molec_alleles_sorted = sorted(molec_alleles, key=lambda x: int(x) if x.isdigit() else float('inf'))  # Sort molecular alleles by length, treating non-numeric as longest
                            # If order of molec_alleles is changed, also change the order of molec_prec_allele in the output
                            if molec_alleles_sorted != molec_alleles:
                                molec_prec_allele1 = metadata[sample_id].get("Molecular_allele2_type")
                                molec_prec_allele2 = metadata[sample_id].get("Molecular_allele1_type")
                            else:
                                molec_prec_allele1 = metadata[sample_id].get("Molecular_allele1_type")
                                molec_prec_allele2 = metadata[sample_id].get("Molecular_allele2_type")
                            gt_allele_lens_sorted = sorted(gt_allele_lens, key=lambda x: int(x) if x != 'NA' else float('inf'))  # Sort GT allele lengths by length, treating 'NA' as longest
                            gt_alleles_sorted = sorted(gt_alleles, key=lambda x: len(x) if x != 'NA' else float('inf'))  # Sort GT alleles by length, treating 'NA' as longest

                            print(f"Sample {sample_id} has genotype {gt_allele_lens} at locus {loc_id} ")
                            out.write("\t".join([tool, sample_id, loc_id, molec_type,
                                            molec_prec_allele1, molec_prec_allele2] +
                                            molec_alleles_sorted +
                                            [str(allele_len) for allele_len in gt_allele_lens_sorted] + gt_alleles_sorted) + "\n")

def parse_straglr_vcf(line):
    # Assumptions:
    # SVLEN is the allele length in the reference in bp
    # RN is the alt allele length in number of motifs
    # RB is the length of the alt allele in bp
    # RUL_REF is the reference repeat unit length in bp
    # RUL is the alt repeat unit length in bp
    # A genotype of 0 corresponds to 0/0 and 1 corresponds to 1/1
    # Example
    # chr1    181146  .       C       <CNV:TR>        .       PASS    END=181463;RUL_REF=29;SVLEN=194;RN=1;RUL=28;RB=381;CIRB=-11,11 GT:DP:AS        1:8:2
    fields = line.strip().split("\t")
    #print(fields)
    GT = fields[9].split(":")[0]
    if len(GT) == 1:
        GT = GT + "/" + GT
    gt_indices = [int(x) for x in GT.replace("|", "/").split("/")]
    # Extract ref len from pos and END
    pos = int(fields[1])
    end = int(fields[7].split("END=")[1].split(";")[0])
    ref_len = end - pos + 1
    # if all GT indices are 0
    if all(i == 0 for i in gt_indices):
        return [ref_len, ref_len]
    alt_lens = [int(x) for x in fields[7].split("RB=")[1].split(";")[0].split(",")]
    allele_lens = [ref_len] + alt_lens
    gt_allele_lens = [allele_lens[i] for i in gt_indices]
    return gt_allele_lens





def parse_metadata_file(metadata_path):
    keep_fields = ["Sample ID", "Known expansion", "Molecular Testing", "Southern", "RP-PCR", "Molecular_allele1_bp", "Molecular_allele2_bp", "Molecular_allele1_type", "Molecular_allele2_type", "Motif_size"]
    metadata = {}
    with open(metadata_path, 'r') as meta:
        header = meta.readline().strip().split("\t")
        field_indices = {field: header.index(field) for field in keep_fields if field in header}
        for line in meta:
            fields = line.strip().split("\t")
            if len(fields) < len(header):
                continue
            sample_id = fields[field_indices["Sample ID"]]
            metadata[sample_id] = {field: fields[field_indices[field]] for field in keep_fields if field in field_indices}
    return metadata

def parse_bed_file(bed_path, gene_col=4):
    pathogenic_loci = []
    with open(bed_path, 'r') as bed:
        for line in bed:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            locus_id = fields[gene_col] if len(fields) > gene_col else None
            pathogenic_loci.append((chrom, start, end, locus_id))
    return pathogenic_loci

def filter_vcf_by_bed(vcf_path, pathogenic_loci):
    filtered_variants = []

    # Check if the file is gzipped    
    if vcf_path.endswith(".gz"):
        open_func = gzip.open
    else:
        open_func = open

    # Process VCF file and filter variants
    with open_func(vcf_path, 'rt') as vcf:
        for line in vcf:
            if line.startswith("#"):
                #print(line.strip())
                continue
            fields = line.strip().split("\t")
            chrom = fields[0]
            pos = int(fields[1])
            print(chrom, pos)

            # Check if ref and alt alleles are present and not symbolic
            if fields[3] == "N" or fields[4] == "<VNTR>":
                raise ValueError(f"Unexpected REF or ALT allele in VCF at {chrom}:{pos}. Expected actual sequences, got REF={fields[3]} and ALT={fields[4]}.")

            # Check if the variant falls within any pathogenic locus
            for locus in pathogenic_loci:
                if chrom == locus[0] and locus[1] -1 <= pos <= locus[2] + 1:  # BED is 0-based, VCF is 1-based + adding some padding
                    # annotate with locus ID if available
                    if locus[3]:
                        fields[7] += f";LOCID={locus[3]}"
                    #print("\t".join(fields))
                    filtered_variants.append("\t".join(fields))

    return filtered_variants

if __name__ == "__main__":
    main()