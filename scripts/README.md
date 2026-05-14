# TR Benchmarking – Utility Scripts & Execution Guide

This document describes helper scripts and execution workflows for STR benchmarking analyses.

---

## 1. Header Fixer

Fix BAM/CRAM headers before running STR genotyping tools.

```bash
sh header_fixer.sh \
  /pl/active/dashnowlab/projects/TR-benchmarking/ont_bam.list \
  /pl/active/dashnowlab/projects/TR-benchmarking/fixed_header_bam/
```
## 2. Mendelian Consistency Calculations

The following commands calculate Mendelian consistency using trio VCFs:

Mother: HG003
Father: HG004
Child: HG002

```bash
python3 mendelian_consistency_calc_ref_alt_debug_test.py \
  --mom ../benchmark-catalog-v2-HG/strdust/HG003.30x.haplotagged.strdust.vcf \
  --dad ../benchmark-catalog-v2-HG/strdust/HG004.30x.haplotagged.strdust.vcf \
  --kids ../benchmark-catalog-v2-HG/strdust/HG002.30x.haplotagged.strdust.vcf \
  --out-prefix ../mendelian_consistency_calc/test/strdust-

  python3 mendelian_consistency_calc_ref_alt_debug_test.py \
  --mom ../benchmark-catalog-v2-HG/atarva/HG003.30x.haplotagged.atarva.vcf \
  --dad ../benchmark-catalog-v2-HG/atarva/HG004.30x.haplotagged.atarva.vcf \
  --kids ../benchmark-catalog-v2-HG/atarva/HG002.30x.haplotagged.atarva.vcf \
  --out-prefix ../mendelian_consistency_calc/test/atarva-

  python3 mendelian_consistency_calc_ref_alt_debug_test.py \
  --mom ../benchmark-catalog-v2-HG/longtr/HG003.30x.haplotagged.longTR.vcf.gz \
  --dad ../benchmark-catalog-v2-HG/longtr/HG004.30x.haplotagged.longTR.vcf.gz \
  --kids ../benchmark-catalog-v2-HG/longtr/HG002.30x.haplotagged.longTR.vcf.gz \
  --out-prefix ../mendelian_consistency_calc/test/longtr- \
  --merge-key id

  python3 mendelian_consistency_calc_ref_alt_debug_test.py \
  --mom ../benchmark-catalog-v2-HG/strkit/HG003.30x.haplotagged.strkit-min_read_1.vcf \
  --dad ../benchmark-catalog-v2-HG/strkit/HG004.30x.haplotagged.strkit-min_read_1.vcf \
  --kids ../benchmark-catalog-v2-HG/strkit/HG002.30x.haplotagged.strkit-min_read_1.vcf \
  --out-prefix ../mendelian_consistency_calc/test/strkit-

  python3 mendelian_consistency_calc_ref_alt_debug_test.py \
  --mom ../benchmark-catalog-v2-HG/vamos/HG003.30x.haplotagged.vamos.fixed.vcf \
  --dad ../benchmark-catalog-v2-HG/vamos/HG004.30x.haplotagged.vamos.fixed.vcf \
  --kids ../benchmark-catalog-v2-HG/vamos/HG002.30x.haplotagged.vamos.fixed.vcf \
  --out-prefix ../mendelian_consistency_calc/test/vamos-

  python3 mendelian_consistency_calc_ref_alt_debug_test.py \
  --mom ../benchmark-catalog-v2-HG/medaka/HG003.30x.haplotagged.medaka.vcf/medaka_to_ref.TR.vcf \
  --dad ../benchmark-catalog-v2-HG/medaka/HG004.30x.haplotagged.medaka.vcf/medaka_to_ref.TR.vcf \
  --kids /pl/active/dashnowlab/projects/TR-benchmarking/benchmark-catalog-v2-HG/medaka/HG002.30x.haplotagged.medaka.vcf/medaka_to_ref.TR.vcf \
  --out-prefix ../mendelian_consistency_calc/test/medaka-

  python3 mendelian_consistency_calc_ref_alt_straglr_test.py \
  --mom ../benchmark-catalog-v2-HG/straglr/HG003.30x.haplotagged.vcf \
  --dad ../benchmark-catalog-v2-HG/straglr/HG004.30x.haplotagged.vcf \
  --kids ../benchmark-catalog-v2-HG/straglr/HG002.30x.haplotagged.vcf \
  --out-prefix ../mendelian_consistency_calc/test/straglr-
  ```

## Extract pathogenic loci from VCFs

Example usage:

```bash
python extract_pathogenic.py --vcfs *.vcf --bed STRchive-disease-loci-v2.15.0.hg38.general.bed --metadata metadata_pathogenic_sciadv.abm5386.tsv --output pathogenic_results.tsv 
```
## 3. Assembly Concordance Workflow

Steps followed for each sample:
1. Download ONT WGS data from HPRC and align using minimap2
2. Download maternal and paternal assemblies from HPRC
3. Align maternal and paternal assemblied to reference GRCh38 Genome using minmap2
4. Genotype the the aligned ONT WGS data using all the seven genotypers
5. Process each tool VCF using `assembly_concordance.py` to generate TSV result file
tool_len	asm_len	len_diff	tool_lev	asm_lev	tool	status

| Column | Description |
|---|---|
| `sample_id` | Sample identifier provided. |
| `catalog_locus` | Locus key {chrom}-{start}-{end} as in the catalog. |
| `motif` | Motif of the TR locus. |
| `motif_length` | Length of the motif. |
| `ref_len` | Length of the TR locus in the reference genome. |
| `vcf_locus` | Locus key {chrom}-{start}-{end} as in the VCF. |
| `haplogroup` | Allele number. Either 1 or 2 for each allele of the TR locus. |
| `asm_vtype` | Variant type as found in the assembly (REF-HOM/REF-HOM/ALT-HOM/ALT-HET). |
| `tool_vtype` | Variant type as found in the tool VCF (REF-HOM/REF-HOM/ALT-HOM/ALT-HET). |
| `tool_len` | Length of the allele as identified by the tool. |
| `asm_len` | Length of the allele as identified by the assembly. |
| `len_diff` | Difference in length between the tool and assembly alleles. |
| `tool_lev` | Levenshtein distance between the tool and assembly alleles. |
| `asm_lev` | Levenshtein distance between the assembly and reference alleles. |
| `status` | Status of the concordance between the tool and assembly alleles (MATCH/ONE-OFF/MOTIF-OFF/MISSING). |

```
Example command for one sample (HG002) and one tool (ATARVA):

python ../../TR-Benchmarking/scripts/assembly_concordance.py -v HG002.atarva.vcf -k XY -a hg002_mat_asm5_hg38.bam hg002_pat_asm5_hg38.bam -bed bechmark_catalog.bed -t atarva -r human_GRCh38_no_alt_analysis_set.fasta -s HG002 -o HG002.atarva.asmcon
```

These assembly concordance files are then processed using R scripts to generate summary statistics and visualizations for each tool and sample. The R scripts can be found in the `visualization` directory.
