# TR Benchmarking – Utility Scripts & Execution Guide

This document describes helper scripts and execution workflows for STR benchmarking analyses.

---

## 1. Header Fixer

Fix BAM/CRAM headers before running STR genotyping tools.

```bash
sh header_fixer.sh \
  /pl/active/dashnowlab/projects/TR-benchmarking/ont_bam.list \
  /pl/active/dashnowlab/projects/TR-benchmarking/fixed_header_bam/

2. Mendelian Consistency Calculations

The following commands calculate Mendelian consistency using trio VCFs:

Mother: HG003
Father: HG004
Child: HG002

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