# TR Benchmarking – Utility Scripts & Execution Guide

This document describes helper scripts and execution workflows for STR benchmarking analyses.

---

# 1️⃣ Header Fixer

Fix BAM/CRAM headers prior to tool execution.

```bash
sh header_fixer.sh \
  /pl/active/dashnowlab/projects/TR-benchmarking/ont_bam.list \
  /pl/active/dashnowlab/projects/TR-benchmarking/fixed_header_bam/
```

---

# 2️⃣ Catalog Builder

## Required BED Resources

Download the following files:

- Alliance Dark Genes Panel  
  https://www.twistbioscience.com/sites/default/files/resources/2023-02/Alliance_Dark_Genes_LR_Pnl_TargetsCaptured_hg38_ann.bed

- GRCh38 Segmental Duplications  
  https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/GRCh38@all/SegmentalDuplications/GRCh38_segdups.bed.gz

- STRchive Disease Loci  
  https://raw.githubusercontent.com/dashnowlab/STRchive/refs/heads/main/data/catalogs/STRchive-disease-loci.hg38.general.bed

## Build Catalog

```bash
bash -x build_catalogs.sh \
  "TR_catalog_for_vamos.shard_*_of_05.*.tsv.gz" \
  /pl/active/dashnowlab/projects/TR-benchmarking/catalogs/hell_catalog/bed_files_filtering/test/ \
  --include \
  /pl/active/dashnowlab/projects/TR-benchmarking/catalogs/hell_catalog/bed_files_filtering/Alliance_Dark_Genes_LR_Pnl_TargetsCaptured_hg38_ann.bed \
  /pl/active/dashnowlab/projects/TR-benchmarking/catalogs/hell_catalog/bed_files_filtering/GRCh38_segdups.bed
```

---

# 3️⃣ Running STR Genotyping Tools

All tools use the same sample list format:

```
CRAM_PATH<TAB>KARYOTYPE
```

Example list file: `test_bam.list`

Move to project directory:

```bash
cd /pl/active/dashnowlab/projects/TR-benchmarking/
```

---

## 🔹 Run Atarva

```bash
while IFS=$'\t' read -r CRAM KARYOTYPE; do
  [[ -z "$CRAM" || "$CRAM" =~ ^[[:space:]]*# ]] && continue
  SAMPLE="$(basename "$CRAM" .haplotagged.cram)"
  SAMPLE="${SAMPLE%.cram}"
  sbatch -J "atarva_${SAMPLE}" run_atarva.sh "$CRAM" "$KARYOTYPE"
done < test_bam.list
```

---

## 🔹 Run LongTR

```bash
while IFS=$'\t' read -r CRAM KARYOTYPE; do
  [[ -z "$CRAM" || "$CRAM" =~ ^[[:space:]]*# ]] && continue
  SAMPLE="$(basename "$CRAM" .haplotagged.cram)"
  SAMPLE="${SAMPLE%.cram}"
  sbatch -J "longtr_${SAMPLE}" run_longtr.sh "$CRAM" "$KARYOTYPE"
done < test_bam.list
```

---

## 🔹 Run Straglr

```bash
while IFS=$'\t' read -r CRAM KARYOTYPE; do
  [[ -z "$CRAM" || "$CRAM" =~ ^[[:space:]]*# ]] && continue
  SAMPLE="$(basename "$CRAM" .haplotagged.cram)"
  SAMPLE="${SAMPLE%.cram}"
  sbatch -J "straglr_${SAMPLE}" run_straglr.sh "$CRAM" "$KARYOTYPE"
done < test_bam.list
```

---

## 🔹 Run TRSV

```bash
while IFS=$'\t' read -r CRAM KARYOTYPE; do
  [[ -z "$CRAM" || "$CRAM" =~ ^[[:space:]]*# ]] && continue
  SAMPLE="$(basename "$CRAM" .haplotagged.cram)"
  SAMPLE="${SAMPLE%.cram}"
  sbatch -J "trsv_${SAMPLE}" run_trsv.sh "$CRAM" "$KARYOTYPE"
done < test_bam.list
```

---

# 4️⃣ Vamos VCF Fixing

Reference genome:

```
/pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta
```

---

## Generic One-Liner (Run Inside Directory)

```bash
for f in *.vamos.vcf; do
  s="${f%.vamos.vcf}"
  python3 fix-vcf.py "$f" "${s}.vamos.fixed.vcf" \
    --ref /pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta
done
```

---

## HG Benchmark Fix

```bash
for f in /pl/active/dashnowlab/projects/TR-benchmarking/benchmark-catalog-v2-HG/vamos/*.vamos.vcf; do
  s="${f%.vamos.vcf}"
  python3 scripts/fix-vcf.py "$f" "${s}.vamos.fixed.vcf" \
    --ref /pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta
done
```

---

## Deveson

```bash
for f in /pl/active/dashnowlab/projects/TR-benchmarking/benchmark-catalog-v2-Deveson/vamos/*URfixed.vamos.vcf; do
  s="${f%.vamos.vcf}"
  python3 scripts/fix-vcf.py "$f" "${s}.vamos.fixed.vcf" \
    --ref /pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta
done
```

---

## Deveson Default

```bash
for f in /pl/active/dashnowlab/projects/TR-benchmarking/benchmark-catalog-V2-Deveson-default/vamos/*URfixed.vamos.vcf; do
  s="${f%.vamos.vcf}"
  python3 scripts/fix-vcf.py "$f" "${s}.vamos.fixed.vcf" \
    --ref /pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta
done
```

---

## HPRC

```bash
for f in /pl/active/dashnowlab/projects/TR-benchmarking/benchmark-catalog-v2-HPRC/vamos/*.vamos.vcf; do
  s="${f%.vamos.vcf}"
  python3 scripts/fix-vcf.py "$f" "${s}.vamos.fixed.vcf" \
    --ref /pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta
done
```

---

# 5️⃣ Mendelian Consistency Calculations

## STRdust

```bash
python3 mendelian_consistency_calc_ref_alt_debug.py \
  --mom ../benchmark-catalog-v2-HG/strdust/HG003.30x.haplotagged.strdust.vcf \
  --dad ../benchmark-catalog-v2-HG/strdust/HG004.30x.haplotagged.strdust.vcf \
  --kids ../benchmark-catalog-v2-HG/strdust/HG002.30x.haplotagged.strdust.vcf \
  --out-prefix ../mendelian_consistency_calc/strdust-
```

---

## Atarva

```bash
python3 mendelian_consistency_calc_ref_alt_debug.py \
  --mom ../benchmark-catalog-v2-HG/atarva/HG003.30x.haplotagged.atarva.vcf \
  --dad ../benchmark-catalog-v2-HG/atarva/HG004.30x.haplotagged.atarva.vcf \
  --kids ../benchmark-catalog-v2-HG/atarva/HG002.30x.haplotagged.atarva.vcf \
  --out-prefix ../mendelian_consistency_calc/atarva-
```

---

## LongTR

```bash
python3 mendelian_consistency_calc_ref_alt_debug.py \
  --mom ../benchmark-catalog-v2-HG/longtr/HG003.30x.haplotagged.longTR.vcf.gz \
  --dad ../benchmark-catalog-v2-HG/longtr/HG004.30x.haplotagged.longTR.vcf.gz \
  --kids ../benchmark-catalog-v2-HG/longtr/HG002.30x.haplotagged.longTR.vcf.gz \
  --out-prefix ../mendelian_consistency_calc/longtr- \
  --merge-key id
```

---

## STRkit

```bash
python3 mendelian_consistency_calc_ref_alt_debug.py \
  --mom ../benchmark-catalog-v2-HG/strkit/HG003.30x.haplotagged.strkit-min_read_1.vcf \
  --dad ../benchmark-catalog-v2-HG/strkit/HG004.30x.haplotagged.strkit-min_read_1.vcf \
  --kids ../benchmark-catalog-v2-HG/strkit/HG002.30x.haplotagged.strkit-min_read_1.vcf \
  --out-prefix ../mendelian_consistency_calc/strkit-
```

---

## Vamos

```bash
python3 mendelian_consistency_calc_ref_alt_debug.py \
  --mom ../benchmark-catalog-v2-HG/vamos/HG003.30x.haplotagged.vamos.fixed.vcf \
  --dad ../benchmark-catalog-v2-HG/vamos/HG004.30x.haplotagged.vamos.fixed.vcf \
  --kids ../benchmark-catalog-v2-HG/vamos/HG002.30x.haplotagged.vamos.fixed.vcf \
  --out-prefix ../mendelian_consistency_calc/vamos-
```

---

## Medaka

```bash
python3 mendelian_consistency_calc_ref_alt_debug.py \
  --mom ../benchmark-catalog-v2-HG/medaka/HG003.30x.haplotagged.medaka.vcf/medaka_to_ref.TR.vcf \
  --dad ../benchmark-catalog-v2-HG/medaka/HG004.30x.haplotagged.medaka.vcf/medaka_to_ref.TR.vcf \
  --kids /pl/active/dashnowlab/projects/TR-benchmarking/benchmark-catalog-v2-HG/medaka/HG002.30x.haplotagged.medaka.vcf/medaka_to_ref.TR.vcf \
  --out-prefix ../mendelian_consistency_calc/medaka-
```

---

## Straglr (Special Script)

```bash
python3 mendelian_consistency_calc_ref_alt_straglr.py \
  --mom ../benchmark-catalog-v2-HG/straglr/HG003.30x.haplotagged.vcf \
  --dad ../benchmark-catalog-v2-HG/straglr/HG004.30x.haplotagged.vcf \
  --kids ../benchmark-catalog-v2-HG/straglr/HG002.30x.haplotagged.vcf \
  --out-prefix ../mendelian_consistency_calc/straglr-
```

---
