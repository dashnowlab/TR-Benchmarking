# TR-Benchmarking

Tandem repeat benchmarking project.

## Run Nextflow

Submit the job with:

```bash
sbatch -J tr_bench   -p amilan --qos=normal   --time=08:00:00 --mem=2G   --output=logs/%x_%j.out --error=logs/%x_%j.err   --wrap="bash run_benchmarking_slurm.sh run_benchmarking.nf --list ont_bam.list --ref /pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta"
```

### What it does
- `-J tr_bench` — sets the Slurm job name.  
- `-p amilan --qos=normal` — selects the partition and QoS.  
- `--time=08:00:00 --mem=2G` — walltime and memory request.  
- `--output/--error` — logs saved to `logs/` with job id.  
- `--wrap="bash …"` — runs the wrapper that launches Nextflow:
  - `run_benchmarking_slurm.sh` invokes `run_benchmarking.nf`
  - `--list ont_bam.list` points to your input list (one BAM/CRAM per line)
  - `--ref ...fasta` provides the GRCh38 reference path


# Separate Shell Scripts

## Run_Atarva
```bash
cd /pl/active/dashnowlab/projects/TR-benchmarking/

while IFS=$'\t' read -r CRAM KARYOTYPE; do
  [[ -z "$CRAM" || "$CRAM" =~ ^[[:space:]]*# ]] && continue
  SAMPLE="$(basename "$CRAM" .haplotagged.cram)"
  SAMPLE="${SAMPLE%.cram}"
  sbatch -J "atarva_${SAMPLE}" run_atarva.sh "$CRAM" "$KARYOTYPE"
done < test_bam.list
```

## Run LongTR
```bash
cd /pl/active/dashnowlab/projects/TR-benchmarking/

while IFS=$'\t' read -r CRAM KARYOTYPE; do
  [[ -z "$CRAM" || "$CRAM" =~ ^[[:space:]]*# ]] && continue
  SAMPLE="$(basename "$CRAM" .haplotagged.cram)"
  SAMPLE="${SAMPLE%.cram}"
  sbatch -J "longtr_${SAMPLE}" run_longtr.sh "$CRAM" "$KARYOTYPE"
done < test_bam.list
```

## Run Straglr
```bash
cd /pl/active/dashnowlab/projects/TR-benchmarking/

while IFS=$'\t' read -r CRAM KARYOTYPE; do
  [[ -z "$CRAM" || "$CRAM" =~ ^[[:space:]]*# ]] && continue
  SAMPLE="$(basename "$CRAM" .haplotagged.cram)"
  SAMPLE="${SAMPLE%.cram}"
  sbatch -J "straglr_${SAMPLE}" run_straglr.sh "$CRAM" "$KARYOTYPE"
done < test_bam.list
```

## Run TRSV
```bash
cd /pl/active/dashnowlab/projects/TR-benchmarking/

while IFS=$'\t' read -r CRAM KARYOTYPE; do
  [[ -z "$CRAM" || "$CRAM" =~ ^[[:space:]]*# ]] && continue
  SAMPLE="$(basename "$CRAM" .haplotagged.cram)"
  SAMPLE="${SAMPLE%.cram}"
  sbatch -J "trsv_${SAMPLE}" run_trsv.sh "$CRAM" "$KARYOTYPE"
done < test_bam.list

```

