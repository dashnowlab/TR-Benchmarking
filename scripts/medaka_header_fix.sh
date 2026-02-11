#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "Usage: $0 <cram_list.txt> <reference.fasta> <out_dir> [partition] [time] [mem]"
  echo "  cram_list.txt: one CRAM path per line (can include blank lines and # comments)"
  echo "  reference.fasta: path to reference fasta (should have .fai рядом)"
  echo "  out_dir: where new CRAMs will be written"
  echo "  partition (optional): default=amilan"
  echo "  time (optional): default=02:00:00"
  echo "  mem (optional): default=16G"
  exit 1
fi

CRAM_LIST="$1"
REF_FASTA="$2"
OUT_DIR="$3"
PARTITION="${4:-amilan}"
TIME="${5:-02:00:00}"
MEM="${6:-16G}"

mkdir -p "$OUT_DIR"

# Basic checks
[[ -f "$CRAM_LIST" ]] || { echo "ERROR: list file not found: $CRAM_LIST"; exit 1; }
[[ -f "$REF_FASTA" ]] || { echo "ERROR: reference fasta not found: $REF_FASTA"; exit 1; }
[[ -f "${REF_FASTA}.fai" ]] || { echo "ERROR: reference fai not found: ${REF_FASTA}.fai (run: samtools faidx $REF_FASTA)"; exit 1; }

# Submit one job per CRAM
while IFS= read -r cram || [[ -n "$cram" ]]; do
  # skip blanks/comments
  [[ -z "$cram" ]] && continue
  [[ "$cram" =~ ^[[:space:]]*# ]] && continue

  if [[ ! -f "$cram" ]]; then
    echo "WARN: CRAM not found, skipping: $cram"
    continue
  fi

  base="$(basename "$cram")"
  base="${base%.cram}"
  out_cram="${OUT_DIR}/${base}.URfixed.cram"
  log="${OUT_DIR}/${base}.URfixed.%j.log"

  sbatch \
    -p "$PARTITION" \
    -t "$TIME" \
    --mem "$MEM" \
    --qos=normal \
    -c 1 \
    -J "URfix_${base}" \
    -o "$log" \
    --export=ALL,IN_CRAM="$cram",REF_FASTA="$REF_FASTA",OUT_CRAM="$out_cram",OUT_DIR="$OUT_DIR" \
<<'SBATCH_SCRIPT'
#!/bin/bash -l
set -euo pipefail

# Ensure module command exists in batch context
if ! command -v module >/dev/null 2>&1; then
  [[ -f /etc/profile.d/lmod.sh ]] && source /etc/profile.d/lmod.sh
  [[ -f /etc/profile.d/modules.sh ]] && source /etc/profile.d/modules.sh
  [[ -f /usr/share/lmod/lmod/init/bash ]] && source /usr/share/lmod/lmod/init/bash
fi

module purge
module load samtools/1.16.1
samtools --version

# echo "IN_CRAM=$IN_CRAM"
# echo "REF_FASTA=$REF_FASTA"
# echo "OUT_DIR=$OUT_DIR"
# echo "OUT_CRAM=$OUT_CRAM"
# echo "HOST=$(hostname)"
# echo "START=$(date -Is)"

# # temp header (INSIDE OUT_DIR; unique per job)
# umask 007
# tmpdir="$(mktemp -d -p "$OUT_DIR" "${SLURM_JOB_ID:-manual}.tmp.URfix.XXXXXXXX")"
# trap 'rm -rf "$tmpdir"' EXIT

# hdr="$tmpdir/header.sam"
# hdr_fixed="$tmpdir/header.fixed.sam"

# # Extract header
# samtools view -H "$IN_CRAM" > "$hdr"

# awk -v ref="$REF_FASTA" 'BEGIN{OFS="\t"}
#   /^@/{
#     for(i=1;i<=NF;i++){
#       if($i ~ /^UR:/){ $i="UR:" ref }
#     }
#   }
#   {print}
# ' "$hdr" > "$hdr_fixed"

# # Reheader -> new CRAM (use reference for CRAM writing)
# samtools reheader "$hdr_fixed" "$IN_CRAM" \
#   | samtools view -C -T "$REF_FASTA" -o "$OUT_CRAM" -

samtools index "$OUT_CRAM"

echo "DONE=$(date -Is)"
SBATCH_SCRIPT

done < "$CRAM_LIST"

echo "Submitted jobs. Output will be in: $OUT_DIR"
