#!/bin/bash

#SBATCH --nodes=4
#SBATCH --ntasks=1
#SBATCH --time=7-00:00:00
#SBATCH --partition=amilan
#SBATCH --qos long
#SBATCH -o ./nf-runner.stdout
#SBATCH -e ./nf-runner.stderr
#SBATCH -J nf-runner
#SBATCH --mail-type END
#SBATCH --mail-type FAIL
#SBATCH --mail-user=akshaykumar.avvaru@cuanschutz.edu

module purge

module load nextflow
module load singularity/3.7.4
module load miniforge/24.11.3-0

# Base directories all tied to out_dir
BASE="/pl/active/dashnowlab/projects/TR-benchmarking/"
WORK_DIR="${BASE}/work_akshay"
ASSETS_DIR="${WORK_DIR}/assets"          # controls projectDir via NXF_ASSETS
NXFHOME_DIR="${BASE}/.nextflow" # controls Nextflow home/cache

# Ensure dirs exist
mkdir -p "$WORK_DIR" "$ASSETS_DIR" "$NXFHOME_DIR" "$WORK_DIR/tmp" "$WORK_DIR/cache" "$WORK_DIR/logs"

# Environment for Singularity & Nextflow
export TMPDIR="${WORK_DIR}/tmp"
export NXF_WORK="${WORK_DIR}"     # workDir
export NXF_ASSETS="${ASSETS_DIR}" # projectDir will be under this
export NXF_HOME="${NXFHOME_DIR}"  # Nextflow home + session/cache/lock files
export SINGULARITY_CACHEDIR=${WORK_DIR}/cache
export NXF_SINGULARITY_CACHEDIR=${WORK_DIR}/cache
export NXF_CONDA_CACHEDIR=${WORK_DIR}/conda_cache

echo "Running sample: $SAMPLE_NAME"
echo "launchDir  -> $BASE"
echo "workDir    -> $WORK_DIR"
echo "projectDir -> (under) $ASSETS_DIR"

echo "== Starting nextflow =="
nextflow run "$@" -profile singularity -w "$WORK_DIR"
echo "== Nextflow complete =="

# Example usage: sbatch run_benchmarking_slurm.sh run_benchmarking.nf --list test_bam.list --ref /pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta -with-overwrite
# Usage with sbatch:
# sbatch -J tr_bench -p amilan --qos=normal --time=08:00:00 --mem=2G \
#        --output=logs/%x_%j.out --error=logs/%x_%j.err \
#         --wrap="bash run_benchmarking_akshay.sh run_benchmarking_akshay.nf --list ont_bam.list --ref /pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta"