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
#SBATCH --mail-user=elbay.aliyev@cuanschutz.edu

module purge

module load nextflow
module load singularity/3.7.4
module load miniforge/24.11.3-0

# Base directories all tied to out_dir
BASE="/pl/active/dashnowlab/projects/TR-benchmarking/"
WORK_DIR="${BASE}/work_elbay"

# Ensure dirs exist
#mkdir -p "$WORK_DIR" "$ASSETS_DIR" "$NXFHOME_DIR" "$WORK_DIR/tmp" "$WORK_DIR/cache" "$WORK_DIR/logs"


echo "== Starting nextflow =="
nextflow run "$@" -profile singularity -w "$WORK_DIR" -resume
echo "== Nextflow complete =="

# Example usage: sbatch run_benchmarking_slurm.sh run_benchmarking.nf --list test_bam.list --ref /pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta -with-overwrite