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

# --- parse -w / --work-dir from nextflow args and set NXF_HOME based on it ---
WORKDIR=""
args=("$@")
for ((i=0; i<${#args[@]}; i++)); do
  case "${args[$i]}" in
    -w|--work-dir)   WORKDIR="${args[$((i+1))]}";;
    --work-dir=*)    WORKDIR="${args[$i]#--work-dir=}";;
    -w=*)            WORKDIR="${args[$i]#-w=}";;
    -w/*)            WORKDIR="${args[$i]#-w}";;
  esac
done

echo "${WORKDIR}"

export NXF_HOME="${WORKDIR}/.nextflow"
export NXF_CACHE_DIR="${WORKDIR}/.nextflow"


echo "== Starting nextflow =="
nextflow run "$@" -profile singularity
echo "== Nextflow complete =="

# Example usage: sbatch run_benchmarking_slurm.sh run_benchmarking.nf --list test_bam.list --ref /pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta -with-overwrite