# TR-Benchmarking

Tandem repeat benchmarking project.

## Run Nextflow

Submit the job with:

```bash
sbatch -J tr_bench_hg -p amilan --qos=normal --time=23:00:00 --mem=2G --output=logs/%x_%j.out --error=logs/%x_%j.err --wrap="bash run_benchmarking_slurm.sh run_benchmarking_hg.nf --list ont_bam.list --ref /pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta -w /pl/active/dashnowlab/projects/TR-benchmarking/benchmark-catalog-v2-HG/work_HG -resume"

sbatch -J tr_bench_hprc -p amilan --qos=long --time=23:00:00 --mem=2G --output=logs/%x_%j.out --error=logs/%x_%j.err --wrap="bash run_benchmarking_slurm.sh run_benchmarking_hprc.nf --list hprc_bam.list --ref /pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta -w /pl/active/dashnowlab/projects/TR-benchmarking/benchmark-catalog-v2-HPRC/work_HPRC -resume"

sbatch -J tr_bench_deveson -p amilan --qos=normal --time=23:00:00 --mem=2G --output=logs/%x_%j.out --error=logs/%x_%j.err --wrap="bash run_benchmarking_slurm.sh run_benchmarking_deveson.nf --list deveson_bam.list --ref /pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta -w /pl/active/dashnowlab/projects/TR-benchmarking/benchmark-catalog-v2-Deveson/work_Deveson -resume"
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




