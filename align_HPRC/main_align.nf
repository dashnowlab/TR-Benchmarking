#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.help = false
if (params.help) {
    log.info """
    Usage:
        nextflow run main.nf --ref <reference genome> --meta <metadata CSV file>

    Options:
        --meta          Sample metadata CSV file
        --ref           Reference genome
        --help          Show this help message
    """
    System.exit(0)
}

// Note: if process.scratch = true is set in nextflow.config. This will cause intermediate files to be deleted after process completion.

// Defaults, can be overridden by setting params in command line
// def sample_metadata = params.meta
def ref = params.ref

// Params defaults are defined in nextflow.config (centralized)

workflow {

    // Get the unaligned BAM files
    // If a metadata CSV is provided (or the default sample_metadata URL), parse it and download files listed in the `path` column.
    // This pipeline expects `params.unaligned_bam` values (no filesystem glob fallback is used).

    // Channel.value(params.meta ?: sample_metadata).set { meta_input }
    Channel.value(params.ref ?: ref).set { ref_input }

    // Check reference genome
    params.ref = '/home/hdashnow@xsede.org/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta'
    ref = file(params.ref, checkIfExists: true)
    fai = file("${params.ref}.fai", checkIfExists: true)

    // fetch metadata and produce tuples (sample_ID, coverage, path, readgroup)
    fetch_meta()

    // Filter metadata inside a process so `filtered_metadata.csv` is produced and published.
    filter_meta(fetch_meta.out)

    // Read the filtered CSV produced by the process and convert to tuples for downstream steps
    filter_meta.out
        // take the first element (path to filtered_metadata.csv)
        .map { t -> t[0] }
        .flatMap { metaFile ->
            def tuples = []
            metaFile.eachLine { line, lineNum ->
                if (lineNum == 1) return // skip header
                // split into at most 3 parts: sample_id, coverage, path, readgroup
                def parts = line.split(',', 4)
                if (parts.size() >= 4) {
                    def sample_id = parts[0].trim()
                    def coverage = parts[1].trim()
                    def url = parts[2].trim()
                    def readgroup = parts[3].trim()
                    // only emit when we have both id and url
                    if (sample_id && url) tuples << tuple(sample_id, coverage, url, readgroup)
                }
            }
            return tuples
        }
        .set { metadata_tuples }

    // Download each path listed in metadata. Emit each downloaded tuple
    // as it completes so downstream processes can start immediately.
    download_file(metadata_tuples)

    // Stream downloaded files directly to `samples` so alignment can run
    // per-sample as soon as a file is available. Include readgroup emitted
    // by the download process so downstream alignment gets it.
    download_file.out
        .map { tup -> tuple(tup[0], tup[1], tup[2]) } // tuple(sample_id, bam, readgroup)
        .set { samples }

    // Align per sample
    align_minimap2(samples, ref, fai)

    index_cram(align_minimap2.out)

    // Calculate sequencing depth
    mosdepth(index_cram.out, ref, fai)

    // Note: workflow-level `emit` removed. Each process uses `publishDir` to persist outputs.
}

/// Workflow stages

// Fetch metadata CSV
process fetch_meta {
    tag { meta }
    publishDir 'metadata', mode: 'copy'
    // input:
    // val meta

    memory { 64.MB }
    time { 10.min }
    output:
    path 'metadata.csv'

    script:
    """
    curl -sL "https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/refs/heads/main/data_tables/sequencing_data/data_ont_pre_release.index.csv" -o metadata.csv
    """
    // To support HTTP(S) URL or local path
    // """
    // if [[ "$meta" =~ ^https?:// ]]; then
        
    // else
    //     cp "$meta" metadata.csv
    // fi
    // """
}

// Download a file given its path (supports http(s), s3, gs and local paths)
process download_file {
    tag { sample_id }
    publishDir 'downloads', mode: 'copy'

    time { 12.h }
    memory { 1.GB }

    input:
    tuple val(sample_id), val(coverage), val(url), val(readgroup)

    output:
        tuple val(sample_id), path("${sample_id}.unaligned.bam"), val(readgroup)

    script:
    """
    set -euo pipefail
    # This pipeline assumes S3 paths only. Use aws CLI to copy.
    if ! command -v aws >/dev/null 2>&1; then
        echo "aws CLI not found; required for S3 downloads" >&2
        exit 2
    fi

    aws s3 --no-sign-request cp ${url} ${sample_id}.unaligned.bam
    """
}

// Filter metadata CSV: enforce coverage minimum and keep highest-coverage row per sample
process filter_meta {
    tag { metadata_file }
    publishDir 'metadata', mode: 'copy'

    memory { 64.MB }
    time { 10.min }

    input:
    path metadata_file

    output:
    tuple path('filtered_metadata.csv'), path('filtered_count.txt')

    script:
    """
    module load python
    python /home/hdashnow@xsede.org/myprojects/git/nf-long-tr/filter_metadata.py \
        --infile ${metadata_file} --outfile filtered_metadata.csv --countfile filtered_count.txt \
        --min_cov 25 --max_cov 30 --chemistry R10 --num_samples 25
    """

}

process align_minimap2 {
    memory { 32.GB * task.attempt }
    cpus 16
    time { 24.h * task.attempt }

    publishDir 'aligned', mode: 'copy'

    input:
    tuple val(sample), path(unaligned_bam), val(readgroup)
    path ref
    path fai

    output:
    // record sample ID and CRAM file
    tuple val("${sample}"), path("${sample}.cram")

    script:
        // Use literal backslash+t sequences ("\\t") so the read-group string
        // contains the two characters '\\' and 't' (escaped tab), not an actual
        // tab character. Minimap2/samtools require escaped tabs in RG lines.

        // Thread allocation logic
        def reserve = (task.cpus > 4) ? 1 : 0
        def sort_threads = Math.max(4, (task.cpus / 4).toInteger())
        def mm_threads = Math.max(1, task.cpus - reserve - sort_threads)
        //println "Thread allocation for sample ${sample}: minimap2=${mm_threads}, samtools_sort=${sort_threads}, reserve=${reserve}"
        // Check that sum of threads is equal to task.cpus
        assert (mm_threads + sort_threads + reserve) == task.cpus : "Thread allocation error"
        // Check that minimap2 has at least 1 thread and samtools sort has at least 1 thread
        assert mm_threads >= 1 : "Minimap2 must have at least 1 thread"
        assert sort_threads >= 1 : "Samtools sort must have at least 1 thread"
    """
    module load samtools
    set -o pipefail # ensure pipeline fails if any command fails
    # Convert input BAM to FASTQ on-the-fly and pipe into minimap2
    samtools fastq \
        -T MM,ML \
        ${unaligned_bam} | \
    minimap2 \
        -a \
        -y \
        -t ${mm_threads} \
        -R "${readgroup}" \
        ${ref} - | \
    samtools sort \
        -@ ${sort_threads} -O CRAM \
        --reference ${ref} \
        -T ${sample}.tmp \
        -o ${sample}.cram

    # Content-based checks using samtools
    if ! samtools quickcheck ${sample}.cram; then
        echo "Error: samtools quickcheck failed for ${sample}.cram" >&2
        exit 1
    else
        echo "samtools quickcheck passed for ${sample}.cram" >&2
    fi

    """
    // Note: the default behavior of minimap2 is ONT reads, so no need to specify the ONT preset with
    //  -x map-ont \

    // Additional code not currently used:
    // # Verify CRAM output: size and basic content checks before declaring success (temporary files are deleted on success)

    // if [ ! -s ${sample}.cram ]; then
    //     echo "Error: ${sample}.cram is empty or missing" >&2
    //     exit 1
    // fi

    // filesize=\$(wc -c < ${sample}.cram || echo 0)
    // min_bytes=1048576
    // if [ "\$filesize" -lt \${min_bytes} ]; then
    //     echo "Error: ${sample}.cram file size \$filesize bytes is smaller than \${min_bytes} bytes (1MB)" >&2
    //     exit 1
    // fi



    // if ! samtools view -H ${sample}.cram >/dev/null 2>&1; then
    //     echo "Error: samtools cannot read header of ${sample}.cram" >&2
    //     exit 1
    // fi
}

process index_cram {
    memory 8.GB
    time 24.h

    publishDir 'aligned', mode: 'copy'

    input:
    tuple val(sample), path(cram)

    output:
    tuple val(sample), path(cram), path("${cram}.crai")

    script:
    """
    module load samtools
    samtools index ${cram}
    """
}

// Calculate sequencing depth
process mosdepth {

    //container 'oras://community.wave.seqera.io/library/mosdepth:0.3.10--21ace4b9c76a055d'

    memory { 8.GB * task.attempt }
    cpus 4
    time 2.h

    publishDir 'aligned', mode: 'copy'

    input:
    tuple val(sample), path(cram), path(crai)
    path ref
    path fai

    output:
    tuple val(sample), path("${sample}.mosdepth.summary.txt")

    """
    mosdepth -n -x \
        -t ${task.cpus} \
        -f ${ref} \
        ${sample} \
        ${cram}
    """
}


