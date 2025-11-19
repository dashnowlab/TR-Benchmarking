#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ----------------- Params & Help -----------------
params.help      = params.help ?: false
params.list      = params.list ?: params.cram_list  // backward compat
params.base_dir  = params.base_dir ?: null
params.ref       = params.ref

if (params.help) {
    log.info """
    Usage:
        nextflow run main.nf --list <paths.list> --ref <reference fasta> [--base_dir <dir>]

    Notes:
        * --list: text file with one path per line (BAM or CRAM). Lines may be absolute,
                  relative (optionally resolved via --base_dir), and may include comments (#).
        * CRAM requires --ref and its .fai alongside (same as for BAM workflows that also use the ref).
        * For BAM, an index of either <file>.bai or <file>.bam.bai must be present.
    """
    System.exit(0)
}

def variantDir = '/pl/active/dashnowlab/projects/TR-benchmarking/'

// ----------------- Helpers -----------------
def resolvePath = { String line ->
    def p = line.trim()
    if (!p) return null
    if (p.startsWith('#')) return null
    if (params.base_dir && !p.startsWith('/'))
        return file("${params.base_dir}/${p}")
    return file(p)
}

// Drop-in replacement
def sampleName = { x ->
    // get leaf filename as String from Path/File/String
    def leaf = (x instanceof java.nio.file.Path) ? x.fileName.toString()
             : (x instanceof java.io.File)       ? x.name
             :                                     x.toString()

    // remove trailing extensions
    def bn = leaf.replaceFirst(/\.(cram|bam)$/, '')
    bn = bn.replaceFirst(/\.merged$/, '')
    return bn
}

def bamIndex = { File bam ->
    def bai1 = file("${bam}.bai")
    def bai2 = file(bam.toString().replaceFirst(/\.bam$/, '.bai'))
    def e1 = new File(bai1.toString()).exists()
    def e2 = new File(bai2.toString()).exists()
    if (e1) return bai1
    if (e2) return bai2
    throw new IllegalArgumentException("Missing BAM index (.bai or .bam.bai) for ${bam}")
}

// ----------------- Workflow -----------------
workflow {

    // Build channel from .list (TAB-separated: <bam|cram><TAB><karyotype>)
    Channel
      .fromPath(params.list, checkIfExists: true)
      .splitText()
      .map { it as String }
      .map { it.replaceAll(/\r$/, '') }                     // strip Windows CRLF if present
      .filter { line -> line.trim() && !line.trim().startsWith('#') }
      .map { line ->
          // Expect exactly 2 columns separated by a TAB
          def parts = line.split('\t', -1)
          if (parts.size() < 2)
              throw new IllegalArgumentException("Expected two TAB-separated columns: <bam|cram>\\t<karyotype>; got: '${line}'")

          def alnPath   = parts[0].trim()
          def karyotype = parts[1].trim()

          def aln = resolvePath(alnPath)
          if (aln == null)
              throw new IllegalArgumentException("Could not resolve path from line: '${line}'")

          def sample = sampleName(aln)

          // Determine index path
          def idx
          if (aln.name.endsWith('.cram')) {
              idx = file("${aln}.crai", checkIfExists: true)
          } else if (aln.name.endsWith('.bam')) {
              idx = bamIndex(aln)
          } else {
              throw new IllegalArgumentException("Unsupported file type for ${aln}; expected .bam or .cram")
          }

          tuple(sample, aln, idx, karyotype)
      }
      .set { aligned_samples }   // emits: [sample, alignment_file, index_file, karyotype]



    // Reference genome (required if any CRAMs are present; harmless otherwise)
    ref = file(params.ref, checkIfExists: true)
    fai = file("${params.ref}.fai", checkIfExists: true)

    // Downstream processes
    // print_aligned_samples(aligned_samples)
    // mosdepth(aligned_samples, ref)
    // strkit(aligned_samples, ref, fai)
    strkitrr(aligned_samples, ref, fai)
    // strdust(aligned_samples, ref, fai)
    // vamos(aligned_samples, ref, fai)
    // longTR(aligned_samples, ref, fai)
    // medaka(aligned_samples, ref, fai)
}
/* -------------------------------------------------------------------------- */
/* Processes                                                                  */
/* -------------------------------------------------------------------------- */
process print_aligned_samples {
  echo true                       // prints stdout immediately
    cpus 1
    memory { 1.GB * task.attempt }
    time { 1.h * task.attempt }
  input:
    tuple val(sample), path(aln), path(idx), val(karyotype)
  script:
    """
    echo "sample: ${sample}"
    echo "alignment: ${aln}"
    echo "index: ${idx}"
    echo "karyotype: ${karyotype}"
    """
}

// Calculate sequencing depth
process mosdepth {

    container 'oras://community.wave.seqera.io/library/mosdepth:0.3.10--21ace4b9c76a055d'

    memory { 8.GB * task.attempt }
    cpus 4
    time { 2.h * task.attempt }

    publishDir variantDir + '/mosdepth', mode: 'copy'

    input:
    tuple val(sample), path(aln), path(idx), val(karyotype)
    path ref

    output:
    tuple val(sample), path("${sample}.mosdepth.summary.txt")

    """
    mosdepth -n -x \
        -t 4 \
        -f ${ref} \
        ${sample} \
        ${aln}
    """
}

process strkit {

    cpus 2
    memory { 2.GB * task.attempt }
    time { 2.h * task.attempt }

    publishDir variantDir + '/strkit', mode: 'copy'

    input:
    tuple val(sample), path(aln), path(idx), val(karyotype)
    path ref
    path fai

    output:
    path "${sample}.strkit.vcf"

    script:
    def strkit_tr_regions = '/pl/active/dashnowlab/projects/TR-benchmarking/catalogs/test-isolated-vc-catalog.strkit.bed'

    """
    strkit call ${aln} --min-reads 10 --ploidy ${karyotype} --realign --ref ${ref} --loc ${strkit_tr_regions} --vcf ${sample}.strkit.vcf
    """
}

process strkitrr {

    cpus 2
    memory { 2.GB * task.attempt }
    time { 2.h * task.attempt }

    publishDir variantDir + '/strkit', mode: 'copy'

    input:
    tuple val(sample), path(aln), path(idx), val(karyotype)
    path ref
    path fai

    output:
    path "${sample}.strkitrr.vcf"

    script:
    def strkit_tr_regions = '/pl/active/dashnowlab/projects/TR-benchmarking/catalogs/test-isolated-vc-catalog.strkit.bed'

    """
    strkit call ${aln} --respect-ref --min-reads 10 --ploidy ${karyotype} --realign --ref ${ref} --loc ${strkit_tr_regions} --vcf ${sample}.strkitrr.vcf
    """
}

process strdust {

    cpus 2
    memory { 2.GB * task.attempt }
    time { 2.h * task.attempt }

    publishDir variantDir + '/strdust', mode: 'copy'

    input:
    tuple val(sample), path(aln), path(idx), val(karyotype)
    path ref
    path fai

    output:
    path "${sample}.strdust.vcf"

    script:
    def strkit_tr_regions = '/pl/active/dashnowlab/projects/TR-benchmarking/catalogs/test-isolated-vc-catalog.strkit.bed'
    def haploid = (karyotype == 'XY') ? "--haploid chrX,chrY" : ""

    """
    /pl/active/dashnowlab/software/STRdust/target/release/STRdust -R ${strkit_tr_regions} --unphased --support 5 --sample ${sample} \\
    ${haploid} ${ref} ${aln} > ${sample}.strdust.vcf
    """
}

process vamos {
    conda 'envs/vamos.yaml'

    cpus 2
    memory { 8.GB * task.attempt }
    time { 2.h * task.attempt }

    publishDir variantDir + '/vamos', mode: 'copy'

    input:
    tuple val(sample), path(aln), path(idx), val(karyotype)
    path ref
    path fai

    output:
    path "${sample}.vamos.vcf"

    script:
    def vamos_tr_regions = '/pl/active/dashnowlab/projects/TR-benchmarking/catalogs/test-isolated-vc-catalog.vamos.bed'

    """
    export LD_LIBRARY_PATH=\${CONDA_PREFIX}/lib:\$LD_LIBRARY_PATH
    /pl/active/dashnowlab/work/aavvaru/vamos/src/vamos --read -b ${aln} -r ${vamos_tr_regions} -s ${sample} -o ${sample}.vamos.vcf
    """
}

process medaka {
    conda 'envs/medaka-2.1.1.yaml'

    cpus 2
    memory { 2.GB * task.attempt }
    time { 2.h * task.attempt }

    publishDir variantDir + '/medaka', mode: 'copy'

    input:
    tuple val(sample), path(aln), path(idx), val(karyotype)
    path ref
    path fai

    output:
    path "${sample}.medaka.vcf"

    script:
    def medaka_tr_regions = '/pl/active/dashnowlab/projects/TR-benchmarking/catalogs/test-isolated-vc-catalog.strkit.bed'
    def sex      = (karyotype == 'XX') ? 'female' : (karyotype == 'XY') ? 'male' : 'unknown'


    """
    medaka tandem --sample_name ${sample} ${aln} ${ref} ${medaka_tr_regions} ${sex} ${sample}.medaka.vcf
    """
}

process longTR {
    conda 'envs/longtr-gcc14-run.yaml'

    cpus 1
    memory { 4.GB * task.attempt }
    time { 4.h * task.attempt }

    publishDir variantDir + '/longtr', mode: 'copy'

    input:
    tuple val(sample), path(aln), path(idx), val(karyotype)
    path ref
    path fai

    output:
    path "${sample}.longTR.vcf.gz"

    script:
    def alignment_params = [-1.0, -0.458675, -1.0, -0.458675, -0.00005800168, -1, -1]
    def tr_regions = '/pl/active/dashnowlab/projects/TR-benchmarking/catalogs/test-isolated-vc-catalog.longtr.bed'

    """
    export LD_LIBRARY_PATH=\${CONDA_PREFIX}/lib:\$LD_LIBRARY_PATH
    MAX_TR_LEN="\$(awk '{print \$3-\$2}' ${tr_regions} | sort -n | tail -n 1)";
    /pl/active/dashnowlab/software/LongTR-1.2/LongTR \\
        --alignment-params ${alignment_params.join(',')} \\
        --fasta ${ref} \\
        --min-reads 1 \\
        --max-tr-len \$MAX_TR_LEN \\
        --regions ${tr_regions} \\
        --bams ${aln} \\
        --tr-vcf ${sample}.longTR.vcf.gz
    """
}

