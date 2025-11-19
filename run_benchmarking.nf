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

          def alnPath  = parts[0].trim()
          def karyo    = parts[1].trim()

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

          tuple(sample, aln, idx, karyo)
      }
      .set { aligned_samples }   // emits: [sample, alignment_file, index_file, karyotype]



    // Reference genome (required if any CRAMs are present; harmless otherwise)
    ref = file(params.ref, checkIfExists: true)
    fai = file("${params.ref}.fai", checkIfExists: true)

    // Downstream processes
    print_aligned_samples(aligned_samples)
    mosdepth(aligned_samples, ref)
    atarva(aligned_samples, ref, fai)
    longTR(aligned_samples, ref, fai)
    medaka(aligned_samples, ref, fai)
    straglr(aligned_samples, ref, fai)
    strdust(aligned_samples, ref, fai)
    strkit(aligned_samples, ref, fai)
    vamos(aligned_samples, ref, fai)
    
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

process atarva {
    //conda 'envs/atarva-3.0.0.yaml'
    container 'dhaksnamoorthy/atarva:v0.3.1'

    cpus 1
    memory { 4.GB * task.attempt }
    time { 4.h * task.attempt }

    publishDir variantDir + '/atarva', mode: 'copy'

    input:
    tuple val(sample), path(aln), path(idx), val(karyotype)
    path ref
    path fai

    output:
    path "${sample}.atarva.vcf"

    script:
    def atarva_tr_regions = '/pl/active/dashnowlab/projects/TR-benchmarking/catalogs/test-isolated-vc-catalog.atarva.bed.gz'

    """
    atarva -t 1 -f ${ref} -b ${aln} -r ${atarva_tr_regions} --format cram --haplotag HP --min-reads 1 --karyotype ${karyotype} -o ${sample}.atarva.vcf
    """
}

process longTR {
    container 'community.wave.seqera.io/library/longtr:1.2--3a7af9434e146eab'

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
    MAX_TR_LEN="\$(awk '{print \$3-\$2}' ${tr_regions} | sort -n | tail -n 1)";
    LongTR --help;
    LongTR \\
        --alignment-params ${alignment_params.join(',')} \\
        --fasta ${ref} \\
        --min-reads 1 \\
        --max-tr-len \$MAX_TR_LEN \\
        --regions ${tr_regions} \\
        --bams ${aln} \\
        --tr-vcf ${sample}.longTR.vcf.gz
    """
}

process straglr {
    container 'community.wave.seqera.io/library/straglr:1.5.5--d8ea229ed1f78ec0'

    cpus 1
    memory { 4.GB * task.attempt }
    time { 4.h * task.attempt }

    publishDir variantDir + '/straglr', mode: 'copy'

    input:
    tuple val(sample), path(aln), path(idx), val(karyotype)
    path ref
    path fai

    output:
    path "${sample}.vcf"

    script:
    def straglr_tr_regions = '/pl/active/dashnowlab/projects/TR-benchmarking/catalogs/test-isolated-vc-catalog.strglr.bed'

    """
    python3 /pl/active/dashnowlab/software/straglr/straglr.py ${aln} ${ref} ${sample} --loci ${straglr_tr_regions} --chroms 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y --min_support 1 --min_cluster_size 1 --max_num_clusters 2 --genotype_in_size --nprocs 1
    """
}

process strkit {
    container 'ghcr.io/davidlougheed/strkit:0.24.2'
    
    cpus 1
    memory { 4.GB * task.attempt }
    time { 4.h * task.attempt }

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

process strdust {

    cpus 1
    memory { 4.GB * task.attempt }
    time { 4.h * task.attempt }

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
    container 'community.wave.seqera.io/library/vamos:2.1.7--45fcc86e7cdc4363'

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
    vamos --read -b ${aln} -r ${vamos_tr_regions} -s ${sample} -o ${sample}.vamos.vcf
    """
}

process medaka {
    conda 'envs/medaka-2.1.1.yaml'

    cpus 1
    memory { 4.GB * task.attempt }
    time { 4.h * task.attempt }

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
