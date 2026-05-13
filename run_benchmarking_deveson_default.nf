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

def variantDir = '/pl/active/dashnowlab/projects/TR-benchmarking/benchmark-catalog-V2-Deveson-default/'

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


    // ch_chroms = Channel.of(
    //     'chr1','chr2','chr3','chr4','chr5',
    //     'chr6','chr7','chr8','chr9','chr10',
    //     'chr11','chr12','chr13','chr14','chr15',
    //     'chr16','chr17','chr18','chr19','chr20',
    //     'chr21','chr22','chrX','chrY'
    // )


    // ch_longtr_in = aligned_samples
    //     .combine(ch_chroms)
    //     .map { sample, aln, idx, karyotype, chrom ->
    //         tuple(sample, chrom, aln, idx, karyotype)
    // }

    // ch_longtr_out = longTR_per_chrom(ch_longtr_in, ref, fai)

    // ch_for_merge = ch_longtr_out
    //     .map { sample, chrom, vcf -> tuple(sample, vcf) }
    //     .groupTuple()

    // merged_longtr = mergeLongTR(ch_for_merge)


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
    cpus { 4 * task.attempt }
    time { 2.h * task.attempt }

    publishDir variantDir + '/mosdepth', mode: 'copy'

    input:
    tuple val(sample), path(aln), path(idx), val(karyotype)
    path ref

    output:
    tuple val(sample), path("${sample}.mosdepth.summary.txt")
    

    """
    mosdepth -n -x \
        -t $task.cpus \
        -f ${ref} \
        ${sample} \
        ${aln}
    """
}

process atarva {
    conda '/pl/active/dashnowlab/projects/TR-benchmarking/envs/atarva-0.5.0.yaml'
    //container 'dhaksnamoorthy/atarva:v0.3.1'

    cpus { 8 * task.attempt }
    memory { 16.GB * task.attempt }
    time { 24.h * task.attempt }

    publishDir variantDir + '/atarva', mode: 'copy'

    input:
    tuple val(sample), path(aln), path(idx), val(karyotype)
    path ref
    path fai

    output:
    path "${sample}.atarva.vcf"

    script:
    //def atarva_tr_regions = '/pl/active/dashnowlab/projects/TR-benchmarking/catalogs/tr_explorer_catalog/TR_catalog_for_atarva.bed.gz'
    def atarva_tr_regions = '/pl/active/dashnowlab/projects/TR-benchmarking/catalogs/benchmark-catalog-v2.atarva.bed.gz'
    

    """
    /projects/ealiyev@xsede.org/software/anaconda/envs/atarva-0.5.0/bin/atarva -t $task.cpus -f ${ref} -b ${aln} -r ${atarva_tr_regions} --format cram --karyotype ${karyotype} -o ${sample}.atarva.vcf
    """
}

process longTR {
    container 'community.wave.seqera.io/library/longtr:1.2--3a7af9434e146eab'

    cpus 1
    memory { 16.GB * task.attempt }
    time { 24.h * task.attempt }

    publishDir variantDir + '/longtr', mode: 'copy'

    input:
    tuple val(sample), path(aln), path(idx), val(karyotype)
    path ref
    path fai

    output:
    path "${sample}.longTR.vcf.gz"

    script:
    def alignment_params = [-1.0, -0.458675, -1.0, -0.458675, -0.00005800168, -1, -1]
    def tr_regions = '/pl/active/dashnowlab/projects/TR-benchmarking/catalogs/benchmark-catalog-v2.longtr.bed'
    def haploid_args = (karyotype == 'XY') ? '--haploid-chrs chrX,chrY' : ''

    """
    MAX_TR_LEN="\$(awk '{print \$3-\$2}' ${tr_regions} | sort -n | tail -n 1)";
    LongTR --help;
    LongTR \\
        --alignment-params ${alignment_params.join(',')} \\
        --fasta ${ref} \\
        --max-tr-len \$MAX_TR_LEN \\
        --regions ${tr_regions} \\
        ${haploid_args} \\
        --bams ${aln} \\
        --bam-samps ${sample} \\
        --bam-libs ${sample} \\
        --tr-vcf ${sample}.longTR.vcf.gz
    """
}

process longTR_per_chrom {
    tag "${sample}:${chrom}"

    container 'community.wave.seqera.io/library/longtr:1.2--3a7af9434e146eab'

    cpus 1
    memory { 16.GB * task.attempt }
    time { 24.h * task.attempt }

    publishDir variantDir + '/longtr_chrom', mode: 'copy'

    input:
    // sample + chrom + alignment tuple
    tuple val(sample), val(chrom), path(aln), path(idx), val(karyotype)
    // reference
    path ref
    path fai

    output:
    // keep chrom for downstream; unique file per chrom
    tuple val(sample), val(chrom), path("${sample}.${chrom}.longTR.vcf.gz")

    script:
    def alignment_params = [-1.0, -0.458675, -1.0, -0.458675, -0.00005800168, -1, -1]
    def tr_regions = '/pl/active/dashnowlab/projects/TR-benchmarking/catalogs/benchmark-catalog-v2.longtr.bed'
    def haploid_args = (karyotype == 'XY') ? '--haploid-chrs chrX,chrY' : ''

    """
    MAX_TR_LEN="\$(awk '{print \$3-\$2}' ${tr_regions} | sort -n | tail -n 1)";
    LongTR \\
        --alignment-params ${alignment_params.join(',')} \\
        --fasta ${ref} \\
        --max-tr-len \$MAX_TR_LEN \\
        --regions ${tr_regions} \\
        ${haploid_args} \\
        --chrom ${chrom} \\
        --bams ${aln} \\
        --tr-vcf ${sample}.${chrom}.longTR.vcf.gz
    """
}

process mergeLongTR {
    container 'community.wave.seqera.io/library/bcftools_samtools:3f506bc690e52c7d'
    
    cpus 1
    memory { 1.GB * task.attempt }
    time { 2.h * task.attempt }

    tag "${sample}"

    publishDir variantDir + '/longtr_chrom', mode: 'copy'

    input:
    // list of VCFs per sample
    tuple val(sample), path(vcfs)

    output:
    path "${sample}.longTR.vcf.gz"

    script:
    """
    for v in $vcfs; do
        if [ ! -f "\${v}.tbi" ] && [ ! -f "\${v}.csi" ]; then
            echo "Indexing \$v"
            tabix -p vcf "\$v"
        fi
    done
    bcftools concat -a -Oz -o ${sample}.longTR.vcf.gz $vcfs
    tabix -p vcf ${sample}.longTR.vcf.gz
    """
}

process straglr {
    container 'community.wave.seqera.io/library/straglr:1.5.5--d8ea229ed1f78ec0'

    cpus { 8 * task.attempt }
    memory { 32.GB * task.attempt }
    time { 24.h * task.attempt }

    publishDir variantDir + '/straglr', mode: 'copy'

    input:
    tuple val(sample), path(aln), path(idx), val(karyotype)
    path ref
    path fai

    output:
    path "${sample}.vcf"

    script:
    def straglr_tr_regions = '/pl/active/dashnowlab/projects/TR-benchmarking/catalogs/benchmark-catalog-v2.strglr.bed'
    def sex      = (karyotype == 'XX') ? 'f' : (karyotype == 'XY') ? 'm' : 'f'

    """
    python3 /pl/active/dashnowlab/software/straglr/straglr.py ${aln} ${ref} ${sample} --loci ${straglr_tr_regions} --chroms 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y --sex ${sex} --genotype_in_size --nprocs $task.cpus
    """
}

process strkit {
    container 'ghcr.io/davidlougheed/strkit:0.24.2'
    
    cpus { 8 * task.attempt }
    memory { 16.GB * task.attempt }
    time { 24.h * task.attempt }

    publishDir variantDir + '/strkit', mode: 'copy'

    input:
    tuple val(sample), path(aln), path(idx), val(karyotype)
    path ref
    path fai

    output:
    path "${sample}.strkit.vcf"

    script:
    //def strkit_tr_regions = '/pl/active/dashnowlab/projects/TR-benchmarking/catalogs/tr_explorer_catalog/TR_catalog_for_strkit.bed'
    def strkit_tr_regions = '/pl/active/dashnowlab/projects/TR-benchmarking/catalogs/benchmark-catalog-v2.strkit.bed'
    

    """
    strkit call ${aln} --min-reads 1 --ploidy ${karyotype} --realign --ref ${ref} --loc ${strkit_tr_regions} --vcf ${sample}.strkit.vcf --processes $task.cpus
    """
}

process strdust {
    cpus { 8 * task.attempt }
    memory { 16.GB * task.attempt }
    time { 24.h }

    publishDir variantDir + '/strdust', mode: 'copy'

    input:
    tuple val(sample), path(aln), path(idx), val(karyotype)
    path ref
    path fai

    output:
    path "${sample}.strdust.vcf"
    path "${sample}.strdust.sorted.vcf"

    script:
    def strdust_tr_regions = '/pl/active/dashnowlab/projects/TR-benchmarking/catalogs/benchmark-catalog-v2.strdust.bed'
    def haploid = (karyotype == 'XY') ? "--haploid chrX,chrY" : ""

    """
    /pl/active/dashnowlab/software/STRdust-0.16.0/STRdust-linux-musl \
      -R ${strdust_tr_regions} --unphased -t ${task.cpus} --sample ${sample} \
      ${haploid} ${ref} ${aln} > ${sample}.strdust.vcf

    { grep '^#' ${sample}.strdust.vcf; \
      grep -v '^#' ${sample}.strdust.vcf | sort -t "\$(printf '\\t')" -k1,1V -k2,2n; \
    } > ${sample}.strdust.sorted.vcf
    """
}

process vamos {
    conda 'envs/vamos-3.0.5.yaml'

    cpus { 8 * task.attempt }
    memory { 16.GB * task.attempt }
    time { 8.h * task.attempt }

    publishDir variantDir + '/vamos', mode: 'copy'

    input:
    tuple val(sample), path(aln), path(idx), val(karyotype)
    path ref
    path fai

    output:
    path "${sample}.vamos.vcf"
    path "${sample}.vamos.fixed.vcf"

    script:
    def vamos_tr_regions = '/pl/active/dashnowlab/projects/TR-benchmarking/catalogs/benchmark-catalog-v2.vamos.bed'
    def fix_vcf         = "${projectDir}/scripts/fix-vcf.py"

    """
    export LD_LIBRARY_PATH=\${CONDA_PREFIX}/lib:\$LD_LIBRARY_PATH
    /pl/active/dashnowlab/software/vamos-3.0.5/vamos/src/vamos --read -b ${aln} -r ${vamos_tr_regions} -s ${sample} -o ${sample}.vamos.vcf -S -Z -t $task.cpus
    python3 ${fix_vcf} ${sample}.vamos.vcf ${sample}.vamos.fixed.vcf --ref ${ref}
    """
}


process medaka {
    conda 'envs/medaka-2.1.1.yaml'

    cpus { 8 * task.attempt }
    memory { 16.GB * task.attempt }
    time { 24.h * task.attempt  }

    publishDir variantDir + '/medaka', mode: 'copy'

    input:
    tuple val(sample), path(aln), path(idx), val(karyotype)
    path ref
    path fai

    output:
    path "${sample}.medaka.vcf"

    script:
    def medaka_tr_regions = '/pl/active/dashnowlab/projects/TR-benchmarking/catalogs/benchmark-catalog-v2.medaka.bed'
    def ref    = '/pl/active/dashnowlab/data/ref-genomes/human_GRCh38_no_alt_analysis_set.fasta'
    def sex      = (karyotype == 'XX') ? 'female' : (karyotype == 'XY') ? 'male' : 'unknown'
    

    """
    medaka tandem ${aln} ${ref} ${medaka_tr_regions} ${sex} ${sample}.medaka.vcf --workers $task.cpus --sample_name ${sample} --debug
    """
}

