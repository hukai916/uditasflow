// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BWA {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename: filename, options: params.options, publish_dir: getSoftwareName(task.process), publish_id: '') }

    container "hukai916/bwa:0.7.17"

    input:
    path ref_genome
    val bam_dir
    path ontarget_read1
    path ontarget_read2
    path offtarget_read1
    path offtarget_read2

    output:
    path 'ontarget/*/ontarget*.both.bam', emit: ontarget_both_bam
    path 'ontarget/*/ontarget*.R1only.bam', emit: ontarget_R1only_bam
    path 'ontarget/*/ontarget*.R2only.bam', emit: ontarget_R2only_bam

    // path 'offtarget/*/offtarget*.both.bam', emit: offtarget_both_bam
    // path 'offtarget/*/offtarget*.R1only.bam', emit: offtarget_R1only_bam
    // path 'offtarget/*/offtarget*.R2only.bam', emit: offtarget_R1only_bam

    script:

    """
    bwa mem -M $ref_genome $ontarget_read1 $ontarget_read2 | samtools view -b -o ontarget/${bam_dir}/${ontarget_read1.simpleName}.both.bam

    samtools sort ontarget/${bam_dir}/${ontarget_read1.simpleName}.both.bam -o ontarget/${bam_dir}/${ontarget_read1.simpleName}.both.bam
    samtools index ontarget/${bam_dir}/${ontarget_read1.simpleName}.both.bam
    samtools stats ontarget/${bam_dir}/${ontarget_read1.simpleName}.both.bam > ontarget/${bam_dir}/${ontarget_read1.simpleName}.both.bam.stat

    multiqc ontarget/${bam_dir}/${ontarget_read1.simpleName}.both.bam.stat --outdir ontarget/${bam_dir}/bam_qc -n ${ontarget_read1.simpleName}.both.bam.stat



    """
}
