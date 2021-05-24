// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CUTADAPTER {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename: filename, options: params.options, publish_dir: getSoftwareName(task.process), publish_id: '') }

    // conda (params.enable_conda ? "dranew:bcl2fastq=2.19.0" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "url_to_singularity_image"
    // } else {
    //     // container "quay.io/biocontainers/multiqc:1.9--pyh9f0ad1d_0"
    //     container "hukai916/bcl2fastq:2.20.0-centos7"
    // }
    container "hukai916/cutadapt_xenial:3.4"

    input:
    path read1
    path read2
    val adapter_read1
    val adapter_read2

    output:
    path 'res_cutadapter/*R1.fastq', emit: cutadapter_read1
    path 'res_cutadapter/*R2.fastq', emit: cutadapter_read2

    script:
    def software = getSoftwareName(task.process)
    """
    mkdir res_cutadapter
    echo TETDONE_${read1}_${read2} > "TEST.txt"
    // cutadapt $options.args \
    //         -a $adapter_read1 \
    //         -A $adapter_read2 \
    //         -o res_cutadapter/${read1.baseName} \
    //         -p res_cutadapter/${read2.baseName} \
    //         $read1 $read2
    """
}
