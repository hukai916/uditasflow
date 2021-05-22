// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process TEST {
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
    container "hukai916/r_demultiplexer:0.4"

    input:
    // tuple path(umi), path(read1), path(read2)
    path(umi)
    path(read1)
    path(read2)

    output:
    path "TEST.txt"

    script:
    def software = getSoftwareName(task.process)
    """
    echo TETDONE_${umi}_${read1}_${read2} > "TEST.txt"
    """
}
