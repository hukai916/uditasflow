// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process DEMULTIPLEX {
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
    container "hukai916/r_demultiplexer:0.3"

    input:
    path index1_file
    path index2_file
    path read1_file
    path read2_file
    path sample_file

    output:
    path 'res_demultiplex/index1_*.fastq.gz', emit: index1
    path 'res_demultiplex/index2_*.fastq.gz', emit: index2
    path 'res_demultiplex/*_R1*.fastq.gz', emit: read1
    path 'res_demultiplex/*_R2*.fastq.gz', emit: read2

    script:
    def software = getSoftwareName(task.process)
    """
    demultiplexer.R --index1_file=$index1_file \
                            --index2_file=$index2_file \
                            --read1_file=$read1_file \
                            --read2_file=$read2_file \
                            --sample_file=$sample_file
    """
}
