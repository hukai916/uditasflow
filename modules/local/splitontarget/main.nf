// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SPLITONTARGET {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename: filename, options: params.options, publish_dir: getSoftwareName(task.process), publish_id: '') }

    container "hukai916/r_demultiplexer:0.4"

    input:
    tuple val(sample_id), val(read_anchor), val(ref_seq)
    path read1_fq
    path read2_fq

    output:
    path './res_split_ontarget/ontarget*R1.fastq.gz', emit: ontarget_read1
    path './res_split_ontarget/ontarget*R2.fastq.gz', emit: ontarget_read2
    path './res_split_ontarget/offtarget*R1.fastq.gz', emit: offtarget_read1
    path './res_split_ontarget/offtarget*R2.fastq.gz', emit: offtarget_read2

    script:
    def software = getSoftwareName(task.process)
    if (read_anchor == "read1") {
      read_file_anchor = read1_fq
      read_file_other  = read2_fq
    } else if (read_anchor == "read2") {
      read_file_anchor = read2_fq
      read_file_other  = read1_fq
    }

    """
    split_ontarget.R $options.args \
                --ref_seq=$ref_seq \
                --read_file_anchor=$read_file_anchor \
                --read_file_other=$read_file_other
    """
}
