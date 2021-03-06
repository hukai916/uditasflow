// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BWA_INDEX {
    label 'process_medium'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename: filename, options: params.options, publish_dir: getSoftwareName(task.process), publish_id: '') }

    container "hukai916/bwa:0.7.17"

    input:
    path ref_genome

    output:
    path 'index', emit: index

    // def software = getSoftwareName(task.process) // task variable not found?
    script:

    """
    mkdir index
    bwa index $options.args $ref_genome -p index/${ref_genome.baseName}
    """
}
