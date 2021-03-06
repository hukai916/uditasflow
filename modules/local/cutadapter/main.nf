// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CUTADAPTER {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename: filename, options: params.options, publish_dir: getSoftwareName(task.process), publish_id: '') }

    container "hukai916/cutadapt_xenial:3.4"

    input:
    path read1
    path read2
    val adapter_read1
    val adapter_read2
    val adapter1_rc
    val adapter2_rc

    output:
    path 'res_cutadapter/*R1.fastq', emit: cutadapter_read1
    path 'res_cutadapter/*R2.fastq', emit: cutadapter_read2

    script:
    def software = getSoftwareName(task.process)
    // Get the complement of a DNA sequence: ref: http://groovyconsole.appspot.com/script/29005
    // Complement table taken from http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html
    String.metaClass.complement = {
    def complements = [ A:'T', T:'A', U:'A', G:'C', C:'G', Y:'R', R:'Y', S:'S', W:'W', K:'M', M:'K', B:'V', D:'H', H:'D', V:'B', N:'N' ]
    delegate.toUpperCase().replaceAll( /./ ) { complements."$it" ?: 'X' }
    }

    if (adapter1_rc) { adapter_read1 = adapter_read1.reverse().complement() }
    if (adapter2_rc) { adapter_read2 = adapter_read2.reverse().complement() }

    """
    mkdir res_cutadapter
    cutadapt $options.args \
          -a $adapter_read1 \
          -A $adapter_read2 \
          -o res_cutadapter/${read1.baseName} \
          -p res_cutadapter/${read2.baseName} \
          $read1 $read2
    """
}
