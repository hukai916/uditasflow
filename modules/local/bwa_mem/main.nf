// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process BWA_MEM {
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
    path 'ontarget/*/bam_qc', emit: ontarget_bam_qc

    path 'offtarget/*/offtarget*.both.bam', emit: offtarget_both_bam
    path 'offtarget/*/offtarget*.R1only.bam', emit: offtarget_R1only_bam
    path 'offtarget/*/offtarget*.R2only.bam', emit: offtarget_R2only_bam
    path 'offtarget/*/bam_qc', emit: offtarget_bam_qc

    script:
    // ontarget_read1.simpleName = ontarget_read1.simpleName[0..-3] // to get rid of _R1/_R2 from simpleName; also can't set readonly property
    ontarget_read_simpleName  = ontarget_read1.simpleName[0..-4]
    offtarget_read_simpleName = offtarget_read1.simpleName[0..-4]

    """
    INDEX=`find -L ./ -name "*.amb" | sed 's/.amb//'`

    targets = ( ontarget offtarget)

    for i in "{targets[@]}"
    do
    ## map both R1 and R2 reads:
    mkdir -p \${i}/${bam_dir}
    bwa mem $options.args \\
            \$INDEX \\
            $\${i}_read1 $\${i}_read2 \\
            | samtools view -b -o \${i}/${bam_dir}/${\${i}_read_simpleName}.both.bam

    samtools sort \${i}/${bam_dir}/${\${i}_read_simpleName}.both.bam -o \${i}/${bam_dir}/${\${i}_read_simpleName}.both.bam
    samtools index \${i}/${bam_dir}/${\${i}_read_simpleName}.both.bam
    samtools stats \${i}/${bam_dir}/${\${i}_read_simpleName}.both.bam > \${i}/${bam_dir}/${\${i}_read_simpleName}.both.bam.stat

    multiqc \${i}/${bam_dir}/${\${i}_read_simpleName}.both.bam.stat --outdir \${i}/${bam_dir}/bam_qc -n ${\${i}_read_simpleName}.both.bam.stat

    ## map R1 only:
    bwa mem $options.args \\
            \$INDEX \\
            $\${i}_read1 \\
            | samtools view -b -o \${i}/${bam_dir}/${\${i}_read_simpleName}.R1only.bam

    samtools sort \${i}/${bam_dir}/${\${i}_read_simpleName}.R1only.bam -o \${i}/${bam_dir}/${\${i}_read_simpleName}.R1only.bam
    samtools index \${i}/${bam_dir}/${\${i}_read_simpleName}.R1only.bam
    samtools stats \${i}/${bam_dir}/${\${i}_read_simpleName}.R1only.bam > \${i}/${bam_dir}/${\${i}_read_simpleName}.R1only.bam.stat

    multiqc \${i}/${bam_dir}/${\${i}_read_simpleName}.R1only.bam.stat --outdir \${i}/${bam_dir}/bam_qc -n ${\${i}_read_simpleName}.R1only.bam.stat

    ## map R2 only:
    bwa mem $options.args \\
            \$INDEX \\
            $\${i}_read2 \\
            | samtools view -b -o \${i}/${bam_dir}/${\${i}_read_simpleName}.R2only.bam

    samtools sort \${i}/${bam_dir}/${\${i}_read_simpleName}.R2only.bam -o \${i}/${bam_dir}/${\${i}_read_simpleName}.R2only.bam
    samtools index \${i}/${bam_dir}/${\${i}_read_simpleName}.R2only.bam
    samtools stats \${i}/${bam_dir}/${\${i}_read_simpleName}.R2only.bam > \${i}/${bam_dir}/${\${i}_read_simpleName}.R2only.bam.stat

    multiqc \${i}/${bam_dir}/${\${i}_read_simpleName}.R2only.bam.stat --outdir \${i}/${bam_dir}/bam_qc -n ${\${i}_read_simpleName}.R2only.bam.stat


    done

  ## ontarget:


  ## offtarget:

    """
}
