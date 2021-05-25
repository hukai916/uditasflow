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

  ## ontarget:
    ## map both R1 and R2 reads:
    mkdir -p ontarget/${bam_dir}
    mkdir -p ontarget/${bam_dir}/bam_qc_${ontarget_read_simpleName}

    bwa mem $options.args \\
            \$INDEX \\
            $ontarget_read1 $ontarget_read2 \\
            | samtools view -b -o ontarget/${bam_dir}/${ontarget_read_simpleName}.both.bam

    samtools sort ontarget/${bam_dir}/${ontarget_read_simpleName}.both.bam -o ontarget/${bam_dir}/${ontarget_read_simpleName}.both.bam
    samtools index ontarget/${bam_dir}/${ontarget_read_simpleName}.both.bam
    samtools stats ontarget/${bam_dir}/${ontarget_read_simpleName}.both.bam > ontarget/${bam_dir}/${ontarget_read_simpleName}.both.bam.stat

    multiqc ontarget/${bam_dir}/${ontarget_read_simpleName}.both.bam.stat --outdir ontarget/${bam_dir}/bam_qc_${ontarget_read_simpleName} -n ${ontarget_read_simpleName}.both.bam.stat

    ## map R1 only:
    bwa mem $options.args \\
            \$INDEX \\
            $ontarget_read1 \\
            | samtools view -b -o ontarget/${bam_dir}/${ontarget_read_simpleName}.R1only.bam

    samtools sort ontarget/${bam_dir}/${ontarget_read_simpleName}.R1only.bam -o ontarget/${bam_dir}/${ontarget_read_simpleName}.R1only.bam
    samtools index ontarget/${bam_dir}/${ontarget_read_simpleName}.R1only.bam
    samtools stats ontarget/${bam_dir}/${ontarget_read_simpleName}.R1only.bam > ontarget/${bam_dir}/${ontarget_read_simpleName}.R1only.bam.stat

    multiqc ontarget/${bam_dir}/${ontarget_read_simpleName}.R1only.bam.stat --outdir ontarget/${bam_dir}/bam_qc_${ontarget_read_simpleName} -n ${ontarget_read_simpleName}.R1only.bam.stat

    ## map R2 only:
    bwa mem $options.args \\
            \$INDEX \\
            $ontarget_read2 \\
            | samtools view -b -o ontarget/${bam_dir}/${ontarget_read_simpleName}.R2only.bam

    samtools sort ontarget/${bam_dir}/${ontarget_read_simpleName}.R2only.bam -o ontarget/${bam_dir}/${ontarget_read_simpleName}.R2only.bam
    samtools index ontarget/${bam_dir}/${ontarget_read_simpleName}.R2only.bam
    samtools stats ontarget/${bam_dir}/${ontarget_read_simpleName}.R2only.bam > ontarget/${bam_dir}/${ontarget_read_simpleName}.R2only.bam.stat

    multiqc ontarget/${bam_dir}/${ontarget_read_simpleName}.R2only.bam.stat --outdir ontarget/${bam_dir}/bam_qc_${ontarget_read_simpleName} -n ${ontarget_read_simpleName}.R2only.bam.stat

  ## offtarget:
    ## map both R1 and R2 reads:
    mkdir -p offtarget/${bam_dir}
    mkdir -p ontarget/${bam_dir}/bam_qc_${ontarget_read_simpleName}

    bwa mem $options.args \\
            \$INDEX \\
            $offtarget_read1 $offtarget_read2 \\
            | samtools view -b -o offtarget/${bam_dir}/${offtarget_read_simpleName}.both.bam

    samtools sort offtarget/${bam_dir}/${offtarget_read_simpleName}.both.bam -o offtarget/${bam_dir}/${offtarget_read_simpleName}.both.bam
    samtools index offtarget/${bam_dir}/${offtarget_read_simpleName}.both.bam
    samtools stats offtarget/${bam_dir}/${offtarget_read_simpleName}.both.bam > offtarget/${bam_dir}/${offtarget_read_simpleName}.both.bam.stat

    multiqc offtarget/${bam_dir}/${offtarget_read_simpleName}.both.bam.stat --outdir offtarget/${bam_dir}/bam_qc_${offtarget_read_simpleName} -n ${offtarget_read_simpleName}.both.bam.stat

    ## map R1 only:
    bwa mem $options.args \\
            \$INDEX \\
            $offtarget_read1 \\
            | samtools view -b -o offtarget/${bam_dir}/${offtarget_read_simpleName}.R1only.bam

    samtools sort offtarget/${bam_dir}/${offtarget_read_simpleName}.R1only.bam -o offtarget/${bam_dir}/${offtarget_read_simpleName}.R1only.bam
    samtools index offtarget/${bam_dir}/${offtarget_read_simpleName}.R1only.bam
    samtools stats offtarget/${bam_dir}/${offtarget_read_simpleName}.R1only.bam > offtarget/${bam_dir}/${offtarget_read_simpleName}.R1only.bam.stat

    multiqc offtarget/${bam_dir}/${offtarget_read_simpleName}.R1only.bam.stat --outdir offtarget/${bam_dir}/bam_qc_${offtarget_read_simpleName} -n ${offtarget_read_simpleName}.R1only.bam.stat

    ## map R2 only:
    bwa mem $options.args \\
            \$INDEX \\
            $offtarget_read2 \\
            | samtools view -b -o offtarget/${bam_dir}/${offtarget_read_simpleName}.R2only.bam

    samtools sort offtarget/${bam_dir}/${offtarget_read_simpleName}.R2only.bam -o offtarget/${bam_dir}/${offtarget_read_simpleName}.R2only.bam
    samtools index offtarget/${bam_dir}/${offtarget_read_simpleName}.R2only.bam
    samtools stats offtarget/${bam_dir}/${offtarget_read_simpleName}.R2only.bam > offtarget/${bam_dir}/${offtarget_read_simpleName}.R2only.bam.stat

    multiqc offtarget/${bam_dir}/${offtarget_read_simpleName}.R2only.bam.stat --outdir offtarget/${bam_dir}/bam_qc_${offtarget_read_simpleName} -n ${offtarget_read_simpleName}.R2only.bam.stat

    """
}
