// Import generic module functions
include { initOptions; saveFiles; getSoftwareName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process ANNOTATE {
    label 'process_low'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename: filename, options: params.options, publish_dir: getSoftwareName(task.process), publish_id: '') }

    container "hukai916/r_annotation:0.1"

    input:
    val genome
    path ontarget_bam_both
    path ontarget_bam_R1only
    path ontarget_bam_R2only

    path offtarget_bam_both
    path offtarget_bam_R1only
    path offtarget_bam_R2only

    output:
    path 'res_fragment_size_dist/*_fragment_size_dist.pdf'
    path 'res_bed/*'
    path 'res_bedstack/*'

    script:

    """
    # bam fragment size distribution
    samtools index $ontarget_bam_both
    samtools index $offtarget_bam_both
    fragment_size_dist.R --bam_file=$ontarget_bam_both
    fragment_size_dist.R --bam_file=$offtarget_bam_both

    # bam2bed
    mkdir res_bed
    bedtools bamtobed -i $ontarget_bam_both > res_bed/${ontarget_bam_both}.bed
    bedtools bamtobed -i $ontarget_bam_R1only > res_bed/${ontarget_bam_R1only}.bed
    bedtools bamtobed -i $ontarget_bam_R2only > res_bed/${ontarget_bam_R2only}.bed

    bedtools bamtobed -i $offtarget_bam_both > res_bed/${offtarget_bam_both}.bed
    bedtools bamtobed -i $offtarget_bam_R1only > res_bed/${offtarget_bam_R1only}.bed
    bedtools bamtobed -i $offtarget_bam_R2only > res_bed/${offtarget_bam_R2only}.bed

    # bed2stack
    mkdir res_bedstack
    cat res_bed/${ontarget_bam_both}.bed | awk -F '\\t' 'BEGIN {OFS="\\t"} { print \$1,\$2,\$3,\$5,\$6 }' | sort | uniq -c | awk 'BEGIN {OFS="\\t"} { print \$2,\$3,\$4,"test",\$5,\$6,\$1 }' > res_bedstack/${ontarget_bam_both}.stacked.bed
    cat res_bed/${ontarget_bam_R1only}.bed | awk -F '\\t' 'BEGIN {OFS="\\t"} { print \$1,\$2,\$3,\$5,\$6}' | sort | uniq -c | awk 'BEGIN {OFS="\\t"} {print \$2,\$3,\$4,"test",\$5,\$6,\$1}' > res_bedstack/${ontarget_bam_R1only}.stacked.bed
    cat res_bed/${ontarget_bam_R2only}.bed | awk -F '\\t' 'BEGIN {OFS="\\t"} { print \$1,\$2,\$3,\$5,\$6}' | sort | uniq -c | awk 'BEGIN {OFS="\\t"} {print \$2,\$3,\$4,"test",\$5,\$6,\$1}' > res_bedstack/${ontarget_bam_R2only}.stacked.bed

    cat res_bed/${offtarget_bam_both}.bed | awk -F '\\t' 'BEGIN {OFS="\\t"} { print \$1,\$2,\$3,\$5,\$6}' | sort | uniq -c | awk 'BEGIN {OFS="\\t"} {print \$2,\$3,\$4,"test",\$5,\$6,\$1}' > res_bedstack/${offtarget_bam_both}.stacked.bed
    cat res_bed/${offtarget_bam_R1only}.bed | awk -F '\\t' 'BEGIN {OFS="\\t"} { print \$1,\$2,\$3,\$5,\$6}' | sort | uniq -c | awk 'BEGIN {OFS="\\t"} {print \$2,\$3,\$4,"test",\$5,\$6,\$1}' > res_bedstack/${offtarget_bam_R1only}.stacked.bed
    cat res_bed/${offtarget_bam_R2only}.bed | awk -F '\\t' 'BEGIN {OFS="\\t"} { print \$1,\$2,\$3,\$5,\$6}' | sort | uniq -c | awk 'BEGIN {OFS="\\t"} {print \$2,\$3,\$4,"test",\$5,\$6,\$1}' > res_bedstack/${offtarget_bam_R2only}.stacked.bed

    # annotate
    annotation.R --bed_file=res_bedstack/${ontarget_bam_both}.stacked.bed --genome=${genome}
    annotation.R --bed_file=res_bedstack/${ontarget_bam_R1only}.stacked.bed --genome=${genome}
    annotation.R --bed_file=res_bedstack/${ontarget_bam_R2only}.stacked.bed --genome=${genome}
    annotation.R --bed_file=res_bedstack/${offtarget_bam_both}.stacked.bed --genome=${genome}
    annotation.R --bed_file=res_bedstack/${offtarget_bam_R1only}.stacked.bed --genome=${genome}
    annotation.R --bed_file=res_bedstack/${offtarget_bam_R2only}.stacked.bed --genome=${genome}

    """
}
