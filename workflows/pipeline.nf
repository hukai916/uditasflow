////////////////////////////////////////////////////
/* --         LOCAL PARAMETER VALUES           -- */
////////////////////////////////////////////////////

params.summary_params = [:]

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// Validate input parameters
// Workflow.validateWorkflowParams(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }


////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

////////////////////////////////////////////////////
/* --       IMPORT MODULES / SUBWORKFLOWS      -- */
////////////////////////////////////////////////////

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? " --title \"$params.multiqc_title\"" : ''

// Modules: local
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions'   addParams( options: [publish_files : ['csv':'']] )

// Modules: nf-core/modules
include { FASTQC                } from '../modules/nf-core/software/fastqc/main'  addParams( options: modules['fastqc']            )

include { MULTIQC               } from '../modules/nf-core/software/multiqc/main' addParams( options: multiqc_options              )

// Subworkflows: local
include { INPUT_CHECK           } from '../subworkflows/local/input_check'        addParams( options: [:]                          )

include { BCL2FASTQ             } from '../modules/local/bcl2fastq/main'          addParams( options: modules['bcl2fastq']         )

include { DEMULTIPLEX           } from '../modules/local/demultiplex/main'        addParams( options: modules['demultiplex']       )

include { PARSEUMI              } from '../modules/local/parseumi/main'           addParams( options: modules['parseumi']          )

include { COLLAPSEUMI           } from '../modules/local/collapseumi/main'        addParams( options: modules['collapseumi']       )

include { CUTADAPTER            } from '../modules/local/cutadapter/main'         addParams( options: modules['cutadapter']        )

include { SPLITONTARGET         } from '../modules/local/splitontarget/main'      addParams( options: modules['splitontarget']     )

include { BWA_INDEX             } from '../modules/local/bwa_index/main'          addParams( options: modules['bwa_index']         )

include { BWA_MEM               } from '../modules/local/bwa_mem/main'            addParams( options: modules['bwa_mem']           )

include { ANNOTATE              } from '../modules/local/annotate/main'            addParams( options: modules['bwa_mem']           )

include { TEST                  } from '../modules/local/test/main'               addParams( options: modules['test']              )


////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

// Info required for completion email and summary
def multiqc_report = []

workflow UDITASFLOW {

    ch_software_versions = Channel.empty()
    ch_bcl_raw           = Channel.fromPath(params.bcl_raw)
    ch_sample_file       = Channel.fromPath(params.sample_file)
    ch_ref_genome        = Channel.fromPath(params.ref_genome)
    // ch_adapter_read1     = Channel.of(params.adapter_read1)
    // ch_adapter_read2     = Channel.of(params.adapter_read2)

    BCL2FASTQ (
      ch_bcl_raw
    )

    DEMULTIPLEX (
      BCL2FASTQ.out.index1,
      BCL2FASTQ.out.index2,
      BCL2FASTQ.out.read1,
      BCL2FASTQ.out.read2,
      ch_sample_file
    )

    if (params.umi_index == "index1") {
      index = DEMULTIPLEX.out.index1.collect().flatten()
    } else {
      index = DEMULTIPLEX.out.index2.collect().flatten()
    }

    PARSEUMI (index,
      params.umi_start,
      params.umi_end
    )

    umi = PARSEUMI.out.umi.toSortedList( { a, b -> a.getName() <=> b.getName() } ).flatten()
    // Note that must sort by filename since path is full path, and files are distributed in random folders under work dir.
    read1 = DEMULTIPLEX.out.read1.toSortedList( { a, b -> a.getName() <=> b.getName() } ).flatten().filter( ~/^(?!.*undetermined).*/ )
    // filter out sample with undetermined in its name: regex ref: https://stackoverflow.com/questions/33159862/regex-match-word-not-containing
    // confirmed that the filtering is working by checking teh results with .filter( ~/^(?!.*S1_R1).*/
    read2 = DEMULTIPLEX.out.read2.toSortedList( { a, b -> a.getName() <=> b.getName() } ).flatten().filter( ~/^(?!.*undetermined).*/ )

    COLLAPSEUMI (
      umi,
      read1,
      read2
    )

    CUTADAPTER (
      COLLAPSEUMI.out.umi_read1,
      COLLAPSEUMI.out.umi_read2,
      params.adapter_read1,
      params.adapter_read2,
      params.adapter1_rc,
      params.adapter2_rc
      // ch_adapter_read2 // note that if using Channel, it will be consumed and will only run for one instance for the process.
    )

    def samples = []
    new File(params.sample_file).splitEachLine(",") {
      fields ->
        if (fields[0] != "Sample_ID") {
          samples.add([fields[0], fields[10], fields[11]])
        }
    }
    sample_csv = Channel.from(samples).toSortedList( {a, b -> a[0] <=> b[0]} ).flatten().collate( 3 )
    // Here, groupTuple should be more robust, but haven't figured out the right syntax.
    cutadapter_read1 = CUTADAPTER.out.cutadapter_read1.toSortedList( { a, b -> a.getName() <=> b.getName() } ).flatten()
    cutadapter_read2 = CUTADAPTER.out.cutadapter_read2.toSortedList( { a, b -> a.getName() <=> b.getName() } ).flatten()

    SPLITONTARGET (
      sample_csv,
      cutadapter_read1,
      cutadapter_read2
    )

    BWA_INDEX (
      ch_ref_genome
    )

    BWA_MEM (
      BWA_INDEX.out.index.collect(),
      // note1: must be a valid path type, otherwise unvalid path error;
      // note2: BWA_INDEX.out.index is a queue channel, it will be consumed, so create a value channel with collect().
      // note2-1: BWA_INDEX.out.index can be used for different processes and will be auto copied to multiple channels, but itself still contains only one element.

      params.bam_dir,
      SPLITONTARGET.out.ontarget_read1,
      SPLITONTARGET.out.ontarget_read2,
      SPLITONTARGET.out.offtarget_read1,
      SPLITONTARGET.out.offtarget_read2
    )

    ANNOTATE (
      params.genome,
      BWA_MEM.out.ontarget_both_bam,
      BWA_MEM.out.ontarget_R1only_bam,
      BWA_MEM.out.ontarget_R2only_bam,
      BWA_MEM.out.offtarget_both_bam,
      BWA_MEM.out.offtarget_R1only_bam,
      BWA_MEM.out.offtarget_R2only_bam
    )
    // TEST (
    //   umi,
    //   read1,
    //   read2
    // )


    // PARSEUMI.out.umi.onComplete {
    //   // tem_umi   = PARSEUMI.out.umi.collect().toSortedList()
    //   // tem_read1 = DEMULTIPLEX.out.read1.collect().toSortedList()
    //   // tem_read2 = DEMULTIPLEX.out.read2.collect().toSortedList()
    //   //
    //   // collapseumi_input_list = []
    //   //
    //   // // for (i in 0..(tem_umi.value.size - 1)) {
    //   // if (tem_umi.value && tem_read1.value && tem_read2.value) {
    //   //   println "tem_umi ready ..."
    //   // for (i in 0..(8 - 1)) {
    //   //   tem = [tem_umi.value[i], tem_read1.value[i], tem_read2.value[i]]
    //   //   collapseumi_input_list.add(tem)
    //   // }
    //   // collapseumi_input_ch = Channel.from(collapseumi_input_list)
    //   // } else {
    //   //   println "wait ..."
    //   // }
    //
    //   process test {
    //     echo true
    //
    //     input:
    //       tuple path(umi), path(read1), path(read2) from collapseumi_input_ch
    //
    //     "echo TEST:$umi, $read1, $read2"
    //   }
    // }

    // appedn tuple (umi_file, read1, read2) to res list and then pass each element of res to next processes
    // so that umi_file, read1, read2 are gunranteed to be matched for each sample.
    // tem_umi   = PARSEUMI.out.umi.collect().toSortedList()
    // tem_read1 = DEMULTIPLEX.out.read1.collect().toSortedList()
    // tem_read2 = DEMULTIPLEX.out.read2.collect().toSortedList()
    //
    // collapseumi_input_list = []
    //
    // // for (i in 0..(tem_umi.value.size - 1)) {
    // if (tem_umi.value && tem_read1.value && tem_read2.value) {
    //   println "tem_umi ready ..."
    // for (i in 0..(8 - 1)) {
    //   tem = [tem_umi.value[i], tem_read1.value[i], tem_read2.value[i]]
    //   collapseumi_input_list.add(tem)
    // }
    // collapseumi_input_ch = Channel.from(collapseumi_input_list)
    // } else {
    //   println "wait ..."
    // }



    // process test {
    //   echo true
    //
    //
    //   input:
    //
    //     tuple path(umi), path(read1), path(read2) from collapseumi_input_ch
    //
    //   when:
    //
    //
    //   "echo $umi, $read1, $read2, testrun"
    // }
    // COLLAPSEUMI (
    //   umi_index
    //   read1_file
    //   read2_file
    // )


    // println res

    // println "For PARSEUMI:"
    // PARSEUMI.out.umi.collect().toSortedList().view()

    // DEMULTIPLEX.out.read1.collect().toSortedList().view()
                // .filter( ~/undetermined*.fastq.gz/ )





    // Below are default:

    /*
     * SUBWORKFLOW: Read in samplesheet, validate and stage input files
     */
    // INPUT_CHECK (
    //     ch_input
    // )
    /*
     * MODULE: Run FastQC
     */
    // FASTQC (
    //     INPUT_CHECK.out.reads
    // )
    // ch_software_versions = ch_software_versions.mix(FASTQC.out.version.first().ifEmpty(null))

    /*
     * MODULE: Run FastQC
     */

    /*
     * MODULE: Pipeline reporting
     */
    // Get unique list of files containing version information
    // ch_software_versions
    //     .map { it -> if (it) [ it.baseName, it ] }
    //     .groupTuple()
    //     .map { it[1][0] }
    //     .flatten()
    //     .collect()
    //     .set { ch_software_versions }
    // GET_SOFTWARE_VERSIONS (
    //     ch_software_versions
    // )

    /*
     * MODULE: MultiQC
     */
    // workflow_summary    = Workflow.paramsSummaryMultiqc(workflow, params.summary_params)
    // ch_workflow_summary = Channel.value(workflow_summary)
    //
    // ch_multiqc_files = Channel.empty()
    // ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    // ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_files = ch_multiqc_files.mix(GET_SOFTWARE_VERSIONS.out.yaml.collect())
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

//     MULTIQC (
//         ch_multiqc_files.collect()
//     )
//     multiqc_report       = MULTIQC.out.report.toList()
//     ch_software_versions = ch_software_versions.mix(MULTIQC.out.version.ifEmpty(null))
}

////////////////////////////////////////////////////
/* --              COMPLETION EMAIL            -- */
////////////////////////////////////////////////////

workflow.onComplete {
    // Completion.email(workflow, params, params.summary_params, projectDir, log, multiqc_report)
    Completion.summary(workflow, params, log)
}

////////////////////////////////////////////////////
/* --                  THE END                 -- */
////////////////////////////////////////////////////
