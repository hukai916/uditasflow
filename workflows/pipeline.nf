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

include { TEST                  } from '../modules/local/test/main'               addParams( options: modules['test']              )

////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

// Info required for completion email and summary
def multiqc_report = []

workflow UDITASFLOW {

    ch_software_versions = Channel.empty()
    ch_bcl_raw = Channel.fromPath(params.bcl_raw)
    ch_sample_file = Channel.fromPath(params.sample_file)


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

    umi = PARSEUMI.out.umi.toSortedList().flatten()
    // read1 = DEMULTIPLEX.out.read1.toSortedList().flatten()
    // read2 = DEMULTIPLEX.out.read2.toSortedList().flatten()

    // // ch1 = Channel.fromPath(PARSEUMI.out.umi.collect()).view()
    // values = PARSEUMI.out.umi.merge(DEMULTIPLEX.out.read1).merge(DEMULTIPLEX.out.read2)
    //
    TEST (
      umi
      // read1,
      // read2
    )


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

    // COLLAPSEUMI (
    //   umi_index
    //   read1_file
    //   read2_file
    // )



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
