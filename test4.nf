nextflow.enable.dsl = 2

process TESTRUN {
    echo true

    input:
      // val y
      tuple val(x), val(y), val(z)

    script:
      """
      echo $x, $z
      """
  }

workflow {

def samples = []

new File("/Users/kaihu/Projects/workflow/test_data/samplesheet_uditas.csv").splitEachLine(",") {
  fields ->
  if (fields[0] != "Sample_ID") {
    samples.add([fields[0], fields[10], fields[11]])
  }

sample_csv = Channel.from(samples).toSortedList( {a, b -> a[0] <=> b[0]} ).flatten().collate( 3 )
  }

TESTRUN(sample_csv)
}
