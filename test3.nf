// res = []
//
// new File("/Users/kaihu/Projects/workflow/test_data/samplesheet_uditas.csv").splitEachLine(",") {
//   fields ->
//     if (fields[0] != "Sample_ID") {
//       res.add([fields[0], fields[10], fields[11]])
//     }
// }
//
//
// Channel.from(res).toSortedList().view()
nextflow.enable.dsl = 2

process splitcsv {
  echo true

  input:
    path sample_sheet

  output:
    file "sample*.txt"

  script:
    """
    echo $res > test.txt
    """
}

ch_sample = Channel.fromPath("/Users/kaihu/Projects/workflow/test_data/samplesheet_uditas.csv")

workflow {
  // sortcsv(ch_sample)

  def samples = []

  new File("/Users/kaihu/Projects/workflow/test_data/samplesheet_uditas.csv").splitEachLine(",") {
    fields ->
      if (fields[0] != "Sample_ID") {
        samples.add([fields[0], fields[10], fields[11]])
      }
  }

  test = Channel.from(samples).toSortedList( {a, b -> a[0] <=> b[0]} ).flatten().collate( 3 ).view()

  sample_file = "/Users/kaihu/Projects/workflow/test_data/samplesheet_uditas.csv"

  // def samples = [1,2,3,4,5,6]
  new File(sample_file).splitEachLine(",") {
    fields ->
      if (fields[0] != "Sample_ID") {
        samples.add([fields[0], fields[10], fields[11]])
      }
  }
  sample_csv = Channel.from(samples).toSortedList( {a, b -> a[0] <=> b[0]} ).flatten().collate( 3 )

  process runtest {
    echo true

    input:
      // val y
      tuple val(x), val(y), val(z)

    script:
      """
      echo $x, $z
      """
  }

  runtest(sample_csv)
}
