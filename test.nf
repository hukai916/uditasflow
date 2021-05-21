index = ["result_test/index1_S1", "result_test/index1_S4", "result_test/index1_S2", "result_test/index1_S3"]
fastq = ["result_test/S3", "result_test/S2", "result_test/S1", "result_test/S4"]

index_ch = Channel.fromPath(index)
fastq_ch = Channel.fromPath(fastq)

// ordered_index = Channel.from(index)
//         .toSortedList()
//
// ordered_fastq = Channel.from(fastq)
//         .toSortedList()

process test1 {
  input:
    path id from index_ch

  output:
    path id1 into index_1



}

process test2 {
  input:
    path fq from fastq_ch

  output:
    path fq1 into fastq_1
}



// println ordered_index.value.size
//

// res = []
//
// for (i in 0..(index.size - 1)) {
//   tem = [ordered_index.value[i], ordered_fastq.value[i]]
//   res.add(tem)
// }
//
//
// if (is.null)
// // println res
//
// intest = Channel.from(res)
//
// process test {
//   echo true
//
//   input:
//     tuple val(x), val(y) from intest
//
//   "echo $x $y"
// }
