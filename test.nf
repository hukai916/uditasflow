index = ["result_test/index1_S1", "result_test/index1_S4", "result_test/index1_S2", "result_test/index1_S3"]
fastq = ["result_test/S3", "result_test/S2", "result_test/S1", "result_test/S4"]
fastq2 = ["result_test/S3", "result_test/S2", "result_test/S1", "result_test/S4"]


index_ch = Channel.fromPath(index)
fastq_ch = Channel.fromPath(fastq)
fastq_ch2 = Channel.fromPath(fastq2)

ordered_index = Channel.fromPath(index).toSortedList()

// Below works fine:
// process runtest1 {
//   echo true
//
//   input:
//   path x from index_ch.toSortedList()
//   path y from fastq_ch.toSortedList()
//
//
//
//   """
//   echo "TEST1 $x-$y"
//   """
// }


ch1 = Channel.from(index_ch.toSortedList())
ch2 = Channel.from(fastq_ch.toSortedList())
ch1.combine(ch2).view()

numbers = Channel.from(1,2,3)
words = Channel.from('hello', 'ciao')
numbers
    .combine(words)
    .view()
// Below won't work: won't output sorted value: you can't use collect together with toSortedLIst, the latter already does collect role
// process runtest2 {
//   echo true
//
//   input:
//   val x from index_ch.collect().toSortedList()
//
//   """
//   echo "TEST2 $x"
//   """
// }

// ordered_index = Channel.from(index)
//         .toSortedList()
//
// ordered_fastq = Channel.from(fastq)
//         .toSortedList()

// process test1 {
//   input:
//     path id from index_ch
//
//   output:
//     path id1 into index_1
//
//   script:
//   "
//   cp $id id1
//   "
//
//
// }
//
// process test2 {
//   input:
//     path fq from fastq_ch
//
//   output:
//     path fq1 into fastq_1
// }
//


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
