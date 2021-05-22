index = ["result_test/index1_S1", "result_test/index1_S4", "result_test/index1_S2", "result_test/index1_S3"]
fastq = ["result_test/S3", "result_test/S2", "result_test/S1", "result_test/S4"]
fastq2 = ["result_test/S5", "result_test/S6", "result_test/S7", "result_test/S8"]


index_ch = Channel.fromPath(index)
fastq_ch = Channel.fromPath(fastq)
fastq_ch2 = Channel.fromPath(fastq2)

// ordered_index = Channel.fromPath(index).toSortedList()
// ordered_fastq = Channel.fromPath(fastq).toSortedList()
// ordered_fastq2 = Channel.fromPath(fastq2).toSortedList()

ordered_index = Channel.fromPath(index)
ordered_fastq = Channel.fromPath(fastq)
ordered_fastq2 = Channel.fromPath(fastq2)
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
// values = Channel.of( [1, 'alpha'], [2, 'beta'], [3, 'delta'] )

// values = index_ch.merge(fastq_ch).merge(fastq_ch2)

values = ordered_index.merge(ordered_fastq).merge(ordered_fastq2)

process tupleExample {
    echo true

    input:
    tuple path(x), path(y), path(z) from values

    """
    echo Processing $x and $y and $z
    """

}

//
// numbers = Channel.from(1,2,3)
// words = Channel.from('hello', 'ciao')
// numbers
//     .combine(words)
//     .view()
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
