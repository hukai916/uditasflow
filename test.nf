index = ["index1_S1", "index1_S4", "index1_S2", "index1_S3"]
fastq = ["S3", "S2", "S1", "S4"]



ordered_index = Channel.from(index)
        .toSortedList()

ordered_fastq = Channel.from(fastq)
        .toSortedList()

// println ordered_index.value.size
//

res = []

for (i in 0..(index.size - 1)) {
  tem = [ordered_index.value[i], ordered_fastq.value[i]]
  res.add(tem)
}

// println res

intest = Channel.from(res)

process test {
  echo true

  input:
    tuple val(x), val(y) from intest

  "echo $x $y"
}
