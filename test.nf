index = ["index1_S1", "index1_S4", "index1_S2", "index1_S3"]
fastq = ["S3", "S2", "S1", "S4"]



ordered_index = Channel.from(index)
        .toSortedList()

ordered_fastq = Channel.from(fastq)
        .toSortedList()

res = []

for (i in 0..(index.size - 1)) {
  tem = [ordered_index.value[i], ordered_fastq.value[i]]
  res.add(tem)
}

println res

Channel.from(res).view()
