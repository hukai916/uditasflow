index = ["index1_S1", "index1_S4", "index1_S2", "index1_S3"]
fastq = ["S3", "S2", "S1", "S4"]

Channel.from(index)
        .toSortedList()
        .view()

Channel.from(fastq)
        .toSortedList()
        .view()
