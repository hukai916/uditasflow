#!/usr/bin/env Rscript
# Collapse UMI reads given umi fastq file, read1 fastq and read2 fastq
library(optparse)

# umi_file <- "/Users/kaihu/Projects/workflow/uditasflow/result_test/umi_index2_S1.fastq.gz"
# read1_file <- "/Users/kaihu/Projects/workflow/uditasflow/result_test/S1_R1.fastq.gz"
# read2_file <- "/Users/kaihu/Projects/workflow/uditasflow/result_test/S1_R2.fastq.gz"
# mismatch_percent <- 0.2 # tolerate this percentage of mismatches
# path_output_fq <- "./res_collapseumi"

option_list = list(
	make_option(c("--umi_file"), type="character", default=NULL, 
							help="full path to file that contains UMI", metavar="character"),
	make_option(c("--read1_file"), type="character", default=NULL, 
							help="full path to index1 file", metavar="character"),
	make_option(c("--read2_file"), type="character", default=NULL, 
							help="full path to index2 file", metavar="character"),
	make_option(c("--mismatch_percent"), type="double", default=0.2, 
							help="percent of mismatch to tolerate for reads [default=%default]", metavar="character"),
	make_option(c("--path_output_fq"), type="character", default="./res_collapseumi", 
							help="Output directory [default=%default]", metavar="character")); 

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if (is.null(opt$umi_file) | is.null(opt$read1_file) | is.null(opt$read2_file)) {
	print_help(opt_parser)
	stop("At least the umi_file, read1_file, read2_file must be provided!")
}

library(ShortRead)
library(purrr)
library(dplyr)
library(stringr)
library(readr)
library(campfin)
library(collections)

collapse_umi <- function(umi_file, read1_file, read2_file, mismatch_percent, path_output_fq) {
	## Check if output folder already exists:
	if (file.exists(path_output_fq)) {
		stop("Result folder already exists! Please delete it first!")
	}
	if (!(dir.create(file.path(path_output_fq)))) {
		stop("can't create output folder, check write permission!")
	}
	
	umi_dict <- dict()
	
	umi_fastq <- readFastq(umi_file)
	read1_fastq <- readFastq(read1_file)
	read2_fastq <- readFastq(read2_file)
	
	umi_seq_vector <- umi_fastq@sread %>% as.vector 
	umi_fastq@id <- umi_fastq@id %>% as.vector %>% str_split("\\s+") %>% sapply(function(x) x[[1]]) %>% BStringSet()
	# by default, the readFastq may/will add " 2:N:0:0" to the end of the original id
	read1_fastq@id <- read1_fastq@id %>% as.vector %>% str_split("\\s+") %>% sapply(function(x) x[[1]]) %>% BStringSet()
	read2_fastq@id <- read2_fastq@id %>% as.vector %>% str_split("\\s+") %>% sapply(function(x) x[[1]]) %>% BStringSet()
	read1_seq_vector <- read1_fastq@sread %>% as.vector
	read2_seq_vector <- read2_fastq@sread %>% as.vector
	
	## data sanity check
	message("Data sanity check ...")
	if (any(umi_fastq@id != read1_fastq@id) | any(read1_fastq@id != read2_fastq@id)) {
		stop("UMI and read fastq files must have matching sequence names!")
	}
	message("Sanity check passed!")
	
	message("Processing starts ...")
	for (i in 1:length(umi_fastq)) {
		umi_seq <- umi_seq_vector[[i]]
		read_id <- umi_fastq@id[[i]] %>% toString
		read1_seq <- read1_seq_vector[[i]]
		read2_seq <- read2_seq_vector[[i]]
		
		if (!(umi_dict$has(umi_seq))) {
			# assign a new dict key/value element:
			tb <- tibble("id" = read_id, "read1" = read1_seq, "read2" = read2_seq, "count" = 1)
			umi_dict$set(umi_seq, tb)
		} else {
			# check if reads matches and process accordingly:
			# each umi element of a tibble consists of the following cols:
			# id	read1	read2	count
			tem <- umi_dict$get(umi_seq)
			tem$dis1 <- str_dist(tem$read1, read1_seq)
			tem$dis2 <- str_dist(tem$read2, read2_seq)
			
			match <- tem[tem$dis1 <= mismatch_percent & tem$dis2 <= mismatch_percent,]
			if (dim(match)[[1]] == 0) {
				# if no matching reads, append the current read to tb as a new row
				updated_tb <- umi_dict$get(umi_seq) %>% add_row("id" = read_id, "read1" = read1_seq, "read2" = read2_seq, "count" = 1)
				umi_dict$set(umi_seq, updated_tb)
			} else {
				# if matching reads, increase the count number of the first match by 1
				count <- tem[tem$dis1 / width(tem$read1) <= mismatch_percent & tem$dis2 / width(tem$read2) <= mismatch_percent,][1,]$count
				tem[tem$dis1 <= mismatch_percent & tem$dis2 <= mismatch_percent,][1,]$count <- count + 1
				umi_dict$set(umi_seq, select(tem, -c(dis1, dis2)))
			}
		}
		
		if ((i %% 5000) == 0) {
			message(paste0("Processed ", i, " (", format(round((i/length(umi_fastq))*100, 2), nsmall = 2), "%) reads ..."))
		}
	}
	message("All reads processed.")
	
	message("Saving output ...")
	id_to_keep <- c(NA)
	length(id_to_keep) <- length(umi_fastq) # to pre-allocate the length
	id_to_keep <- sapply(umi_dict$keys(), function(x) { umi_dict$get(x)$id } ) %>% flatten()
	id_to_keep_index <- which(umi_fastq@id %>% as.vector %in% id_to_keep)
	
	writeFastq(umi_fastq[id_to_keep_index], paste0(path_output_fq, '/', basename(umi_file)))
	writeFastq(read1_fastq[id_to_keep_index], paste0(path_output_fq, '/', basename(read1_file)))
	writeFastq(read2_fastq[id_to_keep_index], paste0(path_output_fq, '/', basename(read2_file)))
	
	# also save a umi summary table
	tem <- map(umi_dict$keys(), function(x) {
		# if (n %% 1000 == 0) {print(n)}
		list(umi = x,
				 count_total = sum(umi_dict$get(x)$count),
				 count_collapsed = length(umi_dict$get(x)$id))
		## Slow: even with pre-allocated length.
		# umi_tb[n,] <<- list(umi = x,
		# 										count_total = sum(umi_dict$get(x)$count),
		# 										count_collapsed = length(umi_dict$get(x)$id))
		# Slowest: add_row() is expensive.
		# umi_tb <<- umi_tb %>% add_row(umi = x,
		# 									 count_total = sum(umi_dict$get(x)$count),
		# 									 count_collapsed = length(umi_dict$get(x)$id))
	})
	
	umi <- map_chr(tem, function(x) {x$umi})
	count_total <- map_int(tem, function(x) {as.integer(x$count_total)})
	count_collapsed <- map_int(tem, function(x) {x$count_collapsed})
	
	umi_tb <- tibble(umi = umi, count_total = count_total, count_collapsed = count_collapsed) # Will there be a better way to convert list of list to tibble though the above one seems fast enough?
	write.table(umi_tb, file = paste0(path_output_fq, "/", "umi_summary.csv"), row.names = FALSE, sep = ",")
	
	# also save a rds for the umi_dict
	saveRDS(umi_dict, file = paste0(path_output_fq, "/", "umi_dict.rds"))
	message("All done.")
}

collapse_umi(umi_file = opt$umi_file, read1_file = opt$read1_file, read2_file = opt$read2_file, mismatch_percent = opt$mismatch_percent, path_output_fq = opt$path_output_fq)

### Development notes:
# map_dfr(tem, ~as_tibble(.)) # convert list of list to tibble
# Ref: https://stackoverflow.com/questions/45452015/how-to-convert-list-of-list-into-a-tibble-dataframe

### Obsolete codes using FastqStreamer that turns out to be 10x slower.
# umi_stream <- FastqStreamer(umi_file, 1)
# read1_stream <- FastqStreamer(read1_file, 1)
# read2_stream <- FastqStreamer(read2_file, 1)
# while (length(umi <- yield(umi_stream))) {
	## Note that yield can be very slow, especially when the original file is large, this is unexpected, I assume the yield should be very fast (O(1)), but seems not.
	# read1 <- yield(read1_stream)
	# read2 <- yield(read2_stream)
	# umi_seq <- umi@sread %>% toString
	# read_id <- umi@id %>% toString
	# read1_seq <- read1@sread %>% toString
	# read2_seq <- read2@sread %>% toString
	# 
	# if (!(umi_dict$has(umi_seq))) {
	# 	# assign a new dict key/value element:
	# 	tb <- tibble("id" = read_id, "read1" = read1_seq, "read2" = read2_seq, "count" = 1)
	# 	umi_dict$set(umi_seq, tb)
	# } else {
	# 	# check if reads matches and process accordingly:
	# 	# each umi element of a tibble consists of the following cols:
	# 	# id	read1	read2	count
	# 	tem <- umi_dict$get(umi_seq)
	# 	tem$dis1 <- str_dist(tem$read1, read1_seq)
	# 	tem$dis2 <- str_dist(tem$read2, read2_seq)
	# 	
	# 	match <- tem[tem$dis1 <= mismatch & tem$dis2 <= mismatch,]
	# 	if (dim(match)[[1]] == 0) {
	# 		# if no matching reads, append the current read to tb as a new row
	# 		updated_tb <- umi_dict$get(umi_seq) %>% add_row("id" = read_id, "read1" = read1_seq, "read2" = read2_seq, "count" = 1)
	# 		umi_dict$set(umi_seq, updated_tb)
	# 	} else {
	# 		# if matching reads, increase the count number of the first match by 1
	# 		count <- tem[tem$dis1 <= mismatch & tem$dis2 <= mismatch,][1,]$count
	# 		tem[tem$dis1 <= mismatch & tem$dis2 <= mismatch,][1,]$count <- count + 1
	# 		umi_dict$set(umi_seq, select(tem, -c(dis1, dis2)))
	# 	}
	# }
# 	n <- n + 1
# 	if ((n %% 100) == 0) {
# 		message(paste0("Processed ", n, " reads ..."))
# 	}
# }
# close(umi_stream)
# close(read1_stream)
# close(read2_stream)