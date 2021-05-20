#!/usr/bin/env Rscript
# Demultiplex samples given index files and read files
library(optparse)
option_list = list(
	make_option(c("--index1_file"), type="character", default=NULL, 
							help="full path to index1 file", metavar="character"),
	make_option(c("--index2_file"), type="character", default=NULL, 
							help="full path to index2 file", metavar="character"),
	make_option(c("--read1_file"), type="character", default=NULL, 
							help="full path to index1 file", metavar="character"),
	make_option(c("--read2_file"), type="character", default=NULL, 
							help="full path to index2 file", metavar="character"),
	make_option(c("--sample_file"), type="character", default=NULL, 
							help="full path to sample sheet file", metavar="character"),
	make_option(c("--index1_mismatch"), type="character", default=1, 
							help="number of mismatch to tolerate for index1 [default=%default]", metavar="character"),
	make_option(c("--index2_mismatch"), type="character", default=1, 
							help="number of mismatch to tolerate for index2 [default=%default]", metavar="character"),
	make_option(c("--reads_per_chunk"), type="integer", default=1000000, 
							help="the number of reads to process in one chunk [default=%default]", metavar="character"),
	make_option(c("--path_output_fq"), type="character", default="./res_demultiplex", 
							help="Output directory [default=%default]", metavar="character")); 

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if (is.null(opt$index1_file) | is.null(opt$index2_file) | is.null(opt$read1_file) | is.null(opt$read2_file) | is.null(opt$sample_file)) {
	print_help(opt_parser)
	stop("At least the index1_file, index2_file, read1_file, read2_file, sample_file must be provided!")
}

library(ShortRead)
library(purrr)
library(dplyr)
library(stringr)
library(campfin)
library(readr)

# index1_file <- "/Users/kaihu/Projects/workflow/test_data/test_Undetermined_S0_L001_I1_001.fastq.gz"
# index2_file <- "/Users/kaihu/Projects/workflow/test_data/test_Undetermined_S0_L001_I2_001.fastq.gz"
# read1_file <- "/Users/kaihu/Projects/workflow/test_data/test_Undetermined_S0_L001_R1_001.fastq.gz"
# read2_file <- "/Users/kaihu/Projects/workflow/test_data/test_Undetermined_S0_L001_R2_001.fastq.gz"
# sample_file <- "/Users/kaihu/Projects/workflow/test_data/samplesheet_uditas.csv"
# path_output_fq <- "/Users/kaihu/Projects/workflow/test_data/demul"
# index1_mismatch <- 1
# index2_mismatch <- 1
# reads_per_chunk <- 10000

demultiplexer <- function(index1_file, index2_file, read1_file, read2_file, sample_file, index1_mismatch, index2_mismatch, reads_per_chunk, path_output_fq) {
	## Check if output folder already exists:
	if (file.exists(path_output_fq)) {
		stop("Result folder already exists! Please delete it first!")
	}
	if (!(dir.create(file.path(path_output_fq)))) {
		stop("can't create output folder, check write permission!")
	}
	
	index1_count <- countLines(index1_file) / 4
	index2_count <- countLines(index2_file) / 4
	read1_count <- countLines(read1_file) / 4
	read2_count <- countLines(read2_file) / 4
	
	## Data sanity check:
	### May also need to make sure the order of each record is the same across index and read files (likely unnecessary), but not implemented here. THis can be added inside the while loop by checking the names of each record for each chunk.
	message("Performing data sanity check ...")
	if (length(unique(c(index1_count, index2_count, read1_count, read2_count))) == 1) {
		message(paste0("Read record numbers are equal...", index1_count))
	} else {
		message("Read record numbers differ across index and read files, fix it before continuing:")
		message(paste0("Index1 records: ", index1_count))
		message(paste0("Index2 records: ", index2_count))
		message(paste0("Read1 records: ", read1_count))
		message(paste0("Read2 records: ", read2_count))
		stop("Read record numbers don't match!")
	}
	## Sample sheet sanity check:
	sample_tb <- read_csv(sample_file)
	if (!(dim(sample_tb)[[1]] == 8)) { stop("Sample sheet must contain 8 columns!") }
	if (!("Sample_Name" %in% colnames(sample_tb))) { stop("Sample sheet must contain 'Sample_Name' column!") }
	if (!("Sample_Name" %in% colnames(sample_tb))) { stop("Sample sheet must contain 'Sample_Name' column!") }
	if (!(length(unique(sample_tb$index1_start)) == 1)) { stop("index1_start must be the same across samples!") }
	if (!(length(unique(sample_tb$index1_end)) == 1)) { stop("index1_start must be the same across samples!") }
	if (!(length(unique(sample_tb$index2_start)) == 1)) { stop("index2_start must be the same across samples!") }
	if (!(length(unique(sample_tb$index2_end)) == 1)) { stop("index2_start must be the same across samples!") }
	
	## Parse into small chunks and do the demultiplexing:
	chunk_number <- ceiling(index1_count / reads_per_chunk)
	message(paste0("To reduce memory usage, will split into ", chunk_number, " chunk(s), each with a max of ", reads_per_chunk, " reads."))
	
	index1_stream <- FastqStreamer(index1_file, reads_per_chunk)
	index2_stream <- FastqStreamer(index2_file, reads_per_chunk)
	read1_stream <- FastqStreamer(read1_file, reads_per_chunk)
	read2_stream <- FastqStreamer(read2_file, reads_per_chunk)
	
	read_chunk <- 1
	
	## Remove pre-existing files:
	walk(c(sample_tb$Sample_Name, "Unknown"), function(x) {
		out_fastq_tem_R1 <- paste0(path_output_fq, "/tem_", x, "_R1.fastq.gz")
		out_fastq_tem_R2 <- paste0(path_output_fq, "/tem_", x, "_R2.fastq.gz")
		out_fastq_tem_index1 <- paste0(path_output_fq, "/tem_index1_", x, ".fastq.gz")
		out_fastq_tem_index2 <- paste0(path_output_fq, "/tem_index2_", x, ".fastq.gz")
		out_fastq_R1 <- paste0(path_output_fq, "/", x, "_R1.fastq.gz")
		out_fastq_R2 <- paste0(path_output_fq, "/", x, "_R2.fastq.gz")
		out_fastq_index1 <- paste0(path_output_fq, "/index1_", x, ".fastq.gz")
		out_fastq_index2 <- paste0(path_output_fq, "/index2_", x, ".fastq.gz")
		
		file_to_remove <- c(out_fastq_tem_R1, out_fastq_tem_R2, out_fastq_tem_index1, out_fastq_tem_index2, out_fastq_R1, out_fastq_R2, out_fastq_index1, out_fastq_index2)
		file.remove(file_to_remove[file.exists(file_to_remove)])
	})
	out_fastq_undetermined_tem_R1 <- paste0(path_output_fq, "/tem_undetermined_R1.fastq.gz")
	out_fastq_undetermined_tem_R2 <- paste0(path_output_fq, "/tem_undetermined_R2.fastq.gz")
	out_fastq_undetermined_R1 <- paste0(path_output_fq, "/undetermined_R1.fastq.gz")
	out_fastq_undetermined_R2 <- paste0(path_output_fq, "/undetermined_R2.fastq.gz")
	out_fastq_undetermined_tem_index1 <- paste0(path_output_fq, "/tem_undetermined_index1.fastq.gz")
	out_fastq_undetermined_tem_index2 <- paste0(path_output_fq, "/tem_undetermined_index2.fastq.gz")
	out_fastq_undetermined_index1 <- paste0(path_output_fq, "/undetermined_index1.fastq.gz")
	out_fastq_undetermined_index2 <- paste0(path_output_fq, "/undetermined_index2.fastq.gz")
	
	file_to_remove <- c(out_fastq_undetermined_tem_R1, out_fastq_undetermined_tem_R2, out_fastq_undetermined_tem_index1, out_fastq_undetermined_tem_index2, out_fastq_undetermined_R1, out_fastq_undetermined_R2, out_fastq_undetermined_index1, out_fastq_undetermined_index2)
	file.remove(file_to_remove[file.exists(file_to_remove)])
	
	## Process each chunk:
	while (1) {df
		if (!(length(index1 <- yield(index1_stream)))) {
			# Remove "tem_" in front of the output name by renaming the tem files:
			tem <- list.files(path_output_fq, pattern = "tem_*")
			file.rename(paste0(path_output_fq, "/", tem), paste0(path_output_fq, "/", str_sub(tem, 5, width(tem))))
			message("Demultiplexing finished!")
			close(index1_stream)
			close(index2_stream)
			close(read1_stream)
			close(read2_stream)
			break
		} else {
			message(paste0("Processing read chunk ", read_chunk, " ..."))
			read_chunk <- read_chunk + 1
			index2 <- yield(index2_stream)
			read1  <- yield(read1_stream)
			read2  <- yield(read2_stream)
		}
		
		# message("Preparing ReadIndex object ...")
		## Create an S4 object with the following slots: slots_from_ShortReadQ, index1_read, index2_read, index1_xxx, ..., index2_xxx, ...)
		index1_set <- sample_tb$index1 %>% str_trim(side = "both") %>% unique()
		index2_set <- sample_tb$index2 %>% str_trim(side = "both") %>% unique()
		slot_content  <- c("BStringSet", "BStringSet") # stores read index for each read
		slot_content  <- c(slot_content, rep("integer", length(index1_set)))
		slot_content  <- c(slot_content, rep("integer", length(index2_set)))
		slot_name     <- c("index1", "index2")
		slot_name     <- c(slot_name, map_chr(index1_set, function(x) { paste0("index1_", x) }))
		slot_name     <- c(slot_name, map_chr(index2_set, function(x) { paste0("index2_", x) }))
		names(slot_content) <- slot_name
		setClass(
			"ReadIndex",
			contains = "ShortReadQ",
			slots = slot_content
		) -> ReadIndex
		read_index <- as(read1, "ReadIndex")
		# message("Preparing ReadIndex object ... done!")
		
		## Fill in read index1 and index2 by parsing index1 and index2 file
		read_index@index1 <- str_sub(index1@sread, sample_tb$index1_start, sample_tb$index1_end) %>% BStringSet()
		read_index@index2 <- str_sub(index2@sread, sample_tb$index2_start, sample_tb$index2_end) %>% BStringSet()
		
		## Fill in index distance
		walk(index1_set, function(x) {
			if (width(x) != mean(width(read_index@index1))) { warn("index1 lengths don't match!") }
			index_slot <- paste0("index1_", x)
			slot(read_index, index_slot) <<- str_dist(read_index@index1, x) %>% as.integer()
		})
		walk(index2_set, function(x) {
			if (width(x) != mean(width(read_index@index2))) { warn("index2 lengths don't match!") }
			index_slot <- paste0("index2_", x)
			slot(read_index, index_slot) <<- str_dist(read_index@index2, x) %>% as.integer()
		})
		
		## Assign to each sample: allow mismatches
		# message("Output fastq files ...")
		pos_determined <- c() # stores assigned reads
		walk(sample_tb$Sample_Name, function(sample) {
			tem_index1 <- paste0("index1_", sample_tb$index1[sample_tb$Sample_ID == sample] %>% str_trim(side = "both"))
			tem_index2 <- paste0("index2_", sample_tb$index2[sample_tb$Sample_ID == sample] %>% str_trim(side = "both"))
			
			pos <- which((slot(read_index, tem_index1) <= index1_mismatch) & (slot(read_index, tem_index2) <= index2_mismatch))
			read_index_sub_R1 <- read_index[pos]
			read_index_sub_R2 <- read2[pos]
			read_index_sub_index1 <- index1[pos]
			read_index_sub_index2 <- index2[pos]
			pos_determined <<- unique(c(pos_determined, pos))
			
			writeFastq(read_index_sub_R1, paste0(path_output_fq, "/tem_", sample, "_R1.fastq.gz"), mode = "a")
			writeFastq(read_index_sub_R2, paste0(path_output_fq, "/tem_", sample, "_R2.fastq.gz"), mode = "a")
			writeFastq(read_index_sub_index1, paste0(path_output_fq, "/tem_index1_", sample, ".fastq.gz"), mode = "a")
			writeFastq(read_index_sub_index2, paste0(path_output_fq, "/tem_index2_", sample, ".fastq.gz"), mode = "a")
		})
		## Assign undetermined reads:
		writeFastq(read1[-pos_determined], paste0(path_output_fq, "/tem_undetermined_R1.fastq.gz"), mode = "a")
		writeFastq(read2[-pos_determined], paste0(path_output_fq, "/tem_undetermined_R2.fastq.gz"), mode = "a")
		writeFastq(index1[-pos_determined], paste0(path_output_fq, "/tem_undetermined_index1.fastq.gz"), mode = "a")
		writeFastq(index2[-pos_determined], paste0(path_output_fq, "/tem_undetermined_index2.fastq.gz"), mode = "a")
	}
}

demultiplexer(index1_file = opt$index1_file, index2_file = opt$index2_file, read1_file = opt$read1_file, read2_file = opt$read2_file, sample_file = opt$sample_file, index1_mismatch = opt$index1_mismatch, index2_mismatch = opt$index2_mismatch, reads_per_chunk = opt$reads_per_chunk, path_output_fq = opt$path_output_fq)
