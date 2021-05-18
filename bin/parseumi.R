#!/usr/bin/env Rscript
# Parse out UMI fastq file given index file and umi position.
library(optparse)

# index_file <- "/Users/kaihu/Projects/workflow/uditasflow/result_test/res_demultiplex/index2_S1.fastq.gz"
# path_output_fq <- "./res_umi"
# umi_start <- 1
# umi_end 	<- 9 
# reads_per_chunk <- 10000

option_list = list(
	make_option(c("--index_file"), type="character", default=NULL, 
							help="full path to index file that contains UMI", metavar="character"),
	make_option(c("--umi_start"), type="character", default=NULL, 
							help="start position of UMI seq", metavar="character"),
	make_option(c("--umi_end"), type="character", default=NULL, 
							help="end position of UMI seq", metavar="character"),
	make_option(c("--reads_per_chunk"), type="integer", default=1000000, 
							help="the number of reads to process in one chunk [default=%default]", metavar="character"),
	make_option(c("--path_output_fq"), type="character", default="./res_umi", 
							help="Output directory [default=%default]", metavar="character")); 

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if (is.null(opt$index_file) | is.null(opt$umi_start) | is.null(opt$umi_end)) {
	print_help(opt_parser)
	stop("At least the index_file, umi_start, umi_end must be provided!")
}

library(ShortRead)
library(purrr)
library(dplyr)
library(stringr)
library(readr)

parse_umi <- function(index_file, umi_start, umi_end, reads_per_chunk, path_output_fq) {
	
	## Check if output folder already exists:
	if (file.exists(path_output_fq)) {
		stop("Result folder already exists! Please delete it first!")
	}
	if (!(dir.create(file.path(opt$path_output_fq)))) {
		stop("can't create output folder, check write permission!")
	}
	read_chunk <- 1
	index_stream <- FastqStreamer(index_file, reads_per_chunk)
	while (1) {
		if (!(length(index <- yield(index_stream)))) {
			# Remove "tem_" in front of the output name by renaming the tem files:
			tem <- list.files(path_output_fq, pattern = "tem_*")
			file.rename(paste0(path_output_fq, "/", tem), paste0(path_output_fq, "/", str_sub(tem, 5, width(tem))))
			message("Parsing UMI finished!")
			close(index_stream)
			break
		} else {
			message(paste0("Processing read chunk ", read_chunk, " ..."))
			read_chunk <- read_chunk + 1
		}
		
		index@sread <- str_sub(index@sread, start = umi_start, end = umi_end) %>% DNAStringSet()
		writeFastq(index, paste0(path_output_fq, "/tem_umi_", basename(index_file)), mode = "a")
	}
}

parse_umi(index_file = opt$index_file, umi_start = opt$umi_start, umi_end = opt$umi_end, reads_per_chunk = opt$reads_per_chunk, path_output_fq = opt$path_output_fq)
