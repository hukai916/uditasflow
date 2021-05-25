#!/usr/bin/env Rscript
# split ontarget/offtarget reads given read1, read2, mismatch_percent, ref_seq
library(optparse)

# ref_seq <- "CTCCGGGGACTGCCGTGCCGGGCGGGAGACCGCCATG"
# read_file_other <- "/Users/kaihu/Projects/workflow/test_data/S1_R1.fastq"
# read_file_anchor <- "/Users/kaihu/Projects/workflow/test_data/S1_R2.fastq"
# mismatch_percent <- 0.25
# path_output_fq <- "./res_split_ontarget"

option_list = list(
	make_option(c("--ref_seq"), type="character", default=NULL, 
							help="the expected reference seq used to compare the starting seq of the reads to in order to determine if reads are on-target or not", metavar="character"),
	make_option(c("--read_file_anchor"), type="character", default=NULL, 
							help="full path to the anchoring read file", metavar="character"),
	make_option(c("--read_file_other"), type="character", default=NULL, 
							help="full path to the other read file", metavar="character"),
	make_option(c("--mismatch_percent"), type="double", default=0.25, 
							help="percent of mismatch (in regard to ref_seq) to tolerate when comparing reads to ref_seq [default=%default]", metavar="character"),
	make_option(c("--path_output_fq"), type="character", default="./res_split_ontarget", 
							help="Output directory [default=%default]", metavar="character"))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if (is.null(opt$ref_seq) | is.null(opt$read_file_anchor) | is.null(opt$read_file_other)) {
	print_help(opt_parser)
	stop("At least the ref_seq, read_file_anchor, read_file_other must be provided!")
}

library(ShortRead)
library(purrr)
library(dplyr)
library(stringr)
library(campfin)

split_ontarget <- function(ref_seq, read_file_anchor, read_file_other, mismatch_percent, path_output_fq) {
	# Check if output folder already exists:
	if (file.exists(path_output_fq)) {
		stop("Result folder already exists! Please delete it first!")
	}
	if (!(dir.create(file.path(path_output_fq)))) {
		stop("can't create output folder, check write permission!")
	}
	
	# Read in sample fastq:
	inputFastqReadAnchor <- readFastq(read_file_anchor)
	inputFastqReadOther <- readFastq(read_file_other)
	
	# Create new slot to store base difference compared to expected ref seq:
	slot_name <- c("differ_base")
	slot_content <- c("integer")
	names(slot_content) <- slot_name
	setClass(
		"ReadFastq",
		contains = "ShortReadQ",
		slots = slot_content
	) -> ReadFastq
	read_anchor_fastq <- as(inputFastqReadAnchor, "ReadFastq")
	
	# Add distance:
	message("Calculating distance ... ")
	slot(read_anchor_fastq, "differ_base") <- str_dist(str_sub(read_anchor_fastq@sread, 1, width(ref_seq)), ref_seq) %>% as.integer()
	
	token_on	<- read_anchor_fastq@differ_base <= width(ref_seq) * mismatch_percent
	token_off <- read_anchor_fastq@differ_base > width(ref_seq) * mismatch_percent
	message("Done.")
	
	# Output:
	message("Writing to output ...")
	writeFastq(read_anchor_fastq[token_on], paste0(path_output_fq, "/ontarget_", sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(read_file_anchor)), ".fastq.gz"), mode = "w")
	writeFastq(inputFastqReadOther[token_on], paste0(path_output_fq, "/ontarget_", sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(read_file_other)), ".fastq.gz"), mode = "w")
	
	writeFastq(inputFastqReadAnchor[token_off], paste0(path_output_fq, "/offtarget_", sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(read_file_anchor)), ".fastq.gz"), mode = "w")
	writeFastq(inputFastqReadOther[token_off], paste0(path_output_fq, "/offtarget_", sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(read_file_other)), ".fastq.gz"), mode = "w")
	message("Done!")	
}

split_ontarget(ref_seq = opt$ref_seq, read_file_anchor = opt$read_file_anchor, read_file_other = opt$read_file_other, mismatch_percent = opt$mismatch_percent, path_output_fq = opt$path_output_fq)