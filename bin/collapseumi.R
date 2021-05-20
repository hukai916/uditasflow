# Collapse UMI reads given umi fastq file, read1 fastq and read2 fastq
setwd("/Users/kaihu/Projects/workflow/uditasflow/bin")
umi_file <- "/Users/kaihu/Projects/workflow/uditasflow/result_test/umi.fastq"

library(ShortRead)
library(purrr)
library(dplyr)
library(stringr)
library(readr)

test <- readFastq(umi_file)

# path_output_fq <- "./res_umi"
# umi_start <- 1
# umi_end 	<- 9 
# reads_per_chunk <- 10000