library(optparse)
option_list = list(
	make_option(c("--path_input_fq"), type="character", default =NULL, 
							help="full path to input fastq file", metavar="character"),
	make_option(c("--path_output_fq"), type="character", default ="./", 
							help="full path to output fastq file [default = %default]", metavar="character"),
	make_option(c("--nt_to_check"), type="character", default ="T", 
							help="which nt to compare to, default to 'T' [default = %default]", metavar="character"),
	make_option(c("--nt_number"), type="integer", default =-1, 
							help="how deep to examine: default to -1 ('max_seq_length') [default = %default]", metavar="character"),
	make_option(c("--error_rate"), type="double", default ="0.1", 
							help="calculated as tatol mismatch / (examined length - 1) [default = %default]", metavar="character"),
	make_option(c("--allow_mismatch_min"), type="integer", default =1, 
							help="tolerate at least this number of mismatch, default to 1 [default = %default]", metavar="character"),
	make_option(c("--min_nt_to_keep"), type="integer", default =0, 
							help="reads must have at least this number of starting poly(nt) to be kept [default = %default]", metavar="character"),
	make_option(c("--min_nt_to_trim"), type="integer", default =0, 
							help="reads must have at least this number of starting poly(nt) to initiate trimming. i.e. if the leading poly(nt) is less than min_nt_to_trim, don't trim that read. [default = %default]", metavar="character"),
	make_option(c("--reads_per_chunk"), type="integer", default =1000000, 
							help="the number of reads to process in one chunk [default = %default]", metavar="character"),
	make_option(c("--flip_end"), type="logical", default =FALSE, 
							help="this is for trimming trailing poly(nt), flip 5' and 3' end, do regular trim leading, then flip back. [default = %default]", metavar="character")
); 

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if (is.null(opt$path_input_fq)) {
	print_help(opt_parser)
	stop("At least the input fastq must be provided!")
}

library(ShortRead)
library(purrr)
library(dplyr)
library(stringr)
library(stringi)
library(campfin)

# Function to trim the leading/trailing poly(nt):
	# path_input_fq: path to fastq file
	# nt_to_check: which nt to compare to, default to "T"
	# nt_number: how deep to examine: default to -1 ("max_seq_length")
	# error_rate: calculated as tatol mismatch / (examined length - 1)
	# allow_mismatch_min: tolerate at least this number of mismatch, default to 1
	# min_nt_to_keep: reads must have at least this number of starting poly(nt) to be kept
	# min_nt_to_trim: reads must have at least this number of starting poly(nt) to initiate trimming. i.e. if the leading poly(nt) is less than min_nt_to_trim, don't trim that read.
	# flip_end: this is for trimming trailing poly(nt), flip 5' and 3' end, do regular trim leading, then flip back.
cut_start <- function(path_input_fq, path_output_fq = "./", nt_to_check = "T", nt_number = -1, error_rate = 0.1, allow_mismatch_min = 1, min_nt_to_keep = 0, min_nt_to_trim = 0, reads_per_chunk = 1000000, flip_end = FALSE) {
	read_total_number <- countLines(path_input_fq) / 4
	print(paste0("Total reads: ", read_total_number, "."))
	read_chunk_number <- ceiling(read_total_number / reads_per_chunk)
	print(paste0("To reduce memory usage, will split into ", read_chunk_number, " chunk(s), each with ", reads_per_chunk, " reads."))
	print("To further reduce memory usage, use smaller number of reads per chunk.")

	param_vector <- c(path_input_fq, path_output_fq, nt_to_check, nt_number, error_rate, allow_mismatch_min, min_nt_to_keep, min_nt_to_trim, reads_per_chunk, flip_end)
	param_vector_name <- c("path_input_fq", "path_output_fq", "nt_to_check", "nt_number", "error_rate", "allow_mismatch_min", "min_nt_to_keep", "min_nt_to_trim", "reads_per_chunk", "flip_end")
	
	cat("\n")
	print("Selected params:")
	walk2(param_vector, param_vector_name, function(.x, .y) {
		print(paste(.y, " : ", .x))
	})
	cat("\n")
		
	f <- FastqStreamer(path_input_fq, reads_per_chunk)
	read_chunk <- 1
	out_fastq <- paste0(path_output_fq, "/cutStart_", nt_to_check, "_", basename(path_input_fq))
	if (file.exists(out_fastq)) {
		file.remove(out_fastq)
	}
	
	while (length(fq <- yield(f))) {
		print(paste0("Processing read chunk ", read_chunk, " ..."))
		
		# Find out the max nt length to inspect
		if (nt_number == -1) {
			max_width <- max(width(fq@sread))
		} else {
			max_width <- as.integer(nt_number)
		}
		
		# Create max_width of new slots to store distance to the target seq (e.g. TTTT)
		slot_name <- c("final", paste0("start:::", rep(1:max_width)))
		slot_content  <- rep("integer", max_width + 1)
		names(slot_content) <- slot_name
		setClass(
			"CutStart",
			contains = "ShortReadQ",
			slots = slot_content
		) -> CutStart
		read_cutstart <- as(fq, "CutStart")
		if (flip_end) {
			read_cutstart@sread <- stri_reverse(read_cutstart@sread) %>% DNAStringSet()
			read_cutstart@quality@quality <- stri_reverse(read_cutstart@quality@quality) %>% BStringSet()
		}
		
		# Add distance to each slot
		print(paste0("Calculating string distance to poly(", nt_to_check, ") ..."))
		walk(slot_name, function(x) {
			if (str_detect(x, ":::")) {
				sub_pos <- as.integer(str_split(x, pattern = ":::")[[1]][2])
				if (sub_pos %% 10 == 1) {
					print(paste0("Processing ", x, " ... "))
				}
				
				tem_vec <- rep(0, max_width)
				pos_more <- which(width(slot(read_cutstart, "sread")) >= sub_pos) 
				# for reads that have the same length as sub_pos
				pos_less <- which(width(slot(read_cutstart, "sread")) < sub_pos)
				
				tem_vec[pos_more] <- str_dist(substr(slot(read_cutstart, "sread")[pos_more], sub_pos, sub_pos), nt_to_check)
				tem_vec[pos_less] <- 2 # assign 2 as distance as it doesn't have enough length
				
				slot(read_cutstart, x) <<- tem_vec %>% as.integer()
			}
		})
		print("Done!")
		
		# Assign the proper position for the final read starts pos
		print("Calculating final start position ...")
		walk(slot_name, function(x) {
			#x <- "start:::10"
			if (str_detect(x, ":::")) {
				sub_pos <- as.integer(str_split(x, pattern = ":::")[[1]][2])
				if (sub_pos %% 10 == 1) {
					print(paste0("Assigning final position ", x, " ... "))
				}
				
				if (sub_pos == 1) {
					# Set default to 0, and if first is specified nucleotide, set it to 1
					default_value <- rep(0L, length(slot(read_cutstart, "sread")))
					# subseq() fails if read length is 0, therefore, add a pre-filter first:
					pos_nonempty <- which(width(slot(read_cutstart, "sread")) > 0)
					default_value[which(subseq(slot(read_cutstart, "sread")[pos_nonempty], 1, 1) == nt_to_check)] <- 1L
					read_cutstart@final <<- default_value
				} else {
					need_update <- read_cutstart@final
					
					# The length must be greater than sub_pos to update:
					need_update_length <- which(width(slot(read_cutstart, "sread")) >= sub_pos)
					
					# The last character must be "T"
					need_update_lastchar <- which(slot(read_cutstart, x) == 0)
					
					# Allow error rate of 0.1 for leading sequences, also at least allow one mismatch
					slot_set <- paste0("start:::", seq(1, sub_pos-1))
					tem <- vector(mode = "integer", length = length(slot(read_cutstart, "sread")))
					walk(slot_set, function(y) {
						tem <<- tem + slot(read_cutstart, y)
					})
					
					tem <- tem - allow_mismatch_min # allow at least this number of mismatch
					tem <- tem / sub_pos
					need_update_error <- which(tem <= error_rate)
					need_update_pos <- intersect(need_update_length, need_update_lastchar) %>% intersect(need_update_error)
					
					need_update[need_update_pos] <- sub_pos
					read_cutstart@final <<- need_update %>% as.integer()
				}
			}
		})
		
		# Collect the final subset of reads to keep
		# Trimmed seq:
		print("Trimming reads ...")
		start_vector <- slot(read_cutstart, "final") + 1
		start_vector[start_vector <= min_nt_to_trim] <- 1
		end_vector <- width(slot(read_cutstart, "sread"))
		# If start > end, change start to end:
		change_start <- which(start_vector >= end_vector)
		start_vector[change_start] <- end_vector[change_start]
		
		# Need to deal with empty reads
		pos_nonempty <- which(width(slot(read_cutstart, "sread")) > 0)
		tem_sread <- read_cutstart@sread
		tem_sread[pos_nonempty] <- subseq(slot(read_cutstart, "sread")[pos_nonempty], start = start_vector[pos_nonempty], end = end_vector[pos_nonempty])
		read_cutstart@sread <- tem_sread
		tem_quality <- read_cutstart@quality@quality
		tem_quality[pos_nonempty] <- subseq(read_cutstart@quality@quality[pos_nonempty], start = start_vector[pos_nonempty], end = end_vector[pos_nonempty])
		read_cutstart@quality@quality <- tem_quality
		print("Trimming done!")
		
		# Must starts with >= 20 "T":
		reads_to_keep <- which(slot(read_cutstart, "final") >= min_nt_to_keep)
		#reads_to_drop <- which(slot(read_cutstart, "final") < min_nt_to_keep)
		read_cutstart@id <- read_cutstart@id[reads_to_keep]
		read_cutstart@sread <- read_cutstart@sread[reads_to_keep]
		read_cutstart@quality <- read_cutstart@quality[reads_to_keep]
		print("Saving output ...")
		out_fastq <- paste0(path_output_fq, "/cutStart_", nt_to_check, "_", basename(path_input_fq))
		if (flip_end) {
			read_cutstart@sread <- stri_reverse(read_cutstart@sread) %>% DNAStringSet()
			read_cutstart@quality@quality <- stri_reverse(read_cutstart@quality@quality) %>% BStringSet()
		}
		writeFastq(read_cutstart, out_fastq, mode = "a")
		
		print(paste0("Read chunk ", read_chunk, "/", read_chunk_number, " done!"))
		read_chunk <- read_chunk + 1
		print(gc())
	}
	close(f)
	
	print("All done!")
}

cut_start(path_input_fq = opt$path_input_fq, path_output_fq = opt$path_output_fq, nt_to_check = opt$nt_to_check, nt_number = opt$nt_number, error_rate = opt$error_rate, allow_mismatch_min = opt$allow_mismatch_min, min_nt_to_keep = opt$min_nt_to_keep, min_nt_to_trim = opt$min_nt_to_trim, reads_per_chunk = opt$reads_per_chunk, flip_end = opt$flip_end)

example <- function() {
	# Examples:
	
	# Trim leading poly(T): drop reads that have less than 20 leading Ts; trim off Ts as long as they appear in the beginning
	reads_per_chunk <- 10000
	nt_to_check <- "T"
	input_fastq <- "/Users/kaihu/Projects/cutendsr/test_R2.fq.gz"
	path_output_fq <- "./"
	cut_start(path_input_fq = input_fastq, path_output_fq = path_output_fq, error_rate = 0.1, reads_per_chunk = reads_per_chunk, min_nt_to_keep = 20, min_nt_to_trim = 0, nt_to_check = nt_to_check, flip_end = FALSE)
	
	# Trim trailing poly(A): keep all reads; only trim if the trailing A is more than 10.
	reads_per_chunk <- 10000
	nt_to_check <- "A"
	input_fastq <- "/Users/kaihu/Projects/cutendsr/test_R1.fq.gz"
	path_output_fq <- "./"
	cut_start(path_input_fq = input_fastq, path_output_fq = path_output_fq, error_rate = 0.1, reads_per_chunk = reads_per_chunk, min_nt_to_keep = 0, min_nt_to_trim = 10, nt_to_check = nt_to_check, flip_end = TRUE)
}