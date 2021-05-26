#!/usr/bin/env Rscript
# Draw bam fragment distribution given bam file
library(optparse)

option_list = list(
	make_option(c("--bam_file"), type="character", default=NULL, 
							help="full path to bam file (from PE reads)", metavar="character"),
	make_option(c("--path_output_dir"), type="character", default="./res_fragment_size_dist", 
							help="Output directory [default=%default]", metavar="character")); 

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

if (is.null(opt$bam_file)) {
	print_help(opt_parser)
	stop("At least the bam file must be provided!")
}

library(ATACseqQC)

get_fragment_size_dist <- function(bam_file, path_output_dir) {
	## Check if output folder already exists:
	if (file.exists(path_output_dir)) {
		warning("Result folder already exists! Will overwrite!")
	}
	if (!(dir.create(file.path(path_output_dir)))) {
		stop("can't create output folder, check write permission!")
	}
	
	bam_file.labels <- gsub(".bam", "", basename(bam_file))
	pdf(file = paste0(path_output_dir, "/", basename(bam_file), "_fragment_size.pdf"),
			width = 4,
			height = 6)
	fragSizeDist(bam_file, bam_file.labels)
	dev.off()
}

get_fragment_size_dist(bam_file = opt$bam_file, path_output_dir = opt$path_output_dir)

