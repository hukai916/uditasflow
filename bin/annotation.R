#!/usr/bin/env Rscript
# Generate peak summary table given stacked bed file with ChIPpeakAnno
library(optparse)

option_list = list(
  make_option(c("--bed_file"), type="character", default=NULL, 
              help="full path to bed file", metavar="character"),
  make_option(c("--genome"), type="character", default=NULL, 
              help="which genome to use to annotate: hg38 or mm10", metavar="character"),
  make_option(c("--path_output_dir"), type="character", default="./res_annotation", 
              help="Output directory [default=%default]", metavar="character")); 

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if (is.null(opt$bed_file) | is.null(opt$genome)) {
  print_help(opt_parser)
  stop("At least the bed file and genome value must be provided!")
}

library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

# bed_file <- "some.stacked.bed"
# genome <- "hg38"
# path_output_dir <- "./res_annotation"

get_annotation <- function(bed_file, genome, path_output_dir) {
  # bed_file <- "/Users/kaihu/Projects/uditas/Yasi/Yasi_fastq_output/bed_stack/S1.bed"
  # genome <- "hg38"
  # path_output_dir <- "."
  
  # Check if output folder already exists:
  if (file.exists(path_output_dir)) {
    warning("Result folder already exists! Will overwrite!")
  }
  if (!(dir.create(file.path(path_output_dir)))) {
    stop("can't create output folder, check write permission!")
  }
  
  # Read in the bed file
  peak <- toGRanges(bed_file, sep = "", format = "BED")
  # Nearest gene symbol
  if (genome == "hg38") {
    knownGene <- "TxDb.Hsapiens.UCSC.hg38.knownGene"
    orgAnn <- "org.Hs.eg.db"
  } else if (genome == "mm10") {
    knownGene <- "TxDb.Mmusculus.UCSC.mm10.knownGene"
    orgAnn <- "org.Mm.eg.db"
  }
  annoDataGene <- toGRanges(get(knownGene))
  gene <- annotatePeakInBatch(peak, 
                              AnnotationData = annoDataGene,
                              output = "both")
  gene <- addGeneIDs(gene, orgAnn = orgAnn,
                     feature_id_type = "entrez_id",
                     IDs2Add = c("symbol"))
  
  df <- data.frame(seqnames=seqnames(gene),
                   starts=start(gene)-1,
                   ends=end(gene),
                   names=names(gene),
                   counts=gene$thickStart,
                   scores=score(gene),
                   strands=strand(gene),
                   feature=gene$feature,
                   symbol=gene$symbol,
                   start_pos=gene$start_position,
                   end_pos=gene$end_position,
                   feature_strand=gene$feature_strand,
                   insideFeature=gene$insideFeature,
                   shortestDistance=gene$shortestDistance,
                   fromOverlappingOrNeareast=gene$fromOverlappingOrNearest
  )
  write.table(df, file=paste0(path_output_dir, "/Annotation_", basename(bed_file), ".txt", sep=""), quote=F, sep="\t", row.names=F, col.names=T)
}

get_annotation(bed_file = opt$bed_file, genome = opt$genome, path_output_dir = opt$path_output_dir)