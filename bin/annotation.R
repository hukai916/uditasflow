#!/usr/bin/env Rscript
# Generate peak summary table given stacked bed file with ChIPpeakAnno

library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)


bed_file <- "some.stacked.bed"
genome <- "hg38"
path_output_dir <- "./"

# Read in the bed file
peak <- toGRanges()

samples = c("S1_proper_pair_stacked", "S2_proper_pair_stacked", "S3_proper_pair_stacked", "S1_unproper_pair_stacked", "S2_unproper_pair_stacked", "S3_unproper_pair_stacked")

for (s in samples) {
  peakS1 <- toGRanges(paste0(s,".bed", sep=""), format="BED")
  
  # Nearest gene symbol
  annoDataGene <- toGRanges(TxDb.Hsapiens.UCSC.hg38.knownGene)
  gene <- annotatePeakInBatch(peakS1, 
                              AnnotationData = annoDataGene,
                              output = "both")
  gene <- addGeneIDs(gene, orgAnn = "org.Hs.eg.db",
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
  
  write.table(df, file=paste0(s,"_summary.txt", sep=""), quote=F, sep="\t", row.names=F, col.names=T)
}