## The effects of gene length on transcription time (Figure 3)
# After examining the effect of transcript stability on mRNA expression dynamics, we want to look
# at the effects of gene length.

### CACNA1C transcription
# We'll first look at one gene in detail, the relatively long calcium channel gene CACNA1C. We start
# by defining a function that accepts a transcript id as input as well as a window_width. The function
# will return the coverage across the specified windows for all introns of the transcript, across all
# time points. Rows of the returned data frame correspond to sequential genomic windows (5' -> 3'),
# columns correspond to time points.

# Load R packages
library(tidyverse)
library(GenomicFeatures)
library(GenomicAlignments)
library(rtracklayer)
library(DESeq2)

if(!file.exists("data/CACNA1C.rds")){
  
  # Begin by loading GTF data to get intron and exon intervals
  TxDb <- makeTxDbFromGFF("data/annotation/GRCh38_spiked.gtf", format="gtf")
  introns <- intronsByTranscript(TxDb, use.names=TRUE)
  
  # Create BamFileList
  fls <- BamFileList(file=paste0("data/alignments/T98G_", stringr::str_pad(seq(0,400,10), width=3, side="left", pad="0"), "_minAligned.sortedByCoord.out.bam"),
                     index=paste0("data/alignments/T98G_", stringr::str_pad(seq(0,400,10), width=3, side="left", pad="0"), "_minAligned.sortedByCoord.out.bam.bai"))
  
  # Get size factors for normalization
  sfs <- estimateSizeFactorsForMatrix(readRDS("data/gene_counts.rds")$counts)
  
  # Define function to calculate overlap with interval, using all reads
  allReadsCounts <- function(i, interval){
    suppressWarnings(countOverlaps(interval,
                                   readGAlignmentPairs(file=fls[[i]], use.names=TRUE, param=ScanBamParam(which = interval, simpleCigar=TRUE), with.which_label=FALSE, strandMode=2),minoverlap=5))
  }
  
  get_intron_cov <- function(tx_id, window_width=1000){
    
    # Get intervals of width ~1Kb, covering the desired transcript
    intervals <-unlist(tile(introns[[tx_id]], width=window_width))
    
    # Count overlaps from BAM files 
    # Load furrr package to speed up running time
    library(furrr);plan(multiprocess)
    # In the resulting data frame, rows correspond to intervals, columns to time points
    allCounts <- data.frame(future_map(1:41, ~ allReadsCounts(., interval=intervals), .progress = TRUE))
    # Give meaningful column names
    colnames(allCounts) <- paste0("T98G_",str_pad(seq(0,400,10),width=3, side="left", pad="0"), "_min")
    
    # normalize 
    allCounts <- allCounts/sfs[col(allCounts)]
    
    # Return the data frame of intron coverage
    return(allCounts)
  }
  
  # Get the intron coverage for selected transcripts and save out, as takes a while
  saveRDS(get_intron_cov("ENST00000399655.5"), file="data/CACNA1C.rds")
}

print("Finished creating data required for Fig_3a")
