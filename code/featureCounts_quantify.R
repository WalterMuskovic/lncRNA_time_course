# Quantify gene expression with featureCounts

# Load R packages
library(tidyverse)
library(Rsubread)

# GRCh38_spiked.gtf -------------------------------------------------------

# Check to see if quantification for genes is already done
if(!file.exists("data/gene_counts.rds")){
  # Get list of all files
  bam_files <- str_glue('data/alignments/T98G_{str_pad(seq(0,400,10), 3, "left",0)}_minAligned.sortedByCoord.out.bam') %>% as.character()
  # Get counts for custom intervals
  gene_counts <- featureCounts(files = bam_files,
                               annot.ext = "data/annotation/GRCh38_spiked.gtf", isGTFAnnotationFile = TRUE, GTF.featureType = "exon", GTF.attrType = "gene_id",
                               isPairedEnd = TRUE,
                               minFragLength = 20,
                               maxFragLength = 5000,
                               strandSpecific = 2,
                               countMultiMappingReads = TRUE,
                               allowMultiOverlap = TRUE,
                               nthreads = 32)
  # Save out counts
  saveRDS(gene_counts, "data/gene_counts.rds")
}



# Stringtie annotation ----------------------------------------------------

# Gene-level
if(!file.exists("data/stringtie_counts_gene.rds")){
  # Get list of all files
  bam_files <- str_glue('data/alignments/T98G_{str_pad(seq(0,400,10), 3, "left",0)}_minAligned.sortedByCoord.out.bam') %>% as.character()
  # Get counts for custom intervals
  stringtie_counts <- featureCounts(files = bam_files,
                                 annot.ext = "data/annotation/T98G_stringtie_stranded.gtf", isGTFAnnotationFile = TRUE, GTF.featureType = "exon", GTF.attrType = "gene_id",
                                 isPairedEnd = TRUE,
                                 minFragLength = 20,
                                 maxFragLength = 5000,
                                 strandSpecific = 2,
                                 countMultiMappingReads = TRUE,
                                 allowMultiOverlap = TRUE,
                                 nthreads = 32)
  # Save out counts
  saveRDS(stringtie_counts, "data/stringtie_counts_gene.rds")
}

# Transcript-level
if(!file.exists("data/stringtie_counts_transcript.rds")){
  # Get list of all files
  bam_files <- str_glue('data/alignments/T98G_{str_pad(seq(0,400,10), 3, "left",0)}_minAligned.sortedByCoord.out.bam') %>% as.character()
  # Get counts for custom intervals
  stringtie_counts <- featureCounts(files = bam_files,
                                    annot.ext = "data/annotation/T98G_stringtie_stranded.gtf", isGTFAnnotationFile = TRUE, GTF.featureType = "exon", GTF.attrType = "transcript_id",
                                    isPairedEnd = TRUE,
                                    minFragLength = 20,
                                    maxFragLength = 5000,
                                    strandSpecific = 2,
                                    countMultiMappingReads = TRUE,
                                    allowMultiOverlap = TRUE,
                                    nthreads = 32)
  # Save out counts
  saveRDS(stringtie_counts, "data/stringtie_counts_transcript.rds")
}


# Full intron length ------------------------------------------------------

# Check if quantification for stringtie full intron intervals is already done
if(!file.exists("data/intron_counts.rds")){
  # Get list of all files
  bam_files <- str_glue('data/alignments/T98G_{str_pad(seq(0,400,10), 3, "left",0)}_minAligned.sortedByCoord.out.bam') %>% as.character()
  # Get counts for custom intervals
  intron_counts <- featureCounts(files = bam_files,
                               annot.ext = "data/annotation/introns.gff3", isGTFAnnotationFile = TRUE, GTF.featureType = "exon", GTF.attrType = "Parent",
                               isPairedEnd = TRUE,
                               minFragLength = 20,
                               maxFragLength = 5000,
                               strandSpecific = 2,
                               countMultiMappingReads = TRUE,
                               allowMultiOverlap = TRUE,
                               nthreads = 32)
  # Save out counts
  saveRDS(intron_counts, "data/intron_counts.rds")
}



# Last 10Kb of intron -----------------------------------------------------

# Check if quantification is already done
if(!file.exists("data/introns_last_10kb_counts.rds")){
  # Get list of all files
  bam_files <- str_glue('data/alignments/T98G_{str_pad(seq(0,400,10), 3, "left",0)}_minAligned.sortedByCoord.out.bam') %>% as.character()
  # Get counts for custom intervals
  intron_counts <- featureCounts(files = bam_files,
                                 annot.ext = "data/annotation/introns_last_10kb.gff3", isGTFAnnotationFile = TRUE, GTF.featureType = "exon", GTF.attrType = "Parent",
                                 isPairedEnd = TRUE,
                                 minFragLength = 20,
                                 maxFragLength = 5000,
                                 strandSpecific = 2,
                                 countMultiMappingReads = TRUE,
                                 allowMultiOverlap = TRUE,
                                 nthreads = 32)
  # Save out counts
  saveRDS(intron_counts, "data/introns_last_10kb_counts.rds")
}



# First 10Kb of intron ----------------------------------------------------

# Check if quantification is already done
if(!file.exists("data/introns_first_10kb_counts.rds")){
  # Get list of all files
  bam_files <- str_glue('data/alignments/T98G_{str_pad(seq(0,400,10), 3, "left",0)}_minAligned.sortedByCoord.out.bam') %>% as.character()
  # Get counts for custom intervals
  intron_counts <- featureCounts(files = bam_files,
                                 annot.ext = "data/annotation/introns_first_10kb.gff3", isGTFAnnotationFile = TRUE, GTF.featureType = "exon", GTF.attrType = "Parent",
                                 isPairedEnd = TRUE,
                                 minFragLength = 20,
                                 maxFragLength = 5000,
                                 strandSpecific = 2,
                                 countMultiMappingReads = TRUE,
                                 allowMultiOverlap = TRUE,
                                 nthreads = 32)
  # Save out counts
  saveRDS(intron_counts, "data/introns_first_10kb_counts.rds")
}



# Whole transcript interval -----------------------------------------------

# Check if quantification is already done
if(!file.exists("data/whole_transcript_counts.rds")){
  # Get list of all files
  bam_files <- str_glue('data/alignments/T98G_{str_pad(seq(0,400,10), 3, "left",0)}_minAligned.sortedByCoord.out.bam') %>% as.character()
  # Get counts for custom intervals
  intron_counts <- featureCounts(files = bam_files,
                                 annot.ext = "data/annotation/whole_tx.gff3",
                                 isGTFAnnotationFile = TRUE,
                                 GTF.featureType = "sequence_feature", GTF.attrType = "tx_name",
                                 isPairedEnd = TRUE,
                                 minFragLength = 20,
                                 maxFragLength = 5000,
                                 strandSpecific = 2,
                                 countMultiMappingReads = TRUE,
                                 allowMultiOverlap = TRUE,
                                 nthreads = 32)
  # Save out counts
  saveRDS(intron_counts, "data/whole_transcript_counts.rds")
}


print("Finished quantifying gene expression with featureCounts")
