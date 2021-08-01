### Expression of three human protein-coding genes and adjacent lncRNAs (Figure 5)

# Load R Packages
library(tidyverse)
library(GenomicFeatures)
library(DESeq2)
library(AnnotationHub)
library(future)
library(furrr)
library(tictoc)
library(Rsamtools)
library(GenomicAlignments)

# Check to see if directory for the current figure exists and if not, create it
if(!dir.exists("figures/Fig_5")){dir.create("figures/Fig_5")}

# Load count data
counts <- readRDS("data/t98g_filtered_coords_clust.rds")

### Gene expression line plots
# Check to see if three data files capturing mRNA and lncRNA expression have already been obtained
if(sum(file.exists(str_glue('data/{c("FOS", "TGFBI", "TGIF1")}_interval_counts.rds')))!=3){
  
  # FOS
  # Get exon and intron counts
  exon_counts <- scale(unlist(counts$exon[[which(counts$gene_name=="FOS")]]))
  intron_counts <- scale(unlist(counts$first_10Kb[[which(counts$gene_name=="FOS")]]))
  # Get counts for lncRNAs
  custom_interval <- data.frame(matrix(nrow=3, ncol=41))
  custom_interval[1,] <- scale(unlist(counts$transcript[[which(counts$transcript_id=="T98G_STRG.41143.1")]]))
  custom_interval[2,] <- scale(unlist(counts$transcript[[which(counts$transcript_id=="T98G_STRG.41159.1")]]))
  custom_interval[3,] <- scale(unlist(counts$transcript[[which(counts$transcript_id=="T98G_STRG.41161.1")]]))
  # save out
  saveRDS(list(exon_counts, intron_counts, custom_interval), "data/FOS_interval_counts.rds")
  
  # TGFBI
  # Get exon and intron counts
  #exon_counts <- scale(unlist(counts$exon[[which(counts$gene_name=="TGFBI")]]))
  #intron_counts <- scale(unlist(counts$first_10Kb[[which(counts$gene_name=="FOS")]]))
  exon_counts <- scale(unlist(data.frame(readRDS("data/whole_transcript_counts.rds")$counts)["T98G_STRG.18105.8",])/readRDS("data/sf.rds"))
  intron_counts <- scale(unlist(data.frame(readRDS("data/introns_first_10kb_counts.rds")$counts)["mRNA10729",])/readRDS("data/sf.rds"))
  # Get counts for other lncRNAs
  custom_interval <- data.frame(matrix(nrow=4, ncol=41))
  custom_interval[1,] <- scale(unlist(counts$transcript[[which(counts$transcript_id=="T98G_STRG.18100.1")]]))
  custom_interval[2,] <- scale(unlist(counts$transcript[[which(counts$transcript_id=="T98G_STRG.18101.1")]]))
  custom_interval[3,] <- scale(unlist(counts$transcript[[which(counts$transcript_id=="T98G_STRG.18102.2")]]))
  custom_interval[4,] <- scale(unlist(counts$transcript[[which(counts$transcript_id=="T98G_STRG.18103.1")]]))
  # save out
  saveRDS(list(exon_counts, intron_counts, custom_interval), "data/TGFBI_interval_counts.rds")
  
  # TGIF1
  # Get exon and intron counts
  exon_counts <- scale(unlist(data.frame(readRDS("data/whole_transcript_counts.rds")$counts)["T98G_STRG.49120.8",])/readRDS("data/sf.rds"))
  intron_counts <- scale(unlist(data.frame(readRDS("data/introns_first_10kb_counts.rds")$counts)["mRNA28656",])/readRDS("data/sf.rds"))
  # Get counts for other lncRNAs
  custom_interval <- data.frame(matrix(nrow=4, ncol=41))
  custom_interval[1,] <- scale(unlist(counts$transcript[[which(counts$transcript_id=="T98G_STRG.49121.1")]]))
  custom_interval[2,] <- scale(unlist(counts$transcript[[which(counts$transcript_id=="T98G_STRG.49127.1")]]))
  custom_interval[3,] <- scale(unlist(counts$transcript[[which(counts$transcript_id=="T98G_STRG.49130.1")]]))
  custom_interval[4,] <- scale(unlist(counts$transcript[[which(counts$transcript_id=="T98G_STRG.49129.1")]]))
  
  # save out
  saveRDS(list(exon_counts, intron_counts, custom_interval), "data/TGIF1_interval_counts.rds")
  
  # Clean up
  rm(exon_counts, intron_counts, custom_interval)
}



### Human gene and lncRNA schematic plots 
#Plot schematics of the genes, lncRNAs and lncRNA splice sites.

# Check to see if schematic data has already been obtained
if(!file.exists("data/human_schematic.rds")){
  
  ## Start by getting coordinates of exons, introns and UTRs for protein-coding 
  # transcripts
  
  # Begin by loading GTF data as a TxDb object
  TxDb <- makeTxDbFromGFF("data/annotation/GRCh38_spiked.gtf", format="gtf")
  cdss <- cdsBy(TxDb, use.names=TRUE)
  exons <- exonsBy(TxDb, use.names=TRUE)
  introns <- intronsByTranscript(TxDb, use.names=TRUE)
  fiveUTRs <- fiveUTRsByTranscript(TxDb, use.names=TRUE)
  threeUTRs <- threeUTRsByTranscript(TxDb, use.names=TRUE)
  
  # Get data for the FOS transcript ENST00000303562.8 
  FOS <- rbind(as.data.frame(cdss[["ENST00000303562.8"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "exon"),
               as.data.frame(introns[["ENST00000303562.8"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "intron"),
               as.data.frame(fiveUTRs[["ENST00000303562.8"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "five_UTR"),
               as.data.frame(threeUTRs[["ENST00000303562.8"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "three_UTR")) %>% 
    mutate(gene="FOS")
  
  # Get data for the TGFBI transcript ENST00000442011.6 
  TGFBI <- rbind(as.data.frame(cdss[["ENST00000442011.6"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "exon"),
                 as.data.frame(introns[["ENST00000442011.6"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "intron"),
                 as.data.frame(fiveUTRs[["ENST00000442011.6"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "five_UTR"),
                 as.data.frame(threeUTRs[["ENST00000442011.6"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "three_UTR")) %>% 
    mutate(gene="TGFBI")
  
  # Get data for the TGIF1 transcript ENST00000343820.9 
  TGIF1 <- rbind(as.data.frame(cdss[["ENST00000343820.9"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "exon"),
                 as.data.frame(introns[["ENST00000343820.9"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "intron"),
                 as.data.frame(fiveUTRs[["ENST00000343820.9"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "five_UTR"),
                 as.data.frame(threeUTRs[["ENST00000343820.9"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "three_UTR")) %>% 
    mutate(gene="TGIF1")
  
    # Clean up
  rm(TxDb, cdss, exons, introns, fiveUTRs, threeUTRs)
  
  ## Add lncRNA data
  TxDb <- makeTxDbFromGFF("data/annotation/T98G_stringtie_stranded.gtf", format="gtf")
  exons <- exonsBy(TxDb, use.names=TRUE)
  introns <- intronsByTranscript(TxDb, use.names=TRUE)

  # Get data for the lncRNAs around FOS 
  FOS_lncRNAs <- rbind(
    rbind(as.data.frame(exons[["T98G_STRG.41143.1"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "exon", gene = "T98G_STRG.41143.1"),
               as.data.frame(introns[["T98G_STRG.41143.1"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "intron", gene = "T98G_STRG.41143.1")),
    rbind(as.data.frame(exons[["T98G_STRG.41159.1"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "exon", gene = "T98G_STRG.41159.1"),
                       as.data.frame(introns[["T98G_STRG.41159.1"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "intron", gene = "T98G_STRG.41159.1")),
    rbind(as.data.frame(exons[["T98G_STRG.41161.1"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "exon", gene = "T98G_STRG.41161.1"),
                       as.data.frame(introns[["T98G_STRG.41161.1"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "intron", gene = "T98G_STRG.41161.1"))
    )
  
  # Get data for the lncRNAs around TGFBI 
  TGFBI_lncRNAs <- rbind(
    rbind(as.data.frame(exons[["T98G_STRG.18100.1"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "exon", gene = "T98G_STRG.18100.1"),
          as.data.frame(introns[["T98G_STRG.18100.1"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "intron", gene = "T98G_STRG.18100.1")),
    rbind(as.data.frame(exons[["T98G_STRG.18101.1"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "exon", gene = "T98G_STRG.18101.1"),
          as.data.frame(introns[["T98G_STRG.18101.1"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "intron", gene = "T98G_STRG.18101.1")),
    rbind(as.data.frame(exons[["T98G_STRG.18102.2"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "exon", gene = "T98G_STRG.18102.2"),
          as.data.frame(introns[["T98G_STRG.18102.2"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "intron", gene = "T98G_STRG.18102.2")),
    rbind(as.data.frame(exons[["T98G_STRG.18103.1"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "exon", gene = "T98G_STRG.18103.1"),
          as.data.frame(introns[["T98G_STRG.18103.1"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "intron", gene = "T98G_STRG.18103.1"))
  )
  
  # Get data for the lncRNAs around TGIF1 
  TGIF1_lncRNAs <- rbind(
    rbind(as.data.frame(exons[["T98G_STRG.49121.1"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "exon", gene = "T98G_STRG.49121.1"),
          as.data.frame(introns[["T98G_STRG.49121.1"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "intron", gene = "T98G_STRG.49121.1")),
    rbind(as.data.frame(exons[["T98G_STRG.49127.1"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "exon", gene = "T98G_STRG.49127.1"),
          as.data.frame(introns[["T98G_STRG.49127.1"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "intron", gene = "T98G_STRG.49127.1")),
    rbind(as.data.frame(exons[["T98G_STRG.49130.1"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "exon", gene = "T98G_STRG.49130.1"),
          as.data.frame(introns[["T98G_STRG.49130.1"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "intron", gene = "T98G_STRG.49130.1")),
    rbind(as.data.frame(exons[["T98G_STRG.49129.1"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "exon", gene = "T98G_STRG.49129.1"),
          as.data.frame(introns[["T98G_STRG.49129.1"]]) %>% dplyr::select(seqnames, start, end, strand) %>% mutate(feature = "intron", gene = "T98G_STRG.49129.1"))
  )

  # Add to protein-coding gene data  
  FOS <- rbind(FOS, FOS_lncRNAs)
  TGFBI <- rbind(TGFBI, TGFBI_lncRNAs)
  TGIF1 <- rbind(TGIF1, TGIF1_lncRNAs)
  
  # Save out
  saveRDS(list(FOS, TGFBI, TGIF1), "data/human_schematic.rds")
  
  # Clean up
  rm(FOS, TGFBI, TGIF1, FOS_lncRNAs, TGFBI_lncRNAs, TGIF1_lncRNAs, exons, introns, TxDb)
}



### Human epigenomics data
#Include human epigenomics data from the Roadmap Epigenomics Project.

# Get annotation hub
ahub <- AnnotationHub()

# Extract BigWig records from EpigenomeRoadMap
roadmap.bw <- subset(ahub, rdataclass == "BigWigFile" & preparerclass == "EpigenomeRoadMapPreparer")
# Create a dataframe containing the AnnotationHub unique identifiers (e.g. AH12345), RoadmapEpigenome name (e.g. E001), target of record the (e.g. H3K27ac), the URL of the resource and the type of BigWigFile (e.g. log10(p-value) or fold enrichment signal tracks)

# Then restrict to these entries;
# BigWigFiles from samples in the uniformly re-processed 111 reference human epigenomes. Note: E060 and E064 aren't present, so we are considering E001-E113.
# DNase-seq data and ChIP-seq data for the following histone modifications; H3K27ac, H3K9ac, H3K4me3, H3K4me1 or H3K36me3

# Imputed pvalue signal tracks for ChIP-seq data
roadmap.bw <- data.frame(id = names(roadmap.bw),
                         epigenome = sapply(strsplit(mcols(roadmap.bw)$title,"-"), `[`, 1),
                         target = str_extract(mcols(roadmap.bw)$title, "DNase|H3K27ac|H3K9ac|H3K4me3|H3K4me1|H3K36me3"),
                         sourceURL = mcols(roadmap.bw)$sourceurl,
                         track_type = str_extract(mcols(roadmap.bw)$title, "fc.signal|imputed.pval.signal|pval.signal|DNase.pval.signal"), stringsAsFactors = FALSE) %>%
  dplyr::filter(epigenome%in%str_glue("E{str_pad(1:113,3,'left',0)}") & target%in%c("DNase", "H3K27ac", "H3K9ac", "H3K4me3", "H3K4me1", "H3K36me3") & track_type%in%c("imputed.pval.signal","DNase.pval.signal"))
head(roadmap.bw)
table(roadmap.bw$target)

# Define ranges to extract data from
hg38 <- suppressWarnings(c(GRanges(str_remove_all("chr14:75,256,000-75,302,000",",")), 
                           GRanges(str_remove_all("chr5:135,989,000-136,080,000",",")),
                           GRanges(str_remove_all("chr18:3,435,000-3,671,000",","))))

# Give names
names(hg38) <- c("FOS", "TGFBI", "TGIF1")

# Get chain file to map hg38 to hg19
chainfiles <- query(ahub , c("hg38", "hg19", "chainfile"))
chainfiles
# Get hg38ToHg19.over.chain.gz
chain <- chainfiles[['AH14108']]
# Lift coordinates for specified range
hg19 <- unlist(liftOver(hg38, chain))
# Width of the region (in Kb)
width(hg19)/1000

# Check to see if Roadmap Epigenomics Project bigwig data has already been imported
if(!file.exists("data/Roadmap_bw_data.rds")){
  # For each Roadmap Epigenomics file, import and store the BigWig data for the specified region - run in parallel
  plan(multiprocess)
  tic()
  bw.region <- future_map(1:nrow(roadmap.bw),  .f = function(x) import(con = roadmap.bw$sourceURL[[x]], format = "bigWig", selection = hg19, as="NumericList"), .progress = TRUE)
  toc()
  x<-import(con = roadmap.bw$sourceURL[[1]], format = "bigWig", selection = hg19, as="NumericList")
  # Save out
  saveRDS(bw.region, "data/Roadmap_bw_data.rds")
  rm(bw.region)
}

# Call in data
roadmap.region <- readRDS("data/Roadmap_bw_data.rds")

# Add data to roadmap.bw data frame
roadmap.bw <- mutate(roadmap.bw,
                     FOS = lapply(roadmap.region, function(x) x[[1]]),
                     TGFBI = lapply(roadmap.region, function(x) x[[5]]),
                     TGIF1 = lapply(roadmap.region, function(x) c(x[[2]], x[[3]], x[[4]])))
saveRDS(roadmap.bw, "data/roadmap_bw.rds")

print("Finished gathering data to plot Fig_5")
