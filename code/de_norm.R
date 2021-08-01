# Goal: identify differentially expressed genes and lncRNAs
library(tidyverse)
library(DESeq2)
library(furrr)
library(tictoc)
library(rtracklayer)
library(GenomicRanges)
library(ComplexHeatmap)
library(circlize)


# 1. Import transcript count  data ----------------------------------------

# Count data is made up of exon + intron counts
t98g <- readRDS("data/whole_transcript_counts.rds")
t98g <- data.frame(t98g$counts)
colnames(t98g) <- paste0("T98G_",str_pad(seq(0,400,10),width=3, side="left", pad="0"), "_min")

# Calculate size factors
sf <- estimateSizeFactorsForMatrix(t98g)



# 2. Filter genes with very low expression --------------------------------

# Require at least one time point has 10 reads detected
t98g <- t98g[rowSums(t98g>10)>0,]

# Require at least three time points in a row have more than 2 counts
keep_tx <- sapply(1:nrow(t98g), function(i) {
  x<-rle(unlist(t98g[i,])>2)
  return(any(x$lengths[x$values]>=3))}
  )
t98g <- t98g[keep_tx,]
rm(keep_tx)



# 3. Define gene types ----------------------------------------------------

# Get exons of all string-tie transcripts
gtx <- readGFF("data/annotation/T98G_stringtie_stranded.gtf") %>%
  filter(type=="exon", transcript_id%in%row.names(t98g))

# Get all protein-coding and lncRNA gene IDs
gtf <- readGFF("data/annotation/GRCh38_spiked.gtf")
pc <- filter(gtf, gene_type=="protein_coding") %>%
  pull(gene_id) %>% 
  unique()
lnc <- filter(gtf, gene_type%in%c("non_coding", "lincRNA", "macro_lncRNA", "antisense", "sense_intronic", "sense_overlapping", "3prime_overlapping_ncRNA", "bidirectional_promoter_lncRNA", "retained_intron")) %>%
  pull(gene_id) %>%
  unique()

# Identify string-tie transcripts that are protein-coding genes
pc_tx <- unique(gtx$transcript_id[gtx$ref_gene_id%in%pc])
gtx <- gtx[!gtx$transcript_id%in%pc_tx,]

# Identify string-tie transcripts that represent annotated lncRNAs
lnc_tx <- unique(gtx$transcript_id[gtx$ref_gene_id%in%lnc])
gtx <- gtx[!gtx$transcript_id%in%lnc_tx,] 

# Include transcripts with partial overlap with annotated lncRNAs
overlaps_lnc <- suppressWarnings(
overlapsAny(query = makeGRangesFromDataFrame(gtx, ignore.strand = TRUE),
            subject = makeGRangesFromDataFrame(filter(gtf, gene_id%in%lnc, type=="exon"), ignore.strand = TRUE))
)
lnc_tx <- c(lnc_tx, unique(gtx$transcript_id[overlaps_lnc]))
gtx <- gtx[!gtx$transcript_id%in%lnc_tx,] 


# Remove transcripts that overlap protein-coding gene exons
overlaps_pc <- suppressWarnings(
overlapsAny(query = makeGRangesFromDataFrame(gtx, ignore.strand = TRUE),
            subject = makeGRangesFromDataFrame(filter(gtf, gene_id%in%pc, type=="exon"), ignore.strand = TRUE))
)
overlaps_pc <- unique(gtx$transcript_id[overlaps_pc])
gtx <- gtx[!gtx$transcript_id%in%overlaps_pc,]

# Remove transcripts that overlap other annotated genes such as pseudogenes
gtx <- gtx[!complete.cases(gtx$ref_gene_id),]

# Define remaining transcripts as novel lncRNAs
novel_lnc_tx <- unique(gtx$transcript_id)

# Create data frame with gene type info
gene_types <- data.frame(
  transcript_id = c(pc_tx, lnc_tx, novel_lnc_tx),
  transcript_type = c(
    rep("protein_coding", length(pc_tx)),
    rep("annotated_lncRNA", length(lnc_tx)),
    rep("novel_lncRNA", length(novel_lnc_tx))
  )
)

# Add gene names where possible
gtx <- readGFF("data/annotation/T98G_stringtie_stranded.gtf") %>%
  filter(type=="transcript")
gene_types$gene_name <- gtx$ref_gene_name[match(gene_types$transcript_id, gtx$transcript_id)]

rm(gtf, gtx, overlaps_lnc, overlaps_pc, pc, lnc, pc_tx, lnc_tx, novel_lnc_tx)



# 4. Filter genes ---------------------------------------------------------

# Restrict counts to protein-coding transcripts and lncRNAs (annotated + novel)
t98g <- t98g[(row.names(t98g)%in%gene_types$transcript_id),]



# 5. Autocorrelation ------------------------------------------------------

# Normalise
t98g_norm <- data.frame(t98g/sf[col(t98g)])
saveRDS(t98g_norm,"data/t98g_norm.rds")
saveRDS(sf, "data/sf.rds")

# Calculate and add -log10 p-values from the Ljung-Box test to counts data frame
BoxTestFun <- function(expression_data) broom::tidy(Box.test(as.numeric(expression_data), lag=1, type = c("Ljung-Box"), fitdf = 0))
t98g_autocor <- t98g_norm %>% rowwise() %>% do(i=BoxTestFun(.)) # apply function
t98g_autocor <- unnest(t98g_autocor[],cols = c(i)) 

# Correct for multiple testing and add transcript ID
t98g_autocor <- mutate(t98g_autocor,
                        p.value.adjusted = p.adjust(p.value, method = "BH"),
                        transcript_id = row.names(t98g)) %>%
  select(transcript_id, p.value.adjusted, p.value)



# 6. Combine and filter ---------------------------------------------------

# Add gene type info to the autocorrelation results
t98g_autocor <- merge(t98g_autocor, gene_types, by="transcript_id")

# Add counts
t98g <- merge(t98g_autocor,
              rownames_to_column(t98g_norm, var = "transcript_id"),
              by = "transcript_id")

# Filter using an adjusted p-value cut-off of 0.01
t98g_filtered <- filter(t98g, p.value.adjusted<0.01)



# 7.  Autcor vs expression ------------------------------------------------

# Perform some quick sanity checks by having a look at the relationship between
# the autocor p-values and the max expression of a gene. We will use 10
# expression bins.

t98g_filtered$max_expression <- log10(rowMax(as.matrix(t98g_filtered[,str_detect(colnames(t98g_filtered),"T98G_")])))
t98g_filtered$mean_expression <- log10(rowMeans(as.matrix(t98g_filtered[,str_detect(colnames(t98g_filtered),"T98G_")])))

library(Hmisc)
t98g_filtered$expression_bin <- as.numeric(cut2(t98g_filtered$max_expression, g=10))
ggplot(t98g_filtered, aes(x=as.factor(expression_bin), y=-log10(p.value))) +
  geom_violin() + facet_wrap(~transcript_type)

boxplot(split(t98g_filtered$max_expression, t98g_filtered$expression_bin), main="log10(max expression) bins")
boxplot(split(-log10(t98g_filtered$p.value), t98g_filtered$expression_bin), main="-log10(autocor p-values) binned by expression")

# For all three transcript classes (protein-coding genes, novel and annotated
# lncRNAs) higher autocor values are observed for higher expression bins  - as
# we would expect.



# 8. Cluster --------------------------------------------------------------

# Get just expression values
gene_expression <- t98g_filtered[,str_detect(colnames(t98g_filtered),"T98G_")]

# Calculate correlation matrix
cor_matrix <- cor(t(gene_expression), method="pearson")

# Get distance matrix
cor_dist <- as.dist(1 - cor_matrix)

# Hierarchical clustering
cor_clust <- fastcluster::hclust(d=cor_dist, method="complete")

# Scale for plotting
mat_scaled = t(scale(t(gene_expression)))
colnames(mat_scaled) <- str_pad(seq(0,400,10),width=3, side="left", pad="0")

# Cluster rows with kmeans
set.seed(1234)
row_clusts <- kmeans(x=mat_scaled, centers = 30, iter.max = 100)
row_clusts <- row_clusts$cluster

# Remove clusters capturing normalisation artefacts and deactivated genes
normalisation_artefacts <- c(3,8,11,16,23,28)
deactivated_genes <- c(2,5,7,10,12,19,29)
keep <- !row_clusts%in%normalisation_artefacts & !row_clusts%in%deactivated_genes
t98g_filtered <- t98g_filtered[keep,]
mat_scaled <- mat_scaled[keep,]

# Create row annotation
row_ha = rowAnnotation(tx_type = t98g_filtered$transcript_type,
                       autocor=-log10(t98g_filtered$p.value.adjusted),
                       col=list(tx_type=c("protein_coding"="red", "novel_lncRNA"="blue", "annotated_lncRNA"="lightblue")))

# Create heatmap
h1 <- Heatmap(matrix = mat_scaled, 
        column_title = "Time points (min)",
        row_title = "Genes",
        col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        name="Z-score",
        cluster_columns = FALSE,
        show_row_names = FALSE,
        show_column_names = TRUE, column_names_side = "bottom",
        row_order = cor_clust$order[keep],
        use_raster=FALSE,
        split=row_clusts[keep],
        right_annotation = row_ha)

# Plot
draw(h1)



# 9. Combine exon and intron counts ---------------------------------------

# Restrict to columns useful for further analysis
t98g_filtered <- select(t98g_filtered,
                        transcript_id,
                        transcript_type,
                        gene_name,
                        colnames(t98g_filtered)[str_detect(colnames(t98g_filtered), "T98G")])

# Replace "T98G" with "transcript"
colnames(t98g_filtered) <- str_replace(colnames(t98g_filtered), "T98G", "transcript")

# Import other count matrices
intron <- data.frame(readRDS("data/intron_counts.rds")$counts)
first_10Kb <- data.frame(readRDS("data/introns_first_10kb_counts.rds")$counts)
last_10Kb <- data.frame(readRDS("data/introns_last_10kb_counts.rds")$counts)
exon <- data.frame(readRDS("data/stringtie_counts_transcript.rds")$counts)

# Change column names
colnames(intron) <- paste0("intron_",str_pad(seq(0,400,10),width=3, side="left", pad="0"), "_min")
colnames(first_10Kb) <- paste0("first_10Kb_",str_pad(seq(0,400,10),width=3, side="left", pad="0"), "_min")
colnames(last_10Kb) <- paste0("last_10Kb_",str_pad(seq(0,400,10),width=3, side="left", pad="0"), "_min")
colnames(exon) <- paste0("exon_",str_pad(seq(0,400,10),width=3, side="left", pad="0"), "_min")

# Change row names to transcript IDs
introns_gtf <- data.frame(readGFF("data/annotation/introns.gff3")) %>%
  select(ID, Name) %>%
  filter(ID %in% row.names(intron))
row.names(intron) <- introns_gtf$Name
row.names(first_10Kb) <- introns_gtf$Name
row.names(last_10Kb) <- introns_gtf$Name
rm(introns_gtf)

# Normalise
intron <- data.frame(intron/sf[col(intron)])
first_10Kb <- data.frame(first_10Kb/sf[col(first_10Kb)])
last_10Kb <- data.frame(last_10Kb/sf[col(last_10Kb)])
exon <- data.frame(exon/sf[col(exon)])

# Add row names as column so we can join with t98g_filtered
intron <- intron %>% rownames_to_column(var="transcript_id")
first_10Kb <- first_10Kb %>% rownames_to_column(var="transcript_id")
last_10Kb <- last_10Kb %>% rownames_to_column(var="transcript_id")
exon <- exon %>% rownames_to_column(var="transcript_id")

# Restrict to genes in t98g_filtered
intron <- intron %>% filter(transcript_id %in% t98g_filtered$transcript_id)
first_10Kb <- first_10Kb %>% filter(transcript_id %in% t98g_filtered$transcript_id)
last_10Kb <- last_10Kb %>% filter(transcript_id %in% t98g_filtered$transcript_id)
exon <- exon %>% filter(transcript_id %in% t98g_filtered$transcript_id)

# Collapse 41 columns containing counts to a single nested column
t98g_filtered <- t98g_filtered %>% nest(transcript = colnames(t98g_filtered[4:44]))
intron <- intron %>% nest(intron = colnames(intron[2:42]))
first_10Kb <- first_10Kb %>% nest(first_10Kb = colnames(first_10Kb[2:42]))
last_10Kb <- last_10Kb %>% nest(last_10Kb = colnames(last_10Kb[2:42]))
exon <- exon %>% nest(exon = colnames(exon[2:42]))

# Join data frames by transcript_id column
t98g_filtered <- list(t98g_filtered, exon, intron, first_10Kb, last_10Kb) %>%
  purrr::reduce(full_join, by = "transcript_id")



# 10. Save ----------------------------------------------------------------

saveRDS(t98g_filtered, file = "data/t98g_filtered.rds")


