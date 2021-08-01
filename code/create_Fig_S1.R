### Fig S1 - Protein-coding gene and lncRNA length

# Load R Packages
library(tidyverse)
library(rtracklayer)

#Import counts
counts <- readRDS("data/t98g_filtered_coords_clust.rds")

# Import FANTOM annotation
FANTOM <- readGFF("data/FANTOM_CAT.lv3_robust.gtf.gz", version = 2) 
# add geneSuperClass to all rows, not just those marked "gene"
super_class_match <- filter(FANTOM, complete.cases(geneSuperClass))
FANTOM$geneSuperClass <- super_class_match$geneSuperClass[match(FANTOM$gene_id, super_class_match$gene_id)]
# Restrict to lncRNAs
FANTOM <- filter(FANTOM, geneSuperClass=="all_lncRNA" & type=="transcript")
# Get transcript lengths
FANTOM_lncRNA <- c(FANTOM$end - FANTOM$start)
rm(super_class_match, FANTOM)

# Import GENCODE data
GENCODE <- readGFF("data/annotation/GRCh38_spiked.gtf", version = 2) 

# Get lengths of protein-coding genes identified in this study
this_study <- readGFF("data/annotation/T98G_stringtie_stranded.gtf", filter = list(type="transcript")) %>%
  dplyr::filter(transcript_id%in%counts$transcript_id[counts$transcript_type=="protein_coding"]) %>%
  mutate(length=end-start) %>%
  pull(length)

# Get lengths of coding genes in GENCODE
all_coding <- GENCODE %>%
  filter(type=="transcript" & gene_type=="protein_coding") %>%
  dplyr::select(start, end, transcript_id, gene_name)
# Get transcript lengths
all_coding <- c(all_coding$end - all_coding$start)

# Get GENCODE lncRNA lengths
GENCODE_lncRNA <- GENCODE %>%
  filter(type=="transcript") %>%
  filter(transcript_type %in% c("non_coding", "lincRNA", "macro_lncRNA", "antisense", "sense_intronic", "sense_overlapping", "3prime_overlapping_ncRNA", "bidirectional_promoter_lncRNA", "retained_intron")) %>%
  dplyr::select(start, end, transcript_id, gene_name)
# Get transcript lengths
GENCODE_lncRNA <- c(GENCODE_lncRNA$end - GENCODE_lncRNA$start)

# Get lengths of lncRNAs identified in this study
study_lncRNA <- readGFF("data/annotation/T98G_stringtie_stranded.gtf", filter = list(type="transcript")) %>%
  dplyr::filter(transcript_id%in%counts$transcript_id[counts$transcript_type!="protein_coding"]) %>%
  mutate(length=end-start) %>%
  pull(length)

# Check to see if figure directory has already been created
if(!dir.exists("figures/Fig_S1")){dir.create("figures/Fig_S1")}

plot_S1 <- function(){
  par(mfcol=c(2,2))
  
  hist(log10(all_coding), ylab="Frequency", xlab="log10(Gene length (Kb))", col="#949494", breaks=20, xlim=c(1.5,6.5), main="All GENCODE protein-coding genes")
  legend("topleft", legend=str_glue('mean: {round(mean(all_coding/1e3),1)} Kb'), bty="n")
  
  hist(log10(this_study), ylab="Frequency", xlab="log10(Gene length (Kb))", col="#558547", breaks=20, xlim=c(1.5,6.5), main="Protein-coding genes identified in this study")
  legend("topleft", legend=str_glue('mean: {round(mean(this_study/1e3),1)} Kb'), bty="n")
  
  hist(log10(FANTOM_lncRNA), ylab="Frequency", xlab="log10(Gene length (Kb))", col="#C1272B", breaks=40, xlim=c(1.5,6.5), main="FANTOM- and GENCODE-annotated lncRNAs")
  hist(log10(GENCODE_lncRNA), ylab="Frequency", xlab="log10(Gene length (Kb))", col="#384471", breaks=20, xlim=c(1.5,6.5), add=T)
  legend("topleft",
         legend=c(str_glue('FANTOM\n mean: {round(mean(FANTOM_lncRNA/1e3),1)} Kb'),
                  str_glue('GENCODE\n mean: {round(mean(GENCODE_lncRNA/1e3),1)} Kb')),
         bty="n", pch=22, fill = c("#C1272B", "#384471"))
  
  hist(log10(study_lncRNA), ylab="Frequency", xlab="log10(Gene length (Kb))", col="#558547", breaks=20, xlim=c(1.5,6.5), main="lncRNAs identified in this study")
  legend("topleft", legend=str_glue('mean: {round(mean(study_lncRNA/1e3),1)} Kb'), bty="n")
  
  par(mfcol=c(1,1))
}

# Histograms
pdf("figures/Fig_S1/Fig_S1_gene_length.pdf", width=10, height=8)
plot_S1()
dev.off()

print("Finished creating Fig_S1")

