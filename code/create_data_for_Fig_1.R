#### Figure 1 - Protein coding and lncRNA expression
#### This figure is made up of:
####  a. Schematic of when the time points were taken - made elsewhere
####  b. lncRNA expression heatmap
####  c. mRNA expression heatmap
####  d. Centroid correlation matric plot

# In this script we will carry out some clustering of:
  #   coding mRNA expression profiles 
  #   lncRNA expression profiles


# 1.  Load packages and data ----------------------------------------------

# Load required R packages
library(tidyverse)
library(rtracklayer)
library(DESeq2)
library(furrr)
library(patchwork)

# Get data
t98g <- readRDS("data/t98g_filtered_coords.rds")



# 2. Get expression data --------------------------------------------------

# Get data frames with just expression data for each transcript type
mrna <- t98g %>%
  filter(transcript_type=="protein_coding") %>%
  select(transcript_id, exon) %>%
  unnest(cols=c(transcript_id, exon))

premrna <- t98g %>%
  filter(transcript_type=="protein_coding") %>%
  select(transcript_id, first_10Kb) %>%
  unnest(cols=c(transcript_id, first_10Kb))

lncrna <- t98g %>%
  filter(transcript_type!="protein_coding") %>%
  select(transcript_id, transcript) %>%
  unnest(cols=c(transcript_id, transcript))



# 3. Determine suitable k -------------------------------------------------

# Define a function to determine a suitable number of clusters for
# visualization by determining at what k the total intra-cluster variation
# is minimized.

get_k <- function(input_expression_data, input_seed = 1234){

  # Remove rows with NA
  input_expression_data <- drop_na(input_expression_data)
  
  # Calculate correlation matrix
  cor_matrix <- cor(t(input_expression_data[,2:42]), method="pearson")
  
  # Define a function to compute the total within-cluster sum of square (wss) -
  # compactness of the clustering, which we seek to minimise -
  wss <- function(k) { kmeans(cor_matrix, k, nstart = 10)$tot.withinss }
  
  # Compute wss for k = 2 to k = 15
  k.values <- 2:15
  
  # extract wss for 2-15 clusters if not already done
  set.seed(input_seed)
  plan(multisession, workers = 20)
  wss_values <- future_map_dbl(k.values, wss, .progress = TRUE, .options = furrr_options(seed = TRUE))
  plan("sequential")
  return(wss_values)
}

# Calculate and save out
if(!file.exists("data/pre_mRNA_wss.rds")){
  saveRDS(get_k(input_expression_data = premrna), "data/pre_mRNA_wss2.rds")
  saveRDS(get_k(input_expression_data = mrna), "data/mRNA_wss2.rds")
  saveRDS(get_k(input_expression_data = lncrna), "data/lncRNA_wss2.rds")
}

# Plot
wss <- tibble(k.values = 2:15,
              mrna_wss = readRDS("data/mRNA_wss2.rds"),
              premrna_wss = readRDS("data/pre_mRNA_wss2.rds"),
              lncrna_wss = readRDS("data/lncRNA_wss2.rds"))

base_p <- ggplot(wss, aes(x=k.values)) +
  geom_point() +
  xlab("Number of clusters K") 

p1 <- base_p + 
  aes(y=mrna_wss) +
  ylab("Total within-clusters sum of squares") +
  ggtitle("mRNA")

p2 <- base_p + 
  aes(y=premrna_wss) +
  ylab("Total within-clusters sum of squares") +
  ggtitle("pre-mRNA")

p3 <- base_p + 
  aes(y=lncrna_wss) +
  ylab("Total within-clusters sum of squares") +
  ggtitle("lncRNA")

ggsave(filename = "figures_track/wss.pdf", device = "pdf",
       plot = p1 + p2 + p3, width = 40, height = 10, units = "cm")

#Following visual inspection, we select 6 clusters (note: for visualisation
#purposes only)
num_clusters <- 6

# Obtain clusters if not done already
if(!file.exists("data/t98g_filtered_coords_clust.rds")){
  
  # pre-mRNA
  premrna_clusters <- drop_na(premrna)
  premrna_clusters$premrna_cluster <- kmeans(cor(t(premrna_clusters[,2:42])),
                             centers = num_clusters,
                             iter.max = 100, nstart = 100)$cluster
  #
  
  # mRNA
  mrna_clusters <- drop_na(mrna)
  mrna_clusters$mrna_cluster <- kmeans(cor(t(mrna_clusters[,2:42])),
                                             centers = num_clusters,
                                             iter.max = 100, nstart = 100)$cluster
  #mrna_clusters <- select(mrna_clusters, transcript_id, mrna_cluster)
  
  # lncRNA
 lncrna_clusters <- drop_na(lncrna)
 lncrna_clusters$lncrna_cluster <- kmeans(cor(t(lncrna_clusters[,2:42])),
                                       centers = num_clusters,
                                       iter.max = 100, nstart = 100)$cluster
 #lncrna_clusters <- select(lncrna_clusters, transcript_id, lncrna_cluster)
  
 
 
 # Rename clusters according to the time point at which the centroid maximum
 # occurs i.e cluster 1 corresponds to an early centroid peak, and the last
 # cluster is the latest centroid peak
 
 # pre-mRNA
 centroid_max <- sapply(1:num_clusters, FUN=function(x) mean(apply(premrna_clusters[premrna_clusters$premrna_cluster==x,],1,which.max))) 
 premrna_clusters$premrna_cluster <- match(premrna_clusters$premrna_cluster, c(1:num_clusters)[order(centroid_max)])
 premrna_clusters <- select(premrna_clusters, transcript_id, premrna_cluster)
 
 # mRNA
 centroid_max <- sapply(1:num_clusters, FUN=function(x) mean(apply(mrna_clusters[mrna_clusters$mrna_cluster==x,],1,which.max))) 
 mrna_clusters$mrna_cluster <- match(mrna_clusters$mrna_cluster, c(1:num_clusters)[order(centroid_max)])
 mrna_clusters <- select(mrna_clusters, transcript_id, mrna_cluster)
 
 # lncRNA
 centroid_max <- sapply(1:num_clusters, FUN=function(x) mean(apply(lncrna_clusters[lncrna_clusters$lncrna_cluster==x,],1,which.max))) 
 lncrna_clusters$lncrna_cluster <- match(lncrna_clusters$lncrna_cluster, c(1:num_clusters)[order(centroid_max)])
 lncrna_clusters <- select(lncrna_clusters, transcript_id, lncrna_cluster)
 
 
 # Add cluster info
 t98g <- list(t98g, mrna_clusters, premrna_clusters, lncrna_clusters) %>%
   purrr::reduce(full_join, by = "transcript_id")
 
 # save  
  saveRDS(t98g, "data/t98g_filtered_coords_clust.rds")
} else {
  t98g <- readRDS("data/t98g_filtered_coords_clust.rds")
}

print("Finished creating clustering data for Fig_1")
