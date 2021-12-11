### Fig.6 explores the relationship between genomic distance, correlation and
### timing of protein-coding gene/lncRNA expression


# Load R packages -------------------------------------------------------
library(tidyverse)
library(future)
library(furrr)
library(mgcv)
library(GenomicRanges)



# Define get_corr.dist ----------------------------------------------------

# Import data
counts <- readRDS("data/t98g_filtered_coords_clust.rds")

# We define a function to obtain a data frame which contains the distance
# between each requested gene type e.g. coding genes and lncRNAs, and the
# Pearson correlation coefficient between their expression profiles. We allow
# the ability to specify for coding genes whether the expression profile from
# the mRNA or first 10 Kb of pre-mRNA should be used.

get_corr.dist <- function(row_vec, col_vec, row_type, col_type, input_df=counts){
  #get_corr.dist <- function(row_vec, col_vec, coding_expression_type="pre_mrna", input_df=counts){
  
  # Correlation
  
  # The following gives us the correlation between the expression profiles of row_vec and col_vec transcripts
  col_expression <- dplyr::select(input_df, col_type)[col_vec,] %>% do.call(rbind, .)
  row_expression <- dplyr::select(input_df, row_type)[row_vec,] %>% do.call(rbind, .)
  corr.results <- data.frame(cor(t(row_expression),t(col_expression), method="pearson"))
  
  
  # Distance
  
  # The following gives us the distance between every coding and lncRNA transcript
  row_dist <- input_df$coords[row_vec]
  col_dist <- input_df$coords[col_vec]
  dist.results <- sapply(col_dist, function(n) abs(row_dist - n)) %>% data.frame
  
  
  # Combine correlation and distance information
  
  # Name columns by col_vec transcript IDs
  colnames(corr.results) <- input_df$transcript_id[col_vec]
  colnames(dist.results ) <- input_df$transcript_id[col_vec]
  # Name rows by row_vec transcript IDs
  rownames(corr.results) <- input_df$transcript_id[row_vec]
  row.names(dist.results ) <- input_df$transcript_id[row_vec]
  
  # Make data long
  corr.results <- rownames_to_column(corr.results, var = "row_vec") %>% 
    pivot_longer(!row_vec, names_to = "col_vec", values_to = "corr")
  dist.results <-  rownames_to_column(dist.results, var = "row_vec") %>%
    pivot_longer(!row_vec, names_to = "col_vec", values_to = "dist")
  
  # Get both togther
  if(!(all(corr.results$row_vec == dist.results$row_vec) & all(corr.results$col_vec == dist.results$col_vec))){print("WARNING: rows don't match")} # sanity check: rows order should be identical
  corr.dist <- mutate(corr.results, dist = dist.results$dist)
  
  
  # Filter unwanted pairs
  
  # In corr.dist we have the correlation and distance between ALL transcripts.
  #This includes transcripts that aren't on the same chromosome - which is meaningless.
  # We will remove these correlations and distances before returning the data frame:
  
  # Add chromosome info
  corr.dist <- mutate(corr.dist,
                      col_chrom = as.character(input_df$seqid)[match(corr.dist$col_vec, input_df$transcript_id)],
                      row_chrom = as.character(input_df$seqid)[match(corr.dist$row_vec, input_df$transcript_id)])
  # Filter out transcripts not on the same chromosome
  corr.dist <- dplyr::filter(corr.dist, col_chrom == row_chrom) %>%
    dplyr::select(-col_chrom, -row_chrom) # remove unneeded columns
  # Filter out self-comparisons
  corr.dist <- dplyr::filter(corr.dist, col_vec != row_vec)
  
  # Sort by distance
  corr.dist <- mutate(corr.dist, log10_dist = log10(dist+1)) %>% # Add a column with log10(genomic distance) 
    arrange(log10_dist)
  
  return(corr.dist)
}


# Determine block width ---------------------------------------------------

print("Determining bootstrap block width")
# Before we look at applying a block bootstrap, we need to estimate the width of
# the blocks that are appropriate for coding genes and lncRNAs. To figure out
# what this distance should be, we will naively permute the genomic distances
# between coding genes and lncRNAs 1000 times respectively, to get a simulation
# envelope. We will then fit a GAM to the unpermuted data and determine at what
# point does the trend exceed the 99% confidence intervals

# Define a function to permute the genomic distances, fit a GAM to the permuted
# data and return the fitted values
naive_permute <- function(input_df) {
  # Permute the genomic distances between transcript types
  input_df <-
    mutate(input_df,
           log10_dist = sample(
             log10_dist,
             size = length(log10_dist),
             replace = FALSE
           )) %>%
    arrange(log10_dist) # Sort by distance
  # Fit GAM
  gam_fit <-
    gam(corr ~ s(log10_dist), data = input_df, method = "REML")
  # Extract fitted values from GAM
  gam_fit <- predict(gam_fit, type = "terms")
  # Return fitted values as a numeric vector
  return(as.numeric(gam_fit))
}

# Use the get_corr.dist function to get a data frame with the distance and 
# correlation between:

# every coding transcript
c_c <- get_corr.dist(row_vec = counts$transcript_type == "protein_coding" & !sapply(counts$first_10Kb, is.null),
                     col_vec = counts$transcript_type == "protein_coding" & !sapply(counts$first_10Kb, is.null),
                     row_type = "first_10Kb",
                     col_type = "first_10Kb",
                     input_df = counts)

# every non-coding transcript
nc_nc <- get_corr.dist(row_vec = counts$transcript_type != "protein_coding" & !sapply(counts$first_10Kb, is.null),
                       col_vec = counts$transcript_type != "protein_coding" & !sapply(counts$first_10Kb, is.null),
                       row_type = "first_10Kb",
                       col_type = "first_10Kb",
                       input_df = counts)

# coding genes and lncRNAs
c_nc <- get_corr.dist(row_vec = counts$transcript_type != "protein_coding" & !sapply(counts$first_10Kb, is.null),
                      col_vec = counts$transcript_type == "protein_coding" & !sapply(counts$first_10Kb, is.null),
                      row_type = "first_10Kb",
                      col_type = "first_10Kb",
                      input_df = counts)

# Save out
saveRDS(c_c, "data/c_c.rds")
saveRDS(nc_nc, "data/nc_nc.rds")
saveRDS(c_nc, "data/c_nc.rds")

# If not already done, randomly permute genomic distance between coding
# transcripts 1000 times
if (!file.exists("data/c_c_percentiles.rds")) {
  # Load furrr package to speed up running time
  plan(multisession, workers = 20)
  # Set how many times the bootstrap should be run
  n = 1000
  # Set the seed for reproducibility
  set.seed(1234)
  # Run function naive_permute n times
  naive_permute_c_c_results <-
    future_map(
      seq_len(n),
      ~ naive_permute(input_df = c_c),
      .progress = TRUE,
      .options = furrr_options(seed = TRUE)
    )
  
  # Get data in data frame instead of a list
  naive_permute_c_c_results <-
    t(bind_cols(naive_permute_c_c_results))
  # Calculate the 1st, 5th, 25th, 75th, 95th and 99th percentiles
  naive_permute_c_c_results <-
    apply(
      naive_permute_c_c_results ,
      2 ,
      quantile ,
      probs = c(0.01, 0.05, 0.25, 0.75, 0.95, 0.99) ,
      na.rm = TRUE
    )
  # Save out, as takes a while
  saveRDS(naive_permute_c_c_results, "data/c_c_percentiles.rds")
}

# If not already done, randomly permute genomic distance between lncRNA transcripts 1000 times
if(!file.exists("data/nc_nc_percentiles.rds")) {
  # Load furrr package to speed up running time
  plan(multisession, workers = 20)
  # Set how many times the bootstrap should be run
  n = 1000
  # Set the seed for reproducibility
  set.seed(1234)
  # Run function naive_permute n times
  naive_permute_nc_nc_results <-
    future_map(
      seq_len(n),
      ~ naive_permute(input_df = nc_nc),
      .progress = TRUE,
      .options = furrr_options(seed = TRUE)
    )
  # Get data in data frame instead of a list
  naive_permute_nc_nc_results <-
    t(bind_cols(naive_permute_nc_nc_results))
  # Calculate the 1st, 5th, 25th, 75th, 95th and 99th percentiles
  naive_permute_nc_nc_results <-
    apply(
      naive_permute_nc_nc_results ,
      2 ,
      quantile ,
      probs = c(0.01, 0.05, 0.25, 0.75, 0.95, 0.99) ,
      na.rm = TRUE
    )
  # Save out, as takes a while
  saveRDS(naive_permute_nc_nc_results, "data/nc_nc_percentiles.rds")
}

# Examine results for lncRNA-lncRNA pairs
# Fit GAM to the unpermuted data
gam_fit <- gam(corr ~ s(log10_dist), data=nc_nc, method = "REML")
# Examine the signifcance of the fit
anova(gam_fit)
#Approximate significance of smooth terms:
#  edf Ref.df     F p-value
#s(log10_dist) 7.918  8.695 114.1  <2e-16
# Get the 1st, 5th, 25th, 75th, 95th and 99th percentiles
quantiles <- readRDS("data/nc_nc_percentiles.rds")
# At what points does the trend exceed the 99% confidence intervals?
nc_nc_block <- 10^(nc_nc$log10_dist[which(as.numeric(predict(gam_fit, type = "terms")) < quantiles["99%",])[1]])
nc_nc_block
#[1] 6397642

# Examine results of the above for coding-coding pairs
# Fit GAM to the unpermuted data
gam_fit <- gam(corr ~ s(log10_dist), data=c_c, method = "REML")
# Examine a summary of the fit
anova(gam_fit)
#Approximate significance of smooth terms:
#  edf Ref.df     F p-value
#s(log10_dist) 8.690  8.972 61.17  <2e-16
# Get the 1st, 5th, 25th, 75th, 95th and 99th percentiles
quantiles <- readRDS("data/c_c_percentiles.rds")
# At what point does the trend exceed the 99% confidence intervals?
c_c_block <- 10^(c_c$log10_dist[which(as.numeric(predict(gam_fit, type = "terms")) < quantiles["99%",])[1]])
c_c_block
#[1] 4283435

# Examine results for lncRNA-coding pairs
# Fit GAM to the unpermuted data
gam_fit <- gam(corr ~ s(log10_dist), data=c_nc, method = "REML")
# Examine a summary of the fit
anova(gam_fit)
#Approximate significance of smooth terms:
#  edf Ref.df     F p-value
#s(log10_dist) 8.197  8.822 80.99  <2e-16

# Clean up
rm(quantiles, gam_fit, c_c, nc_nc, c_nc, naive_permute)




# Apply the block bootstrap -----------------------------------------------

print("Applying the block bootstrap")
# Now that the block sizes have been determined, create a function that will
# carry out the block bootstrap. The function takes as input the counts data
# frame and returns another counts data frame generated by sampling genomic
# blocks with replacement.

## Start by tiling chromosomes with blocks of appropriate size and allocating
## coding genes and lncRNAs (seperately) into a numbered bin

# Define vector of chromosome lengths for GRCh38
chromosome_length <- c(chr1=248956422, chr10=133797422, chr11=135086622,
                       chr12=133275309, chr13=114364328, chr14=107043718,
                       chr15=101991189, chr16=90338345, chr17=83257441,
                       chr18=80373285, chr19=58617616, chr2=242193529,
                       chr20=64444167, chr21=46709983, chr22=50818468,
                       chr3=198295559, chr4=190214555, chr5=181538259,
                       chr6=170805979, chr7=159345973, chr8=145138636,
                       chr9=138394717, chrX=156040895, chrY=57227415)

# Get tiles across the genome, separately for coding genes and lncRNAs
coding_blocks <- tileGenome(chromosome_length,
                            tilewidth = c_c_block, cut.last.tile.in.chrom = TRUE)
non_coding_blocks <- tileGenome(chromosome_length,
                                tilewidth = nc_nc_block, cut.last.tile.in.chrom = TRUE)

# Mark which bin each gene falls into and the distance of each gene from the
# start of the genomic bin (coords_within_block)
coding_genes <- filter(counts, transcript_type == "protein_coding") %>%
  mutate(block = findOverlaps(
    GRanges(
      seqnames = seqid,
      ranges = IRanges(start = coords, end = coords),
      strand = "*"
    ),
    coding_blocks,
    select = "first"
  )) %>%
  mutate(coords_within_block = sapply(1:n(), function(x)
    coords[x] - start(coding_blocks[block[x]])))

noncoding_genes <- filter(counts, transcript_type != "protein_coding") %>%
  mutate(block = findOverlaps(
    GRanges(
      seqnames = seqid,
      ranges = IRanges(start = coords, end = coords),
      strand = "*"
    ),
    non_coding_blocks,
    select = "first"
  )) %>%
  mutate(coords_within_block = sapply(1:n(), function(x)
    coords[x] - start(non_coding_blocks[block[x]])))

# Bind together
counts <- rbind(coding_genes, noncoding_genes) 

# Clean up
rm(coding_genes, noncoding_genes)

## Define `bbstrap` function to apply the block bootstrap. Generate new
#chromosomes by sampling (with replacement) an appropriate number of blocks
#(enough to cover each chromosome). Blocks from different chromosomes can be
#selected to create new chromosomes.

bbstrap <- function(input_df = counts,
                    cb = coding_blocks, ncb = non_coding_blocks,
                    ccb = c_c_block, ncncb = nc_nc_block){
  # Sample tiles - separately for coding and non-coding - with replacement to
  # give new shuffled chromosomes
  shuff_c <- lapply(table(seqnames(cb)),
                    function(n) sample(cb,size = n, replace = TRUE))
  shuff_nc <- lapply(table(seqnames(ncb)),
                     function(n) sample(ncb,size = n, replace = TRUE))
  
  # Add new chromosome name and tile start positions
  for(i in seq_along(shuff_c)) { 
    shuff_c[[i]]$new_seqid <- names(shuff_c[i])
    shuff_c[[i]]$new_start <- (seq_along(shuff_c[[i]])-1)*ccb 
  }
  for(i in seq_along(shuff_nc)) {
    shuff_nc[[i]]$new_seqid <- names(shuff_nc[i])
    shuff_nc[[i]]$new_start <- (seq_along(shuff_nc[[i]])-1)*ncncb
  }
  
  # Collapse tiles to GRanges objects
  shuff_c <- unlist(as(shuff_c, "GRangesList"))
  shuff_nc <- unlist(as(shuff_nc, "GRangesList"))
  
  # split counts into coding and non-coding
  coding <- dplyr::filter(counts, transcript_type=="protein_coding")
  noncoding <- dplyr::filter(counts, transcript_type!="protein_coding")
  
  # Create ranges
  coding_ranges <- dplyr::mutate(coding, seqid=seqid, start=coords, end=coords) %>%
    makeGRangesFromDataFrame()
  noncoding_ranges <- dplyr::mutate(noncoding, seqid=seqid, start=coords, end=coords) %>%
    makeGRangesFromDataFrame()
  
  # Match genes with tiles
  # coding
  hits <- findOverlaps(query = coding_ranges, subject = shuff_c) %>% data.frame()
  coding <- coding[hits$queryHits,] %>%
    mutate(new_seqid = shuff_c$new_seqid[hits$subjectHits],
           new_start = shuff_c$new_start[hits$subjectHits])
  #non-coding
  hits <- findOverlaps(query = noncoding_ranges, subject = shuff_nc) %>% data.frame()
  noncoding <- noncoding[hits$queryHits,] %>%
    mutate(new_seqid = shuff_nc$new_seqid[hits$subjectHits],
           new_start = shuff_nc$new_start[hits$subjectHits])
  
  # Update chromosome positions
  coding <- mutate(coding, coords = new_start + coords_within_block, seqid = new_seqid)
  noncoding <- mutate(noncoding, coords = new_start + coords_within_block, seqid = new_seqid)
  
  # Remove any duplicated transcript IDs
  bbstrap_counts <- rbind(coding, noncoding)
  bbstrap_counts <- mutate(bbstrap_counts, transcript_id=paste0(transcript_id,"###",row_number()))
  
  return(bbstrap_counts)
}

# The block bootstrap is carried out by the function `bbstrap`, which we just
# defined above. We will write another function that calls `bbstrap` and returns
# a GAM fit to the bootstrapped data. From these values we can take the 5th and
# 95th percentiles to get our confidence intervals. We'll call this new function
# `bbstrap_iterate`. Get data frame with correlation and log10(distance) for the
# un-bootstrapped data
c_nc <- readRDS("data/c_nc.rds")
# Define function to apply bootstrap multiple times
bbstrap_iterate <- function(df_in=counts){
  boot_data <- bbstrap()
  # Get corr.dist data frame using the block bootsrap data
  boot_data <- get_corr.dist(row_vec = boot_data$transcript_type!="protein_coding" & !sapply(boot_data$first_10Kb, is.null), col_vec = boot_data$transcript_type=="protein_coding" & !sapply(boot_data$first_10Kb, is.null), row_type = "first_10Kb", col_type = "first_10Kb", input_df = boot_data)
  # Fit GAM to bootstrapped data
  gam_fit <- gam(corr ~ s(log10(dist+1)), data=boot_data, method = "REML")
  # Return fitted GAM values for vector of log10_dist from unbooted data
  return(as.numeric(predict(gam_fit, c_nc)))
}

# If not already done, run the bootsrap n times to obtain confidence intervals -
# using ALL lncRNAs (Fig. 6)
if(!file.exists("data/block_bootstrap_quantiles.rds")){
  # Load furrr package to speed up running the block bootstrap
  plan(multisession, workers = 20)
  # Set how many times the bootstrap should be run
  n=1000
  # Set the seed for reproducibility
  set.seed(1234)
  # Run block bootsrap n times
  block_bootstrap_results <- future_map(seq_len(n), ~ bbstrap_iterate(), .progress = TRUE, .options = furrr_options(seed = TRUE))
  # Get data in data frame instead of a list
  block_bootstrap_results <- t(bind_cols(block_bootstrap_results))
  # Calculate percentiles
  quantiles <- apply(block_bootstrap_results , 2, quantile, probs = c(0.01,0.05,0.25,0.75,0.95,0.99) , na.rm = TRUE )
  # Save out, as takes a while
  saveRDS(quantiles, "data/block_bootstrap_quantiles.rds")
  # Clean up
  rm(n, quantiles, block_bootstrap_results)
}



# Prepare data for Supplementary Fig. 4 -----------------------------------

# Run the bootstrap using only lncRNAs that overlap annotated lncRNAs
# (Supplementary Fig. 4 pt I)
subset_counts <- filter(counts, transcript_type=="annotated_lncRNA" | transcript_type=="protein_coding")
c_nc <- get_corr.dist(row_vec = subset_counts$transcript_type != "protein_coding" & !sapply(subset_counts$first_10Kb, is.null),
                      col_vec = subset_counts$transcript_type == "protein_coding" & !sapply(subset_counts$first_10Kb, is.null),
                      row_type = "first_10Kb",
                      col_type = "first_10Kb",
                      input_df = subset_counts)
saveRDS(c_nc, "data/c_nc_anno.rds")
if(!file.exists("data/S4_overlapping_block_bootstrap_quantiles.rds")){
  # Load furrr package to speed up running the block bootstrap
  plan(multisession, workers = 20)
  # Set how many times the bootstrap should be run
  n=1000
  # Set the seed for reproducibility
  set.seed(1234)
  # Run block bootsrap n times
  block_bootstrap_results <- future_map(seq_len(n), ~ bbstrap_iterate(df_in = subset_counts), .progress = TRUE, .options = furrr_options(seed = TRUE))
  # Get data in data frame instead of a list
  block_bootstrap_results <- t(bind_cols(block_bootstrap_results))
  # Calculate percentiles
  quantiles <- apply(block_bootstrap_results , 2, quantile, probs = c(0.01,0.05,0.25,0.75,0.95,0.99) , na.rm = TRUE )
  # Save out, as takes a while
  saveRDS(quantiles, "data/S4_overlapping_block_bootstrap_quantiles.rds")
  # Clean up
  rm(n, quantiles, block_bootstrap_results)
}

# Run the bootsrap using only lncRNAs that do NOT overlap annotated lncRNAs
# (Supplementary Fig. 4 pt II)
subset_counts <- filter(counts, transcript_type=="novel_lncRNA" | transcript_type=="protein_coding")
c_nc <- get_corr.dist(row_vec = subset_counts$transcript_type != "protein_coding" & !sapply(subset_counts$first_10Kb, is.null),
                      col_vec = subset_counts$transcript_type == "protein_coding" & !sapply(subset_counts$first_10Kb, is.null),
                      row_type = "first_10Kb",
                      col_type = "first_10Kb",
                      input_df = subset_counts)
if(!file.exists("data/S4_not_overlapping_block_bootstrap_quantiles.rds")){
  # Load furrr package to speed up running the block bootstrap
  plan(multisession, workers = 20)
  # Set how many times the bootstrap should be run
  n=1000
  # Set the seed for reproducibility
  set.seed(1234)
  # Run block bootsrap n times
  block_bootstrap_results <- future_map(seq_len(n), ~ bbstrap_iterate(df_in = subset_counts), .progress = TRUE, .options = furrr_options(seed = TRUE))
  # Get data in data frame instead of a list
  block_bootstrap_results <- t(bind_cols(block_bootstrap_results))
  # Calculate percentiles
  quantiles <- apply(block_bootstrap_results , 2, quantile, probs = c(0.01,0.05,0.25,0.75,0.95,0.99) , na.rm = TRUE )
  # Save out, as takes a while
  saveRDS(quantiles, "data/S4_not_overlapping_block_bootstrap_quantiles.rds")
  # Clean up
  rm(n, quantiles, block_bootstrap_results, subset_counts)
}



# Prepare data for Supplementary Fig. 5 -----------------------------------

# Run the bootstrap using only lncRNAs with TSS near a promoter element
# (Supplementary Fig. 5 pt I)
h.lnc <- readRDS("data/h.lnc.rds")
counts <- counts[counts$transcript_id%in%h.lnc$transcript_id,]
h.lnc <- h.lnc[h.lnc$transcript_id%in%counts$transcript_id,]
h.lnc <- arrange(h.lnc,transcript_id)
counts <- arrange(counts,transcript_id)
table(h.lnc$transcript_id==counts$transcript_id)
counts$reg <- h.lnc$reg
subset_counts <- filter(counts, (transcript_type!="protein_coding" & reg=="promoter") | transcript_type=="protein_coding")
c_nc <- get_corr.dist(row_vec = subset_counts$transcript_type != "protein_coding" & !sapply(subset_counts$first_10Kb, is.null),
                      col_vec = subset_counts$transcript_type == "protein_coding" & !sapply(subset_counts$first_10Kb, is.null),
                      row_type = "first_10Kb",
                      col_type = "first_10Kb",
                      input_df = subset_counts)
saveRDS(c_nc, "data/c_nc_promoter.rds")
if(!file.exists("data/S5_promoter_block_bootstrap_quantiles.rds")){
  # Load furrr package to speed up running the block bootstrap
  plan(multisession, workers = 20)
  # Set how many times the bootstrap should be run
  n=1000
  # Set the seed for reproducibility
  set.seed(1234)
  # Run block bootsrap n times
  block_bootstrap_results <- future_map(seq_len(n), ~ bbstrap_iterate(df_in = subset_counts), .progress = TRUE, .options = furrr_options(seed = TRUE))
  # Get data in data frame instead of a list
  block_bootstrap_results <- t(bind_cols(block_bootstrap_results))
  # Calculate percentiles
  quantiles <- apply(block_bootstrap_results , 2, quantile, probs = c(0.01,0.05,0.25,0.75,0.95,0.99) , na.rm = TRUE )
  # Save out, as takes a while
  saveRDS(quantiles, "data/S5_promoter_block_bootstrap_quantiles.rds")
  # Clean up
  rm(n, quantiles, block_bootstrap_results)
}

# Run the bootstrap using only lncRNAs with TSS near an enhancer element
# (Supplementary Fig. 5 pt II)
subset_counts <- filter(counts, (transcript_type!="protein_coding" & reg=="enhancer") | transcript_type=="protein_coding")
c_nc <- get_corr.dist(row_vec = subset_counts$transcript_type != "protein_coding" & !sapply(subset_counts$first_10Kb, is.null),
                      col_vec = subset_counts$transcript_type == "protein_coding" & !sapply(subset_counts$first_10Kb, is.null),
                      row_type = "first_10Kb",
                      col_type = "first_10Kb",
                      input_df = subset_counts)
saveRDS(c_nc, "data/c_nc_enhancer.rds")
if(!file.exists("data/S5_enhancer_block_bootstrap_quantiles.rds")){
  # Load furrr package to speed up running the block bootstrap
  plan(multisession, workers = 20)
  # Set how many times the bootstrap should be run
  n=1000
  # Set the seed for reproducibility
  set.seed(1234)
  # Run block bootsrap n times
  block_bootstrap_results <- future_map(seq_len(n), ~ bbstrap_iterate(df_in = subset_counts), .progress = TRUE, .options = furrr_options(seed = TRUE))
  # Get data in data frame instead of a list
  block_bootstrap_results <- t(bind_cols(block_bootstrap_results))
  # Calculate percentiles
  quantiles <- apply(block_bootstrap_results , 2, quantile, probs = c(0.01,0.05,0.25,0.75,0.95,0.99) , na.rm = TRUE )
  # Save out, as takes a while
  saveRDS(quantiles, "data/S5_enhancer_block_bootstrap_quantiles.rds")
  # Clean up
  rm(n, quantiles, block_bootstrap_results, subset_counts)
}

# Clean up
rm(bbstrap_iterate, c_nc, subset_counts)



# Examine time lags of expression profiles --------------------------------

#Below we compare coding gene expression profiles with lags of the lncRNA
#expression profiles. We use the
#[ccf](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/acf.html)
#function to reveal how the correlation between the two time series changes as
#their separation in time changes.

# Check to see if the file already exists and if not, create it
if(!file.exists("data/c_nc_lagged_corr.rds")){
  
  # Get data frame with ALL coding/lncRNA pairs, including pairs that are on
  # different chromosomes, as the block bootstrap mixes these.
  c_nc <- readRDS("data/c_nc.rds") %>% dplyr::select(c("row_vec", "col_vec"))
  c_nc <- expand.grid(unique(c_nc$row_vec), unique(c_nc$col_vec))
  colnames(c_nc) <- c("row_vec", "col_vec")
  
  # Modify the data frame to include the expression data
  c_nc <- mutate(c_nc,
                 ncRNA = counts$first_10Kb[match(c_nc$row_vec, counts$transcript_id)],
                 pre_mRNA = counts$first_10Kb[match(c_nc$col_vec, counts$transcript_id)],
                 mRNA = counts$exon[match(c_nc$col_vec, counts$transcript_id)])

  # Define a function to calculate the lagged correlation between each pair in
  # the created data frame
  get_lag_cor <- function(x, coding_expression_type = "pre_mRNA"){
    if(x%%1000==0){print(str_glue('{round(100*(x/nrow(c_nc)),3)}% complete'))} # Update progress periodically
    if(coding_expression_type=="pre_mRNA"){
      return(as.numeric(ccf(x = unlist(c_nc$ncRNA[[x]]), y = unlist(c_nc$pre_mRNA[[x]]), lag.max = 20, type="correlation", na.action = na.omit, plot = FALSE)$acf))
    } else {
      return(as.numeric(ccf(x = unlist(c_nc$ncRNA[[x]]), y = unlist(c_nc$mRNA[[x]]), lag.max = 20, type="correlation", na.action = na.omit, plot = FALSE)$acf))
    }
  }

  # Calculate and return the required lagged correlation values using the
  # pre-mRNA expression data, convert to a data frame and give meaningful column
  # names
  lagged_cor_pre_mRNA <- lapply(1:nrow(c_nc), get_lag_cor)
  lagged_cor_pre_mRNA <- do.call(rbind, lagged_cor_pre_mRNA)
  colnames(lagged_cor_pre_mRNA) <- str_glue('premRNA_lag_{seq(-200,200,10)}')
  
  # Calculate and return the required lagged correlation values using the mRNA
  # expression data, Convert to a data frame and give meaningful column names
  lagged_cor_mRNA <- lapply(1:nrow(c_nc), get_lag_cor, coding_expression_type="mRNA")
  lagged_cor_mRNA <- do.call(rbind, lagged_cor_mRNA)
  colnames(lagged_cor_mRNA) <- str_glue('lag_mRNA_{seq(-200,200,10)}')
  
  # Add to c_nc data frame and save out
  c_nc <- cbind(c_nc, lagged_cor_pre_mRNA)
  c_nc <- cbind(c_nc, lagged_cor_mRNA)
  saveRDS(dplyr::select(c_nc, -ncRNA, -pre_mRNA, -mRNA), file = "data/c_nc_lagged_corr.rds")
  
  # Clean up
  rm(c_nc, get_lag_cor, lagged_cor_mRNA, lagged_cor_pre_mRNA)
}

# Check to see if block bootstrap data has already been obtained, if not get it
if(!file.exists("data/lagged_cor_bbstrap.rds")){
  #import lagged correlation data
  lagged_cor <- readRDS("data/c_nc_lagged_corr.rds")
  lagged_cor <- lagged_cor[rowSums(is.na(lagged_cor))==0,]
  
  # Define a function to apply the block bootstrap and return the mean ACF
  # values
  lag_bbstrap <-function(iteration_number){
    # Create new counts data frame using the block bootstrap
    counts_bbstrap <- bbstrap()
    # Get corr dist from bootstrap data frame
    counts_bbstrap <- get_corr.dist(row_vec = counts_bbstrap$transcript_type != "protein_coding" & !sapply(counts_bbstrap$first_10Kb, is.null),
                          col_vec = counts_bbstrap$transcript_type == "protein_coding" & !sapply(counts_bbstrap$first_10Kb, is.null),
                          row_type = "first_10Kb",
                          col_type = "first_10Kb",
                          input_df = counts_bbstrap)
    ## No need to re-calculate lags - take from lagged_cor data frame
    
    # Create ids for coding-ncRNA pairs from counts_bbstrap data frame
    counts_bbstrap$row_vec <- str_split_fixed(counts_bbstrap$row_vec, "###",2)[,1]
    counts_bbstrap$col_vec <- str_split_fixed(counts_bbstrap$col_vec, "###",2)[,1]
    bb_id <- as.character(str_glue('{counts_bbstrap$row_vec}#{counts_bbstrap$col_vec}'))
    
    # Create ids for coding-ncRNA pairs from lagged_cor data frame
    lc_id <- as.character(str_glue('{lagged_cor$row_vec}#{lagged_cor$col_vec}'))
    
    # Join lagged corr data to counts_bbstrap data frame
    counts_bbstrap <- cbind(counts_bbstrap, lagged_cor[match(bb_id, lc_id), str_detect(colnames(lagged_cor), "lag_")])
    counts_bbstrap <- counts_bbstrap[rowSums(is.na(counts_bbstrap))==0,]
    
    # Update progress
    print(paste0("Completed iteration ", iteration_number))
    # Calculate mean ACF value across all lags for the following distances (between coding genes and ncRNAs)
    return(list(pre_mRNA_10Kb = colMeans(counts_bbstrap[counts_bbstrap$dist<1E4, str_detect(colnames(counts_bbstrap), "premRNA_lag_")]),
                pre_mRNA_100Kb = colMeans(counts_bbstrap[counts_bbstrap$dist<1E5, str_detect(colnames(counts_bbstrap), "premRNA_lag_")]),
                pre_mRNA_250Kb = colMeans(counts_bbstrap[counts_bbstrap$dist<2.5E5, str_detect(colnames(counts_bbstrap), "premRNA_lag_")]),
                pre_mRNA_1Mb = colMeans(counts_bbstrap[counts_bbstrap$dist<1E6, str_detect(colnames(counts_bbstrap), "premRNA_lag_")]),
                pre_mRNA_g1Mb = colMeans(counts_bbstrap[counts_bbstrap$dist>1E6, str_detect(colnames(counts_bbstrap), "premRNA_lag_")]),
                mRNA_10Kb = colMeans(counts_bbstrap[counts_bbstrap$dist<1E4, str_detect(colnames(counts_bbstrap), "lag_mRNA_")]),
                mRNA_100Kb = colMeans(counts_bbstrap[counts_bbstrap$dist<1E5, str_detect(colnames(counts_bbstrap), "lag_mRNA_")]),
                mRNA_250Kb = colMeans(counts_bbstrap[counts_bbstrap$dist<2.5E5, str_detect(colnames(counts_bbstrap), "lag_mRNA_")]),
                mRNA_1Mb = colMeans(counts_bbstrap[counts_bbstrap$dist<1E6, str_detect(colnames(counts_bbstrap), "lag_mRNA_")]),
                mRNA_g1Mb = colMeans(counts_bbstrap[counts_bbstrap$dist>1E6, str_detect(colnames(counts_bbstrap), "lag_mRNA_")])))
  }
  # Set how many times lag_bbstrap should be run
  n=1000
  # Set the seed for reproducibility
  set.seed(1234)
  # Run lag_bbstrap n times
  lagged_cor_bbstrap <- lapply(1:n, lag_bbstrap)
  # Calculate the 5th, 25th, 75th and 95th percentiles
  lagged_cor_bbstrap <- list(pre_mRNA_10Kb = apply(do.call(rbind, lapply(lagged_cor_bbstrap, function(i) i[[1]])) , 2 , quantile , probs = c(0.01,0.05,0.25,0.75,0.95,0.99) , na.rm = TRUE),
                             pre_mRNA_100Kb = apply(do.call(rbind, lapply(lagged_cor_bbstrap, function(i) i[[2]])) , 2 , quantile , probs = c(0.01,0.05,0.25,0.75,0.95,0.99) , na.rm = TRUE),
                             pre_mRNA_250Kb = apply(do.call(rbind, lapply(lagged_cor_bbstrap, function(i) i[[3]])) , 2 , quantile , probs = c(0.01,0.05,0.25,0.75,0.95,0.99) , na.rm = TRUE),
                             pre_mRNA_1Mb = apply(do.call(rbind, lapply(lagged_cor_bbstrap, function(i) i[[4]])) , 2 , quantile , probs = c(0.01,0.05,0.25,0.75,0.95,0.99) , na.rm = TRUE),
                             pre_mRNA_g1Mb = apply(do.call(rbind, lapply(lagged_cor_bbstrap, function(i) i[[5]])) , 2 , quantile , probs = c(0.01,0.05,0.25,0.75,0.95,0.99) , na.rm = TRUE),
                             mRNA_10Kb = apply(do.call(rbind, lapply(lagged_cor_bbstrap, function(i) i[[6]])) , 2 , quantile , probs = c(0.01,0.05,0.25,0.75,0.95,0.99) , na.rm = TRUE),
                             mRNA_100Kb = apply(do.call(rbind, lapply(lagged_cor_bbstrap, function(i) i[[7]])) , 2 , quantile , probs = c(0.01,0.05,0.25,0.75,0.95,0.99) , na.rm = TRUE),
                             mRNA_250Kb = apply(do.call(rbind, lapply(lagged_cor_bbstrap, function(i) i[[8]])) , 2 , quantile , probs = c(0.01,0.05,0.25,0.75,0.95,0.99) , na.rm = TRUE),
                             mRNA_1Mb = apply(do.call(rbind, lapply(lagged_cor_bbstrap, function(i) i[[9]])) , 2 , quantile , probs = c(0.01,0.05,0.25,0.75,0.95,0.99) , na.rm = TRUE),
                             mRNA_g1Mb = apply(do.call(rbind, lapply(lagged_cor_bbstrap, function(i) i[[10]])) , 2 , quantile , probs = c(0.01,0.05,0.25,0.75,0.95,0.99) , na.rm = TRUE))
  # Save out
  saveRDS(lagged_cor_bbstrap, "data/lagged_cor_bbstrap.rds")

  # Clean up
  rm(lagged_cor, lag_bbstrap, lagged_cor_bbstrap, n)
}

print("Finished creating data for Fig. 6, Fig. S3 and Fig. S4")
