### Fig S2 - Estimate mean Pol II transcription elongation rate

# Load R Packages
library(tidyverse)
library(rtracklayer)
library(GenomicFeatures)
library(DESeq2)
library(Rsamtools)
library(GenomicAlignments)
library(minpack.lm)
library(RColorBrewer)
library(glue)

# Check to see if data has already been obtained
if(!dir.exists("figures/Fig_S2")){dir.create("figures/Fig_S2")}
if(!file.exists("data/elongation_rate_data.rds")){
  
  # Import counts data
  counts <- readRDS("data/custom_counts_norm_overlap_dist_clust.rds")
  
  # Create data frame to hold the data we will need - ~100 suitable genes
  elongation_rate_data <- data.frame(tx_id=as.character(counts$id[counts$transcript_type=="coding"][c(12,27,39,62,92,171,255,259,290,291,317,325,327,344,348,407,457,484,487,503,558,674,696,697,717,731,770,872,894,906,918,941,1006,1027,1053,1071,1120,1165,1183,1188,1260,1261,1265,1274,1275,1298,1299,1328,1332,1346,1369,1371,1384,1385,1400,1433,1470,1472,1477,1526,1535,1547,1557,1599,1612,1616,1685,1751,1775,1787,1838,1841,1875,1885,1912,1969,2023,2047,2052,2055,2073,2077,2093,2236,2239,2247,2293,2315,2448,2504,2573,2616,2673,2696,2726,2740,2775,2823,2845,2921,2930,2937)]),
                                     gene_name=NA,
                                     distance=NA,
                                     first_10Kb_expression=NA,
                                     last_10Kb_expression=NA,
                                     first_10Kb_fit=NA,
                                     last_10Kb_fit=NA,
                                     time=NA)
  
  # Get GENCODE gtf so that we can add gene names
  GENCODE <- readGFF("data/annotation/GRCh38_spiked.gtf", version = 2) %>%
    filter(type=="transcript") %>%
    dplyr::select(transcript_id, gene_name)
  elongation_rate_data$gene_name <- GENCODE$gene_name[match(elongation_rate_data$tx_id, GENCODE$transcript_id)]
  rm(GENCODE)
  
  # Load GTF data as a TxDb class object so we can get intron intervals
  TxDb <- makeTxDbFromGFF("data/annotation/GRCh38_spiked.gtf", format="gtf")
  introns <- intronsByTranscript(TxDb, use.names=TRUE)
  # Restrict to txs of interest
  introns <- introns[names(introns)%in%elongation_rate_data$tx_id]
  # same order as elongation rate data frame
  introns <- introns[match(as.character(elongation_rate_data$tx_id), names(introns))]
  rm(TxDb)
  
  # Get size factors for normalization
  # Import gene counts
  gene_counts <- readRDS("data/gene_counts.rds")
  sfs <- estimateSizeFactorsForMatrix(gene_counts$counts)
  rm(gene_counts)
  
  # Create BamFileList
  fls <- BamFileList(file=paste0("data/alignments/T98G_", stringr::str_pad(seq(0,400,10), width=3, side="left", pad="0"), "_minAligned.sortedByCoord.out.bam"),
                     index=paste0("data/alignments/T98G_", stringr::str_pad(seq(0,400,10), width=3, side="left", pad="0"), "_minAligned.sortedByCoord.out.bam.bai"))
  
  # Define function to calculate overlap with an interval, using RNA-seq reads from all time points
  allReadsCounts <- function(i, interval){
    suppressWarnings(countOverlaps(interval,
                                   readGAlignmentPairs(file=fls[[i]], use.names=TRUE, param=ScanBamParam(which = interval, simpleCigar=TRUE), with.which_label=FALSE, strandMode=2),minoverlap=5))
  }
  
  # Define a function that accepts a transcript ID as input and finds the expression profile of the first/last 10Kb of intronic regions
  get_plot_data <- function(tx_id){
    # Get GRanges that define the first and last 10Kb
    # Check if tx is on the forward strand
    if(as.logical(strand(introns[[tx_id]])[1]=="+")){
      # Tile introns with 1bp ranges
      intron_tiles <- unlist(tile(x = introns[[tx_id]], width = 1))
      first_10Kb <- GenomicRanges::reduce(intron_tiles[1:1e4])
      last_10Kb <- GenomicRanges::reduce(intron_tiles[(length(intron_tiles)-1e4+1):length(intron_tiles)])
      # if tx is on the reverse strand:
    } else {
      # Tile introns with 1bp ranges
      intron_tiles <- unlist(tile(x = introns[[tx_id]], width = 1))
      last_10Kb <- GenomicRanges::reduce(intron_tiles[1:1e4])
      first_10Kb <- GenomicRanges::reduce(intron_tiles[(length(intron_tiles)-1e4+1):length(intron_tiles)])
    }
    
    # Get required data which we will later add to the elongation_rate_data dataframe
    
    # Distance
    transcription_dist <- start(intron_tiles[length(intron_tiles)-5e3])-start(intron_tiles[5e3])
    # first_10Kb
    first_10Kb <- colSums(data.frame(lapply(1:41, allReadsCounts, interval=first_10Kb)))/sfs
    # last_10Kb
    last_10Kb <- colSums(data.frame(lapply(1:41, allReadsCounts, interval=last_10Kb)))/sfs
    
    return(list(transcription_dist, first_10Kb, last_10Kb))
  }
  
  # Get data
  dist_expr <- lapply(elongation_rate_data$tx_id, get_plot_data)
  # Add to elongation_rate_data dataframe
  elongation_rate_data$distance <- sapply(dist_expr, function(x) x[[1]])
  elongation_rate_data$first_10Kb_expression <- lapply(dist_expr, function(x) x[[2]])
  elongation_rate_data$last_10Kb_expression <- lapply(dist_expr, function(x) x[[3]])
  rm(get_plot_data, allReadsCounts, fls, sfs, dist_expr, introns)
  
  # Next we will fit a smooth function to the first/last 10Kb expression profiles. We will use the 9 parameter impulse model that allows up to three transitions
  impulse3 <- function(h0, h1, h2, h3, lambda1, lambda2, t, t1, t2, t3) { (1/h1)*((h0+(h1-h0)*(1/(1+exp(-lambda1*(t-t1)))))*(h2+(h1-h2)*(1/(1+exp(-lambda1*(t-t2)))))*(h3+(h2-h3)*(1/(1+exp(-lambda2*(t-t3)))))) }
  # Define time
  t<-seq(0,400,10)
  # Define function that will fit the impulse3 model and return the best fit after 100 iterations
  fit_impulse3 <- function(input_time_series, t){
    temp <- data.frame(t=t, M=input_time_series)
    best_NRMSD <- 1E6
    best_fit <- NA
    for(i in 1:100){
      # Set fit to 0 to start with
      impulse_fit <- 0
      # Continue trying to the model to the data until a fit is found
      while(is.numeric(impulse_fit)){
        impulse_fit <- tryCatch({
          nlsLM(M ~ impulse3(h0, h1, h2, h3, lambda1,lambda2, t, t1, t2, t3), 
                data = temp, start = c(h0=runif(1, min=min(input_time_series), max = max(input_time_series)),
                                       h1=runif(1, min=min(input_time_series), max = max(input_time_series)),
                                       h2=runif(1, min=min(input_time_series), max = max(input_time_series)),
                                       h3=runif(1, min=min(input_time_series), max = max(input_time_series)),
                                       lambda1=runif(n=1, min=0, max=1),
                                       lambda2=runif(n=1, min=0, max=1),
                                       t1=sample(1:max(t),1),
                                       t2=sample(1:max(t),1),
                                       t3=sample(1:max(t),1)),
                control = list(maxiter = 1024,warnOnly=TRUE), trace = F)
        }, error = function(e) { # if nlsLM can't converge, ouput zero
          c(0)
        }, finally = {})
      }
      # calculate the NRMSD
      current_NRMSD <- sqrt(sum((residuals(impulse_fit))^2)/length(residuals(impulse_fit)))/diff(range(temp$M))
      # Check if current fit is more optimal
      if(current_NRMSD < best_NRMSD){ 
        best_NRMSD <- current_NRMSD
        best_fit <- impulse_fit
      }
    }
    return(best_fit)
  }
  # Fit model to the first and last 10Kb expression profiles
  set.seed(1234)
  for(i in 1:nrow(elongation_rate_data)){
    elongation_rate_data$first_10Kb_fit[i] <- list(fit_impulse3(input_time_series = elongation_rate_data$first_10Kb_expression[[i]], t=t))
    elongation_rate_data$last_10Kb_fit[i] <- list(fit_impulse3(input_time_series = elongation_rate_data$last_10Kb_expression[[i]], t=t))
  }
  
  # Determine time lag (in minutes) that optimizes the fit between the first and last 10Kb
  for(i in 1:nrow(elongation_rate_data)){
    # Get expression profiles at 1-min resolution
    first_10Kb <- predict(elongation_rate_data$first_10Kb_fit[[i]], data.frame(t=seq(0,400,1)))
    last_10Kb <- predict(elongation_rate_data$last_10Kb_fit[[i]], data.frame(t=seq(0,400,1)))
    # Find lag
    best_lag <- ccf(x= last_10Kb, y = first_10Kb, lag.max = 300, type="correlation", na.action = na.omit, plot = FALSE)
    # Save out lag that gives the best correlation
    elongation_rate_data$time[i] <- best_lag$lag[which.max(best_lag$acf)]
  }
  
  rm(impulse3, fit_impulse3, t, i, first_10Kb, last_10Kb, best_lag)
  
  #save
  saveRDS(elongation_rate_data, "data/elongation_rate_data.rds")
}

# Import data
elongation_rate_data <- readRDS("data/elongation_rate_data.rds")

# Re-define impulse3 function to use in plotting
impulse3 <- function(h0, h1, h2, h3, lambda1, lambda2, t, t1, t2, t3) { (1/h1)*((h0+(h1-h0)*(1/(1+exp(-lambda1*(t-t1)))))*(h2+(h1-h2)*(1/(1+exp(-lambda1*(t-t2)))))*(h3+(h2-h3)*(1/(1+exp(-lambda2*(t-t3)))))) }
# Define function to create some explanatory lag plots
lag_plots <- function(gene_name, gene_cols){
  par(mfcol=c(3,1))
  ## Plot 1: first and last 10Kb expression profiles
  # Get scaled expression for the first 10Kb + fit
  t<-seq(0,400,10)
  first_10Kb <- elongation_rate_data$first_10Kb_expression[[which(elongation_rate_data$gene_name==gene_name)]]
  first_10Kb_fit <- predict(elongation_rate_data$first_10Kb_fit[[which(elongation_rate_data$gene_name==gene_name)]], data.frame(t=seq(0,400,1)))
  first_10Kb_fit <- first_10Kb_fit-min(first_10Kb)
  first_10Kb <- first_10Kb-min(first_10Kb)
  first_10Kb_fit <- first_10Kb_fit/max(first_10Kb)
  first_10Kb <- first_10Kb/max(first_10Kb)
  # Get scaled expression for the last 10Kb + fit
  last_10Kb <- elongation_rate_data$last_10Kb_expression[[which(elongation_rate_data$gene_name==gene_name)]]
  last_10Kb_fit <- predict(elongation_rate_data$last_10Kb_fit[[which(elongation_rate_data$gene_name==gene_name)]], data.frame(t=seq(0,400,1)))
  last_10Kb_fit <- last_10Kb_fit-min(last_10Kb)
  last_10Kb <- last_10Kb-min(last_10Kb)
  last_10Kb_fit <- last_10Kb_fit/max(last_10Kb)
  last_10Kb <- last_10Kb/max(last_10Kb)
  #plot
  plot(t, first_10Kb, col=scales::alpha("#4575b4", 0.6), pch=19, ylab="Scaled expression", xlab="Time (min)", main=gene_name, ylim=range(first_10Kb, first_10Kb_fit, last_10Kb, last_10Kb_fit))
  lines(seq(0,400,1), first_10Kb_fit, col="#4575b4", lwd=2)
  points(t, last_10Kb, col=scales::alpha("#d73027", 0.6), pch=19)
  lines(seq(0,400,1), last_10Kb_fit, col="#d73027", lwd=2)
  legend("topright", bty="n", legend=c("First 10 Kb intron", "Last 10 Kb intron"), col=c("#4575b4", "#d73027"), lwd=2)
  
  ## Plot 2: lagged versions of first 10Kb vs last 10Kb
  plot(seq(0,400,1), last_10Kb_fit, col="#d73027", type="l", lwd=2, lty=2, ylab="Scaled expression", xlab="Time (min)", main=gene_name, ylim=range(first_10Kb, first_10Kb_fit, last_10Kb, last_10Kb_fit))
  lines(seq(0,400,1), first_10Kb_fit, col="#4575b4", lwd=2, lty=2)
  lines(seq(0,400,1)[51:400], first_10Kb_fit[1:350], lwd=2, col=gene_cols[4])
  lines(seq(0,400,1)[101:400], first_10Kb_fit[1:300], lwd=2, col=gene_cols[5])
  lines(seq(0,400,1)[151:400], first_10Kb_fit[1:250], lwd=2, col=gene_cols[6])
  lines(seq(0,400,1)[201:400], first_10Kb_fit[1:200], lwd=2, col=gene_cols[7])
  lines(seq(0,400,1)[251:400], first_10Kb_fit[1:150], lwd=2, col=gene_cols[8])
  lines(seq(0,400,1)[301:400], first_10Kb_fit[1:100], lwd=2, col=gene_cols[9])
  legend("topright", bty="n", legend=c("First 10 Kb, 0 min lag", "First 10 Kb, 50 min lag", "First 10 Kb, 100 min lag", "First 100 Kb, 150 min lag", "First 10 Kb, 200 min lag", "First 10 Kb, 250 min lag", "First 10 Kb, 300 min lag", "Last 10 Kb intron"), col=c("#4575b4", gene_cols[4:9],"#d73027"), lwd=2, lty=c(2,1,1,1,1,1,1,2))
  
  ## Plot 3: Cross-correlation vs lag
  first_10Kb <- predict(elongation_rate_data$first_10Kb_fit[[which(elongation_rate_data$gene_name==gene_name)]], data.frame(t=seq(0,400,1)))
  last_10Kb <- predict(elongation_rate_data$last_10Kb_fit[[which(elongation_rate_data$gene_name==gene_name)]], data.frame(t=seq(0,400,1)))
  lags <- ccf(x= last_10Kb, y = first_10Kb, lag.max = 300, type="correlation", na.action = na.omit, plot = FALSE)
  # Plot line for cross correlation as a function of lag
  plot(as.numeric(lags$lag)[301:601], as.numeric(lags$acf)[301:601], type="l", ylab="Cross-correlation", xlab="Time lag (min)", lwd=2, main=gene_name, ylim=c(min((as.numeric(lags$acf)[301:601]))-0.1, max((as.numeric(lags$acf)[301:601]))+0.1))
  # Add vertical line at time that achieves max cross correlation
  segments(x0=lags$lag[which.max(lags$acf)], y0=(min(lags$acf[301:601])-0.1), x1=lags$lag[which.max(lags$acf)], y1=lags$acf[which.max(lags$acf)], col=gene_cols[8], lty=4, lwd=2)
  # Add points for 0, 50, 100, 150, 200, 250 and 300 min lags
  points(as.numeric(lags$lag)[301+seq(0,300,50)], as.numeric(lags$acf)[301+seq(0,300,50)], pch=21, cex=3, col="black", lwd=2, bg=c("white",gene_cols[4:9]))
  # Add legend
  legend("topright", legend = c(glue('Optimum lag = {lags$lag[which.max(lags$acf)]} min')), bty="n", lty=4, lwd=2, col=gene_cols[8])
  par(mfcol=c(1,1))
}

# Save out
pdf("figures/Fig_S2/FigS2a.pdf", width=5, height=10);lag_plots(gene_name = "LDLRAD4", gene_cols = brewer.pal(9, "Greens"));dev.off()
pdf("figures/Fig_S2/FigS2b.pdf", width=5, height=10);lag_plots(gene_name = "PARVA", gene_cols = brewer.pal(9, "Purples"));dev.off()
pdf("figures/Fig_S2/FigS2c.pdf", width=5, height=10);lag_plots(gene_name = "TIMP3", gene_cols = brewer.pal(9, "Oranges"));dev.off()
# include small PNG
png("figures_track/LDLRAD4_lag_plot.png", width=500, height=700);lag_plots(gene_name = "LDLRAD4", gene_cols = brewer.pal(9, "Greens"));dev.off()

# Define function to plot distance vs transcription time
plot_dist_time <- function(){
  plot(x = elongation_rate_data$time, y=elongation_rate_data$distance,
       xlab="Transcription time (min)", ylab=" Distance (Kb)", pch=19, xlim=c(0,250), ylim=c(0,7e5), xaxt="n", yaxt="n", col=rgb(0,0,0,0.8))
  # Add x-axis
  axis(side=1, at=seq(0,250,50), cex.axis=0.6)
  # Add y-axis
  axis(side=2, at=seq(0,6.5e5,1e5), labels = c("0", "100,000", "200,000", "300,000", "400,000", "500,000", "600,000"), las=1, cex.axis=0.6)
  # Fit linear model
  data_fit <- lm(distance ~ time, data=elongation_rate_data)
  # Add dashed red line for model fit
  abline(data_fit, col="black", lty=1, lwd=2)
  # examine fit
  summary(data_fit)
  # Add some of this information to the plot in a legend
  legend("top", legend = c(glue('Adj. R2 = {format(summary(data_fit)$adj.r.squared, digits=3)}')), bty = "n", cex=0.6, lty=1, lwd=2)
}
# Add lines and points for specfic genes
add_pts <- function(gene_name, gene_col){
  gene_time <- elongation_rate_data$time[which(elongation_rate_data$gene_name==gene_name)]
  gene_length <- elongation_rate_data$distance[which(elongation_rate_data$gene_name==gene_name)]
  segments(x0=gene_time, y0=0, x1=gene_time, y1=elongation_rate_data$distance[which(elongation_rate_data$gene_name==gene_name)], col=gene_col, lty=4, lwd=2)
  # tx length
  segments(x0=0, y0=elongation_rate_data$distance[which(elongation_rate_data$gene_name==gene_name)], x1=gene_time, y1=elongation_rate_data$distance[which(elongation_rate_data$gene_name==gene_name)], col=gene_col, lty=2, lwd=2)
  # Point
  points(elongation_rate_data$time[which(elongation_rate_data$gene_name==gene_name)],
         elongation_rate_data$distance[which(elongation_rate_data$gene_name==gene_name)],col="black",bg=gene_col, pch=21, cex=1.5)
}

# Save out
pdf("figures/Fig_S2/FigS2d.pdf", width=5, height=10)
par(mfcol=c(2,1))
plot_dist_time()
# Add points
add_pts(gene_name="LDLRAD4",gene_col = brewer.pal(9, "Greens")[8])
add_pts(gene_name="PARVA",gene_col = brewer.pal(9, "Purples")[8])
add_pts(gene_name="TIMP3",gene_col = brewer.pal(9, "Oranges")[8])
# Histogram of transcription elongation rates
hist((elongation_rate_data$distance/1e3)/elongation_rate_data$time, breaks = 10, xlim=c(0,5), main="",
     ylab="Number of genes", xlab="Transcription elongation rate (Kb/min)", col="#386cb0")
legend("topleft", legend=glue('mean: {round(mean((elongation_rate_data$distance/1e3)/elongation_rate_data$time),1)} Kb/min'), bty="n")
# Histogram of gene lengths
par(mfcol=c(1,1))
dev.off()

# Include a small png
png("figures_track/distance_vs_transcription_time.png", width=500, height=800)
par(mfcol=c(2,1))
plot_dist_time()
# Add points
add_pts(gene_name="LDLRAD4",gene_col = brewer.pal(9, "Greens")[8])
add_pts(gene_name="PARVA",gene_col = brewer.pal(9, "Purples")[8])
add_pts(gene_name="TIMP3",gene_col = brewer.pal(9, "Oranges")[8])
# Histogram of transcription elongation rates
hist((elongation_rate_data$distance/1e3)/elongation_rate_data$time, breaks = 10, xlim=c(0,5), main="",
     ylab="Number of genes", xlab="Transcription elongation rate (Kb/min)", col="#386cb0")
legend("topleft", legend=glue('mean: {round(mean((elongation_rate_data$distance/1e3)/elongation_rate_data$time),1)} Kb/min'), bty="n")
par(mfcol=c(1,1))
dev.off()

# It would also be handy to have plots indicating the gene schematics and how distance is calculated

# Begin by loading GTF data to get intron and exon intervals
TxDb <- makeTxDbFromGFF("data/annotation/GRCh38_spiked.gtf", format="gtf")
exons <- cdsBy(TxDb, use.names=TRUE)
introns <- intronsByTranscript(TxDb, use.names=TRUE)
fiveUTRs <- fiveUTRsByTranscript(TxDb, use.names=TRUE)
threeUTRs <- threeUTRsByTranscript(TxDb, use.names=TRUE)

plot_gene <- function(tx_id, gene_col){
  # Get feature coordinates for current gene
  current_exons <- exons[[tx_id]]
  current_introns <- introns[[tx_id]]
  current_fiveUTRs <- fiveUTRs[[tx_id]]
  current_threeUTRs <- threeUTRs[[tx_id]]
  
  par(lend=2, mar=c(0,0,0,0))
  # Create empty plot to add polygons too
  plot("","", xlim=c(min(start(c(current_exons, current_introns, current_fiveUTRs, current_threeUTRs))),
                     max(end(c(current_exons, current_introns, current_fiveUTRs, current_threeUTRs)))),
       type="n", ylim=c(0.4,0.7), xlab="", ylab="", xaxt="none", yaxt="none", bty="n")
  
  # Add boxes for introns
  for(i in 1:length(current_introns)){
    polygon(x=c(start(current_introns)[i], start(current_introns)[i], end(current_introns)[i], end(current_introns)[i]), y=c(0.495, 0.505, 0.505, 0.495), col="black", border = "black", lwd=1)
  }
  
  # Add boxes for 5'UTR
  for(i in 1:length(current_fiveUTRs)){
    polygon(x=c(start(current_fiveUTRs)[i], start(current_fiveUTRs)[i], end(current_fiveUTRs)[i], end(current_fiveUTRs)[i]), y=c(0.45, 0.55, 0.55, 0.45), col="black", border = "black", lwd=1)
  }
  
  # Add boxes for 3'UTR
  for(i in 1:length(current_threeUTRs)){
    polygon(x=c(start(current_threeUTRs)[i], start(current_threeUTRs)[i], end(current_threeUTRs)[i], end(current_threeUTRs)[i]), y=c(0.45, 0.55, 0.55, 0.45), col="black", border = "black", lwd=1)
  }
  
  # Add boxes for exons
  for(i in 1:length(current_exons)){
    polygon(x=c(start(current_exons)[i], start(current_exons)[i], end(current_exons)[i], end(current_exons)[i]), y=c(0.4, 0.6, 0.6, 0.4), col="black", border = "black", lwd=1)
  }
  
  # Add boxes for first 10Kb of intron
  first_10Kb <- GenomicRanges::reduce(unlist(tile(x = current_introns, width = 1))[1:1e4])
  for(i in 1:length(first_10Kb)){
    if(as.logical(strand(current_exons)[1]=="+")){
      polygon(x=c(start(first_10Kb)[i], start(first_10Kb)[i], end(first_10Kb)[i], end(first_10Kb)[i]), y=c(0.65, 0.7, 0.7, 0.65), col="#4575b4", border = NA)
    } else {
      polygon(x=c(start(first_10Kb)[i], start(first_10Kb)[i], end(first_10Kb)[i], end(first_10Kb)[i]), y=c(0.65, 0.7, 0.7, 0.65), col="#d73027", border = NA)
    }
  }
  
  # Add boxes for last 10Kb of intron
  last_10Kb <- unlist(tile(x = current_introns, width = 1))
  last_10Kb <- GenomicRanges::reduce(last_10Kb[(length(last_10Kb)-1e4+1):length(last_10Kb)])
  for(i in 1:length(last_10Kb)){
    if(as.logical(strand(current_exons)[1]=="+")){
      polygon(x=c(start(last_10Kb)[i], start(last_10Kb)[i], end(last_10Kb)[i], end(last_10Kb)[i]), y=c(0.65, 0.7, 0.7, 0.65), col="#d73027", border = NA)
    } else {
      polygon(x=c(start(last_10Kb)[i], start(last_10Kb)[i], end(last_10Kb)[i], end(last_10Kb)[i]), y=c(0.65, 0.7, 0.7, 0.65), col="#4575b4", border = NA)
    }
  }
  
  # Add line indicating distance
  intron_tiles <- unlist(tile(x = introns[[tx_id]], width = 1))
  segments(x0=start(intron_tiles[5e3]), y0=0.71, x1=start(intron_tiles[length(intron_tiles)-5e3]), y1=0.71, col = gene_col, lwd=2, lty=2)
}

# Save out the gene structures
pdf("figures/Fig_S2/FigS2e.pdf", width=5, height=10)
par(mfcol=c(3,1))
plot_gene(elongation_rate_data$tx_id[elongation_rate_data$gene_name=="LDLRAD4"], gene_col = brewer.pal(9, "Greens")[8])
plot_gene(elongation_rate_data$tx_id[elongation_rate_data$gene_name=="PARVA"], gene_col = brewer.pal(9, "Purples")[8])
plot_gene(elongation_rate_data$tx_id[elongation_rate_data$gene_name=="TIMP3"], gene_col = brewer.pal(9, "Oranges")[8])
par(mfcol=c(1,1))
dev.off()

print("Finished creating panels for Fig_S2")
