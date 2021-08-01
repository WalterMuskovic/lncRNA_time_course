### Plot spatial correlation trend with lncRNAs that do and do not overlap annotated lncRNAs (Supplementary Figure 4)

library(tidyverse)
library(mgcv)
library(vioplot)

# Get data frames with distance and correlation (data was created in `create_data_for_Fig_6.R` script):
# Protein-coding genes and lncRNAs that overlap annotated lncRNAs
c_nc_anno <- readRDS("data/c_nc_anno.rds")
# Protein-coding genes and lncRNAs with no overlap with annotated lncRNAs
c_nc_no_anno <- readRDS("data/c_nc_no_anno.rds")

# Fit GAM to the coding-lncRNAs with overlap
gam_fit_c_nc_anno <- gam(corr ~ s(log10_dist), data=c_nc_anno, method = "REML")
# Examine a summary of the fit
anova(gam_fit_c_nc_anno)
#Approximate significance of smooth terms:
#  edf Ref.df     F p-value
#s(log10_dist) 7.183  8.214 35.04  <2e-16
# Get the percentiles
quantiles_c_nc_anno <- readRDS("data/S4_overlapping_block_bootstrap_quantiles.rds")

# Fit GAM to the coding-lncRNAs with no overlap
gam_fit_c_nc_no_anno <- gam(corr ~ s(log10_dist), data=c_nc_no_anno, method = "REML")
# Examine a summary of the fit
anova(gam_fit_c_nc_no_anno)
#Approximate significance of smooth terms:
#  edf Ref.df     F p-value
#s(log10_dist) 7.769  8.614 50.64  <2e-16
# Get the percentiles
quantiles_c_nc_no_anno <- readRDS("data/S4_not_overlapping_block_bootstrap_quantiles.rds")

#Define functions to create the required plot. 
# The first function will produce violin plots of Pearson's correlation values between genes at binned distances.
# The second function will overlay the GAM fit over the binned data. 
# The third function will plot the confidence intervals obtained from the block bootstrap over the binned data.
# Finally a wrapper function will combine everything in one plot.

## 1. Define plotting function to plot relationship between correlation and distance using violin plots
vioplot_corr_dist <- function(bounds = c(0, 1E4, 2.5E4, 5E4, 1E5, 2.5E5, 5E5, 1E6, 5E6, 1E7, 2.5E7),
                              bin_labels = c("<10Kb", "10-25Kb", "25-50Kb", "50-100Kb", "100-250Kb", "250-500Kb", "500Kb-1Mb", "1-5Mb", "5-10Mb", "10-25Mb", ">25Mb"),
                              input_df,
                              plot_title,
                              x_axis_label="Distance between coding and lncRNA transcripts"){
  # Add the maximum distance to the bin bounds
  bounds <- c(bounds, max(input_df$dist))
  
  # Create bins
  bins<-list()
  for(i in 1:c(length(bounds)-1)){
    bins <- append(bins, list(input_df$corr[which(input_df$dist>=bounds[i] & input_df$dist<=bounds[i+1])]))
  }
  
  # Create violin plots
  par(mar=c(7,4,3,1)) # bottom, left, top and right margins
  vioplot(bins,
          names=NA, 
          las=2, cex.lab=1.3,cex.axis=0.8,
          xlab="", ylab="Pearson correlation coefficient", main=plot_title, h=0.1, wex=0.8,
          col=rgb(0.7568627, 0.8392157, 0.9176471), rectCol=rgb(0.2, 0.4, 0.6), lineCol = rgb(0.2, 0.4, 0.6),
          pchMed=19, colMed="white")
  axis(side=1, at=1:length(bin_labels), labels = FALSE, tck=-0.025)
  text(1:length(bin_labels), par("usr")[3] - 0.2, labels = bin_labels, srt = 45, pos = 1, xpd = TRUE)
  mtext(x_axis_label, side=1, line=5, cex=1)
}

## 2. Define another function to overlay the GAM fit over the binned data
# Define function to return mean value for data binned by bounds vector, so that we can plot the continuous GAM fit over the boxplots
boxplot_line <- function(bounds=c(0, 1E4, 2.5E4, 5E4, 1E5, 2.5E5, 5E5, 1E6, 5E6, 1E7, 2.5E7), gam_vector, distance_vector){
  bin_means <- c()
  for(i in 1:(length(bounds)-1)){
    bin_means[i] <- mean(gam_vector[distance_vector>=bounds[i] & distance_vector<bounds[i+1]])
  }
  bin_means[length(bounds)] <- mean(gam_vector[distance_vector>bounds[length(bounds)]])
  return(bin_means)  
}

## 3. Define a third function to add the confidence intervals obtained using the block bootsrap 
# Write function that, similarly to boxplot_line() above, will plot the average of the continuous confidence intervals for each binned distance

# Import block bootstrap quantile data
confInt_lines <- function(bounds=c(0, 1E4, 2.5E4, 5E4, 1E5, 2.5E5, 5E5, 1E6, 5E6, 1E7, 2.5E7), quantiles_vectors, distance_vector){
  quantile_means <- data.frame(matrix(nrow=6, ncol=length(bounds)))
  row.names(quantile_means) <- c("1%","5%","25%","75%","95%","99%")
  for(n in 1:nrow(quantile_means)){
    for(i in 1:(length(bounds)-1)){
      quantile_means[n,i] <- mean(quantiles_vectors[n,][distance_vector>=bounds[i] & distance_vector<bounds[i+1]])
    }
    quantile_means[n,length(bounds)] <- mean(quantiles_vectors[n,][distance_vector>bounds[length(bounds)]])
  }
  return(quantile_means)
}

# Create function to combine all of the above elements in one plot
fig_S4 <- function(){
  par(mfcol=c(1,2))
  # Violin plot 1 - with overlap
  vioplot_corr_dist(input_df = c_nc_anno, plot_title = "lncRNAs that overlap annotated lncRNAs",x_axis_label="Distance between coding gene and lncRNA")
  # Add line with the mean of the GAM fit for each distance bin
  lines(1:11, boxplot_line(gam_vector = as.numeric(predict(gam_fit_c_nc_anno), type="terms"), distance_vector=c_nc_anno$dist),col="red", lwd=2)
  # Add confidence intervals
  quantiles <- confInt_lines(quantiles_vectors =quantiles_c_nc_anno, distance_vector = c_nc_anno$dist)
  polygon(c(1:11, rev(1:11)), c(quantiles["1%",], rev(quantiles["99%",])), col = adjustcolor("#de2d26", alpha=0.4), border = NA)
  
  # Violin plot 2 - without overlap
  vioplot_corr_dist(input_df = c_nc_no_anno, plot_title = "lncRNAs that do not overlap annotated lncRNAs", x_axis_label="Distance between coding gene and lncRNA")
  # Add line with the mean of the GAM fit for each distance bin
  lines(1:11, boxplot_line(gam_vector = as.numeric(predict(gam_fit_c_nc_no_anno), type="terms"), distance_vector=c_nc_no_anno$dist),col="red", lwd=2)
  # Add confidence intervals
  quantiles <- confInt_lines(quantiles_vectors =quantiles_c_nc_no_anno, distance_vector = c_nc_no_anno$dist)
  polygon(c(1:11, rev(1:11)), c(quantiles["1%",], rev(quantiles["99%",])), col = adjustcolor("#de2d26", alpha=0.4), border = NA)
}

# Check to see if figure directory has already been created
if(!dir.exists("figures/Fig_S4")){dir.create("figures/Fig_S4")}
# Save out figure and plot a copy in the markdown document
pdf("figures/Fig_S4/lncRNA_overlap_spatial_correlation.pdf", width=14, height=5)
fig_S4()
dev.off()

print("Finished creating Fig. S4")
