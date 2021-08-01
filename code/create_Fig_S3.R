### Plot spatial correlation trend among protein-coding genes and lncRNAs
### (Supplementary Figure 3)

library(tidyverse)
library(mgcv)
library(vioplot)

# Get data frames with distance and correlation between (data was created in `create_data_for_Fig_6.R` script):
# every coding transcript 
c_c <- readRDS("data/c_c.rds")
# every non-coding transcript
nc_nc <- readRDS("data/nc_nc.rds")

# Fit GAM to the coding-coding data
gam_fit_c_c <- gam(corr ~ s(log10_dist), data=c_c, method = "REML")

# Fit GAM to the lncRNA-lncRNA data
gam_fit_nc_nc <- gam(corr ~ s(log10_dist), data=nc_nc, method = "REML")

# Define functions to create the required plot.
# The first function will produce violin plots of Pearson's correlation values between genes at binned distances.
# The second function will overlay the GAM fit over the binned data. 
# Finally the third wrapper function will combine everything in one plot.

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

## 3. Create function to combine all of the above elements in one plot
fig_S3 <- function(){
  par(mfcol=c(1,2))
  # Violin plot 1 - coding/coding
  bounds <- c(0, 1E4, 2.5E4, 1E5, 2.5E5, 5E5, 1E6, 5E6, 1E7, 2.5E7)
  bin_labels = c("<10Kb", "10-25Kb", "25-100Kb", "100-250Kb", "250-500Kb", "500Kb-1Mb", "1-5Mb", "5-10Mb", "10-25Mb", ">25Mb")
  vioplot_corr_dist(input_df = c_c, plot_title = "Coding/coding genomic distance \nvs. expression correlation", bounds = bounds, bin_labels=bin_labels,
                    x_axis_label="Distance between coding transcripts")
  # Add line with the mean of the GAM fit for each distance bin
  lines(1:10, boxplot_line(gam_vector = as.numeric(predict(gam_fit_c_c), type="terms"), distance_vector=c_c$dist, bounds = bounds),col="red", lwd=2)
  
  # Violin plot 2 - lncRNA/lncRNA
  vioplot_corr_dist(input_df = nc_nc, plot_title = "lncRNA/lncRNA genomic distance \nvs. expression correlation", bounds = bounds, bin_labels=bin_labels,
                    x_axis_label="Distance between lncRNA transcripts")
  # Add line with the mean of the GAM fit for each distance bin
  lines(1:10, boxplot_line(gam_vector = as.numeric(predict(gam_fit_nc_nc), type="terms"), distance_vector=nc_nc$dist, bounds = bounds),col="red", lwd=2, bounds = bounds)
}

# Check to see if figure directory has already been created
if(!dir.exists("figures/Fig_S3")){dir.create("figures/Fig_S3")}

# Save out figure and plot a copy in the markdown document
pdf("figures/Fig_S3/lncRNA_coding_spatial_correlation.pdf", width=14, height=5)
fig_S3()
dev.off()

# Clean up
rm(c_c, nc_nc, gam_fit_c_c, gam_fit_nc_nc, fig_S3, boxplot_line, vioplot_corr_dist)

print("Finished creating Fig. S3")
