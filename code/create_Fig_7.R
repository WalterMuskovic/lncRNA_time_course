### Plot relationship between genomic distance and coding/lncRNA expression profiles

library(tidyverse)
library(vioplot)
library(mgcv)

#Define functions to create the required plot. The first function will produce violin plots of Pearson's correlation values between genes at binned distances. The second function will overlay the GAM fit over the binned data. The third function will plot the confidence intervals obtained from the block bootstrap over the binned data. Finally a wrapper function will combine everything in one plot.
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
confInt_lines <- function(bounds=c(0, 1E4, 2.5E4, 5E4, 1E5, 2.5E5, 5E5, 1E6, 5E6, 1E7, 2.5E7), quantiles_vectors = readRDS("data/block_bootstrap_quantiles_mouse.rds"), distance_vector){
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
c_nc <- readRDS("data/c_nc_mouse.rds")
gam_fit <- gam(corr ~ s(log10_dist), data=c_nc, method = "REML")
fig7a <- function(){
  # Violin plot of the relationship between genomic distance and expression profile correlation
  vioplot_corr_dist(input_df = c_nc, plot_title = "Coding/lncRNA genomic distance \nvs. expression correlation")
  # Add line with the mean of the GAM fit for each distance bin
  lines(1:11, boxplot_line(gam_vector = as.numeric(predict(gam_fit), type="terms"), distance_vector=c_nc$dist),col="red", lwd=2)
  # Add confidence intervals
  quantiles <- confInt_lines(quantiles_vectors = readRDS("data/block_bootstrap_quantiles_mouse.rds"), distance_vector = c_nc$dist)
  polygon(c(1:11, rev(1:11)), c(quantiles["1%",], rev(quantiles["99%",])), col = adjustcolor("#de2d26", alpha=0.4), border = NA)
  # Add legend
  par(xpd=TRUE)
  legend(x = -1, y = 1.4, legend="1st-99th percentile", bty="n", col=adjustcolor("#de2d26", alpha=0.4), pch=15)
}

# Save out figure
pdf("figures/Fig_7/corr_dist_plot_mouse.pdf", width=7, height=5)
fig7a()
dev.off()



### Plot time-lagged data
# Read data in
lagged_cor <- readRDS("data/c_nc_lagged_corr_mouse.rds")
lagged_cor_bbstrap <- readRDS("data/lagged_cor_bbstrap_mouse.rds")

# Add distance info
c_nc <- readRDS("data/c_nc_mouse.rds")
lagged_cor$dist <- c_nc$dist[match(str_glue('{lagged_cor$row_vec}###{lagged_cor$col_vec}'),
                                   str_glue('{c_nc$row_vec}###{c_nc$col_vec}'))]
lagged_cor <- lagged_cor[!is.na(lagged_cor$dist),]

# Define function to plot lags
fig7bc <- function(){
  par(mfcol=c(5,2))
  par(mar=c(2.5,4,0.5,1))
  for(i in 1:10){
    # Define variables for plotting
    distances <- list(lagged_cor$dist<1E4, lagged_cor$dist<1E5, lagged_cor$dist<2.5E5, lagged_cor$dist<1E6, lagged_cor$dist>1E6, lagged_cor$dist<1E4, lagged_cor$dist<1E5,  lagged_cor$dist<2.5E5, lagged_cor$dist<1E6, lagged_cor$dist>1E6)[[i]]
    distances_text <- str_glue('mRNA/lncRNA pairs \n{c("< 10 Kb", "< 100 Kb", "< 250 Kb", "< 1 Mb", "> 1 Mb", "< 10 Kb", "< 100 Kb", "< 250 Kb", "< 1 Mb", "> 1 Mb")} apart')
    cols <- c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00")[i]
    mRNA_measurement_type <- c(rep("premRNA_lag_",5),rep("lag_mRNA_",5))[i]
    # Create a blank plot to add lines and polygons to
    plot("", xlab="Time lag (min)", ylab="Cross correlation", ylim=c(-0.4,0.7), xlim=c(-90,90), cex.axis=0.8, cex.lab=0.9,tck=-0.05, xaxt='n')
    # Add dashes for every time point on the x-axis
    axis(side=1, at=seq(-90, 90, 15), labels=c(-90, NA, -60, NA, -30, NA, 0, NA, 30, NA, 60, NA, 90), tck=-0.025)
    # Add a verical line to indicate 0 min
    abline(v=0, lty=2, lwd=1, col="darkgrey")
    # Add polygon for 99% confidence interval obtained using the block bootstrap
    polygon(c(seq(-90,90,15), rev(seq(-90,90,15))), c(lagged_cor_bbstrap[[i]]["1%",], rev(lagged_cor_bbstrap[[i]]["99%",])),col=adjustcolor(cols, alpha=0.3), border = NA)
    # Add horizontal line at rp=0
    abline(h=0)
    # add lines
    lines(seq(-90,90,15), colMeans(lagged_cor[distances, str_detect(colnames(lagged_cor), mRNA_measurement_type)], na.rm = TRUE), col=cols, lwd=2)
    # Add legend 1
    legend("topright", legend=distances_text[i], bty="n", lty=1, lwd=2, col=cols)
    # Add legend 2
    legend("topleft", legend=c("Simulation envelope",expression('1'^"st"*'-99'^"th"*' percentile')), bty="n", fill=c(NA,adjustcolor(cols, alpha=0.3)), border=NA)
  }
  par(mfcol=c(1,1))
  par(mar=c(5.1, 4.1, 4.1, 2.1))
}

# Save out figure
pdf("figures/Fig_7/lag_plot_mouse.pdf", width=7, height=10)
fig7bc()
dev.off()

print("Finished creating all panels for Fig. 7")
