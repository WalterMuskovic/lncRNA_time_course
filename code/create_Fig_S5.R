### Plot spatial correlation trend with promoter-/enhancer-associated lncRNAs (Supplementary Figure 5)

library(tidyverse)
library(mgcv)
library(vioplot)

# Get data frames with distance and correlation (data was created in `create_data_for_Fig_6.R` script):
# Promoter-associated
c_nc_promoter <- readRDS("data/c_nc_promoter.rds")
# Enhancer-associated
c_nc_enhancer <- readRDS("data/c_nc_enhancer.rds")

# Fit GAM - promoter
gam_fit_c_nc_promoter <- gam(corr ~ s(log10_dist), data=c_nc_promoter, method = "REML")
# Examine a summary of the fit
anova(gam_fit_c_nc_promoter)
#Approximate significance of smooth terms:
#  edf Ref.df     F p-value
#s(log10_dist) 6.492  7.636 34.15  <2e-16
# Get the percentiles
quantiles_c_nc_promoter <- readRDS("data/S5_promoter_block_bootstrap_quantiles.rds")

# Fit GAM - enhancer
gam_fit_c_nc_enhancer <- gam(corr ~ s(log10_dist), data=c_nc_enhancer, method = "REML")
# Examine a summary of the fit
anova(gam_fit_c_nc_enhancer)
#Approximate significance of smooth terms:
#  edf Ref.df     F p-value
#s(log10_dist) 6.733  7.861 50.47  <2e-16
# Get the percentiles
quantiles_c_nc_enhancer <- readRDS("data/S5_enhancer_block_bootstrap_quantiles.rds")

#Define functions to create the required plot. 
# The first function will produce violin plots of Pearson's correlation values between genes at binned distances.
# The second function will overlay the GAM fit over the binned data. 
# The third function will plot the confidence intervals obtained from the block bootstrap over the binned data.
# Finally a wrapper function will combine everything in one plot.



### Plot relationship between genomic distance and coding/lncRNA expression profiles

#Define functions to create the required plot. The first function will produce violin plots of Pearson's correlation values between genes at binned distances. The second function will overlay the GAM fit over the binned data. The third function will plot the confidence intervals obtained from the block bootstrap over the binned data. Finally a wrapper function will combine everything in one plot.
## 1. Define plotting function to plot relationship between correlation and distance using violin plots
vioplot_corr_dist <- function(bounds = c(0, 1E4, 5E4, 1E5, 2.5E5, 5E5, 1E6, 5E6, 1E7, 2.5E7),
                              bin_labels = c("<10Kb", "10-50Kb", "50-100Kb", "100-250Kb", "250-500Kb", "500Kb-1Mb", "1-5Mb", "5-10Mb", "10-25Mb"),
                              input_df,
                              plot_title,
                              x_axis_label="Distance between coding and lncRNA transcripts"){
  # Add the maximum distance to the bin bounds
  bounds <- c(bounds, max(input_df$dist))
  
  # Create bins
  bins<-list()
  for(i in 1:c(length(bounds)-2)){
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
boxplot_line <- function(bounds=c(0, 1E4, 5E4, 1E5, 2.5E5, 5E5, 1E6, 5E6, 1E7, 2.5E7), gam_vector, distance_vector){
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
confInt_lines <- function(bounds=c(0, 1E4, 5E4, 1E5, 2.5E5, 5E5, 1E6, 5E6, 1E7, 2.5E7), quantiles_vectors = readRDS("data/XXX_quantiles.rds"), distance_vector){
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

# Create function to [combine all of the above elements in one plot fig S5a
c_nc <- readRDS("data/c_nc_promoter.rds")
gam_fit <- gam(corr ~ s(log10_dist), data=c_nc,method = "REML")
figS5a <- function(){
  # Violin plot of the relationship between genomic distance and expression profile correlation
  vioplot_corr_dist(input_df = c_nc, plot_title = "Promoter-associated lncRNAs")
  # Add line with the mean of the GAM fit for each distance bin
  lines(1:10, boxplot_line(gam_vector = as.numeric(predict(gam_fit), type="terms"), distance_vector=c_nc$dist),col="red", lwd=2)
  # Add confidence intervals
  quantiles <- confInt_lines(quantiles_vectors = quantiles_c_nc_promoter, distance_vector = c_nc$dist)
  polygon(c(1:9, rev(1:9)), c(quantiles["1%",1:9], rev(quantiles["99%",1:9])), col = adjustcolor("#de2d26", alpha=0.4), border = NA)
  # Add legend
  par(xpd=TRUE)
  legend(x = -1, y = 1.4, legend="1st-99th percentile", bty="n", col=adjustcolor("#de2d26", alpha=0.4), pch=15)
}

# Save out figure
pdf("figures/Fig_S5/Fig_S5a.pdf", width=7, height=5)
figS5a()
dev.off()

# Create function to [combine all of the above elements in one plot fig S5b
c_nc <- readRDS("data/c_nc_enhancer.rds")
gam_fit <- gam(corr ~ s(log10_dist), data=c_nc,method = "REML")
figS5b <- function(){
  # Violin plot of the relationship between genomic distance and expression profile correlation
  vioplot_corr_dist(input_df = c_nc, plot_title = "Enhancer-associated lncRNAs")
  # Add line with the mean of the GAM fit for each distance bin
  lines(1:10, boxplot_line(gam_vector = as.numeric(predict(gam_fit), type="terms"), distance_vector=c_nc$dist),col="red", lwd=2)
  # Add confidence intervals
  quantiles <- confInt_lines(quantiles_vectors = quantiles_c_nc_enhancer, distance_vector = c_nc$dist)
  polygon(c(1:9, rev(1:9)), c(quantiles["1%",1:9], rev(quantiles["99%",1:9])), col = adjustcolor("#de2d26", alpha=0.4), border = NA)
  # Add legend
  par(xpd=TRUE)
  legend(x = -1, y = 1.4, legend="1st-99th percentile", bty="n", col=adjustcolor("#de2d26", alpha=0.4), pch=15)
}

# Save out figure
pdf("figures/Fig_S5/Fig_S5b.pdf", width=7, height=5)
figS5b()
dev.off()


### In response to a reviewer question:  for each lncRNA's TSS, plot the
### distance to the nearest protein-coding TSS

# Get data
t98g <- readRDS("data/t98g_filtered_coords.rds") %>% 
  select(transcript_id, seqid, coords, transcript_type) %>%
  mutate(transcript_type = case_when(
    transcript_type == "novel_lncRNA" ~ "lncRNA",
    transcript_type == "annotated_lncRNA" ~ "lncRNA",
    TRUE ~ "protein_coding"
  ))

# Define function to return distance of nearest transcript of opposite type
dist_fun <- function(i){
  subset.transcript_type <-t98g$transcript_type[i]
  subset.seqid <- t98g$seqid[i]
  t98g.subset <- filter(t98g, seqid==subset.seqid,
                        transcript_type != subset.transcript_type)
  return(min(abs(t98g$coords[i] - t98g.subset$coords)))
}

# Apply function
t98g$nearest_pc <- sapply(1:nrow(t98g), dist_fun)

# Restrict to lncRNAs
t98g <- filter(t98g, transcript_type=="lncRNA")

# Check how many are within 10Kb
table(t98g$nearest_pc<10000)/nrow(t98g)*100
#FALSE     TRUE 
#88.33393 11.66607 

# What is the median distance
median(t98g$nearest_pc)
# [1] 139263

# Add reg info so we can split into enhancer-associated and promoter-associated
t98g <- merge(t98g, select(readRDS("data/h.lnc.rds"), transcript_id, reg), by="transcript_id") %>%
  filter(reg%in%c("promoter", "enhancer")) %>%
  mutate(reg = case_when(
    reg == "promoter" ~ "promoter-associated",
    reg == "enhancer" ~ "enhancer-associated"
  ))

# Plot
p1 <- ggplot(t98g, aes(x=log10(nearest_pc))) +
  geom_histogram(aes(y=..density..),bins=30) + 
  geom_density(alpha=.2, aes(fill=reg)) +
  facet_wrap(~reg, ) +
  xlab("log10(distance between lncRNA TSS and nearest protein-coding TSS (bp))") +
  guides(fill=guide_legend(title="lncRNA")) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

# Save out figure
ggsave(filename = "figures/Fig_S5/Fig_S5c.pdf",
       plot = p1,
       device = "pdf", width = 30, height = 12, units = "cm")
