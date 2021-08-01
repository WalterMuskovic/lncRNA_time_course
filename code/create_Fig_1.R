### Heatmap summarizing lncRNA  expression (Fig. 1b)
# Having generated all of our clustering info with 'get_data_for_Fig_1.R', we can create the heatmaps now

# Load required R packages
library(tidyverse)
library(ggsci)
library(RColorBrewer)
library(gplots)

# Load required data
counts <- readRDS("data/t98g_filtered_coords_clust.rds")

# Check to see if figure directories already exist and if not, create them
if(!dir.exists("figures")){dir.create("figures")}
if(!dir.exists("figures/Fig_1")){dir.create("figures/Fig_1")}

# Get just non-coding transcript data and sort by cluster number
noncoding <-  filter(counts, transcript_type!="protein_coding") %>%
  arrange(lncrna_cluster, coords)

# Get colours for sidebars which indicate cluster membership
cols <- ggsci::pal_npg("nrc", alpha = 1)(length(unique(noncoding$lncrna_cluster)))
cols <-  cols[noncoding$lncrna_cluster]

# Restrict to just expression data and scale rows to have mean zero and standard
# deviation one
plot_data <- noncoding %>%
  dplyr::select(transcript) %>%
  unnest(cols = transcript)
# Scale the data
plot_data <- t(apply(plot_data, 1, scale))

# Note that we scale the input expression matrix on rows (genes) to Z-scores to
# facilitate meaningful visual comparison. Some genes have extreme values (a
# very rapid induction - and correspondingly rapid decrease) and they tend to
# make the rest of the heatmap a little faint. What percentage of values are
# outside the z-score range of -3 to 3?
paste0(round(100*sum(as.vector(plot_data)>2 |  as.vector(plot_data)< -2)/length(as.vector(plot_data)),2),"%")
# To get around this we'll use the breaks argument
breaks <- seq(-2, 2, by=0.2)
# add outliers
breaks <- append(breaks, max(as.vector(plot_data)))
breaks <- append(breaks, min(as.vector(plot_data)),0)

# Get colours for the heatmap
heatmap_cols <- colorRampPalette(brewer.pal(9,"PuBuGn"))(length(breaks)-1)

# Before plotting, a quick note on controlling panels in the heatmap:
# The graphics device is divided into a grid where each element of the heatmap
# will be plotted: color key, dendrograms (not plotted here) and the heatmap. By
# default four components will be displayed in the plot. At the top left is the
# color key, top right is the column dendrogram, bottom left is the row
# dendrogram, bottom right is the image plot.  When you add the colsideColors or
# rowsidecolors the grid expands by 1 in the appropriate dimension (e.g. 1 more
# row when adding colsidecolors). This layout can be overriden by specifiying
# appropriate values for lmat, lwid, and lhei. lhei is the relative heights of
# the rows in the plot. lwid = column width. lmat controls the relative
# postition of each element. The device is divided up into as many rows and
# columns as there are in the lmat matrix. Each value in the matrix must be 0 or
# a positive integer

# Define function to plot fig1b
fig1b <- function(){
  heatmap.2(plot_data,
            Rowv=NA, Colv=NA,
            col=heatmap_cols,
            scale="none", # already done
            density.info="none",
            trace="none",
            main="",
            RowSideColors=cols,
            dendrogram="none",
            key=FALSE,
            labRow="",
            labCol=c("0","","","","","50","","","","","100","","","","","150","","","","","200","","","","","250","","","","","300","","","","","350","","","","","400"), srtCol = 0, cexCol=1.1,
            adjCol=c(NA,0.2), # shift x-axis labels slightly left or right
            offsetCol=0.5, # move x axis labels up and down
            breaks=breaks,
            lmat=rbind(c(9,5,6), c(4,1,2),c(7,8,3)),
            lwid = c(1,1,9), # 1st column - no dendrogram, 2nd column - coloured bars, 3rd column -heatmap
            lhei = c(1,10,1), # 1st row - key, 2nd row heatmap
            margins = c(0.5, 5),
            useRaster = TRUE)
}

# Plot heatmap and save out as pdf - later to be assembled as multipanel pdf
print("plotting Fig. 1b")
pdf("figures/Fig_1/Fig1b.pdf", width=8, height=10)
fig1b()
dev.off()

# save a small png as well
png("figures_track/Fig1b.png", width=600, height=1200)
fig1b()
dev.off()

# Function to plot row z-score key
color.bar <- function(lut=heatmap_cols, min=-2, max=2, nticks=5, ticks=seq(min, max, len=nticks), title='', rectangle_edges=breaks) {
  rectangle_edges[1] <- min
  rectangle_edges[length(rectangle_edges)] <- max
  
  plot(c(min,max),c(0,10), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(1, ticks, las=1, labels=FALSE, lwd=2)
  for (i in 1:(length(lut))) {
    x1=rectangle_edges[i]
    x2=rectangle_edges[i+1]
    rect(x1,0,x2,10, col=lut[i], border=NA)
  }
}

# Save out the row z-score key
pdf("figures/Fig_1/Fig1b_zscore_key.pdf", width=3, height=0.5)
par(mar=c(1,0.2,0.2,0.2)) # bottom, left, top and right margins
color.bar()
dev.off()

# Clean up 
rm(breaks, heatmap_cols, cols, plot_data, noncoding, fig1b, color.bar)

### Heatmap summarizing protein-coding (mRNA) expression (Fig. 1c)
# Get just coding data and sort by cluster number
coding <-  filter(counts, transcript_type=="protein_coding" & complete.cases(mrna_cluster)) %>%
  arrange(mrna_cluster, coords)

# Get colours for sidebars which indicate cluster membership
cols <- ggsci::pal_npg("nrc", alpha = 1)(length(unique(coding$mrna_cluster)))
cols <-  cols[coding$mrna_cluster]

# Restrict to just expression data and scale rows to have mean zero and standard deviation one
plot_data <- coding %>%
  dplyr::select(exon) %>%
  unnest(cols=exon)
# Scale the data
plot_data <- t(apply(plot_data, 1, scale))

# Note that we again scale the input expression matrix on rows (genes) to Z-scores to facilitate meaningful visual comparison. What percentage of values are outside the z-score range of -2 to 2?
paste0(round(100*sum(as.vector(plot_data)>2 |  as.vector(plot_data)< -2)/length(as.vector(plot_data)),2),"%")
# To get around this we'll use the breaks argument
breaks <- seq(-2, 2, by=0.2)
# add outliers
breaks <- append(breaks, max(as.vector(plot_data)))
breaks <- append(breaks, min(as.vector(plot_data)),0)

# Get colours for the heatmap
heatmap_cols <- colorRampPalette(brewer.pal(9,"PuBuGn"))(length(breaks)-1)

fig1c <- function(){
  library(gplots)
  heatmap.2(plot_data,
            Rowv=NA, Colv=NA,
            col=heatmap_cols,
            scale="none", # already done
            density.info="none",
            trace="none",
            main="",
            RowSideColors=cols,
            dendrogram="none",
            key=FALSE,
            labRow="",
            labCol=c("0","","","","","50","","","","","100","","","","","150","","","","","200","","","","","250","","","","","300","","","","","350","","","","","400"), srtCol = 0, cexCol=1.1,
            adjCol=c(NA,0.2), # shift x-axis labels slightly left or right
            offsetCol=0.5, # move x axis labels up and down
            breaks=breaks,
            lmat=rbind(c(9,5,6), c(4,1,2),c(7,8,3)),
            lwid = c(1,1,9), # 1st column - no dendrogram, 2nd column - coloured bars, 3rd column -heatmap
            lhei = c(1,10,1), # 1st row - key, 2nd row heatmap
            margins = c(0.5, 5),
            useRaster = TRUE)
}

# Plot heatmap and save out as pdf - later to be assembled as multipanel pdf
print("plotting Fig. 1c")
pdf("figures/Fig_1/Fig1c.pdf", width=8, height=10)
fig1c()
dev.off()

# save a small png as well
png("figures_track/Fig1c.png", width=600, height=1200)
fig1c()
dev.off()

# Function to plot row z-score key
color.bar <- function(lut=heatmap_cols, min=-2, max=2, nticks=5, ticks=seq(min, max, len=nticks), title='', rectangle_edges=breaks) {
  rectangle_edges[1] <- min
  rectangle_edges[length(rectangle_edges)] <- max
  
  plot(c(min,max),c(0,10), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(1, ticks, las=1, labels=FALSE, lwd=2)
  for (i in 1:(length(lut))) {
    x1=rectangle_edges[i]
    x2=rectangle_edges[i+1]
    rect(x1,0,x2,10, col=lut[i], border=NA)
  }
}

# Save out the row z-score key
pdf("figures/Fig_1/Fig1c_zscore_key.pdf", width=3, height=0.5)
par(mar=c(1,0.2,0.2,0.2)) # bottom, left, top and right margins
color.bar()
dev.off()

# Clean up 
rm(breaks, heatmap_cols, cols, plot_data, coding, fig1c, color.bar)


### Plot matrix of cluster centroid correlations (lncRNA/mRNA - Fig. 1d)
#Comparison of lncRNA and mRNA centroids.

# Define function to find the mRNA cluster centroids, 5th and 95th percentiles of expression values
get_mRNA_centroids <- function(input_df){
  # Get number of clusters
  num_clusters <- length(unique(input_df$mrna_cluster))
  # Restrict to just expression data and scale rows to have mean zero and standard deviation one
  expr_only <- input_df %>%
    dplyr::select(exon) %>%
    unnest(cols = c(exon))
  # Scale the data
  expr_only <- t(apply(expr_only, 1, scale))
  
  # define data frame to hold centroids
  centroids <- data.frame(matrix(nrow=num_clusters, ncol=41))
  # define data frame to hold 5th percentile
  cent5 <- data.frame(matrix(nrow=num_clusters, ncol=41))
  # define data frame to hold 95th percentile
  cent95 <- data.frame(matrix(nrow=num_clusters, ncol=41))
  for(n in 1:num_clusters){
    # Select one of the clusters
    cluster_x <- as.matrix(expr_only[input_df$mrna_cluster==n,])
    centroid<-colMeans(cluster_x)
    # Save cluster centroid
    centroids[n,] <- centroid
    # Save 5th percentile of the cluster
    cent5[n,]<-apply(cluster_x , 2 , quantile , probs = c(0.05) , na.rm = TRUE )
    # Save 95th percentile of the cluster
    cent95[n,]<-apply(cluster_x , 2 , quantile , probs = c(0.95) , na.rm = TRUE )
  }
  return(list(centroids, cent5, cent95))
}

# Define function to find the lncRNA cluster centroids, 5th and 95th percentiles of expression values
get_noncoding_centroids <- function(input_df){
  # Get number of clusters
  num_clusters <- length(unique(input_df$lncrna_cluster))
  # Restrict to just expression data and scale rows to have mean zero and standard deviation one
  expr_only <- input_df %>%
     dplyr::select(transcript) %>%
    unnest(cols = c(transcript))
  # Scale the data
  expr_only <- t(apply(expr_only, 1, scale))
  
  # define data frame to hold centroids
  centroids <- data.frame(matrix(nrow=num_clusters, ncol=41))
  # define data frame to hold 5th percentile
  cent5 <- data.frame(matrix(nrow=num_clusters, ncol=41))
  # define data frame to hold 95th percentile
  cent95 <- data.frame(matrix(nrow=num_clusters, ncol=41))
  for(n in 1:num_clusters){
    # Select one of the clusters
    cluster_x <- as.matrix(expr_only[input_df$lncrna_cluster==n,])
    centroid<-colMeans(cluster_x)
    # Save cluster centroid
    centroids[n,] <- centroid
    # Save 5th percentile of the cluster
    cent5[n,]<-apply(cluster_x , 2 , quantile , probs = c(0.05) , na.rm = TRUE )
    # Save 95th percentile of the cluster
    cent95[n,]<-apply(cluster_x , 2 , quantile , probs = c(0.95) , na.rm = TRUE )
  }
  return(list(centroids, cent5, cent95))
}

# Define function to plot centroids plus 95th and 5th percentiles
plot_centroid <- function(centroid_num, centroids_df){
  # define time
  t <- seq(0,400,10)
  # Get cumber of clusters
  num_clusters <- nrow(centroids_df[[1]])
  # Get colours for plotting:
  # light coloured-background representing 5% and 95% precentiles
  percentile_cols <- ggsci::pal_npg("nrc", alpha = 0.3)(num_clusters)
  # Dark coloured centroids
  centroid_cols <- ggsci::pal_npg("nrc", alpha = 1)(num_clusters)
  
  plot(t,
       unlist(centroids_df[[1]][centroid_num,]), 
       ylim=range(c(centroids_df[[1]][centroid_num,],
                    centroids_df[[2]][centroid_num,],
                    centroids_df[[3]][centroid_num,])),
       type="l",
       col=centroid_cols[centroid_num],
       lwd=2,
       xlab="", ylab="", xaxt="n", yaxt="n")
  polygon(c(t, rev(t)), c(centroids_df[[2]][centroid_num,], rev(centroids_df[[3]][centroid_num,])),
          col = percentile_cols[centroid_num], border = NA)
}

# Get centroid info
coding_centroids <- get_mRNA_centroids(filter(counts, transcript_type=="protein_coding"))
noncoding_centroids <- get_noncoding_centroids(filter(counts, transcript_type!="protein_coding"))

fig1d <- function(){
  # Set up a 7x7 grid
  par(mfcol=c(7,7))
  # Set thin plot margins
  par(mar=c(0.2,0.2,0.2,0.2)) # bottom, left, top and right margins
  
  # Create a blank plot for the top left corner
  plot(0,type='n',axes=FALSE,ann=FALSE)
  
  # Plot all of the non-coding centroids in the left-most column of the grid
  plot_centroid(centroid_num = 1,centroids_df = noncoding_centroids)
  plot_centroid(centroid_num = 2,centroids_df = noncoding_centroids)
  plot_centroid(centroid_num = 3,centroids_df = noncoding_centroids)
  plot_centroid(centroid_num = 4,centroids_df = noncoding_centroids)
  plot_centroid(centroid_num = 5,centroids_df = noncoding_centroids)
  plot_centroid(centroid_num = 6,centroids_df = noncoding_centroids)
  
  # Define colours for background of boxes with correlation values
  #bg_cols <- scales::alpha(colorRampPalette(brewer.pal(9,"RdBu"))(21), 0.7)
  bg_cols <- scales::alpha(colorRampPalette(brewer.pal(9,"Blues"))(11), 0.7)
  bg_cols <- c(rep(bg_cols[1],10), bg_cols)
  
  # Plot coding centroids
  for(i in 1:nrow(coding_centroids[[1]])){
    # Plot coding gene centroid in the top row of the grid
    plot_centroid(centroid_num = i, centroids_df = coding_centroids)
    
    # Plot coding/lncRNA centroid correlation value
    for(n in 1:nrow(noncoding_centroids[[1]])){
      # Determine the Pearson correlation between the coding and lncRNA centroids
      cor_value <-round(cor(unlist(coding_centroids[[1]][i,]), unlist(noncoding_centroids[[1]][n,]), method = "pearson"),2)
      # Determine the colour of the box
      box_col <-  bg_cols[which(cor_value < seq(-0.9,1,0.1))[1]]
      # plot
      plot(0,axes=TRUE, xlab="", ylab="", xaxt="n", yaxt="n", col="white")
      rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = box_col)
      legend(x="center", legend=paste0(cor_value,"     "), bty="n", cex=2.5)  
    }
  }
  par(mfcol=c(1,1))
}

# Plot figure 1d, saving out as svg - later to be assembled in multipanel pdf
print("plotting Fig. 1d")
pdf("figures/Fig_1/Fig1d.pdf", width=12, height=8)
fig1d()
dev.off()

# save a small png as well
png("figures_track/Fig1d.png", width=400, height=350)
fig1d()
dev.off()

# Clean up
rm(coding_centroids, noncoding_centroids, fig1d, get_mRNA_centroids, get_noncoding_centroids, plot_centroid)

print("Finished plotting Fig. 1")
