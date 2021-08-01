### Block bootstrap - schematic

#library(tidyverse)
library(RColorBrewer)

# Get 5 colours for plotting using the RColorBrewer package
cols <-brewer.pal(5, name="Dark2")

# Define function to make clusters
makeClusters <- function(start, end, num.features) {
  sort(sample(end - start, num.features) + start)
}

# Define parameters for making clusters
total.length <- 1000
cluster_size <- 200
num.features <- 5

# Vectors to hold pch and colours
# Colour is a stand in for an archetypal expression pattern. Two points that are green have a similar expression pattern.
col_vec <- rep(NA, total.length)
# The two point shapes, triangle and circle, are a stand-in for transcript type (coding or non-coding)
pch_vec <- rep(NA, total.length)

# Set the seed for reproducibility
set.seed(12345)
# Generate 5 loose clusters
for(i in 1:5){
  cluster_start <- sample(1:c(total.length-cluster_size),1)
  coords <- makeClusters(start = cluster_start, end = (cluster_start+cluster_size), num.features = num.features)
  # Save colour and pch
  col_vec[coords] <- cols[i]
  pch_vec[coords] <- sample(c(21,24), num.features, replace = TRUE)
}

# Plot the results
png(filename = "figures_track/block_bootstrap_schematic.png", width = 8, height = 6, units = "in", res=300)

# Set the grid layout
layout(matrix(c(1,2,2,3,4,4,5,5,6,7,7,8), 
              3, # 3 rows
              4, # 2 columns
              byrow = TRUE),
       widths=c(1,1,1),
       heights=c(1,2,1))

# Set margins
par(mar=c(3,4,2,2)) # bottom, left, top, right

# Blank plot
plot(0,type='n',axes=FALSE,ann=FALSE)
legend(x = "left", legend = c("Transcript type","coding", "non-coding"), bty="n", pch=c(NA,21,24), pt.bg=c(NA,"grey", "grey"), col=c("black"), pt.cex=2, cex=0.9)

# First plot shows the original data
plot(which(complete.cases(col_vec)),
     rep(1, sum(complete.cases(col_vec))),
     bg=col_vec[complete.cases(col_vec)],
     pch=pch_vec[complete.cases(col_vec)],
     cex=2, col="black", ylab="", xlab="", yaxt="n", main="",
     ylim=c(0.8, 1.2), xlim=c(0,1000))
axis(2,at=1,label=c("data"),las=1)
# Blank plot
plot(0,type='n',axes=FALSE,ann=FALSE)
legend(x = "left", legend = c("Expression pattern", "1","2","3","4","5"), bty="n", pch=c(NA,22,22,22,22,22), pt.bg=c(NA, cols), col=c("black"), pt.cex=2, cex=0.9)

# Second plot (circles only)
# Define function for getting blocks for the bootstrap
blockBoot <- function(f, length.total, length.block) {
  num.block <- total.length / length.block
  block.starts <- sample(length.total - length.block, num.block, replace=TRUE)
  boot.idx <- do.call(c,lapply(block.starts, function(x) x:(x+length.block - 1)))
  f[boot.idx]
}
length.block <- 200
num.block <- total.length / length.block
# Plot circles
plot(which(pch_vec==21),
     rep(2, sum(pch_vec==21, na.rm = TRUE)),
     bg=col_vec[which(pch_vec==21)],
     pch=pch_vec[which(pch_vec==21)],
     cex=2, col="black", ylab="", xlab="" ,yaxt="n", main="",
     ylim=c(0.8, 2.2), xlim=c(0,1000))
axis(2,at=2:1,label=c("data", "boot"),las=1)
# Add circles for block bootstrap
# Get blocks
set.seed(1234)
block.starts <- sample(total.length - length.block, num.block, replace=TRUE)
boot.idx <- do.call(c,lapply(block.starts, function(x) x:(x+length.block - 1)))
# Plot points
points(which(pch_vec[boot.idx]==21),
       rep(1, sum(pch_vec[boot.idx]==21, na.rm = TRUE)),
       bg=col_vec[boot.idx][which(pch_vec[boot.idx]==21)],
       pch=pch_vec[boot.idx][which(pch_vec[boot.idx]==21)],
       cex=2, col="black")
# Add polygons
for (i in 1:num.block) {
  polygon(c(i*length.block,block.starts[i] + length.block,block.starts[i],(i-1)*length.block),c(1,2,2,1),col=rgb(0,0,0,.3),border=NA)
}

# Third plot (triangles only)
plot(which(pch_vec==24),
     rep(2, sum(pch_vec==24, na.rm = TRUE)),
     bg=col_vec[which(pch_vec==24)],
     pch=pch_vec[which(pch_vec==24)],
     cex=2, col="black", ylab="", xlab="", yaxt="n", main="",
     ylim=c(0.8, 2.2), xlim=c(0,1000))
axis(2,at=2:1,label=c("data", "boot"),las=1)
# Add triangles for block bootstrap
# Get blocks
set.seed(123)
block.starts <- sample(total.length - length.block, num.block, replace=TRUE)
boot.idx <- do.call(c,lapply(block.starts, function(x) x:(x+length.block - 1)))
# Plot triangles
points(which(pch_vec[boot.idx]==24),
       rep(1, sum(pch_vec[boot.idx]==24, na.rm = TRUE)),
       bg=col_vec[boot.idx][which(pch_vec[boot.idx]==24)],
       pch=pch_vec[boot.idx][which(pch_vec[boot.idx]==24)],
       cex=2, col="black")
# Add polygons
for (i in 1:num.block) {
  polygon(c(i*length.block,block.starts[i] + length.block,block.starts[i],(i-1)*length.block),c(1,2,2,1),col=rgb(0,0,0,.3),border=NA)
}

# Fourth plot (rearranged circles and triangles together)
# Blank plot
plot(0,type='n',axes=FALSE,ann=FALSE)
# Plot triangles
plot(which(pch_vec[boot.idx]==24),
     rep(1, sum(pch_vec[boot.idx]==24, na.rm = TRUE)),
     bg=col_vec[boot.idx][which(pch_vec[boot.idx]==24)],
     pch=pch_vec[boot.idx][which(pch_vec[boot.idx]==24)],
     cex=2, col="black",
     ylab="", xlab="", yaxt="n",main="",
     ylim=c(0.8, 1.2), xlim=c(0,1000))
axis(2,at=1,label=c("combined boot"),las=1)
# Get rearranged circles again and add to the plot
# Get blocks
set.seed(1234)
block.starts <- sample(total.length - length.block, num.block, replace=TRUE)
boot.idx <- do.call(c,lapply(block.starts, function(x) x:(x+length.block - 1)))
# Plot points
points(which(pch_vec[boot.idx]==21),
       rep(1, sum(pch_vec[boot.idx]==21, na.rm = TRUE)),
       bg=col_vec[boot.idx][which(pch_vec[boot.idx]==21)],
       pch=pch_vec[boot.idx][which(pch_vec[boot.idx]==21)],
       cex=2, col="black")
# Blank plot
plot(0,type='n',axes=FALSE,ann=FALSE)

dev.off()

# Clean up
rm(cols, makeClusters, total.length, cluster_size, num.features, col_vec, pch_vec, cluster_start, coords, blockBoot, length.block, num.block, block.starts, boot.idx, i)

