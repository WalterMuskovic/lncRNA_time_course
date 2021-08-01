# Summarise the read depth for the T98G time point samples

# Load R packages
library(tidyverse)
library(RColorBrewer)
library(glue)

# Import T98G RNA-seq read depth data 
read_depth <- read_csv("data_track/T98G_counts.txt", col_names = FALSE) %>% unlist()
read_depth <- read_depth/4E6
# Get plot colours (red=lower read depth, blue=higher read depth)
rbPal <- colorRampPalette(brewer.pal(5,"RdBu"))
plot_cols <- rbPal(10)[as.numeric(cut(read_depth, breaks = 10))]

png("figures_track/T98G_read_depth.png", width=800, height=400)
# Plot
plot(seq(0,400,10), read_depth, xlab="Time point sample (min)", ylab="Read depth (millions of reads)", pch=21, col="black", bg=plot_cols, cex=2.5, main= "T98G Time Point Sample Read Depth", lwd=2, ylim=c(10,60))
# Add line indicating the mean read depth
abline(h=mean(read_depth), col="grey", lty=2, lwd=2)
# Label max and min points
# Max
text(x=seq(0,400,10)[which.max(read_depth)], y=read_depth[which.max(read_depth)],labels=glue('max = {round(max(read_depth),1)} milion reads'), pos=4, offset=0.65, cex=0.8)
# Min
text(x=seq(0,400,10)[which.min(read_depth)], y=read_depth[which.min(read_depth)],labels=glue('min = {round(min(read_depth),1)} milion reads'), pos=4, offset=0.65, cex=0.8)
# Also add label for the mean read depth
text(x=20, y=mean(read_depth),labels=glue('mean = \n{round(mean(read_depth),1)} milion\nreads'), pos=3, offset=0.65, cex=0.8)
dev.off()

# Clean up
rm(plot_cols, read_depth, rbPal)

print("Completed figure illustrating read depth for the T98G time series samples")