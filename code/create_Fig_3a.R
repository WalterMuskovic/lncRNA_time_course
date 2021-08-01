## The effects of gene length on transcription time (Figure 3)
# Now that the data for CACNA1C has been obtained (create_data_for_Fig_3a.R), we call the data back in, process it and
# plot the intron coverage over the gene body for each time point.

# Load R packages
library(tidyverse)
library(RColorBrewer)

# Check to see if directory for the current figure exists and if not, create it
if(!dir.exists("figures/Fig_3")){dir.create("figures/Fig_3")}

# Get the CACNA1C data
CACNA1C <- readRDS("data/CACNA1C.rds")

# Smooth the data to get rid of jagged peaks
CACNA1C_smooth <- data.frame(matrix(ncol=41, nrow=length(zoo::rollmean(CACNA1C[,1], k=20))))
for(i in 1:ncol(CACNA1C_smooth)){ CACNA1C_smooth[, i] <- zoo::rollmean(CACNA1C[, i], k=20) }

# Create a data frame to hold the colours for each polygon
polygon_cols <- data.frame(matrix(nrow=nrow(CACNA1C_smooth)-1, ncol=ncol(CACNA1C_smooth)))

# Define function that will take a numeric vector, scale it to be between zero and one, and return a vector of colours
get_cols <- function(input_vector){
  scaled_data <- input_vector
  scaled_data<-scaled_data-min(scaled_data)
  scaled_data<-scaled_data/max(scaled_data)
  return(rgb(colorRamp(RColorBrewer::brewer.pal(9,"Blues"))(scaled_data), maxColorValue = 256))
}

# For each define a colour to use: dark colours correspond to the time when expression is maximal, light colours when expression is minimal
for (i in 1:nrow(polygon_cols)){
  # Get data
  input_df <- data.frame(time=seq(0,400,10), expression=colMeans(CACNA1C_smooth[i:(i+1),]))
  # Fit loess
  lw1 <- loess(expression ~ time, data=input_df, span = 0.5)
  # Get colours
  polygon_cols[i,] <- get_cols(lw1$fitted)
}

# Scale the data for plotting
CACNA1C_smooth <- CACNA1C_smooth - min(CACNA1C_smooth)
CACNA1C_smooth <- CACNA1C_smooth/max(CACNA1C_smooth)
CACNA1C_smooth <- CACNA1C_smooth*50
for(i in 1:ncol(CACNA1C_smooth)){ CACNA1C_smooth[,i] <-  CACNA1C_smooth[,i]+seq(0,400,10)[i] }

# Define function to plot CACNA1C intron expression across the gene body
CACNA1C_fig <- function(){
  # Create a blank plot
  par(mar=c(1,1,1,1))
  plot("","",ylim=c(10,410), xlim=c(0,nrow(CACNA1C_smooth)), xlab="", type="n", ylab="", main="",xaxs="i", bty="n", yaxt="n", xaxt="n")
  # add y-axis ticks
  axis(side = 2, at=c(0,100,200,300,400), labels = FALSE, tick = TRUE, lwd=0, lwd.ticks=0.5, tck=-0.0075)
  # add y-axis ticks
  axis(side = 1, at=c(0,100,200,300,400, 500, 600), labels = FALSE, tick = TRUE, lwd=0, lwd.ticks=0.5, tck=-0.0075)
  # Draw a box around the plot
  box(lwd=0.5)
  
  # The following will draw a large number of polygons: one for each genomic interval and time point.
  # We could just draw one polygon for each time point, but this method allows us to adjust the colour
  # for each genomic interval and time point
  
  # Cycle through time points
  for(tp in 1:41){
    # cycle through genomic bins
    for(genomic_bin in 1:nrow(polygon_cols)){ 
      # Get the current time
      current_time <- seq(0,400,10)[tp]
      # Get the y-axis values for the top verices of the polygon
      top_left <- CACNA1C_smooth[genomic_bin,tp]
      top_right <- CACNA1C_smooth[genomic_bin+1,tp]
      # Get the x-axis values for the polygon
      left_val <- genomic_bin
      right_val <- genomic_bin+1
      
      # Add the shaded area.
      polygon(x=c(left_val, left_val, right_val, right_val), y=c(current_time, top_left, top_right,current_time), col=polygon_cols[genomic_bin,tp], border = NA)
    }
  }
  # Add some black lines to the ridges
  for(tp in 1:41){
    # Add a line
    lines(CACNA1C_smooth[,tp], lwd=0.5)
  }
}

pdf("figures/Fig_3/CACNA1C_ridge_plot.pdf", width=4, height=2)
CACNA1C_fig()
dev.off()

# Plot a small png as well
png("figures_track/CACNA1C_ridge_plot.png", width=400, height=200)
CACNA1C_fig()
dev.off()

# Function to plot color bar for the peak colours
color.bar <- function(lut, min=0, max=1, nticks=2, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(4, ticks, las=1, labels=FALSE)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

# Save out the colour bar to be added to the figure later
pdf("figures/Fig_3/CACNA1C_color_scale.pdf", width=0.8, height=1)
par(mar=c(0,0,0,3))
color.bar(scales::alpha(colorRampPalette(RColorBrewer::brewer.pal(9,"Blues"))(21), 0.7), -1)
dev.off()

print("Finished plotting Fig_3a")
