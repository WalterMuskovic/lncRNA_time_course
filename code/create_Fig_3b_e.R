### pre-mRNA expression dynamics for genes of different length
# We now have the data we need (obtained with create_data_for_Fig_3b_e.R) to create our plots for individual
# genes showing the effects of gene length on transcription time.

# Load R Packages
library(scales)
library(tidyverse)

# Get the gene plot data
fig3 <- readRDS("data/Fig_3_gene_plot_data.rds")

# Define function to plot gene
length_plot <- function(raw_data, fit_data, legend_x_coord, legend_y_coord, arrow_height){
  # Plot mRNA
  par(mar=c(1,1,0,1))
  plot(seq(0,400,10), raw_data$mRNA, col=scales::alpha("#bababa", 0.8), pch=16, xlab="", ylab="", ylim=c(-0.02,1.1), xlim=c(-5,405), axes = FALSE, yaxs = "i", xaxs = "i")
  lines(seq(0,400,1), fit_data$mRNA, col=scales::alpha("#1a1a1a", 0.8), lwd=2)
  # Add lines and points for the first 10 Kb of intron 
  points(seq(0,400,10), raw_data$first_10Kb, col=scales::alpha("#abd9e9", 0.8), pch=16)
  lines(seq(0,400,1), fit_data$first_10Kb, col=scales::alpha("#4575b4", 0.8), lwd=2)
  # Add lines and points for the last 10 Kb of intron 
  points(seq(0,400,10), raw_data$last_10Kb, col=scales::alpha("#f46d43", 0.8), pch=16)
  lines(seq(0,400,1), fit_data$last_10Kb, col=scales::alpha("#d73027", 0.8), lwd=2)
  # Add axes
  # x-axis
  axis(side = 1, at=c(0,100,200,300,400), labels = FALSE, tck=-0.0075)
  # y-axis
  axis(side = 2, at=c(0,0.2,0.4,0.6,0.8,1), labels = FALSE, tck=-0.0075)
  # arrow between first_10Kb and last_10Kb peak
  arrows(x0 = seq(0,400,1)[which.max(fit_data$first_10Kb)], y0 = arrow_height, # coordinates of points FROM which to draw.
         x1 = seq(0,400,1)[which.max(fit_data$last_10Kb)], y1 = arrow_height, #coordinates of points TO which to draw.
         code=3, # head is drawn at both ends of the arrow
         length = 0.05) # length of the edges of the arrow head (in inches).
  # Add legend
  legend(x = legend_x_coord, y = legend_y_coord, legend=c("","",""), lty=c(1,1,1), lwd=c(2,2,2), col=c(scales::alpha("#4575b4", 0.8),scales::alpha("#d73027", 0.8), scales::alpha("#1a1a1a", 0.8)),bty="n", seg.len=1)
}

# Create plots
# CACNA1C
pdf("figures/Fig_3/CACNA1C_line_plot.pdf", width=4, height=2.5)
length_plot(raw_data = fig3[["CACNA1C_raw"]], fit_data = fig3[["CACNA1C_fit"]], legend_x_coord = 110, legend_y_coord = 0.95, arrow_height = 1.05)
dev.off()
# Include a small png
png("figures_track/CACNA1C_line_plot.png", width=400, height=250)
length_plot(raw_data = fig3[["CACNA1C_raw"]], fit_data = fig3[["CACNA1C_fit"]], legend_x_coord = 110, legend_y_coord = 0.95, arrow_height = 1.05)
dev.off()

# LDLRAD4
pdf("figures/Fig_3/LDLRAD4_line_plot.pdf", width=4, height=2.5)
length_plot(raw_data = fig3[["LDLRAD4_raw"]], fit_data = fig3[["LDLRAD4_fit"]], legend_x_coord = 280, legend_y_coord = 0.5, arrow_height = 0.95)
dev.off()

# VCL
pdf("figures/Fig_3/VCL_line_plot.pdf", width=4, height=2.5)
length_plot(raw_data = fig3[["VCL_raw"]], fit_data = fig3[["VCL_fit"]], legend_x_coord = 200, legend_y_coord = 0.75, arrow_height = 1.05)
dev.off()

# LIMA1
pdf("figures/Fig_3/LIMA1_line_plot.pdf", width=4, height=2.5)
length_plot(raw_data = fig3[["LIMA1_raw"]], fit_data = fig3[["LIMA1_fit"]], legend_x_coord = 200, legend_y_coord = 0.6, arrow_height = 1.05)
dev.off()

# What is the delay between the first 10Kb and last 10Kb peaks?
# CACNA1C
paste0(seq(0,400,1)[which.max(fig3[["CACNA1C_fit"]]$last_10Kb)] - seq(0,400,1)[which.max(fig3[["CACNA1C_fit"]]$first_10Kb)], " minutes")
#[1] "221 minutes"

# LDLRAD4
paste0(seq(0,400,1)[which.max(fig3[["LDLRAD4_fit"]]$last_10Kb)] - seq(0,400,1)[which.max(fig3[["LDLRAD4_fit"]]$first_10Kb)], " minutes")
#[1] "146 minutes"

# VCL
paste0(seq(0,400,1)[which.max(fig3[["VCL_fit"]]$last_10Kb)] - seq(0,400,1)[which.max(fig3[["VCL_fit"]]$first_10Kb)], " minutes")
#[1] "39 minutes"

# LIMA1
paste0(seq(0,400,1)[which.max(fig3[["LIMA1_fit"]]$last_10Kb)] - seq(0,400,1)[which.max(fig3[["LIMA1_fit"]]$first_10Kb)], " minutes")
#[1] "16 minutes"

print("Finished plotting Fig_3b-e")
