## Figure 2 - The effects of transcript stability
# Having generated all of the required data with 'get_data_for_Fig_2.R', we can create the plot now

# Load required R packages
library(tidyverse)
library(viridis)

# Check to see if directory for the current figure exists and if not, create it
if(!dir.exists("figures/Fig_2")){dir.create("figures/Fig_2")}

# Get data required for plotting
fig_2_data <- readRDS("data_track/fig_2_data.rds")

# Define time (in minutes)
t<-seq(0,400,10)

# Get some colours for plotting each gene
cols_lines <- viridis(nrow(fig_2_data)+1, alpha=1, begin=0, end=1)
cols_points <-viridis(nrow(fig_2_data)+1, alpha=0.7, begin=0, end=1)

# Define function to plot Fig. 2
fig2_plot <- function(){
  
  # Set up plot grid
  par(mfrow=c(4,4))
  # Set outer margins; bottom, left, top right
  par(oma=c(2,0,1,0))
  
  # Set individual plot margins; bottom, left, top right
  par(mar=c(0.8,4,2,1))
  
  ## FOS pre-mRNA
  plot(t, fig_2_data$last_10Kb[[1]], xlab="", ylab="", xaxt="n", yaxt="n", pch=16, col=cols_points[1],
       ylim=range(c(fig_2_data$last_10Kb[[1]], fig_2_data$impulse_model_fit[[1]])))
  lines(t, fig_2_data$impulse_model_fit[[1]], col=cols_lines[1], lwd=2)
  # x-axis
  axis(side = 1, at = c(0, 100, 200, 300, 400), labels = FALSE, tick = TRUE, lwd.ticks=0.5, tck=-0.03)
  # y-axis
  axis(side = 2, at = c(0,400,800,1200), labels = c("0","400","800","1,200"), tick = TRUE, lwd.ticks=0.5, tck=-0.03, las=1, mgp=c(3,0.35,0))
  # add legend
  legend("topright", legend=c("pre-mRNA", ""), lty=c(NA, 1), pch=c(16, NA), col=c(cols_points[1], cols_lines[1]), bty="n", seg.len=1)
  
  ## HES1 pre-mRNA
  plot(t, fig_2_data$last_10Kb[[2]], xlab="", ylab="", xaxt="n", yaxt="n", pch=16, col=cols_points[2],
       ylim=range(c(fig_2_data$last_10Kb[[2]], fig_2_data$impulse_model_fit[[2]])))
  lines(t, fig_2_data$impulse_model_fit[[2]], col=cols_lines[2], lwd=2)
  # x-axis
  axis(side = 1, at = c(0, 100, 200, 300, 400), labels = FALSE, tick = TRUE, lwd.ticks=0.5, tck=-0.03)
  # y-axis
  axis(side = 2, at = c(0,150,300,450), labels = TRUE, tick = TRUE, lwd.ticks=0.5, tck=-0.03, las=1, mgp=c(3,0.35,0))
  # add legend
  legend("topright", legend=c("pre-mRNA", ""), lty=c(NA, 1), pch=c(16, NA), col=c(cols_points[2], cols_lines[2]), bty="n", seg.len=1)
  
  ## FOSB pre-mRNA
  plot(t, fig_2_data$last_10Kb[[3]], xlab="", ylab="", xaxt="n", yaxt="n", pch=16, col=cols_points[3],
       ylim=range(c(fig_2_data$last_10Kb[[3]], fig_2_data$impulse_model_fit[[3]])))
  lines(t, fig_2_data$impulse_model_fit[[3]], col=cols_lines[3], lwd=2)
  # x-axis
  axis(side = 1, at = c(0, 100, 200, 300, 400), labels = FALSE, tick = TRUE, lwd.ticks=0.5, tck=-0.03)
  # y-axis
  axis(side = 2, at = c(0,800,1600,2400), labels = c("0","800","1,600","2,400"), tick = TRUE, lwd.ticks=0.5, tck=-0.03, las=1, mgp=c(3,0.35,0))
  # add legend
  legend("topright", legend=c("pre-mRNA", ""), lty=c(NA, 1), pch=c(16, NA), col=c(cols_points[3], cols_lines[3]), bty="n", seg.len=1)
  
  ## CTGF pre-mRNA
  plot(t, fig_2_data$last_10Kb[[4]], xlab="", ylab="", xaxt="n", yaxt="n", pch=16, col=cols_points[4],
       ylim=range(c(fig_2_data$last_10Kb[[4]], fig_2_data$impulse_model_fit[[4]])))
  lines(t, fig_2_data$impulse_model_fit[[4]], col=cols_lines[4], lwd=2)
  # x-axis
  axis(side = 1, at = c(0, 100, 200, 300, 400), labels = FALSE, tick = TRUE, lwd.ticks=0.5, tck=-0.03)
  # y-axis
  axis(side = 2, at = c(0,250,500,750), labels = c("0","250","500","750"), tick = TRUE, lwd.ticks=0.5, tck=-0.03, las=1, mgp=c(3,0.35,0))
  # add legend
  legend("topright", legend=c("pre-mRNA", ""), lty=c(NA, 1), pch=c(16, NA), col=c(cols_points[4], cols_lines[4]), bty="n", seg.len=1)
  
  # Set individual plot margins; bottom, left, top right
  par(mar=c(3,4,0,1))
  
  ## FOS mRNA
  plot(t, fig_2_data$mRNA_expression[[1]], xlab="", ylab="", xaxt="n", yaxt="n", pch=16, col=cols_points[1],
       ylim=range(c(fig_2_data$mRNA_expression[[1]], fig_2_data$transcription_model_fit[[1]])))
  lines(t, fig_2_data$transcription_model_fit[[1]], col=cols_lines[1], lwd=2)
  # x-axis
  axis(side = 1, at = c(0, 100, 200, 300, 400), labels = TRUE, tick = TRUE, lwd.ticks=0.5, tck=-0.03, mgp=c(3,0.35,0))
  # y-axis
  axis(side = 2, at = c(0,5000,10000,15000), labels = c("0","5,000","10,000","15,000"), tick = TRUE, lwd.ticks=0.5, tck=-0.03, las=1, mgp=c(3,0.35,0))
  # add legend
  legend("topright", legend=c("mRNA", "model fit"), lty=c(NA, 1), pch=c(16, NA), col=c(cols_points[1], cols_lines[1]), bty="n", seg.len=1)
  
  ## HES1 mRNA
  plot(t, fig_2_data$mRNA_expression[[2]], xlab="", ylab="", xaxt="n", yaxt="n", pch=16, col=cols_points[2],
       ylim=range(c(fig_2_data$mRNA_expression[[2]], fig_2_data$transcription_model_fit[[2]])))
  lines(t, fig_2_data$transcription_model_fit[[2]], col=cols_lines[2], lwd=2)
  # x-axis
  axis(side = 1, at = c(0, 100, 200, 300, 400), labels = TRUE, tick = TRUE, lwd.ticks=0.5, tck=-0.03, mgp=c(3,0.35,0))
  # y-axis
  axis(side = 2, at = c(0,1500,3000,4500), labels = c("0","1,500","3,000","4,500"), tick = TRUE, lwd.ticks=0.5, tck=-0.03, las=1, mgp=c(3,0.35,0))
  # add legend
  legend("topright", legend=c("mRNA", "model fit"), lty=c(NA, 1), pch=c(16, NA), col=c(cols_points[2], cols_lines[2]), bty="n", seg.len=1)
  
  ## FOSB mRNA
  plot(t, fig_2_data$mRNA_expression[[3]], xlab="", ylab="", xaxt="n", yaxt="n", pch=16, col=cols_points[3],
       ylim=range(c(fig_2_data$mRNA_expression[[3]], fig_2_data$transcription_model_fit[[3]])))
  lines(t, fig_2_data$transcription_model_fit[[3]], col=cols_lines[3], lwd=2)
  # x-axis
  axis(side = 1, at = c(0, 100, 200, 300, 400), labels = TRUE, tick = TRUE, lwd.ticks=0.5, tck=-0.03, mgp=c(3,0.35,0))
  # y-axis
  axis(side = 2, at = c(0,2000,4000,6000), labels = c("0","2,000","4,000","6,000"), tick = TRUE, lwd.ticks=0.5, tck=-0.03, las=1, mgp=c(3,0.35,0))
  # add legend
  legend("topright", legend=c("mRNA", "model fit"), lty=c(NA, 1), pch=c(16, NA), col=c(cols_points[3], cols_lines[3]), bty="n", seg.len=1)
  
  ## CTGF mRNA
  plot(t, fig_2_data$mRNA_expression[[4]], xlab="", ylab="", xaxt="n", yaxt="n", pch=16, col=cols_points[4],
       ylim=range(c(fig_2_data$mRNA_expression[[4]], fig_2_data$transcription_model_fit[[4]])))
  lines(t, fig_2_data$transcription_model_fit[[4]], col=cols_lines[4], lwd=2)
  # x-axis
  axis(side = 1, at = c(0, 100, 200, 300, 400), labels = TRUE, tick = TRUE, lwd.ticks=0.5, tck=-0.03, mgp=c(3,0.35,0))
  # y-axis
  axis(side = 2, at = c(0,10000,20000,30000), labels = c("0","10,000","20,000","30,000"), tick = TRUE, lwd.ticks=0.5, tck=-0.03, las=1, mgp=c(3,0.35,0))
  # add legend
  legend("bottom", legend=c("mRNA", "model fit"), lty=c(NA, 1), pch=c(16, NA), col=c(cols_points[4], cols_lines[4]), bty="n", seg.len=1)
  
  # Set individual plot margins; bottom, left, top right
  par(mar=c(0.8,4,2,1))
  
  ## SRF pre-mRNA
  plot(t, fig_2_data$last_10Kb[[5]], xlab="", ylab="", xaxt="n", yaxt="n", pch=16, col=cols_points[5],
       ylim=range(c(fig_2_data$last_10Kb[[5]], fig_2_data$impulse_model_fit[[5]])))
  lines(t, fig_2_data$impulse_model_fit[[5]], col=cols_lines[5], lwd=2)
  # x-axis
  axis(side = 1, at = c(0, 100, 200, 300, 400), labels = FALSE, tick = TRUE, lwd.ticks=0.5, tck=-0.03)
  # y-axis
  axis(side = 2, at = c(100,450,800,1150), labels = c("100","450","800","1,150"), tick = TRUE, lwd.ticks=0.5, tck=-0.03, las=1, mgp=c(3,0.35,0))
  # add legend
  legend("topright", legend=c("pre-mRNA", ""), lty=c(NA, 1), pch=c(16, NA), col=c(cols_points[5], cols_lines[5]), bty="n", seg.len=1)
  
  ## DSTN pre-mRNA
  plot(t, fig_2_data$last_10Kb[[6]], xlab="", ylab="", xaxt="n", yaxt="n", pch=16, col=cols_points[6],
       ylim=range(c(fig_2_data$last_10Kb[[6]], fig_2_data$impulse_model_fit[[6]])))
  lines(t, fig_2_data$impulse_model_fit[[6]], col=cols_lines[6], lwd=2)
  # x-axis
  axis(side = 1, at = c(0, 100, 200, 300, 400), labels = FALSE, tick = TRUE, lwd.ticks=0.5, tck=-0.03)
  # y-axis
  axis(side = 2, at = c(300,1100,1900,2700), labels = c("300","1,100","1,900","2,700"), tick = TRUE, lwd.ticks=0.5, tck=-0.03, las=1, mgp=c(3,0.35,0))
  # add legend
  legend("topright", legend=c("pre-mRNA", ""), lty=c(NA, 1), pch=c(16, NA), col=c(cols_points[6], cols_lines[6]), bty="n", seg.len=1)
  
  ## TPM4 pre-mRNA
  plot(t, fig_2_data$last_10Kb[[7]], xlab="", ylab="", xaxt="n", yaxt="n", pch=16, col=cols_points[7],
       ylim=range(c(fig_2_data$last_10Kb[[7]], fig_2_data$impulse_model_fit[[7]])))
  lines(t, fig_2_data$impulse_model_fit[[7]], col=cols_lines[7], lwd=2)
  # x-axis
  axis(side = 1, at = c(0, 100, 200, 300, 400), labels = FALSE, tick = TRUE, lwd.ticks=0.5, tck=-0.03)
  # y-axis
  axis(side = 2, at = c(100, 400, 700, 1000), labels = c("100", "400", "700", "1,000"), tick = TRUE, lwd.ticks=0.5, tck=-0.03, las=1, mgp=c(3,0.35,0))
  # add 150
  legend("topright", legend=c("pre-mRNA", ""), lty=c(NA, 1), pch=c(16, NA), col=c(cols_points[7], cols_lines[7]), bty="n", seg.len=1)
  
  # Add one blank plot
  plot(0,type='n',axes=FALSE,ann=FALSE)
  
  # Set individual plot margins; bottom, left, top right
  par(mar=c(3,4,0,1))
  
  ## SRF mRNA
  plot(t, fig_2_data$mRNA_expression[[5]], xlab="", ylab="", xaxt="n", yaxt="n", pch=16, col=cols_points[5],
       ylim=range(c(fig_2_data$mRNA_expression[[5]], fig_2_data$transcription_model_fit[[5]])))
  lines(t, fig_2_data$transcription_model_fit[[5]], col=cols_lines[5], lwd=2)
  # x-axis
  axis(side = 1, at = c(0, 100, 200, 300, 400), labels = TRUE, tick = TRUE, lwd.ticks=0.5, tck=-0.03, mgp=c(3,0.35,0))
  # y-axis
  axis(side = 2, at = c(500, 1700, 2900, 4100), labels = c("500", "1,700", "2,900", "4,100"), tick = TRUE, lwd.ticks=0.5, tck=-0.03, las=1, mgp=c(3,0.35,0))
  # add legend
  legend("bottom", legend=c("mRNA", "model fit"), lty=c(NA, 1), pch=c(16, NA), col=c(cols_points[5], cols_lines[5]), bty="n", seg.len=1)
  
  ## DSTN mRNA
  plot(t, fig_2_data$mRNA_expression[[6]], xlab="", ylab="", xaxt="n", yaxt="n", pch=16, col=cols_points[6],
       ylim=range(c(fig_2_data$mRNA_expression[[6]], fig_2_data$transcription_model_fit[[6]])))
  lines(t, fig_2_data$transcription_model_fit[[6]], col=cols_lines[6], lwd=2)
  # x-axis
  axis(side = 1, at = c(0, 100, 200, 300, 400), labels = TRUE, tick = TRUE, lwd.ticks=0.5, tck=-0.03, mgp=c(3,0.35,0))
  # y-axis
  axis(side = 2, at = c(15000,20000,25000,30000), labels = c("15,000","20,000","25,000","30,000"), tick = TRUE, lwd.ticks=0.5, tck=-0.03, las=1, mgp=c(3,0.35,0))
  # add legend
  legend("bottomright", legend=c("mRNA", "model fit"), lty=c(NA, 1), pch=c(16, NA), col=c(cols_points[6], cols_lines[6]), bty="n", seg.len=1)
  
  ## TPM4 mRNA
  plot(t, fig_2_data$mRNA_expression[[7]], xlab="", ylab="", xaxt="n", yaxt="n", pch=16, col=cols_points[7],
       ylim=range(c(fig_2_data$mRNA_expression[[7]], fig_2_data$transcription_model_fit[[7]])))
  lines(t, fig_2_data$transcription_model_fit[[7]], col=cols_lines[7], lwd=2)
  # x-axis
  axis(side = 1, at = c(0, 100, 200, 300, 400), labels = TRUE, tick = TRUE, lwd.ticks=0.5, tck=-0.03, mgp=c(3,0.35,0))
  # y-axis
  axis(side = 2, at = c(13000, 21000, 29000, 37000), labels = c("13,000", "21,000", "29,000", "37,000"), tick = TRUE, lwd.ticks=0.5, tck=-0.03, las=1, mgp=c(3,0.35,0))
  # add legend
  legend("bottomright", legend=c("mRNA", "model fit"), lty=c(NA, 1), pch=c(16, NA), col=c(cols_points[7], cols_lines[7]), bty="n", seg.len=1)
  par(mfcol=c(1,1))
}

# Save out
pdf("figures/Fig_2/Fig_2_line_plots.pdf", width=10, height=6)
fig2_plot()
dev.off()

# save a small png as well
png("figures_track/Fig_2_line_plots.png", width=600, height=350)
fig2_plot()
dev.off()

# We also want to include a barplot with the half-life of each transcript, as inferred using our simple model of transcription.

# Rather than plotting the degradation rates, which are perhaps less intuitive, we will plot the mRNA half lives
fig_2_data <- mutate(fig_2_data, mRNA_half_life = log(2)/alpha)

# Define function to plot Fig. 2 barplot
fig2_barplot <- function(){
  
  # Set individual plot margins; bottom, left, top right
  par(mar=c(2.5,2,1,0))
  
  x<-barplot(fig_2_data$mRNA_half_life, ylim=c(0,200), col=cols_lines, names.arg = "", yaxt="n")
  
  # y-axis
  axis(side = 2, at = c(0, 60, 120, 180),
       labels = c("0", "1", "2", "3"),
       tick = TRUE, lwd.ticks=0.5, tck=-0.03, las=1, mgp=c(3,0.3,0), cex.axis=0.7)  
  
  # x-axis
  text(cex=0.7, x=x, y=-0.002, as.character(fig_2_data$gene_name), xpd=TRUE, srt=45, adj=1)
}

# Save out
pdf("figures/Fig_2/Fig_2_barplot.pdf", width=2.1, height=3)
fig2_barplot()
dev.off()

print("Finished plotting Fig. 2")
