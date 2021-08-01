### Create panels for Fig. 5 : three examples of human protein-coding gene and
### adjacent lncRNA expression; FOS, TGFBI & TGIF1

# Load R packages
library(tidyverse)
library(RColorBrewer)


### Gene expression line plots

# Define function to plot gene and lncRNA line plots
RNAseq_lineplots <- function(input_data, gene_name, xlimits=c(0, 400), ylimits){
  pdf(str_glue('figures/Fig_5/{gene_name}_line_plot.pdf'), width=8, height=4)
  # Create empty plot
  plot("", xlim=xlimits, ylim=ylimits, xlab="Time (min)", type="n", ylab="Scaled normalized counts", main=gene_name, bty="l", xaxt="n")
  # Add x-axis
  axis(side=1, at=seq(0,400,50), labels = c("0","50","100","150","200","250","300","350","400"))
  axis(side=1, at=seq(0,400,10), labels = NA)
  # Add pale green shaded area that contains the range of all lncRNA values
  polygon(c(seq(0,400,10), rev(seq(0,400,10))), c(sapply(input_data[[3]], max), rev(sapply(input_data[[3]], min))),
          col="#bae4b3", border = NA)
  # Add line for the lncRNA mean - dark green line
  lines(seq(0,400,10), colMeans(input_data[[3]]), col="#238b45", lwd=2)
  # exons - black line
  lines(seq(0,400,10), scale(input_data[[1]]), xlab="Time (min)", ylab="Scaled normalized counts", main=gene_name, type="l", ylim=ylimits, lwd=2)
  # introns - red line
  lines(seq(0,400,10), scale(input_data[[2]]), col="#d73027", lwd=2);
  
  # Add legend
  legend("topright", legend=c("spliced mRNA", "unspliced precursor", "lncRNA"), bty="n", lty=1, lwd=2, col=c("black", "#d73027", "#238b45"), fill = c(NA, NA, "#bae4b3"), border=NA, merge = FALSE)
  
  dev.off()
}

# Save out plots if not already done
if(!file.exists("figures/Fig_5/FOS_line_plot.pdf")){
  RNAseq_lineplots(input_data=readRDS("data/FOS_interval_counts.rds"), gene_name="FOS", xlimits = c(0,200), ylimits=c(-0.5,4.5))
  RNAseq_lineplots(input_data=readRDS("data/TGFBI_interval_counts.rds"), gene_name="TGFBI", ylimits=c(-2.5,3))
  RNAseq_lineplots(input_data=readRDS("data/TGIF1_interval_counts.rds"), gene_name="TGIF1", ylimits=c(-2,4.5))
}
# Clean up
rm(RNAseq_lineplots)



### Plot a schematic of the genes/lncRNAs

# Define function to plot; a scale for the plot, genes and lncRNAs
plot_gene_ncRNA <- function(input_df, gene_name, region_start, region_end, tick_distance = 1E4){
  
  ## Plot genomic scale with 10Kb intervals
  # Create vector of where to put tick marks. We want the ticks at tick_distance intervals
  axis_tickmarks <- seq(region_start, region_end, tick_distance)
  # Remove end tick
  axis_tickmarks <- axis_tickmarks[1:length(axis_tickmarks)-1] # remove last
  # Remove first tick
  axis_tickmarks <- axis_tickmarks[-1]
  #plot
  pdf(file=str_glue("figures/Fig_5/{gene_name}_schematic.pdf"), width = 8, height = 3)
  
  par(lend=2, ljoin=1, mar=c(0,0,3,0))
  plot("","", xlim=c(region_start, region_end),
       type="n", ylim=c(-2,3), xlab="", ylab="", xaxs="i", bty="n", yaxt="n", xaxt="n", main=str_glue('{gene_name} {tick_distance/1E3} Kb genomic intervals'))
  axis(side=3, at = c(region_start, region_end),lwd.ticks=0, labels = NA)
  axis(side=3, at = c(axis_tickmarks), col=NA,col.ticks = 1, labels = NA)
  
  ## Plot gene
  # add black boxes for introns
  introns <- filter(input_df, feature=="intron")
  for(i in 1:nrow(introns)){
    polygon(x=c(introns$start[i], introns$start[i], introns$end[i], introns$end[i]), y=c(2+(6.5/14), 2+(7.5/14), 2+(7.5/14), 2+(6.5/14)), col="black", border = "black", lwd=1)
  }
  # add black boxes for exons
  exons <- filter(input_df, feature=="exon")
  for(i in 1:nrow(exons)){
    polygon(x=c(exons$start[i], exons$start[i], exons$end[i], exons$end[i]), y=c(2, 3, 3, 2), col="black", border = "black", lwd=1)
  }
  # add black boxes for utrs
  utrs <- filter(input_df, feature=="five_UTR" | feature=="three_UTR")
  for(i in 1:nrow(utrs)){
    polygon(x=c(utrs$start[i], utrs$start[i], utrs$end[i], utrs$end[i]), y=c(2+(3/14), 2+(11/14), 2+(11/14), 2+(3/14)), col="black", border = "black", lwd=1)
  }
  
  # Finish off plot
  dev.off()
}

# Save out plots
plot_gene_ncRNA(input_df=readRDS("data/human_schematic.rds")[[1]], gene_name="FOS", tick_distance = 5E3, region_start=75256000, region_end=75302000)
plot_gene_ncRNA(input_df=readRDS("data/human_schematic.rds")[[2]], gene_name="TGFBI", tick_distance = 1E4,region_start=135989000, region_end=136080000)
plot_gene_ncRNA(input_df=readRDS("data/human_schematic.rds")[[3]], gene_name="TGIF1", tick_distance = 2E4, region_start=3435000, region_end=3671000)

# Clean up
rm(plot_gene_ncRNA)



### Plot human epigenomics data 
roadmap.bw <- readRDS("data/roadmap_bw.rds")

# Define function to plot epigenome line plots
plot_epi_lineplots <- function(input_df, target_gene, region_start, region_end, bp_res, smooth_width){
  pdf(file=str_glue("figures/Fig_5/{target_gene}_line_plots.pdf"), width = 8, height = 3)
  par(oma=c(2,0,1,1))
  par(mar=c(1,3,0,0))
  par(mfrow=c(3,1))
  for(i in 1:3){
    target_data <- c("DNase" ,"H3K4me1","H3K4me3")[i]
    colour_pallette <- c("Blues","Oranges", "Reds")[i]
    
    print(str_glue("Plotting {target_data} line plot for {target_gene} - epi mark {i} of 3"))
    
    # Get target data as a numeric matrix
    region_data <- input_df %>% 
      filter(target==target_data) %>%
      dplyr::select(target_gene)
    region_data <- do.call(rbind, region_data[,target_gene])
    
    # Reduce to n bp resolution
    region_data.nbp <- t(aggregate(t(region_data),list(rep(1:(ncol(region_data)%/%bp_res+1),each=bp_res,len=ncol(region_data))),mean)[-1])
    
    # Get quantiles
    region_quantiles <- apply(region_data.nbp , 2 , quantile , probs = seq(0,1,0.05))
    
    # Get mean and sd
    region_mean <- colMeans(region_data.nbp)
    region_sd <- apply(region_data.nbp, 2, sd)
    
    # Smooth quantiles, mean and sd
    region_quantiles.smoothed <- lapply(1:nrow(region_quantiles), function(i) zoo::rollmean(x = region_quantiles[i,], k=smooth_width)) %>% do.call(rbind,.)
    row.names(region_quantiles.smoothed) <- row.names(region_quantiles)
    region_mean.smoothed <- zoo::rollmean(x = region_mean, k=50)
    region_sd.smoothed <- zoo::rollmean(x = region_sd, k=50)
    
    # Get plot colours
    plot_colours <- RColorBrewer::brewer.pal(9, colour_pallette)
    
    # Get x-axis values
    x_vals <- seq(from = region_start, to = region_end, length.out = ncol(region_quantiles.smoothed))
    
    ## Plot
    plot_range <- range(c(region_quantiles.smoothed["10%",], region_quantiles.smoothed["90%",]))
    if(i!=3){
      plot(x = x_vals, y = region_quantiles.smoothed["50%",], type="l", ylim = plot_range, col = plot_colours[9], lwd=2, xlab="",ylab="", xaxt='n')
    } else {
      plot(x = x_vals, y = region_quantiles.smoothed["50%",], type="l", ylim = plot_range, col = plot_colours[9], lwd=2, xlab="",ylab="")
    }
      # Shaded region from 90% to 10%
    polygon(c(x_vals, rev(x_vals)),
            c(region_quantiles.smoothed["10%",], rev(region_quantiles.smoothed["90%",])),col=plot_colours[3], border = NA)
    # Shaded region from 25% to 75%
    polygon(c(x_vals, rev(x_vals)),
            c(region_quantiles.smoothed["25%",], rev(region_quantiles.smoothed["75%",])),col=plot_colours[6], border = NA)
    # Add lines around shaded region    
    lines(x_vals, region_quantiles.smoothed["10%",], col = plot_colours[5], lwd=0.5)
    lines(x_vals, region_quantiles.smoothed["90%",], col = plot_colours[5], lwd=0.5)
    # Re-plot median line
    lines(x_vals, region_quantiles.smoothed["50%",], col = plot_colours[9], lwd=2)
    legend("topright",legend = str_glue('{target_gene} - {target_data}'), bty="n")
  }
  dev.off()
  print(str_glue("Finished plotting epigenomics heatmaps for sample {target_gene}"))
}

plot_epi_lineplots(input_df=roadmap.bw, target_gene="FOS", region_start=75256000, region_end=75302000, bp_res = 10, smooth_width = 50)
plot_epi_lineplots(input_df=roadmap.bw, target_gene="TGFBI", region_start=135989000, region_end=136080000, bp_res = 10, smooth_width = 100)
plot_epi_lineplots(input_df=roadmap.bw, target_gene="TGIF1", region_start=3435000, region_end=3671000, bp_res = 20, smooth_width = 100)

print("Finished creating panels for Fig_5")
