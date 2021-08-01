### Gene schematics
# It would also be useful to have some small schematic of the genes included in the plot.
# We code this in R so that everything is in the correct proportion. 
# As well as CACNA1c, ee'll also create schematics for; LDLRAD4, VCL and LIMA1.

# Load R packages
library(GenomicFeatures)
library(GenomicRanges)
library(rtracklayer)

# Begin by loading GTF data to get intron and exon intervals
TxDb <- makeTxDbFromGFF("data/annotation/GRCh38_spiked.gtf", format="gtf")
exons <- cdsBy(TxDb, use.names=TRUE)
introns <- intronsByTranscript(TxDb, use.names=TRUE)
fiveUTRs <- fiveUTRsByTranscript(TxDb, use.names=TRUE)
threeUTRs <- threeUTRsByTranscript(TxDb, use.names=TRUE)

# Define function to plot the gene schematic
plot_gene <- function(tx_id){
  # Get feature coordinates for current gene
  current_exons <- exons[[tx_id]]
  current_introns <- introns[[tx_id]]
  current_fiveUTRs <- fiveUTRs[[tx_id]]
  current_threeUTRs <- threeUTRs[[tx_id]]
  
  par(lend=2, mar=c(0,0,0,0))
  # Create empty plot to add polygons too
  plot("","", xlim=c(min(start(c(current_exons, current_introns, current_fiveUTRs, current_threeUTRs))),
                     max(end(c(current_exons, current_introns, current_fiveUTRs, current_threeUTRs)))),
       type="n", ylim=c(0.4,0.7), xlab="", ylab="", xaxt="none", yaxt="none", bty="n")
  
  # Add boxes for introns
  for(i in 1:length(current_introns)){
    polygon(x=c(start(current_introns)[i], start(current_introns)[i], end(current_introns)[i], end(current_introns)[i]), y=c(0.495, 0.505, 0.505, 0.495), col="black", border = "black", lwd=1)
  }
  
  # Add boxes for 5'UTR
  for(i in 1:length(current_fiveUTRs)){
    polygon(x=c(start(current_fiveUTRs)[i], start(current_fiveUTRs)[i], end(current_fiveUTRs)[i], end(current_fiveUTRs)[i]), y=c(0.45, 0.55, 0.55, 0.45), col="black", border = "black", lwd=1)
  }
  
  # Add boxes for 3'UTR
  for(i in 1:length(current_threeUTRs)){
    polygon(x=c(start(current_threeUTRs)[i], start(current_threeUTRs)[i], end(current_threeUTRs)[i], end(current_threeUTRs)[i]), y=c(0.45, 0.55, 0.55, 0.45), col="black", border = "black", lwd=1)
  }
  
  # Add boxes for exons
  for(i in 1:length(current_exons)){
    polygon(x=c(start(current_exons)[i], start(current_exons)[i], end(current_exons)[i], end(current_exons)[i]), y=c(0.4, 0.6, 0.6, 0.4), col="black", border = "black", lwd=1)
  }
  
  # Add boxes for first 10Kb of intron
  first_10Kb <- GenomicRanges::reduce(unlist(tile(x = current_introns, width = 1))[1:1e4])
  for(i in 1:length(first_10Kb)){
    if(as.logical(strand(current_exons)[1]=="+")){
      polygon(x=c(start(first_10Kb)[i], start(first_10Kb)[i], end(first_10Kb)[i], end(first_10Kb)[i]), y=c(0.65, 0.7, 0.7, 0.65), col="#4575b4", border = NA)
    } else {
      polygon(x=c(start(first_10Kb)[i], start(first_10Kb)[i], end(first_10Kb)[i], end(first_10Kb)[i]), y=c(0.65, 0.7, 0.7, 0.65), col="#d73027", border = NA)
    }
  }
  
  # Add boxes for last 10Kb of intron
  last_10Kb <- unlist(tile(x = current_introns, width = 1))
  last_10Kb <- GenomicRanges::reduce(last_10Kb[(length(last_10Kb)-1e4+1):length(last_10Kb)])
  for(i in 1:length(last_10Kb)){
    if(as.logical(strand(current_exons)[1]=="+")){
      polygon(x=c(start(last_10Kb)[i], start(last_10Kb)[i], end(last_10Kb)[i], end(last_10Kb)[i]), y=c(0.65, 0.7, 0.7, 0.65), col="#d73027", border = NA)
    } else {
      polygon(x=c(start(last_10Kb)[i], start(last_10Kb)[i], end(last_10Kb)[i], end(last_10Kb)[i]), y=c(0.65, 0.7, 0.7, 0.65), col="#4575b4", border = NA)
    }
  }
}

# Save out the gene structure
# CACNA1C
pdf("figures/Fig_3/CACNA1C_gene_structure.pdf", width=4, height=0.5)
plot_gene("ENST00000399655.5")
dev.off()

# LDLRAD4
pdf("figures/Fig_3/LDLRAD4_gene_structure.pdf", width=4, height=0.5)
plot_gene("ENST00000399848.7")
dev.off()

# VCL
pdf("figures/Fig_3/VCL_gene_structure.pdf", width=4, height=0.5)
plot_gene("ENST00000372755.7")
dev.off()

# LIMA1
pdf("figures/Fig_3/LIMA1_gene_structure.pdf", width=4, height=0.5)
plot_gene("ENST00000552783.5")
dev.off()
# Create a small PNG as well
png("figures_track/LIMA1_gene_structure.png", width=500, height=80)
plot_gene("ENST00000552783.5")
dev.off()

print("Finished plotting gene schematics for Fig_3")
