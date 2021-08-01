### pre-mRNA expression dynamics for genes of different length
# Now we will plot the expression profile for the first and last 10Kb of each gene's introns.
# We'll also include the mRNA expression on the plots.
# To make it easier to see the trend in the expression profile for each gene, we will also include impulse model fits.

# Check to see if plot data has already been obtained
if(!file.exists("data/Fig_3_gene_plot_data.rds")){
  
  # Load R Packages
  library(DESeq2)
  library(tidyverse)
  library(rtracklayer)
  library(GenomicFeatures)
  library(GenomicAlignments)
  library(Rsamtools)
  library(minpack.lm)
  
  # Import gene counts
  gene_counts <- readRDS("data/gene_counts.rds")
  
  # Get size factors for normalization
  sfs <- estimateSizeFactorsForMatrix(gene_counts$counts)
  
  # normalize the counts
  gene_counts$counts <- gene_counts$counts/sfs[col(gene_counts$counts)]
  
  # Get just the counts
  gene_counts <- data.frame(cbind(id=row.names(gene_counts$counts), as.data.frame(gene_counts$counts)))
  # Give sensible column names
  colnames(gene_counts) <- c("gene_id", paste0("T98G_", str_pad(seq(0,400,10),3,"left", 0),"_min"))
  
  # We want to match the gene_counts expression profiles with specified transcript IDs
  # Get GENCOE gtf so that we can match gene_id with the transcript IDs
  GENCODE <- readGFF("data/annotation/GRCh38_spiked.gtf", version = 2) %>%
    filter(type=="transcript") %>%
    dplyr::select(gene_id, transcript_id)
  
  # Load GTF data as a TxDb class object to get intron, exon and UTR intervals
  TxDb <- makeTxDbFromGFF("data/annotation/GRCh38_spiked.gtf", format="gtf")
  introns <- intronsByTranscript(TxDb, use.names=TRUE)
  
  # Create BamFileList
  fls <- BamFileList(file=paste0("/data/alignments/T98G_", stringr::str_pad(seq(0,400,10), width=3, side="left", pad="0"), "_minAligned.sortedByCoord.out.bam"),
                     index=paste0("/data/alignments/T98G_", stringr::str_pad(seq(0,400,10), width=3, side="left", pad="0"), "_minAligned.sortedByCoord.out.bam.bai"))
  
  # Define function to calculate overlap with an interval, using RNA-seq reads from all time points
  allReadsCounts <- function(i, interval){
    suppressWarnings(countOverlaps(interval,
                                   readGAlignmentPairs(file=fls[[i]], use.names=TRUE, param=ScanBamParam(which = interval, simpleCigar=TRUE), with.which_label=FALSE, strandMode=2),minoverlap=5))
  }
  
  # Define a function that accepts a transcript ID as input and returns the expression profile of the first/last 10Kb of intronic regions as well as the mRNA expression profile
  get_plot_data <- function(tx_id){
    # Get GRanges that define the first and last 10Kb
    if(as.logical(strand(introns[[tx_id]])[1]=="+")){
      first_10Kb <- GenomicRanges::reduce(unlist(tile(x = introns[[tx_id]], width = 1))[1:1e4])
      last_10Kb <- unlist(tile(x = introns[[tx_id]], width = 1))
      last_10Kb <- GenomicRanges::reduce(last_10Kb[(length(last_10Kb)-1e4+1):length(last_10Kb)])
    } else {
      first_10Kb <- unlist(tile(x = introns[[tx_id]], width = 1))
      first_10Kb <- GenomicRanges::reduce(first_10Kb[(length(first_10Kb)-1e4+1):length(first_10Kb)])
      last_10Kb <- GenomicRanges::reduce(unlist(tile(x = introns[[tx_id]], width = 1))[1:1e4])
    }
    
    # Get expression values for first and last 10Kb, normalize and scale
    # first_10Kb
    first_10Kb <- colSums(data.frame(lapply(1:41, allReadsCounts, interval=first_10Kb)))/sfs
    first_10Kb <- first_10Kb-min(first_10Kb); first_10Kb <- first_10Kb/max(first_10Kb)
    # last_10Kb
    last_10Kb <- colSums(data.frame(lapply(1:41, allReadsCounts, interval=last_10Kb)))/sfs
    last_10Kb <- last_10Kb-min(last_10Kb); last_10Kb <- last_10Kb/max(last_10Kb)
    
    # Get mRNA and scale
    mRNA <- unlist(gene_counts[gene_counts$gene_id==GENCODE$gene_id[which(GENCODE$transcript_id==tx_id)],2:42])
    mRNA <- mRNA-min(mRNA); mRNA <- mRNA/max(mRNA)
    
    return(data.frame(mRNA, first_10Kb, last_10Kb))
  }
  
  # Get expression data needed for plotting
  CACNA1C_raw <- get_plot_data("ENST00000399655.5")
  LDLRAD4_raw <- get_plot_data("ENST00000399848.7")
  VCL_raw <- get_plot_data("ENST00000372755.7")
  LIMA1_raw <- get_plot_data("ENST00000552783.5")
  
  ## One more thing we'd like for plotting is impulse model fits to the raw data.
  # Define time
  t<-seq(0,400,10)
  
  # Define 7-parameter impulse model function that allows two transitions, where the transitions can have different slopes, defined by lambda1 and lambda2.
  impulse2b <- function(h0, h1, h2, lambda1, lambda2, t, t1, t2) { (1/h1)*(h0+(h1-h0)*(1/(1+exp(-lambda1*(t-t1)))))*(h2+(h1-h2)*(1/(1+exp(lambda2*(t-t2))))) }
  
  # Define a function that will fit the impulse2b model and return the best fit after num_it iterations
  fit_impulse2b <- function(input_time_series, t=seq(0,400,10), num_it=100){
    temp <- data.frame(t=t, M=input_time_series)
    best_NRMSD <- 1E6
    best_fit <- NA
    for(i in 1:num_it){
      # Set fit to 0 to start with
      impulse_fit <- 0
      # Continue trying to the model to the data until a fit is found
      while(is.numeric(impulse_fit)){
        impulse_fit <- tryCatch({
          nlsLM(M ~ impulse2b(h0, h1, h2, lambda1,lambda2, t, t1, t2), 
                data = temp, start = c(h0=runif(1, min=min(input_time_series), max = max(input_time_series)), # the initial expression value
                                       h1=runif(1, min=min(input_time_series), max = max(input_time_series)), # the peak expression value
                                       h2=runif(1, min=min(input_time_series), max = max(input_time_series)), # the steady-state expression value
                                       lambda1=runif(n=1, min=0, max=1), # Slope of the first transition
                                       lambda2=runif(n=1, min=0, max=1), # Slope of the second transition
                                       t1=sample(1:max(t),1), # onset time is the time of the first transition (where rise or fall is maximal)
                                       t2=sample(1:max(t),1)), # offset time is the time of the second transition
                control = list(maxiter = 1024,warnOnly=TRUE), trace = F)
        }, error = function(e) { # if nlsLM can't converge, ouput zero
          c(0)
        }, finally = {
          
        })
      }
      # calculate the NRMSD
      current_NRMSD <- sqrt(sum((residuals(impulse_fit))^2)/length(residuals(impulse_fit)))/diff(range(temp$M))
      # Check if current fit is more optimal
      if(current_NRMSD < best_NRMSD){ 
        best_NRMSD <- current_NRMSD
        best_fit <- impulse_fit
      }
    }
    return(best_fit)
  }
  
  # Get impulse model fits to expression data for plotting purposes
  CACNA1C_fit <- data.frame(mRNA = predict(fit_impulse2b(CACNA1C_raw$mRNA), data.frame(t=seq(0,400,1))),
                            first_10Kb = predict(fit_impulse2b(CACNA1C_raw$first_10Kb), data.frame(t=seq(0,400,1))),
                            last_10Kb = predict(fit_impulse2b(CACNA1C_raw$last_10Kb), data.frame(t=seq(0,400,1))))
  
  LDLRAD4_fit <- data.frame(mRNA = predict(fit_impulse2b(LDLRAD4_raw$mRNA), data.frame(t=seq(0,400,1))),
                            first_10Kb = predict(fit_impulse2b(LDLRAD4_raw$first_10Kb), data.frame(t=seq(0,400,1))),
                            last_10Kb = predict(fit_impulse2b(LDLRAD4_raw$last_10Kb), data.frame(t=seq(0,400,1))))
  
  VCL_fit <- data.frame(mRNA = predict(fit_impulse2b(VCL_raw$mRNA), data.frame(t=seq(0,400,1))),
                        first_10Kb = predict(fit_impulse2b(VCL_raw$first_10Kb), data.frame(t=seq(0,400,1))),
                        last_10Kb = predict(fit_impulse2b(VCL_raw$last_10Kb), data.frame(t=seq(0,400,1))))
  
  LIMA1_fit <- data.frame(mRNA = predict(fit_impulse2b(LIMA1_raw$mRNA), data.frame(t=seq(0,400,1))),
                          first_10Kb = predict(fit_impulse2b(LIMA1_raw$first_10Kb), data.frame(t=seq(0,400,1))),
                          last_10Kb = predict(fit_impulse2b(LIMA1_raw$last_10Kb), data.frame(t=seq(0,400,1))))
  
  # Save out
  saveRDS(list(CACNA1C_raw=CACNA1C_raw, LDLRAD4_raw=LDLRAD4_raw, VCL_raw=VCL_raw, LIMA1_raw=LIMA1_raw,
               CACNA1C_fit=CACNA1C_fit, LDLRAD4_fit=LDLRAD4_fit, VCL_fit=VCL_fit, LIMA1_fit=LIMA1_fit), file = "data/Fig_3_gene_plot_data.rds")
}
  
print("Finished creating data required for Fig_3b-e")
