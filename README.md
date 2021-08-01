# lncRNA Time Series Analysis
This repository was created to assist with the reproducibility of all the analyses and figures presented in the manuscript entitled “No evidence for lncRNA cis-regulatory roles from high temporal resolution RNA-seq time-course data”. 

To run the analysis from start to finish the script [pipeline.sh](https://github.com/WalterMuskovic/lncRNA_time_course/blob/master/code/pipeline.sh) can be used. To provide a high-level understanding of the important steps in the analysis, a description of each step is included below along with some examples of the techniques used. If you want to reproduce a particular figure, I recommend reading the relevant text below then following the link to the associated script.

All of the data and software used in this analysis are freely available, with links to software and datasets provided in the sections below. This document is split into the following parts:

1. [Software overview](#software-overview) - links to R packages and other software used
2. [Data pre-processing](#data-pre-processing) - download and process raw data to produce analysis-ready files
3. [Figure 1](#figure-1) - mRNA and lncRNA expression
4. [Figure 2](#figure-2) - the effect of transcript stability
5. [Figure 3](#figure-3) - the effect of gene length
6. [Figure 4](#figure-4) - compare lncRNA and coding gene pre-mRNA expression dynamics
7. [Figure 5](#figure-5) - examples of human gene and adjacent lncRNA expression
8. [Figure 6](#figure-6) - relationship between genomic distance, correlation and timing of human protein-coding gene/lncRNA expression
9. [Figure 7](#figure-7) - relationship between genomic distance, correlation and timing of mouse protein-coding gene/lncRNA expression
10. [Figure S1](#figure-s1) - protein-coding gene and lncRNA length
11. [Figure S2](#figure-s2)  - estimate mean Pol II transcription elongation rate
12. [Figure S3](#figure-s3)  - spatial correlation trend among protein-coding genes and among lncRNAs
13. [Figure S4](#figure-s4)  - spatial correlation trend between protein-coding genes and lncRNAs that do, and do not, overlap annotated lncRNAs

## Software overview
Most of the statistical analysis and plotting was done with the R programming language
* [R](https://www.r-project.org/)

Data manipulation and helper functions
* [tidyverse](https://www.tidyverse.org/)
* [glue](https://cran.r-project.org/package=glue)
* [tictoc](https://cran.r-project.org/package=tictoc)
* [furrr](https://cran.r-project.org/package=furrr)

Manipulation of genomics data files
* [Rsamtools](https://bioconductor.org/packages/release/bioc/html/Rsamtools.html)
* [rtracklayer](https://www.bioconductor.org/packages/release/bioc/html/rtracklayer.html)
* [GenomicAlignments](https://bioconductor.org/packages/release/bioc/html/GenomicAlignments.html)
* [GenomicFeatures](https://bioconductor.org/packages/release/bioc/html/GenomicFeatures.html)
* [samtools](http://www.htslib.org/)

Access to publicly available genomics data files
* [AnnotationHub](http://www.bioconductor.org/packages/AnnotationHub)

Plotting
* [vioplot](https://cran.r-project.org/package=vioplot)
* [gplots](https://cran.r-project.org/package=gplots)
* [RColorBrewer](https://cran.r-project.org/package=RColorBrewer)
* [ggsci](https://cran.r-project.org/package=ggsci)
* [viridis](https://cran.r-project.org/package=viridis)

Combining, arranging and annotating plot panels
* [Inkscape](https://inkscape.org/)

Solving differential equations and nonlinear regression
* [deSolve](https://cran.r-project.org/package=deSolve)
* [minpack.lm](https://cran.r-project.org/package=minpack.lm)

Applying generalized additive models
* [mgcv](https://cran.r-project.org/package=mgcv)

Quantifying and normalizing gene expression data
* [Rsubread](https://bioconductor.org/packages/release/bioc/html/Rsubread.html)
* [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)

RNA-seq data alignment and processing
* [cutadapt](https://github.com/marcelm/cutadapt)
* [STAR](https://github.com/alexdobin/STAR)

lncRNA annotation
* [StringTie](https://ccb.jhu.edu/software/stringtie/)

## Data pre-processing
Data preparation can be broken down into the following main steps
* [Download sequence and annotation files](#download-and-pre-process-external-sequence-and-annotation-files)
* [Prepare FASTQ files for alignment](#prepare-t98g-fastq-files-for-alignment)
* [Prepare STAR genome index](#prepare-grch38-star-index)
* [Align human RNA-seq reads](#align-t98g-rna-seq-reads)
* [Download and process the mouse LPS-response time course data](#download-and-pre-process-gse56977)
* [Prepare lncRNA annotations](#prepare-lncrna-annotation-files)
* [Quantify gene expression with featureCounts](#quantify-gene-expression-using-featurecounts)

Intermediate annotation and bed files produced can be viewed in a genome browser such as the Broad Institute's excellent [IGV](https://software.broadinstitute.org/software/igv/) browser.

### Download and pre-process external sequence and annotation files
The human genome (GRCh38) and RNA [sequin](https://doi.org/10.1038/nmeth.3958) spike-in sequences files used were:
* [GRCh38.primary_assembly.genome.fa.gz](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh38.primary_assembly.genome.fa.gz)
* [chrIS.v2.fa](https://s3.amazonaws.com/sequins/annotations/chrIS.v2.fa)

The sequence files were combined with
```
cat GRCh38.primary_assembly.genome.fa chrIS.v2.fa > GRCh38_spiked.fa
```

GRCh38 annotation and spike-in annotation files used were:
* [gencode.v29.annotation.gtf.gz](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz)
* [RNAsequins.v2.2.gtf](https://s3.amazonaws.com/sequins/annotations/RNAsequins.v2.2.gtf)

Annotations were combined with
```
cat gencode.v29.annotation.gtf RNAsequins.v2.2.gtf > GRCh38_spiked.gtf
```

GRCm38 mouse genome and gene annotation files used were:
* [gencode.vM21.annotation.gtf.gz](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/gencode.vM21.annotation.gtf.gz)
* [GRCm38.primary_assembly.genome.fa.gz](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M21/GRCm38.primary_assembly.genome.fa.gz)

### Prepare T98G FASTQ files for alignment
Prior to alignment, Illumina adapters were trimmed from the T98G stranded paired-end RNA-seq reads with cutadapt as follows:
```
for time_point in {000..400..10}
do 
  cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
  -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
  --minimum-length 15 \
  -o "./trimmed_FASTQ/trimmed_T98G_${time_point}_min_R1.FASTQ.gz" \
  -p "./trimmed_FASTQ/trimmed_T98G_${time_point}_min_R2.FASTQ.gz" \
  "T98G_${time_point}_min_R1.FASTQ.gz" "T98G_${time_point}_min_R2.FASTQ.gz"
done
```

### Prepare GRCh38 STAR index
A GRCh38 STAR genome index was prepared for 2-pass alignment in two steps, starting with first-pass mapping of all T98G totalRNA samples:
```
# Create index
STAR --runMode genomeGenerate --genomeFastaFiles .data/genome/GRCh38_spiked.fa --sjdbGTFfile .data/annotation/GRCh38_spiked.gtf --sjdbOverhang 124 --runThreadN 10 --genomeDir .data/STAR_index --outFileNamePrefix STAR_GRCh38

# Align all samples
for time_point in {000..400..10}
do
STAR --genomeDir .data/STAR_index/ --readFilesIn ".data/trimmed_FASTQ/trimmed_T98G_${time_point}_min_R1.FASTQ.gz" ".data/trimmed_FASTQ/trimmed_T98G_${time_point}_min_R2.FASTQ.gz" \
    --readFilesCommand zcat --runThreadN 12 --genomeLoad NoSharedMemory \
    --alignSJoverhangMin 8 --alignSJDBoverhangMin 3 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
    --outSAMheaderCommentFile ,data/temp_STAR/COfile.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate \
    --outFilterType BySJout --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --outFilterMultimapNmax 20 \
    --quantMode TranscriptomeSAM GeneCounts --sjdbScore 1 \
    --limitBAMsortRAM 60000000000 --outFileNamePrefix ".data/temp_STAR/T98G_${time_point}_min" \
    --sjdbGTFfile .data/annotation/GRCh38_spiked.gtf \
    --alignEndsType Local
done
```
In the second step, all *SJ.out.tab files were combined and used to rebuild the STAR genome index with novel junctions included
```
# concatenate SJ.out.tab files produced by all alignments
cat .data/temp_STAR/*SJ.out.tab > .data/STAR_index/combined.SJ.out.tab 

# Re-build STAR genome index, supplying combined junctions file combined.SJ.out.tab from first pass of all time points
STAR --runMode genomeGenerate --genomeFastaFiles .data/genome/GRCh38_spiked.fa --sjdbGTFfile .data/annotation/GRCh38_spiked.gtf --sjdbFileChrStartEnd .data/STAR_index/combined.SJ.out.tab --limitSjdbInsertNsj 1500000 --sjdbOverhang 124 --runThreadN 10 --genomeDir .data/STAR_index --outFileNamePrefix STAR_T98G
```

### Align T98G RNA-seq reads
Using the newly created genome index, human glioblastoma T98G time series adapter-trimmed fastq files were aligned with the following:
```
for time_point in {000..400..10}
do
# Carry out STAR 2nd-pass paired-end read mapping
STAR --genomeDir .data/STAR_index/ --readFilesIn ".data/trimmed_FASTQ/trimmed_T98G_${time_point}_min_R1.FASTQ.gz"  ".data/trimmed_FASTQ/trimmed_T98G_${time_point}_min_R2.FASTQ.gz"  \
    --readFilesCommand zcat --runThreadN 12 --genomeLoad NoSharedMemory \
    --alignSJoverhangMin 8 --alignSJDBoverhangMin 3 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
    --outSAMheaderCommentFile .data/alignments/COfile.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate \
    --outFilterType BySJout --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --outFilterMultimapNmax 20 \
    --quantMode TranscriptomeSAM GeneCounts --sjdbScore 1 --limitSjdbInsertNsj 1500000 \
    --limitBAMsortRAM 60000000000 --outFileNamePrefix ".data/alignments/T98G_${time_point}_min" \
    --sjdbGTFfile .data/annotation/GRCh38_spiked.gtf \
    --alignEndsType Local

# Create BAM index file
samtools-1.3.1/bin/samtools index "T98G_${time_point}_minAligned.sortedByCoord.out.bam"
done
```
The read depth obtained for each time point is summarised below:
![alt text](https://github.com/WalterMuskovic/lncRNA_time_course/raw/master/figures_track/T98G_read_depth.png "T98G sample read depth")

### Download and pre-process GSE56977
In addition to the human glioblastoma T98G time course data we access a publicly available mouse time course dataset described in [Rabani et al, Cell,  2014](https://doi.org/10.1016/j.cell.2014.11.015). The data is available using this accession [GSE56977](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE56977). The dataset describes mouse dendritic cells responding to LPS at 0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165 and 180 min intervals post-treatment, with one total RNA-seq sample per time point. The following code will retrieve and process the GSE56977 mouse LPS-response time series data (SRA > FASTQ > (QC + trimming) > BAM).
```
# Prefetch GSE56977 SRA files
for accession in {1258347..1258359}
do
  prefetch -O ./data/mouse_temp "SRR${accession}" &
done
wait

# Extract fastq files
for accession in {1258347..1258359}
do
  fastq-dump.2 --gzip --split-3 -O ./data/mouse_temp ./data/mouse_temp/"SRR${accession}.sra" &
done
wait

# Trim Illumina TruSeq adapter sequences with cutadapt
cd ./data/mouse_temp
for accession in {1258347..1258359}
do 
  cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
  --minimum-length 15 \
  -o "trimmed_SRR${accession}_1.FASTQ.gz" \
  -p "trimmed_SRR${accession}_2.FASTQ.gz" \
  "SRR${accession}_1.FASTQ.gz" "SRR${accession}_2.FASTQ.gz"
done

# Give files meaningful names and move to the fastq directory

# Create STAR genome index for GRCm38
mkdir ./data/STAR_index_mouse_101bp
STAR --runMode genomeGenerate --genomeFastaFiles ./data/genome_mouse/GRCm38.primary_assembly.genome.fa --sjdbGTFfile ./data/annotation_mouse/gencode.vM21.annotation.gtf --sjdbOverhang 100 --runThreadN 10 --genomeDir ./data/STAR_index_mouse_101bp --outFileNamePrefix STAR_GRCm38_101bp

# Align samples using the newly created genome index
for time_point in total_000 total_015 total_030 total_045 total_060 total_075 total_090 total_105 total_120 total_135 total_150 total_165 total_180
do
  # Align
  cd ./data/alignments
  STAR --genomeDir ./data/STAR_index_mouse_101bp --readFilesIn "./data/FASTQ/Rabani_${time_point}_min_R1.FASTQ.gz" "./data/FASTQ/Rabani_${time_point}_min_R2.FASTQ.gz" \
      --readFilesCommand zcat --runThreadN 10 --genomeLoad NoSharedMemory \
      --alignSJoverhangMin 8 --alignSJDBoverhangMin 3 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
      --outSAMheaderHD @HD VN:1.4 SO:coordinate --outSAMunmapped Within --outSAMtype BAM SortedByCoordinate \
      --outFilterType BySJout --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --outFilterMultimapNmax 20 \
      --quantMode GeneCounts --sjdbScore 1 \
      --limitBAMsortRAM 60000000000 --outFileNamePrefix "./data/alignments/Rabani_${time_point}_min" \
      --sjdbGTFfile ./data/annotation_mouse/gencode.vM21.annotation.gtf \
      --alignEndsType Local

  # Create BAM index file
  samtools index "Rabani_${time_point}_minAligned.sortedByCoord.out.bam"
done
```
The read depth for each time point in GSE56977 is summarised below:
![alt text](https://github.com/WalterMuskovic/lncRNA_time_course/raw/master/figures_track/GSE56977_read_depth.png "GSE56977 sample read depth")

### Prepare lncRNA annotation files

A *de novo* transcriptome annotation was created with StringTie using the following:
```
# Merge BAM files
samtools merge --threads 64 merged.bam -b bam_file_list.txt 
samtools merge --threads 64 Rabani_merged.bam -b $data_dir/Rabani_bam_file_list.txt 
samtools index merged.bam
samtools index Rabani_merged.bam

# Run StringTie
stringtie merged.bam -G GRCh38_spiked.gtf -o T98G_stringtie.gtf -j 1 -a 10 -A T98G_gene_abundance.tab -m 200 -l T98G_STRG -p 64 -v --rf -x 'chrM,chrIS'
stringtie Rabani_merged.bam -G gencode.vM21.annotation.gtf -o Rabani_stringtie.gtf -j 1 -a 10 -A Rabani_gene_abundance.tab -m 200 -l Rabani_STRG -p 64 -v --rf -x 'chrM'
```


### Quantify gene expression using featureCounts
Using the BAM files generated above with STAR, we quantify the expression using Rsubread's featureCounts function:

```
gene_counts <- featureCounts(files = bam_files,
                   annot.ext = anno.gtf, isGTFAnnotationFile = TRUE, GTF.featureType = "exon", GTF.attrType = "gene_id",
                   isPairedEnd = TRUE,
                   minFragLength = 20,
                   maxFragLength = 5000,
                   strandSpecific = 2,
                   countMultiMappingReads = TRUE,
                   allowMultiOverlap = TRUE,
                   nthreads = 10)
```

## Figure 1
This figure is made up of:

 * a. Schematic of when the time points were taken - made elsewhere
 * b. lncRNA expression heatmap
 * c. mRNA expression heatmap
 * d. Centroid correlation matrix plot

The code used to generate the panels for the figure can be found [here](https://github.com/WalterMuskovic/lncRNA_time_course/blob/master/code/create_Fig_1.R). We use K-means to cluster genes. To determine a suitable number of clusters for visualisation we determine at what k the total intra-cluster variation is minimized as follows:
```
# Calculate correlation matrix
cor_matrix <- cor(t(gene_expression), method="pearson")
# Define a function to compute the total within-cluster sum of square (wss) - compactness of the clustering, which we seek to minimise -
wss <- function(k) { kmeans(cor_matrix, k, nstart = 10)$tot.withinss }
# Compute and plot wss for k = 2 to k = 15
k.values <- 2:15
plan(multiprocess)
wss_values <- future_map_dbl(k.values, wss)
```
![alt text](https://github.com/WalterMuskovic/lncRNA_time_course/blob/master/figures_track/mRNA_wss.png "mRNA WSS")

To plot the heatmaps we can use the `heatmap.2` function from [gplots](https://cran.r-project.org/package=gplots):
```
# Get colours for sidebars which indicate cluster membership
cols <- ggsci::pal_npg("nrc", alpha = 1)(length(unique(cluster_info)))
cols <-  cols[cluster_info]

# Scale rows to have mean zero and standard deviation one
plot_data <- t(apply(plot_data, 1, scale))

# Note that we've scaled the input expression matrix on rows (genes) to Z-scores to facilitate meaningful visual comparison. Some genes have extreme values (a very rapid induction - and correspondingly rapid decrease) and tend to make the rest of the heatmap a little faint. To facilitate visualisation we can use the breaks argument to reduce the effect of the small numer of genes with expression outside a z-score range of -3 to 3.
breaks <- seq(-3, 3, by=0.2)
# add outliers
breaks <- append(breaks, max(as.vector(plot_data)))
breaks <- append(breaks, min(as.vector(plot_data)),0)

# Get colours for the heatmap
heatmap_cols <- colorRampPalette(brewer.pal(9,"PuBuGn"))(length(breaks)-1)

# Plot the heatmap
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
```
![alt text](https://github.com/WalterMuskovic/lncRNA_time_course/blob/master/figures_track/Fig1b.png "mRNA WSS")

## Figure 2
This figure illustrates the effect of transcript stability. We use a simple model of transcription:

![alt text](https://github.com/WalterMuskovic/lncRNA_time_course/blob/master/figures_track/transcription_eq.png "transcription model")

to infer the half-lives of these transcripts. We fit the transcription model using the following:
```
# Define some example expression data
pre_mRNA <- c(82,231,851,1210,1257,1099,972,823,792,612,531,452,459,383,393,348,386,332,296,313,272,257,269,258,302,219,219,197,233,194,180,177,182,184,159,206,167,183,173,157,184)
mRNA <- c(493,548,1119,1952,2155,2724,3273,3358,4205,3688,3898,4349,4191,3939,4135,4152,4200,3840,3740,3969,3568,3729,3392,3594,3549,2892,2910,2732,2740,2433,2788,2556,2353,2016,2100,2004,2085,1996,2021,1841,1683)


# state variable
state <- c(m=mRNA[1])

# Time specification.
times <- seq(0,400,10)

# Define P(t) as an external variable
p <- approxfun(as.data.frame(list(times = times, p = pre_mRNA)), rule = 2)

# Model equations
get_m <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # Define rates of change
    P <- p(t)
    dm <- beta*P - alpha*m
    
    # Return the rate of change
    list(dm)
  })
}

# Define function that accepts values for alpha & beta then returns the values for M(t) this gives
get_fit_m <- function(alpha, beta){
  parameters <- c(alpha=alpha, beta=beta)
  out <- data.frame(ode(y = state, times = times, func = get_m, parms = parameters, method = "euler"))
  out$m
}

# Fit the model 
NRMSD <- 1E10
best_model <- NA
set.seed(1234)

for(n in 1:10){
  transcription_model_fit <- 0
  while(is.numeric(transcription_model_fit)){
    transcription_model_fit <- tryCatch({
      nlsLM(m ~ get_fit_m(alpha, beta), data = data.frame(m=mRNA), start = list(alpha=runif(1,0,1), beta=runif(1,0,1)), control = list(maxiter = 1024, warnOnly=TRUE), trace = F)
    }, error = function(e) { # if nlsLM can't converge, ouput zero
      c(0)
    }, finally = {
      
    })
  }
  
  current_NRMSD <- sqrt(sum((residuals(transcription_model_fit))^2)/length(residuals(transcription_model_fit)))/diff(range(mRNA))
  
  if(current_NRMSD < NRMSD){
    NRMSD <- current_NRMSD
    best_model <- transcription_model_fit
  }
}

# Get the transcription model fit
best_model_fit <- predict(best_model)

# Get the inferred degradation rate
alpha <- coefficients(best_model)["alpha"]

#Plot
par(mfcol=c(1,2))
plot(t, pre_mRNA, pch=16, xlab="time (min)", ylab="pre-mRNA expression")
plot(t, mRNA, pch=16, xlab="time (min)", ylab="mRNA expression"); lines(t, best_model_fit)
```
![alt text](https://github.com/WalterMuskovic/lncRNA_time_course/blob/master/figures_track/transcription_fit_example.png "transcription model fit example")

To assist in visualizing the pre-mRNA expression data we also fit impulse models to the expression estimates obtained from RNA-seq. The impulse model is described by [Chechik and Koller](https://doi.org/10.1089/cmb.2008.13TT). They describe a parametric model that is a product of two sigmoid functions. Their original model has six free parameters; h0, h1, h2, t1, t2, and lambda. Three amplitude parameters determine the initial amplitude h0, the peak amplitude h1, and the steady state amplitude h2 of the response. The onset time t1 is the time of first transition (where rise or fall is maximal) and the offset t2 is the time of second transition. Lambda is the slope parameter. As the authors point out, this model is easily generalisable and we modify it slightly here to include lambda1 and lambda2, such that the two transitions may have different slopes. So we get:

![alt text](https://github.com/WalterMuskovic/lncRNA_time_course/blob/master/figures_track/impulse2b_eq.png "impulse 2b")

To estimate the seven free parameters we use the `nlsLM()` function from the [minpack.lm](https://cran.r-project.org/package=minpack.lm) R package. Note that we also use the six-parameter version as well as a nine-parameter model that allows up to three transitions. You can see the code used to apply these in [create_Fig_2.R](https://github.com/WalterMuskovic/lncRNA_time_course/blob/master/code/create_Fig_2.R), but here we will just look at an example of fitting the 7-parameter model:
```
# Define 7-parameter impulse model function that allows two transitions, where the transitions can have different slopes, defined by lambda1 and lambda2
impulse2b <- function(h0, h1, h2, lambda1, lambda2, t, t1, t2) { (1/h1)*(h0+(h1-h0)*(1/(1+exp(-lambda1*(t-t1)))))*(h2+(h1-h2)*(1/(1+exp(lambda2*(t-t2))))) }

# Define time
t<-seq(0,400,10)

# Define some example expression data
expression_data <- c(82,231,851,1210,1257,1099,972,823,792,612,531,452,459,383,393,348,386,332,296,313,272,257,269,258,302,219,219,197,233,194,180,177,182,184,159,206,167,183,173,157,184)

# Define function that will fit the impulse2b model and return the best fit after 100 iterations
fit_impulse2b <- function(input_time_series, t){
  temp <- data.frame(t=t, M=input_time_series)
  best_NRMSD <- 1E6
  best_fit <- NA
  for(i in 1:100){
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

# Fit the model
set.seed(1234)
model_fit <- predict(fit_impulse2b(input_time_series = expression_data, t=t), data.frame(t=t))

#Plot
plot(t, expression_data, pch=16, xlab="time (min)", ylab="expression"); lines(t, model_fit)
```
![alt text](https://github.com/WalterMuskovic/lncRNA_time_course/blob/master/figures_track/impulse2b_example.png "impulse 2b example fit")

## Figure 3
The previous figure illustrated the effect of transcript stability on mRNA expression dynamics. Here, in Fig. 3 we illustrate the effect of gene length. Panel A captures in detail transcription across the relatively long calcium channel gene CACNA1C. We capture this data using a function that accepts a transcript id as input as well as a window_width. The function will return the coverage across the specified windows for all introns of the transcript, across all time points. Rows of the returned data frame correspond to sequential genomic windows (5' -> 3'), columns correspond to time points:
```
# Load GTF data to get intron and exon intervals
TxDb <- makeTxDbFromGFF("data/annotation/GRCh38_spiked.gtf", format="gtf")
introns <- intronsByTranscript(TxDb, use.names=TRUE)

# Create BamFileList
fls <- BamFileList(file=bam_files, index=bam_file_indexes)

# Get size factors for normalization
sfs <- estimateSizeFactorsForMatrix(readRDS("data/gene_counts.rds")$counts)

# Define function to calculate overlap with interval, using all reads
allReadsCounts <- function(i, interval){
  countOverlaps(interval,
  readGAlignmentPairs(file=fls[[i]],
  use.names=TRUE,
  param=ScanBamParam(which = interval, simpleCigar=TRUE),
  with.which_label=FALSE, strandMode=2),
  minoverlap=5)
}

# Define function to pull coverage for a specific transcript
get_intron_cov <- function(tx_id, window_width=1000){
  # Get intervals of width ~1Kb, covering the desired transcript
  intervals <-unlist(tile(introns[[tx_id]], width=window_width))
  
  # Count overlaps from BAM files 
  # Load furrr package to speed up running time (not required, but helpful)
  library(furrr);plan(multiprocess)
  # In the resulting data frame, rows correspond to intervals, columns to time points
  allCounts <- data.frame(future_map(1:41, ~ allReadsCounts(., interval=intervals), .progress = TRUE))
  # Give meaningful column names
  colnames(allCounts) <- paste0("T98G_",str_pad(seq(0,400,10),width=3, side="left", pad="0"), "_min")
  
  # normalize 
  allCounts <- allCounts/sfs[col(allCounts)]
  
  # Return the data frame of intron coverage
  return(allCounts)
}

# Get the intron coverage for CACNA1C
get_intron_cov("ENST00000399655.5")
```
The data retrieved with `get_intron_cov` is used to make the following figure:
![alt text](https://github.com/WalterMuskovic/lncRNA_time_course/blob/master/figures_track/CACNA1C_ridge_plot.png "CACNA1C ridge plot")

To help understand how the expression is being quantified (using the first and last 10Kb of intron), we also include schematics of the genes presented in the figure. These were created using this code:
```
# Load GTF data to get intron and exon intervals
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

plot_gene("ENST00000552783.5")
```
![alt text](https://github.com/WalterMuskovic/lncRNA_time_course/blob/master/figures_track/LIMA1_gene_structure.png "LIMA1 gene schematic")

To get the data presented in the panels of Fig. 3b-e we use the two functions `allReadsCounts` and `get_plot_data` defined below:
```
# Define function to calculate overlap with an interval, using RNA-seq reads from all time points
  allReadsCounts <- function(i, interval){
    countOverlaps(interval,
    readGAlignmentPairs(file=fls[[i]],
    use.names=TRUE,
    param=ScanBamParam(which = interval, simpleCigar=TRUE),
    with.which_label=FALSE, strandMode=2),
    minoverlap=5))
  }
  
# Define a function that accepts a transcript ID as input and returns the expression profile of
# the first/last 10Kb of intronic regions as well as the mRNA expression profile
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
```
To facilitate visualisation we also include impulse model fits, obtained as described for Fig. 2 above. Together these make the the line plots presented in panels b-e, like this:

![alt text](https://github.com/WalterMuskovic/lncRNA_time_course/blob/master/figures_track/CACNA1C_line_plot.png "CACNA1C first and last 10Kb")

## Figure 4
Here we compare the expression dynamics of the first 10 Kb of intron of coding genes with the same gene's mRNA dynamics. We then compare lncRNA expression dynamics with coding gene pre-mRNA dynamics. 

The code used to produce these figures is very similar to the code used to create Fig. 1. In fact, the expression data in Fig. 4b is exactly the same as in Fig. 1b, the rows (genes) are just arranged differently. Genes are first ordered by the approximate peak in mRNA expression which we obtain with the following simple function:
```
# Define function to find the approximate time when mRNA expression is at its peak 
get_mRNA_peak_time <- function(input_mRNA){
  # Get mRNA expression
  mRNA_expression <- data.frame(time = seq(0,400,10), mRNA=unlist(input_mRNA))
  # Fit loess
  lw1 <- loess(mRNA ~ time, data = mRNA_expression, span = 0.5)
  # Return the time when mRNA expression is ~ maximal
  return(seq(0,400,10)[which.max(lw1$fitted)])
}
```
We can see how this function works with this simple example:
```
# Create some noisy gene expression data
set.seed(123)
test_data <- data.frame(time = seq(0,400,10),
                        true_gene_expression=sin(seq(0,4,0.1)),
                        gene_expression=sin(seq(0,4,0.1)) + rnorm(41, mean = 0, sd=0.1))

# Plot and add vertical line for our estimate of when gene expression peaks
ggplot(test_data, aes(x=time, y=gene_expression)) + 
  geom_point() +
  geom_line(aes(y = true_gene_expression)) +
  geom_smooth(span = 0.5, se = FALSE) +
  geom_vline(xintercept=get_mRNA_peak_time(test_data$gene_expression))
```
![alt text](https://github.com/WalterMuskovic/lncRNA_time_course/blob/master/figures_track/get_mRNA_peak_time.png "get mRNA peak time example")

Then genes are reordered by the pre-mRNA cluster number to give the final ordering. Note that the gene order is the same in Fig. 4a and Fig. 4b, we are just visualising the first 10 Kb of intron expression and mature mRNA expression respectively. The full code used to produce Fig.4 can be found in [create_Fig_4.R](https://github.com/WalterMuskovic/lncRNA_time_course/blob/master/code/create_Fig_4.R).

## Figure 5
This figure highlights the expression patterns of three protein-coding genes and adjacent lncRNAs. The data used to generate this figure is created using [create_data_for_Fig_5.R](https://github.com/WalterMuskovic/lncRNA_time_course/blob/master/code/create_data_for_Fig_5.R) and the figures are plotted with [create_Fig_5.R](https://github.com/WalterMuskovic/lncRNA_time_course/blob/master/code/create_Fig_5.R). The gene schematics at the top of panels a-c were created using the same approach from Fig. 3. 
The human epigenomics data were pulled from the [Roadmap Epigenomics Project](http://www.roadmapepigenomics.org/) using the [AnnotationHub R package](https://bioconductor.org/packages/release/bioc/html/AnnotationHub.html):

```
# Get annotation hub
ahub <- AnnotationHub()

# Extract BigWig records from the Roadmap Epigenomics Project
roadmap.bw <- subset(ahub, rdataclass == "BigWigFile" & preparerclass == "EpigenomeRoadMapPreparer")
# Create a dataframe containing the AnnotationHub unique identifiers (e.g. AH12345), RoadmapEpigenome name (e.g. E001), target of record the (e.g. H3K27ac), the URL of the resource and the type of BigWigFile (e.g. log10(p-value) or fold enrichment signal tracks)

# Then restrict to these entries;
# BigWigFiles from samples in the uniformly re-processed 111 reference human epigenomes. Note: E060 and E064 aren't present, so we are considering E001-E113.
# DNase-seq data and ChIP-seq data for the following histone modifications; H3K4me3 & H3K4me1

# Imputed pvalue signal tracks for ChIP-seq data
roadmap.bw <- data.frame(id = names(roadmap.bw),
                         epigenome = sapply(strsplit(mcols(roadmap.bw)$title,"-"), `[`, 1),
                         target = str_extract(mcols(roadmap.bw)$title, "H3K4me3|H3K4me1"),
                         sourceURL = mcols(roadmap.bw)$sourceurl,
                         track_type = str_extract(mcols(roadmap.bw)$title, "fc.signal|imputed.pval.signal|pval.signal|DNase.pval.signal"), stringsAsFactors = FALSE) %>%
  dplyr::filter(epigenome%in%str_glue("E{str_pad(1:113,3,'left',0)}") & target%in%c("DNase", "H3K4me3", "H3K4me1") & track_type%in%c("imputed.pval.signal","DNase.pval.signal"))

# Define ranges to extract data from
hg38 <- suppressWarnings(c(GRanges(str_remove_all("chr14:75,256,000-75,302,000",",")), 
                           GRanges(str_remove_all("chr5:135,989,000-136,080,000",",")),
                           GRanges(str_remove_all("chr18:3,435,000-3,671,000",","))))

# Give names
names(hg38) <- c("FOS", "TGFBI", "TGIF1")

# Get chain file to map hg38 to hg19
chainfiles <- query(ahub , c("hg38", "hg19", "chainfile"))
chainfiles
# Get hg38ToHg19.over.chain.gz
chain <- chainfiles[['AH14108']]
# Lift coordinates for specified range
hg19 <- unlist(liftOver(hg38, chain))

  # For each Roadmap Epigenomics file, import and store the BigWig data for the specified region  bw.region <- lapply(1:nrow(roadmap.bw), function(x) import(con = roadmap.bw$sourceURL[[x]], format = "bigWig", selection = hg19, as="NumericList"))
```
The figure panels were wrangled into shape using Inkscape. Panels d-f are simple line plots.

## Figure 6
Fig.6 explores the relationship between genomic distance, correlation and timing of protein-coding gene/lncRNA expression. To capture the relationship between genomic distance and the correlation between the expression profiles of coding and lncRNA transcripts we fit a generalized additive model (GAM) using the `gam()` function from the [mgcv R package](https://cran.r-project.org/web/packages/mgcv/index.html). To obtain confidence intervals around this fit we use a block bootstrap approach. To understand why we need this, you can read Mike Love's [excellent blog post](https://mikelove.wordpress.com/2012/07/28/block-bootstrap/) which gives a great overview. I adapted some of Mike's code to generate the following schematic plot of how the block bootstrap is employed for our data:

![alt text](https://github.com/WalterMuskovic/lncRNA_time_course/blob/master/figures_track/block_bootstrap_schematic.png "block bootstrap schematic")

In the top panel of the schematic we illustrate the spatial distribution along the genome of some hypothetical protein-coding genes (circles) and lncRNAs (triangles). The colours represent expression patterns. Pattern 1 may represent genes that are rapidly upregulated then rapidly downregulated (such as immediate early transcription factors) while genes from expression pattern 2 may represent genes with more gradually increasing expression. Genes or lncRNAs that are close to one another may tend to share the same expression pattern, for example clusters of lncRNAs with similar expression profiles. If we naively shuffled gene positions we would ignore this dependence. Instead we separate coding genes and lncRNAs (middle panels) and bin them into blocks that take into account this spatial dependence. The shuffled blocks are then combined to produce the bootstrap data seen in the bottom panel. The schematic is created using [bootstrap_schematic.R](https://github.com/WalterMuskovic/lncRNA_time_course/blob/master/code/bootstrap_schematic.R).

In panels b & c we compare coding gene expression profiles with lags of the lncRNA expression profiles. To do this efficiently we can use the [stats::ccf() R function](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/acf.html) to investigate how the correlation between two time series changes as their separation in time changes. This simple example script illustrates what `ccf` is calculating:

```
library(ggplot2)
library(patchwork)
library(forecast)

# Define function to generate expression data using an impulse model and adding some noise
impulse_plus_noise <- function(h0=0, h1=10, h2=0, lambda=0.1, t=seq(0,400,10), t1) { 
  t2=t1+120 # set offset time to 120 min after onset
  (1/h1)*((h0+(h1-h0)*(1/(1+exp(-lambda*(t-t1)))))*(h2+(h1-h2)*(1/(1+exp(lambda*(t-t2)))))) + rnorm(length(t), mean = 1, sd = 0.5)
}

plot_lags <- function(shift){
  # Create some gene expression data
  test_data <- data.frame(time = seq(0,400,10),
                          gene_1_expression = impulse_plus_noise(t1 = 100),
                          gene_2_expression = impulse_plus_noise(t1 = 100 + shift))
  
  # Plot expression of the two genes + ccf plot
  g1 <- ggplot(test_data, aes(x=time, y=gene_1_expression)) + geom_point()
  g2 <- ggplot(test_data, aes(x=time, y=gene_2_expression)) + geom_point()
  g3 <- ggCcf(x = test_data$gene_1_expression, y = test_data$gene_2_expression, lag.max = 20, type="correlation", plot = TRUE) + ggtitle("lagged expression")
  
  return(g1 | g2 | g3)
}

# Shift gene 2's expression forward and backwards in time relative to gene 1
set.seed(1234)
plot_lags(-60) / plot_lags(0) / plot_lags(60)
```
![alt text](https://github.com/WalterMuskovic/lncRNA_time_course/blob/master/figures_track/lag_example.png "lag example")


## Figure 7
The code used to produce this figure is very similar to that used for Fig. 5, but using the Rabani et al mouse LPS-response dataset (GSE56977). The code can be found in [create_data_for_Fig_7.R](https://github.com/WalterMuskovic/lncRNA_time_course/blob/master/code/create_Fig_5.R) and [create_Fig_7.R](https://github.com/WalterMuskovic/lncRNA_time_course/blob/master/code/create_Fig_7.R). 

## Figure S1
This figure illustrates the differences between protein-coding gene and lncRNA lengths. This one is pretty straightforward, the code used to generate it can be found in [create_Fig_S1.R](https://github.com/WalterMuskovic/lncRNA_time_course/blob/master/code/create_Fig_S1.R).

## Figure S2
This figure illustrates how we can use the time series data to come up with a rough estimate of the mean Pol II transcription elongation rate in our stimulated T98G cells. The following three-panel figure illustrates the principle using the gene LDLRAD4 as an example.

![alt text](https://github.com/WalterMuskovic/lncRNA_time_course/blob/master/figures_track/LDLRAD4_lag_plot.png "LDLRAD4 lag plots")

The top panel depicts the gene expression profile of the first (blue) and last (red) 10 Kb of intron of the gene. The solid lines represent impulse models fit to the expression data using non-linear least squares via the `nlsLM()` function, as was described for Fig. 2. Panel 2 illustrates how we then 'slide' the impulse model fit to the first 10 Kb forwards in time to find the closest match to the last 10 Kb. Panel 3 illustrates that the best fit is obtained using a lag of 140 min. This is then repeated and plotted against distance to produce the following figure:

![alt text](https://github.com/WalterMuskovic/lncRNA_time_course/blob/master/figures_track/distance_vs_transcription_time.png "distance vs transcription time")

The histogram in the bottom panel summarises the distribution of calculated elongation rates (distance in Kb / transcription time in min).

## Figure S3
This figure demonstrates the spatial correlation trend between protein-coding genes and between lncRNAs. It demonstrates that, similar to the relationship between lncRNAs and protein-coding genes illustrated in Fig. 6, protein-coding genes and lncRNAs that are near each other tend to have more similar gene expression profiles. The code used to generate these figures is essentially the same as in Fig. 6 and can be found in [create_Fig_S3.R](https://github.com/WalterMuskovic/lncRNA_time_course/blob/master/code/create_Fig_S3.R).

## Figure S4
This figure splits the lncRNAs identified in this study into those which are known (overlap annotated lncRNAs) and those which are novel (no overlap with annotated lncRNAs) and examines the spatial correlation trend between their expression and that of protein-coding genes. It demonstrates that both known and novel lncRNAs display the same trend: mirroring the expression of adjacent protein-coding genes. The code used to generate these two figures is very similar to Fig. 6 and can be found in [create_Fig_S4.R](https://github.com/WalterMuskovic/lncRNA_time_course/blob/master/code/create_Fig_S4.R).
