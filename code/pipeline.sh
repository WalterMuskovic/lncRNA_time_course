#!/bin/bash

#*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/
#         File: pipeline.sh
#               Execute lncRNA time series analysis
#       Author: Walter Muskovic
#  Modified on: 2021/01/02
#      Version: 1.0.0
#      Example: ./pipeline.sh
#               This script will reproduce the analysis and figures 
#               presented in the manuscript entitled "No evidence for lncRNA 
#               cis-regulatory roles from high temporal resolution RNA-seq 
#               time-course data”.
# -*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*/



# Define directories
export code_dir=lncRNA_time_course/code
export data_dir=lncRNA_time_course/data
export data_track_dir=lncRNA_time_course/data_track
export figure_dir=lncRNA_time_course/figures



### Data pre-processing

## Get human and mouse sequence/annotation data
# We start by downloading the human genome sequence (GRCh38) and RNA 
# [sequin](https://doi.org/10.1038/nmeth.3958) spike-in sequences
bash $code_dir/download_human_seq.sh
# We also need GRCh38 annotation and the RNA spike-in sequins annotation files
bash $code_dir/download_human_annotation.sh
# We will also need the mouse (GRCm38) genome sequence and gene annotation
bash $code_dir/download_mouse_seq_annotation.sh

## Prepare and align T98G FASTQ files
# Check md5sums and read depth of FASTQ files (stored in $data_dir/fastq - 
# available from GSE138662)
bash $code_dir/check_T98G_checksums.sh
# Summarise the read depth for the T98G time point samples
Rscript $code_dir/T98G_read_depth.R
# Trim Illumina adapters with cutadapt
bash $code_dir/trim_T98G_adapters.sh
# Prepare STAR GRCh38 genome index
bash $code_dir/prepare_human_STAR_index.sh
# Align T98G RNA-seq reads
bash $code_dir/align_T98G.sh

## Download and process mouse LPS-response time course data
bash $code_dir/process_GSE56977.sh
# Summarise the read depth for the mouse LPS time point samples
Rscript $code_dir/GSE56977_read_depth.R

## Create annotation
bash $code_dir/stringtie.sh
Rscript $code_dir/intron_annotations.R

## Quantify gene expression
# Using the BAM files created with STAR, quantify gene expression using 
# featureCounts from the Rsubread R package
Rscript $code_dir/featureCounts_quantify.R
Rscript $code_dir/featureCounts_quantify_mouse.R

## Combine count data with transcript coordinates
# Identify DE genes and normalize count data
Rscript $code_dir/de_norm.R
Rscript $code_dir/de_norm_mouse.R
# Add genomic coordinates
Rscript $code_dir/get_genomic_coords.R
Rscript $code_dir/get_genomic_coords_mouse.R

## Download cis reg element data from ENCODE
Rscript $code_dir/encode_cis_reg.R



### Fig. 1 - mRNA and lncRNA expression
Rscript $code_dir/create_data_for_Fig_1.R
Rscript $code_dir/create_Fig_1.R



### Figure 2 - The effects of transcript stability
Rscript $code_dir/create_data_for_Fig_2.R
Rscript $code_dir/create_Fig_2.R



### Figure 3 - The effects of gene length on transcription time
# Panel a - transcription acrosss CACNA1C
Rscript $code_dir/create_data_for_Fig_3a.R
Rscript $code_dir/create_Fig_3a.R
# Gene schematics
Rscript $code_dir/create_Fig_3_gene_schematics.R
# Panels b-e
Rscript $code_dir/create_data_for_Fig_3b_e.R
Rscript $code_dir/create_Fig_3b_e.R



### Figure 4 - Comparison of lncRNA dynamics with protein-coding gene pre-mRNA 
# dynamics
Rscript $code_dir/create_Fig_4.R



### Figure 5 - Examples of human coding gene and adjacent lncRNA expression
Rscript $code_dir/create_data_for_Fig_5.R
Rscript $code_dir/create_Fig_5.R



### Figure 6
Rscript $code_dir/create_data_for_Fig_6.R
Rscript $code_dir/bootstrap_schematic.R
Rscript $code_dir/lag_example.R
Rscript $code_dir/create_Fig_6.R



### Figure 7
Rscript $code_dir/create_data_for_Fig_7_updated.R
Rscript $code_dir/create_Fig_7_updated.R


### Supplementary figures
## Figure S1
Rscript $code_dir/create_Fig_S1.R
## Figure S2
Rscript $code_dir/create_Fig_S2.R
## Figure S3
Rscript $code_dir/create_Fig_S3.R
## Figure S4
Rscript $code_dir/create_Fig_S4.R
## Figure S5
Rscript $code_dir/create_Fig_S5.R
