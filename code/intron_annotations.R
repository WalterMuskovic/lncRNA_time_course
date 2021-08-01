# Load libraries ----------------------------------------------------------
library(tidyverse)
library(rtracklayer)
library(GenomicFeatures)



# Stringtie annotation information ----------------------------------------

## Some information on the annotation format, provided byt he stringtie authors

#The primary output of StringTie is a Gene Transfer Format (GTF) file that
#contains details of the transcripts that StringTie assembles from RNA-Seq data.

# Description of each column's values:

#seqname: Denotes the chromosome, contig, or scaffold for this transcript. 

#source: The source of the GTF file. This column simply shows 'StringTie'.

#feature: Feature type; e.g., exon, transcript, mRNA, 5'UTR.

#start: Start position of the feature (exon, transcript, etc), using a 1-based
#index.

#end: End position of the feature, using a 1-based index.

#score: A confidence score for the assembled transcript. A constant value of
#1000 iis used by stringtie.

#strand: If the transcript resides on the forward strand, '+'. If the transcript
#resides on the reverse strand, '-'.

#frame: Frame or phase of CDS features. StringTie does not use this field and
#simply records a ".".

#attributes: A semicolon-separated list of tag-value pairs, providing additional
#information about each feature. Depending on whether an instance is a
#transcript or an exon and on whether the transcript matches the reference
#annotation file provided by the user, the content of the attributes field will
#differ. The following list describes the possible attributes shown in this
#column:

  #gene_id: A unique identifier for a single gene and its child transcript and
  #exons based on the alignments' file name.

  #transcript_id: A unique identifier for a single transcript and its child
  #exons based on the alignments' file name.

  #exon_number: A unique identifier for a single exon, starting from 1, within a
  #given transcript.

  #reference_id: The transcript_id in the reference annotation (optional) that
  #the instance matched.
  
  #ref_gene_id: The gene_id in the reference annotation (optional) that the
  #instance matched.

  #ref_gene_name: The gene_name in the reference annotation (optional) that the
  #instance matched.

  #cov: The average per-base coverage for the transcript or exon.

  #FPKM: Fragments per kilobase of transcript per million read pairs. This is
  #the number of pairs of reads aligning to this feature, normalized by the
  #total number of fragments sequenced (in millions) and the length of the
  #transcript (in kilobases).

  #TPM: Transcripts per million. This is the number of transcripts from this
  #particular gene normalized first by gene length, and then by sequencing depth
  #(in millions) in the sample.



# Select transcript -------------------------------------------------------

## Choose one transcript per gene - the tx with highest TPM. This produces a
## very small number of ties, in which case, we take the longest tx.
gtf <- readGFF("data/annotation/T98G_stringtie_stranded.gtf")
gtf_mouse <- readGFF("data/annotation/Rabani_stringtie_stranded.gtf")

tx <- gtf %>%
  mutate(TPM = as.numeric(TPM)) %>%
  filter(type=="transcript") %>%
  group_by(gene_id) %>%
  filter(TPM==max(TPM)) %>%
  mutate(tx_length = end-start) %>%
  filter(tx_length==max(tx_length)) %>%
  pull(transcript_id)

tx_mouse <- gtf_mouse %>%
  mutate(TPM = as.numeric(TPM)) %>%
  filter(type=="transcript") %>%
  group_by(gene_id) %>%
  filter(TPM==max(TPM)) %>%
  mutate(tx_length = end-start) %>%
  filter(tx_length==max(tx_length)) %>%
  pull(transcript_id)



# Intron annotation files -------------------------------------------------

## Create three intron annotation files for the selected transcripts containing;
  # All introns
  # The first 10 Kb of intron
  # The last 10 Kb of intron

txdb <- makeTxDbFromGFF("data/annotation/T98G_stringtie_stranded.gtf")
introns <- intronsByTranscript(txdb,
                               use.names = TRUE)

txdb_mouse <- makeTxDbFromGFF("data/annotation/Rabani_stringtie_stranded.gtf")
introns_mouse <- intronsByTranscript(txdb_mouse,
                               use.names = TRUE)

# Reduce to selected transcripts
introns <- introns[tx]
introns_mouse <- introns_mouse[tx_mouse]

# Remove transcripts that have no introns
introns <- introns[lengths(introns)>0]
introns_mouse <- introns_mouse[lengths(introns_mouse)>0]

# Define function to return 1bp GRanges tiles covering an input GRanges
tile_tx <- function(input_tx){
  # Get GRangesList with 1bp tiles covering each of the input ranges
  intron_tiles <- GenomicRanges::tile(x = input_tx, width = 1)
  # GRangesList >> GRanges
  intron_tiles <- unlist(intron_tiles)
  return(intron_tiles)
  }

# Define function to return the first 10 Kb of input GRanges
  # If the first 10 Kb is needed, set strand="+"
  # If the first 10 Kb is needed, set strand="-"
first_last_10Kb <- function(input_ranges, input_strand){
  # If number of 1bp intron ranges is less than 10Kb, just return what was input
  if(length(input_ranges)<=1E4){return(reduce(input_ranges))}
  # Deal with strand
  if(as.character(strand((input_ranges[1])))==input_strand){
    return(reduce(input_ranges[1:1E4]))
  } else {
    length_ranges <- length(input_ranges)
    return(reduce(input_ranges[(length_ranges-9999):length_ranges]))
  }
  }

# Get first 10Kb of intron intervals
introns_first_10kb <-
  lapply(seq_along(introns), function(x)
    first_last_10Kb(tile_tx(introns[[x]]), input_strand="+"))
names(introns_first_10kb) <- names(introns)
introns_first_10kb <- GRangesList(introns_first_10kb)

introns_first_10kb_mouse <-
  lapply(seq_along(introns_mouse), function(x)
    first_last_10Kb(tile_tx(introns_mouse[[x]]), input_strand="+"))
names(introns_first_10kb_mouse) <- names(introns_mouse)
introns_first_10kb_mouse <- GRangesList(introns_first_10kb_mouse)

# Get last 10Kb of intron intervals
introns_last_10kb <-
  lapply(seq_along(introns), function(x)
    first_last_10Kb(tile_tx(introns[[x]]), input_strand="-"))
names(introns_last_10kb) <- names(introns)
introns_last_10kb <- GRangesList(introns_last_10kb)

introns_last_10kb_mouse <-
  lapply(seq_along(introns_mouse), function(x)
    first_last_10Kb(tile_tx(introns_mouse[[x]]), input_strand="-"))
names(introns_last_10kb_mouse) <- names(introns_mouse)
introns_last_10kb_mouse <- GRangesList(introns_last_10kb_mouse)

# Export the three gtf files with the appropriate intron ranges
export.gff3(object = introns,
           con = "data/annotation/introns.gff3")
export.gff3(object = introns_first_10kb,
            con = "data/annotation/introns_first_10kb.gff3")
export.gff3(object = introns_last_10kb,
            con = "data/annotation/introns_last_10kb.gff3")

export.gff3(object = introns_mouse,
            con = "data/annotation/introns_mouse.gff3")
export.gff3(object = introns_first_10kb_mouse,
            con = "data/annotation/introns_first_10kb_mouse.gff3")
export.gff3(object = introns_last_10kb_mouse,
            con = "data/annotation/introns_last_10kb_mouse.gff3")


# Whole transcript interval (exons + introns) -----------------------------

whole_tx <- unlist(transcriptsBy(txdb))
whole_tx <- whole_tx[whole_tx$tx_name%in%tx]
whole_tx <- sort(whole_tx)
export.gff3(whole_tx, "data/annotation/whole_tx.gff3")

whole_tx_mouse <- unlist(transcriptsBy(txdb_mouse))
whole_tx_mouse <- whole_tx_mouse[whole_tx_mouse$tx_name%in%tx_mouse]
whole_tx_mouse <- sort(whole_tx_mouse)
export.gff3(whole_tx_mouse, "data/annotation/whole_tx_mouse.gff3")
