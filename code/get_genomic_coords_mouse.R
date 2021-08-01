# For our analysis we will need to get the distances between each coding gene
# and lncRNA. The appropriate coordinates to use would be the transcription
# start sites for all of these genes. For genes on the forward strand this will
# be the most 5' coordinate, while for genes on the reverse stand it will be the
# most 3' coordinate.

# Load R packages
library(rtracklayer)
library(tidyverse)


# 1.  Get data ------------------------------------------------------------

# Get a quick summary of the data we have
mouse <- readRDS("data/mouse_filtered.rds")

# Get annotation data
gtf <- readGFF("data/annotation/Rabani_stringtie_stranded.gtf",
               filter=(list(type="transcript"))) %>%
  filter(transcript_id %in% mouse$transcript_id)




# 2. Get coordinate -------------------------------------------------------

# Calculate appropriate coordinate according to the strand
coords <- rep(NA, nrow(gtf))
# Transcripts on the positive strand
coords[gtf$strand=="+"] <- gtf$start[gtf$strand=="+"]
# Transcripts on the negative strand
coords[gtf$strand=="-"] <- gtf$end[gtf$strand=="-"]

# Change to data frame
coords <- data.frame(seqid = gtf$seqid,
                     coords,
                     transcript_id = gtf$transcript_id)
# What does this look like?
head(coords)



# 3. Combine and save -----------------------------------------------------

# Add chromosome and coordinate to ccounts 
mouse <- merge(x = coords, y = mouse, by = "transcript_id") %>%
  arrange(seqid, coords)

# Restrict to reference chromosomes
mouse <- filter(mouse, mouse$seqid%in%c(str_glue('chr{1:19}'), "chrY", "chrX"))

# Save out
saveRDS(mouse, "data/mouse_filtered_coords.rds")

# Clean up
rm(coords, gtf, mouse)
