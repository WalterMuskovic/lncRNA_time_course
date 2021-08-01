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
t98g <- readRDS("data/t98g_filtered.rds")

# Get annotation data
gtf <- readGFF("data/annotation/T98G_stringtie_stranded.gtf",
               filter=(list(type="transcript"))) %>%
  filter(transcript_id %in% t98g$transcript_id)




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
t98g <- merge(x = coords, y = t98g, by = "transcript_id") %>%
  arrange(seqid, coords)

# Restrict to reference chromosomes
t98g <- filter(t98g, t98g$seqid%in%c(str_glue('chr{1:22}'), "chrY", "chrX"))

# Save out
saveRDS(t98g, "data/t98g_filtered_coords.rds")

# Clean up
rm(coords, gtf, t98g)
