
# Load R packages ---------------------------------------------------------

library(tidyverse)
library(AnnotationHub)
library(GenomicFeatures)
library(rtracklayer)
library(vioplot)
library(mgcv)


# Load enhancer/promoter data ---------------------------------------------

if(!file.exists("data/mm10-PLS.bed")){
# Download human and mouse candidate promoter elements from ENCODE from the 
# SCREEN server http://screen.encodeproject.org/

  # Promoters
  download.file(url = "https://api.wenglab.org/screen_v13/fdownloads/cCREs/GRCh38-PLS.bed",
                destfile = "data/GRCh38-PLS.bed")
  download.file(url = "https://api.wenglab.org/screen_v13/fdownloads/cCREs/mm10-PLS.bed",
                destfile = "data/mm10-PLS.bed")
  
  # Enhancers
  download.file(url = "https://api.wenglab.org/screen_v13/fdownloads/cCREs/GRCh38-ELS.bed",
                destfile = "data/GRCh38-ELS.bed")
  download.file(url = "https://api.wenglab.org/screen_v13/fdownloads/cCREs/mm10-ELS.bed",
                destfile = "data/mm10-ELS.bed")
  }

# Import as point coordinate granges
getgr <- function(file.name){
  bed.cols <- c("seqid", "start", "end", "X1", "X2", "type")
  input.ranges <- read_tsv(str_glue("data/{file.name}"), col_names = bed.cols)
  new.coord <- as.integer(round(input.ranges$start + (input.ranges$end - input.ranges$start)/2))
  input.ranges$start <- new.coord
  input.ranges$end <- new.coord
  input.ranges <- makeGRangesFromDataFrame(input.ranges)
  return(input.ranges)
  }

he <- getgr("GRCh38-ELS.bed")
hp <- getgr("GRCh38-PLS.bed")
me <- getgr("mm10-ELS.bed")
mp <- getgr("mm10-PLS.bed")

rm(getgr)



# Load lncRNA TSS data ----------------------------------------------------

h.lnc <- readRDS("data/t98g_filtered_coords_clust.rds")
h.lnc.gr  <- h.lnc %>%
  mutate(start=coords, end=coords) %>%
  makeGRangesFromDataFrame()

m.lnc <- readRDS("data/mouse_filtered_coords.rds")
m.lnc.gr  <- m.lnc %>%
  mutate(start=coords, end=coords) %>%
  makeGRangesFromDataFrame()



# Classify as enhancer- or promoter-associated lncRNA ---------------------

## Human
# Get distance to closest element, label according to closest element, if >1Kb mark as none
h.lnc <- dplyr::select(h.lnc, transcript_id, transcript_type) %>%
  mutate(enhancer_dist = distanceToNearest(x = h.lnc.gr, subject = he)@elementMetadata$distance,
         promoter_dist = distanceToNearest(x = h.lnc.gr, subject = hp)@elementMetadata$distance) %>%
  mutate(reg = ifelse(test = enhancer_dist < promoter_dist,
                      yes = "enhancer",
                      no = "promoter")) %>%
  mutate(reg_dist = ifelse(test = reg=="enhancer",
                           yes = enhancer_dist,
                           no = promoter_dist)) %>%
  mutate(nearby = enhancer_dist<=300 | promoter_dist<=300) %>%
  mutate(reg = ifelse(test = nearby == FALSE,
                      yes = "none",
                      no = reg)) %>%
  dplyr::select(-nearby)
  
## Mouse
m.lnc <- dplyr::select(m.lnc, transcript_id, transcript_type) %>%
  mutate(enhancer_dist = distanceToNearest(x = m.lnc.gr, subject = me)@elementMetadata$distance,
         promoter_dist = distanceToNearest(x = m.lnc.gr, subject = mp)@elementMetadata$distance) %>%
  mutate(reg = ifelse(test = enhancer_dist < promoter_dist,
                      yes = "enhancer",
                      no = "promoter")) %>%
  mutate(reg_dist = ifelse(test = reg=="enhancer",
                           yes = enhancer_dist,
                           no = promoter_dist)) %>%
  mutate(nearby = enhancer_dist<=300 | promoter_dist<=300) %>%
  mutate(reg = ifelse(test = nearby == FALSE,
                      yes = "none",
                      no = reg)) %>%
  dplyr::select(-nearby)

rm(h.lnc.gr, m.lnc.gr, he, hp, me, mp)



# Save out ----------------------------------------------------------------

saveRDS(h.lnc, "data/h.lnc.rds")
saveRDS(m.lnc, "data/m.lnc.rds")
