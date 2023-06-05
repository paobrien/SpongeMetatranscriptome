## 03 pathway pre-processing

# Subset and format KO and CAZY annotations by host

setwd("~/Documents/R/Metagenomics/Metabolic_pathways/")

library(tidyverse)
#library(reshape2)

## import and format data ----

# MAGs with metadata
mags <- read.table(
  "~/Documents/R/Metagenomics/Bin_stats/95_ani/final_bins_stats_coverage_das_95ANI_short.tsv", 
  sep = '\t',
  row.names = 1, 
  header = T, 
  strip.white = T)

# remove unmapped data
mags <- mags[-416,]

# get taxonomy
mags <- unite(
  mags, "Taxonomy", Domain:Phylum:Class:Order:Family:Genus:Species, sep = ";", remove = F
  )
tax <- select(mags, "Taxonomy")

# get ko table
ko <- read.table(
  "~/Documents/R/Metagenomics/Data/Enrichm/KO/ko_frequency_table_das_95.tsv", 
  sep = "\t",
  row.names = 1, 
  header = T, 
  strip.white = T, 
  check.names = F)

ko <- as.data.frame(t(ko))


# get cazy table
cazy <- read.table(
  "~/Documents/R/Metagenomics/Data/Enrichm/CAZY/cazy_frequency_table_allMags.tsv", 
  sep = "\t", 
  row.names = 1, 
  header = T, 
  strip.white = T, 
  check.names = F)

cazy <- as.data.frame(t(cazy))

# subset for each host species
car <- select(
  mags, "Taxonomy", "Car1", "Car2", "Car3", "Car4", "Car5"
  )
car <- car[apply(car[,-1], 1, function(x) !all(x==0)),] # keep only rows that have values above zero

irb <- select(
  mags, "Taxonomy", "IrB1", "IrB2", "IrB3", "IrB4", "IrB5"
  )
irb <- irb[apply(irb[,-1], 1, function(x) !all(x==0)),]

irc <- select(
  mags, "Taxonomy", "Irc1", "Irc2", "Irc3", "Irc4", "Irc5", "CS69", "CS70", "CS70b", "CS71", "CS72", "CS73"
  )
irc <- irc[apply(irc[,-1], 1, function(x) !all(x==0)),]

# make ko and cazy data set for each host
# Fist get common mag_id column
ko$mag_id <- row.names(ko)
cazy$mag_id <- row.names(cazy)
mags$mag_id <- row.names(mags)
car$mag_id <- row.names(car)
irb$mag_id <- row.names(irb)
irc$mag_id <- row.names(irc)

# join dataframes
#car
ko_car <- left_join(car, ko, by = "mag_id")
cazy_car <- left_join(car, cazy, by = "mag_id")

#irb
ko_irb <- left_join(irb, ko, by = "mag_id")
cazy_irb <- left_join(irb, cazy, by = "mag_id")

#irc
ko_irc <- left_join(irc, ko, by = "mag_id")
cazy_irc <- left_join(irc, cazy, by = "mag_id")

#all
ko_all <- cbind(tax, ko)
cazy_all <- cbind(tax, cazy)

#combine to list
ko_list <- list(ko_car, ko_irb, ko_irc, ko_all)
cazy_list <- list(cazy_car, cazy_irb, cazy_irc, cazy_all)
names(ko_list) <- c("Car", "Irb", "Irc", "All")
names(cazy_list) <- c("Car", "Irb", "Irc", "All")

## Save ----
saveRDS(ko_list, file = "ko_by_host_metagenome_list")
saveRDS(cazy_list, file = "cazy_by_host_metagenome_list")
ko_list <- readRDS("ko_by_host_metagenome_list")
cazy_list <- readRDS("cazy_by_host_metagenome_list")
names(ko_list)
names(cazt_list)











