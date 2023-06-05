## 01 Format metagenomics data 

# create dataframe of taxonomy, coverage, genome stats

setwd("~/Documents/R/Metagenomics") 

library(tidyverse)

## Import and format metagenome data ----

## import GTDB classification
# bacteria
tax_bac <- read.table(
  "Data/GTDB/dastool_95.bac120.summary.tsv", 
  sep = '\t', 
  row.names = 1, 
  header = T, 
  strip.white = T)

# archaea
tax_arc <- read.table(
  "Data/GTDB/dastool_95.ar122.summary.tsv", 
  sep = '\t', 
  row.names = 1, 
  header = T, 
  strip.white = T)

# join dfs
tax <- rbind(tax_arc, tax_bac)
dim(tax)

## import checkM results
checkm <- read.table(
  "Data/Checkm/checkm_derep_das_95_long.tsv", 
  sep = '\t', 
  row.names = 1, 
  header = T, 
  strip.white = T)


## import coverM table
coverm <- read.table(
  "Data/Coverm/all_samples_genome_coverage_dastool_95ANI.txt", 
  sep = '\t', 
  row.names = 1, 
  header = T, 
  strip.white = T)

## combine GTDB with CheckM
bin_stats <- cbind(tax, checkm)
bin_stats <- tax # if not cbinding with checkm

# separate tax column
bin_stats <- separate(
  bin_stats, 
  classification, 
  c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), 
  sep = ";", 
  remove = TRUE, 
  convert = FALSE, 
  extra = "warn", 
  fill = "warn")

# save
write.table(
  bin_stats, 
  "/Bin_stats/95_ani/final_bins_stats_das_95ANI.tsv", 
  sep = "\t")

# merge bin stats with coverage
coverm$Bin_id <- row.names(coverm)
bin_stats$Bin_id <- row.names(bin_stats)
df <- right_join(bin_stats, coverm)

# move bin_id column to the front
df <- df %>% select(Bin_id, everything())

# save
write.table(
  df, 
  "/Bin_stats/95_ani/final_bins_stats_coverage_das_95ANI.tsv", 
  sep = "\t")





