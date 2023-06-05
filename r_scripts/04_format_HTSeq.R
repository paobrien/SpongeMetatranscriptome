## 03 format HTseq

# Format results from HTseq metatranscrtipome count data

setwd("~/Documents/R/Metatranscriptomics/HTSeq")

library(tidyverse)
library(edgeR)

## Import files ----

## Import data
file_list_c <- list.files(path = "Data", pattern = "^car")
file_list_i <- list.files(path = "Data", pattern = "^irb")
file_list_r <- list.files(path = "Data", pattern = "^irc")

# HTseq file
car <- list()
for(i in seq_along(file_list_c)) {
  car[[i]] <- read.table(paste0("Data/", file_list_c[[i]]), 
                         sep = "\t", header = TRUE)
}

irb <- list()
for(i in seq_along(file_list_i)) {
  irb[[i]] <- read.table(paste0("Data/", file_list_i[[i]]), 
                         sep = "\t", header = TRUE)
}

irc <- list()
for(i in seq_along(file_list_r)) {
  irc[[i]] <- read.table(paste0("Data/", file_list_r[[i]]),
                         sep = "\t", header = TRUE)
}

# GFF files
file_list_gff <- list.files(path = "Data", pattern = "^combined")
sponge <- c("car", "irb", "irc")

gff <- list()
for (i in seq_along(file_list_gff)) {
  gff[[i]] <- read.table(paste0("Data/", file_list_gff[i]), 
                         header = F, 
                         strip.white = T)
}
names(gff) <- sponge

# MAG -> Contig IDs
mag_ids_c <- read.table(
  "Data/mag_contig_carterio.tsv", sep = "\t", header = TRUE
  )

mag_ids_i <- read.table(
  "Data/mag_contig_iblue.tsv", sep = "\t", header = TRUE
  )

mag_ids_r <- read.table(
  "Data/mag_contig_iramosa.tsv", sep = "\t", header = TRUE
  )


# MAG data
mag_data <- read.table(
  "~/Documents/R/Metagenomics/Bin_stats/All_mags/final_bins_stats_all_mags.tsv",
  sep = "\t", 
  header = TRUE)

## format data ----

# get gene lengths (from contig coordinates)
for(i in seq_along(gff)) {
  gff[[i]]$gene_length <- gff[[i]]$V5 - gff[[i]]$V4 
}

# get contig id
# first parse columns
for(i in seq_along(gff)) {
  gff[[i]] <- separate(gff[[i]], 
                       col = V9, 
                       into = c("V10", "V11", "V12", "V13",
                                "V14", "V15", "V16", "V17"),
                       sep = ";",
                       extra = "warn",
                       fill = "warn")
}

# create id column
for(i in seq_along(gff)) {
  gff[[i]]$Contig <- gsub("seq_id=", "", gff[[i]]$V10)
}

# subset to contig id and gene length
for(i in seq_along(gff)) {
  gff[[i]] <- select(gff[[i]], Contig, gene_length)
}

# add gene length to HTseq file
for(i in seq_along(car)) {
  car[[i]] <- left_join(car[[i]], gff$car, by = "Contig")
  irb[[i]] <- left_join(irb[[i]], gff$irb, by = "Contig")
  irc[[i]] <- left_join(irc[[i]], gff$irc, by = "Contig")
}

# join dataframes
#car
car.df <- full_join(car[[1]], car[[2]], 
                     by = c("Contig", "Annotation", "gene_length"), 
                     suffix = c("_car1", "_car2"))

car.df <- full_join(car.df, car[[3]], 
                     by = c("Contig", "Annotation", "gene_length"))

colnames(car.df)[6] <- "Reads_aligned_car3"
car.df <- select(car.df, Contig, Annotation, gene_length, everything())

#irb
irb.df <- full_join(irb[[1]], irb[[2]], 
                    by = c("Contig", "Annotation", "gene_length"), 
                    suffix = c("_irb1", "_irb2"))

irb.df <- full_join(irb.df, irb[[3]], 
                    by = c("Contig", "Annotation", "gene_length"))

colnames(irb.df)[6] <- "Reads_aligned_irb3"
irb.df <- select(irb.df, Contig, Annotation, gene_length, everything())

#irc
irc.df <- full_join(irc[[1]], irc[[2]], 
                    by = c("Contig", "Annotation", "gene_length"), 
                    suffix = c("_irc1", "_irc2"))

irc.df <- full_join(irc.df, irc[[3]], 
                    by = c("Contig", "Annotation", "gene_length"))

colnames(irc.df)[6] <- "Reads_aligned_irc3"
irc.df <- select(irc.df, Contig, Annotation, gene_length, everything())


## normalise read counts ----

#note: since we aren't looking at DE between samples/libraries, it is
#appropriate to normalise by gene length. Doing this before filtering will 
#give equal chance for a gene to be removed due to low expression levels

#normalise
car.rpkm <- data.frame(sapply(select(car.df, starts_with("Reads_aligned")), 
                              function(column) 10^9 * column / car.df$gene_length / sum(column)))

irb.rpkm <- data.frame(sapply(select(irb.df, starts_with("Reads_aligned")), 
                              function(column) 10^9 * column / irb.df$gene_length / sum(column)))

irc.rpkm <- data.frame(sapply(select(irc.df, starts_with("Reads_aligned")), 
                              function(column) 10^9 * column / irc.df$gene_length / sum(column)))

#add contig to join df
car.rpkm <- cbind(select(car.df, Contig), car.rpkm)
irb.rpkm <- cbind(select(irb.df, Contig), irb.rpkm)
irc.rpkm <- cbind(select(irc.df, Contig), irc.rpkm)

# join with raw counts
car.all <- full_join(car.df, car.rpkm, 
                     by ="Contig", 
                     suffix = c(".raw", ".rpkm"))

irb.all <- full_join(irb.df, irb.rpkm, 
                     by ="Contig", 
                     suffix = c(".raw", ".rpkm"))

irc.all <- full_join(irc.df, irc.rpkm, 
                     by ="Contig", 
                     suffix = c(".raw", ".rpkm"))

#get count totals
car.all$total <- rowSums(select(car.all, ends_with(".raw")))
irb.all$total <- rowSums(select(irb.all, ends_with(".raw")))
irc.all$total <- rowSums(select(irc.all, ends_with(".raw")))

## filter low abundnace reads ----

# sort by total aligned reads
car.all <- car.all[order(car.all$total, decreasing = T),]
irb.all <- irb.all[order(irb.all$total, decreasing = T),]
irc.all <- irc.all[order(irc.all$total, decreasing = T),]

# check count distribution
tmp <- car.all%>% 
  filter(total < 100) %>% 
  filter(total > 1)

ggplot(tmp, aes(x = reorder(Contig, -Reads_aligned_car2.rpkm), 
                y = Reads_aligned_car2.rpkm)) +
  geom_bar(stat = 'identity', colour = "black", fill = "black") +
  scale_y_continuous(breaks = seq(0, 0.9, 0.05)) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

# remove genes with <0.1 rpkm
# first change <0.1 to 0 values
car.tmp <- data.frame(sapply(select(car.all, ends_with(".rpkm")),
                        function(column) ifelse(column < 0.1, 0, column)))

irb.tmp <- data.frame(sapply(select(irb.all, ends_with(".rpkm")),
                        function(column) ifelse(column < 0.1, 0, column)))

irc.tmp <- data.frame(sapply(select(irc.all, ends_with(".rpkm")),
                        function(column) ifelse(column < 0.1, 0, column)))

#add contig column
car.tmp <- cbind(select(car.all, Contig), car.tmp)
irb.tmp <- cbind(select(irb.all, Contig), irb.tmp)
irc.tmp <- cbind(select(irc.all, Contig), irc.tmp)

# filter zero values
# note: also removes unaligned, ambigous etc reads (which we want to do anyway)
car.tmp <- filter(car.tmp, rowSums(select(car.tmp, ends_with(".rpkm"))) > 0)
irb.tmp <- filter(irb.tmp, rowSums(select(irb.tmp, ends_with(".rpkm"))) > 0)
irc.tmp <- filter(irc.tmp, rowSums(select(irc.tmp, ends_with(".rpkm"))) > 0)

# join with dataframe
car.all <- right_join(select(
  car.all, !c(Reads_aligned_car1.rpkm, Reads_aligned_car2.rpkm, Reads_aligned_car3.rpkm)),
  car.tmp, by = "Contig")

irb.all <- right_join(select(
  irb.all, !c(Reads_aligned_irb1.rpkm, Reads_aligned_irb2.rpkm, Reads_aligned_irb3.rpkm)),
  irb.tmp, by = "Contig")

irc.all <- right_join(select(
  irc.all, !c(Reads_aligned_irc1.rpkm, Reads_aligned_irc2.rpkm, Reads_aligned_irc3.rpkm)),
  irc.tmp, by = "Contig")

## Add mag IDs
# Car
# Edit contig column
car.all$Contig <- sub("_[^_]+$", "", car.all$Contig) # remove _* at the end of the string
irb.all$Contig <- sub("_[^_]+$", "", irb.all$Contig)
irc.all$Contig <- sub("_[^_]+$", "", irc.all$Contig)

# Join dataframes
car.all <- left_join(car.all, mag_ids_c, by = "Contig")
irb.all <- left_join(irb.all, mag_ids_i, by = "Contig")
irc.all <- left_join(irc.all, mag_ids_r, by = "Contig")

## Add mag taxonomy

# select taxa columns
mag_data <- select(
  mag_data, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  )
mag_data$MAG_id <- rownames(mag_data)

# Car
car.all <- left_join(car.all, mag_data, by = "MAG_id")

# Irb
irb.all <- left_join(irb.all, mag_data, by = "MAG_id")

# Irc
irc.all <- left_join(irc.all, mag_data, by = "MAG_id")


## save data ----

combined_data <- list(car.all, irb.all, irc.all)
names(combined_data) <- c("car", "irb", "irc")

saveRDS(combined_data, file = "sponge_combined_expression_data_filtered_20220923")

# load data
combined_data <- readRDS("sponge_combined_expression_data_filtered_20220923")


## sort data and get alignment values ----

# sort by reads aligned
for (i in seq_along(car)) {
  car[[i]] <- car[[i]][order(car[[i]]$Reads_aligned, decreasing = T),]
  print(car[[i]] %>% head)
}

for (i in seq_along(irb)) {
  irb[[i]] <- irb[[i]][order(irb[[i]]$Reads_aligned, decreasing = T),]
  print(irb[[i]] %>% head)
}

for (i in seq_along(irc)) {
  irc[[i]] <- irc[[i]][order(irc[[i]]$Reads_aligned, decreasing = T),]
  print(irc[[i]] %>% head)
}


# get proportion of reads aligned
for (i in seq_along(car)) {
  car[[i]]$Proportion <- car[[i]]$Reads_aligned / sum(car[[i]]$Reads_aligned)
  print(car[[i]][1:10,])
}

for (i in seq_along(irb)) {
  irb[[i]]$Proportion <- irb[[i]]$Reads_aligned / sum(irb[[i]]$Reads_aligned)
  print(irb[[i]][1:10,])
}

for (i in seq_along(irc)) {
  irc[[i]]$Proportion <- irc[[i]]$Reads_aligned / sum(irc[[i]]$Reads_aligned)
  print(irc[[i]][1:10,])
}

# check mean proportion of unaligned reads
tmp_c <- list()
tmp_i <- list()
tmp_r <- list()
for (i in seq_along(car)) {
  tmp_c[[i]] <- car[[i]]$Proportion[1]
  tmp_i[[i]] <- irb[[i]]$Proportion[1]
  tmp_r[[i]] <- irc[[i]]$Proportion[1]
}

tmp_c <- do.call("c", tmp_c)
tmp_i <- do.call("c", tmp_i)
tmp_r <- do.call("c", tmp_r)

stats_c <- c()
stats_c$mean <- mean(tmp_c)
stats_c$std_er <- sd(tmp_c)/sqrt(length(tmp_c))
stats_c

stats_i <- c()
stats_i$mean <- mean(tmp_i)
stats_i$std_er <- sd(tmp_i)/sqrt(length(tmp_i))
stats_i

stats_r <- c()
stats_r$mean <- mean(tmp_r)
stats_r$std_er <- sd(tmp_r)/sqrt(length(tmp_r))
stats_r

# check mean proportion of ambiguous aligned reads
tmp_c <- list()
tmp_i <- list()
tmp_r <- list()
for (i in seq_along(car)) {
  tmp_c[[i]] <- car[[i]]$Proportion[2]
  tmp_r[[i]] <- irc[[i]]$Proportion[2]
  tmp_i[[i]] <- irb[[i]]$Proportion[2]
}

tmp_c <- do.call("c", tmp_c)
tmp_i <- do.call("c", tmp_i)
tmp_r <- do.call("c", tmp_r)

stats_c <- c()
stats_c$mean <- mean(tmp_c)
stats_c$std_er <- sd(tmp_c)/sqrt(length(tmp_c))
stats_c

stats_i <- c()
stats_i$mean <- mean(tmp_i)
stats_i$std_er <- sd(tmp_i)/sqrt(length(tmp_i))
stats_i

stats_r <- c()
stats_r$mean <- mean(tmp_r)
stats_r$std_er <- sd(tmp_r)/sqrt(length(tmp_r))
stats_r

# check number of aligned reads
for (i in seq_along(car)) {
  print(sum(car[[i]]$Reads_aligned) - car[[i]]$Reads_aligned[1])
}

for (i in seq_along(irb)) {
  print(sum(irb[[i]]$Reads_aligned) - irb[[i]]$Reads_aligned[1])
}

for (i in seq_along(irc)) {
  print(sum(irc[[i]]$Reads_aligned) - irc[[i]]$Reads_aligned[1])
}
