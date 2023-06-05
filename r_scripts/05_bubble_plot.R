## 04_bubble_plot

# Plot metatranscriptome expression by taxonomy and compare to metagenome relative abundnace

setwd("~/Documents/R/Metatranscriptomics")

library(tidyverse)
library(reshape2)
library(ggpubr)

## import data ----

# metatranscriptome
meta.t <- readRDS(
  "HTSeq/Output/sponge_combined_expression_data_filtered_20220923"
  )

# metagenome
meta.g <- read.table(
  "../Metagenomics/Bin_stats/95_ani/final_bins_stats_coverage_das_95ANI_short.tsv", 
  sep = '\t',
  row.names = 1, 
  header = T, 
  strip.white = T)

## format data ----

## metatranscriptome

# subset data
car_df <- select(meta.t[[1]], 
                 Reads_aligned_car1.raw,
                 Reads_aligned_car2.raw, 
                 Reads_aligned_car3.raw, 
                 Phylum)
irb_df <- select(meta.t[[2]], 
                 Reads_aligned_irb1.raw,
                 Reads_aligned_irb2.raw,
                 Reads_aligned_irb3.raw,
                 Phylum)
irc_df <- select(meta.t[[3]], 
                 Reads_aligned_irc1.raw, 
                 Reads_aligned_irc2.raw, 
                 Reads_aligned_irc3.raw, 
                 Phylum)

# get proportion
# car
reads_cols <- names(car_df)[1:3]
car_df[reads_cols] <- lapply(car_df[reads_cols], function(x) x / sum(x))

# irb
reads_cols <- names(irb_df)[1:3]
irb_df[reads_cols] <- lapply(irb_df[reads_cols], function(x) x / sum(x))

# irc
reads_cols <- names(irc_df)[1:3]
irc_df[reads_cols] <- lapply(irc_df[reads_cols], function(x) x / sum(x))

# change to long format
car_df <- melt(car_df)
irb_df <- melt(irb_df)
irc_df <- melt(irc_df)

# add sponge group
car_df$Sponge <- rep("P. foliascens", nrow(car_df))
irb_df$Sponge <- rep("I. microconulosa", nrow(irb_df))
irc_df$Sponge <- rep("I. ramosa", nrow(irc_df))

# join together
df <- rbind(car_df, irb_df, irc_df)

# clean up column values
df$Phylum <- gsub("p__", "", df$Phylum)
df$variable <- gsub("Reads_aligned_", "", df$variable)
df$variable <- gsub("\\.raw", "", df$variable)

# change column names
names(df)[c(2,3)] <- c("Sample", "Proportion_aligned")

# edit prop
df$Proportion_aligned <- df$Proportion_aligned * 100

# group by phylum
df_sum <- group_by(df, Sample, Phylum, Sponge) %>% 
  summarise(Proportion_aligned = sum(Proportion_aligned, na.rm = T))

## metagenome

# remove unmapped data
meta.g <- meta.g[-416,]

# subset data
meta.g <- select(
  meta.g, Phylum, starts_with(c("car", "cs", "ir"), ignore.case = T)
  )

# edt colnames
names(meta.g)[2:6] <- paste0("Pfol", rep(1:5))
names(meta.g)[7:12] <- paste0("Irc", rep(6:11))

# change to long format
meta.g <- melt(meta.g)

# add sponge group
meta.g$Sponge <- ifelse(grepl("Pfol.*", meta.g$variable), "P. foliascens",
                        ifelse(grepl("IrB.*", meta.g$variable), "I. microconulosa",
                               ifelse(grepl("Irc.*", meta.g$variable), "I. ramosa", "no_match")))

# clean up column values
meta.g$Phylum <- gsub("p__", "", meta.g$Phylum)

# change column names
names(meta.g)[c(2,3)] <- c("Sample", "Relative abundance")

# group by phylum
meta.g_sum <- group_by(meta.g, Sample, Phylum, Sponge) %>% 
  summarise(`Relative abundance` = sum(`Relative abundance`, na.rm = T))


## plot ----

# edit sample labels
# metatranscriptome
df_sum$Sample <- gsub("car", "Pfol", df_sum$Sample)
df_sum$Sample <- gsub("irb", "Imic", df_sum$Sample)
df_sum$Sample <- gsub("irc", "Iram", df_sum$Sample)

# metagenome
meta.g_sum$Sample <- gsub("Irc", "Iram", meta.g_sum$Sample)
meta.g_sum$Sample <- gsub("IrB", "Imic", meta.g_sum$Sample)

# relevel sampple order (metagenome)
meta.g_sum$Sample <- factor(meta.g_sum$Sample, levels = c("Pfol1", "Pfol2",  "Pfol3",
                                                  "Pfol4", "Pfol5", "Iram1",
                                                  "Iram2", "Iram3", "Iram4",
                                                  "Iram5", "Iram6", "Iram7", 
                                                  "Iram8", "Iram9", "Iram10", 
                                                  "Iram11", "Imic1", "Imic2",
                                                  "Imic3",  "Imic4", "Imic5"))
meta.g_sum$Sample %>% levels()

# metatranscriptome
p <- ggplot(subset(df_sum, Proportion_aligned > 0), aes(x = Sample, y = Phylum)) +
  geom_point(aes(size = Proportion_aligned, fill = Sponge), shape = 21) + 
  scale_size(range = range(df_sum$Proportion_aligned)) +
  labs(size = "Relative expression \n(Metatranscriptome)") +
  scale_size_area(max_size = 15) +
  scale_fill_manual(values = c("#0079d9", "#dcba44", "grey")) +
  guides(fill = "none") +
  theme_bw() +
  theme(#axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text = element_text(size=13), 
        legend.title=element_text(size=15))
p <- p + facet_wrap("Sponge", ncol = 3, scales = "free_x") +
  theme(strip.text = element_text(size = 15, face = "italic"))

# metagenome
p2 <- ggplot(subset(meta.g_sum, `Relative abundance` > 0), aes(x = Sample, y = Phylum)) +
  geom_point(aes(size = `Relative abundance`, fill = Sponge), 
             shape = 21) + 
  scale_size(range = range(meta.g_sum$`Relative abundance`)) +
  labs(size = "Relative abundance \n(Metagenome)") +
  scale_size_area(max_size = 15) +
  scale_fill_manual(values = c("#0079d9", "#dcba44", "grey")) + 
  guides(fill = "none") +
  theme_bw() +
  theme(#axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 15),
        legend.text = element_text(size=13), 
        legend.title = element_text(size=15))
p2 <- p2 + facet_wrap("Sponge", ncol = 3, scales = "free_x") +
  theme(strip.text = element_text(size = 15, face = "italic"))


# join into one plot
ggarrange(p + 
            theme(axis.title.x = element_blank(),
                  title = element_text(size = 15)), 
          p2 + 
            theme(title = element_text(size = 15)), 
          heights = c(0.82, 0.9, 1.0),
          nrow = 2, ncol = 1,  common.legend = F, legend = "right")

# save
ggsave("metagenome_metatranscriptome_abundnace.svg", device = "svg", 
       width = 30, height = 35, units = "cm", dpi = 300)


