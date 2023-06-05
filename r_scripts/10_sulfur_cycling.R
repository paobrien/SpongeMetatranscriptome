## 10 Sulfur cycling

# Genome-centric metatranscriptomic analysis of sulfur cycling

setwd("~/Documents/R/Metatranscriptomics/Metabolic_pathways")

library(tidyverse)
library(ape)
library(plyr)
library(reshape2)

## Load functions
source("~/Documents/R/Metatranscriptomics/Scripts/R_functions_metatranscriptomics.R")

## Import files ----

# KO
ko_list <- readRDS(
  "~/Documents/R/Metagenomics/Metabolic_pathways/ko_by_host_metagenome_list"
  )
ko_list <- ko_list[c("Car", "Irb", "Irc")]

# Expression data - HTseq
expression_data <- readRDS(
  "~/Documents/R/Metatranscriptomics/HTSeq/Output/sponge_combined_expression_data_filtered_20220923"
  )

# Tree
bact <- ape::read.tree(
  "~/Documents/R/Metagenomics/Data/GTDB/all_mags.bac120.classify.tree"
  )
arch <- ape::read.tree(
  "~/Documents/R/Metagenomics/Data/GTDB/all_mags.ar122.classify.tree"
  )
full_tree <- bind.tree(bact, arch)

## Make KO dataframe ----

# subset to sulfur cycling KOs
sulf <- lapply(ko_list, function(x) select(x, "mag_id", "Taxonomy", "K17222", 
                                           "K17223", "K17224", "K17225", 
                                           "K17226", "K17227", "K00958", 
                                           "K00394", "K00395", "K11180", 
                                           "K11181", "K15551", "K15552", 
                                           "K10831", "K03119"))
# get ko list
ko <- get_ko(sulf$Car)

# create vector of sulfur cycling process
process <- c(rep("Thiosulfate_oxidation", 6), 
             rep("Sulfate_Reduction", 5), 
             rep("Taurine_metabolsim", 4))

# combine data
sulf_cyc <- cbind(ko, process) %>% as.data.frame

# change to presence absence
sulf <- lapply(sulf, pres_abs)


## Format expression data ----

# Subset data and combine reads column
#names(expression_data) <- c("Car", "Irb", "Irb")
#expression_data <- lapply(expression_data, subset_combine)

# Subset expression data to sulfur KOs
exp_sulf <- lapply(expression_data, function(x) exp_pathway(ko, x))

# Clean up annotation column
for(i in seq_along(exp_sulf)) {
  exp_sulf[[i]]$Annotation <- sub("\\,.*", "", exp_sulf[[i]]$Annotation)
  exp_sulf[[i]]$Annotation <- sub("K02050&", "", exp_sulf[[i]]$Annotation)
}

# Combine duplicate rows (after removing trailing annotations)
for(i in seq_along(exp_sulf)) {
  exp_sulf[[i]] <- ddply(exp_sulf[[i]], 
                         c("MAG_id", "Annotation"),numcolwise(sum))
}

# Convert to wide format
for(i in seq_along(exp_sulf)) {
  exp_sulf[[i]] <- dcast(exp_sulf[[i]], MAG_id ~ Annotation, 
                         value.var ="total")
  exp_sulf[[i]][is.na(exp_sulf[[i]])] <- 0
}

# Join the two dataframes
tmp <- list()
for(i in seq_along(exp_sulf)) {
  names(exp_sulf[[i]])[1] <- "mag_id"
  tmp[[i]] <- left_join(sulf[[i]], exp_sulf[[i]], by = "mag_id")
  tmp[[i]][is.na(tmp[[i]])] <- 0
}

# clean up data frame 
df_list <- lapply(tmp, expression_df)
names(df_list) <- c("car", "irb", "irc")

# Our final dataframe consists of presence/absence matrix, 
# where 0 is no gene found, 1 is gene found in metagenome but not expressed, 
# and 2 is gene in metagenome and expressed 

## Save data ----

saveRDS(df_list, "Formatted_data/sulfur_exp_data_filtered_formatted_20220923")
df_list <- readRDS("Formatted_data/sulfur_exp_data_filtered_formatted_20220923")

## Plot heatmap ----

# get mag phylogeny
car_tree <- keep.tip(full_tree, df_list$car$mag_id)
irb_tree <- keep.tip(full_tree, df_list$irb$mag_id)
irc_tree <- keep.tip(full_tree, df_list$irc$mag_id)

# clean tax strings in df
df_list <- lapply(df_list, clean_tax_strings)

# reduce tax strings to fam (if necessary to group plot by fam)
df_list <- lapply(df_list, reduce_tax_strings)

# Format dataframe
# reorder by phylogeny
#car
df_list$car <- df_list$car[match(car_tree$tip.label, df_list$car$mag_id),] 

# irb
df_list$irb <- df_list$irb[match(irb_tree$tip.label, df_list$irb$mag_id),] 

# irc
df_list$irc <- df_list$irc[match(irc_tree$tip.label, df_list$irc$mag_id),] 

# relevel factors
# mag id
#car
df_list$car$mag_id <- factor(df_list$car$mag_id, levels = car_tree$tip.label)

#irb
df_list$irb$mag_id <- factor(df_list$irb$mag_id, levels = irb_tree$tip.label)

#irc
df_list$irc$mag_id <- factor(df_list$irc$mag_id, levels = irb_tree$tip.label)

# taxonomy
for(i in seq_along(df_list)) {
  df_list[[i]]$Tonomy <- factor(df_list[[i]]$Tonomy, 
                                levels = unique(df_list[[i]]$Tonomy))
  names(df_list[[i]])[2] <- "Taxonomy"
}

# melt data
for (i in seq_along(df_list)) {
  df_list[[i]] <- melt(df_list[[i]], 
                       id.vars = c("mag_id", "Taxonomy"), 
                       value.name = "Expression")
  colnames(df_list[[i]])[3] <- c("KO")
}

# add sulfur cycling steps
for(i in seq_along(df_list)) {
  df_list[[i]]$Process <- sulf_cyc$process[match(df_list[[i]]$KO, sulf_cyc$ko)]
}
head(df_list$car)

# relevel steps
# process
for(i in seq_along(df_list)) {
  df_list[[i]]$Process <- factor(df_list[[i]]$Process, 
                                 levels = c("Taurine_metabolsim", 
                                            "Sulfate_Reduction",
                                            "Thiosulfate_oxidation"))
}

# Group by taxonomy
df_list2 <- list()
for(i in seq_along(df_list)) {
  df_list2[[i]] <- group_by(df_list[[i]], Process, Taxonomy, KO) %>% 
    dplyr::summarise(Expression = max(Expression))
  df_list2[[i]]$Expression <- df_list2[[i]]$Expression %>% as.factor
}

names(df_list2) <- c("car", "irb", "irc")
head(df_list2$car)

# edit x-axis labels
x_labels <- c("K17222" = "soxA (K17222)",
              "K17223" = "soxX (K17223)",
              "K17224" = "soxB (K17224)",
              "K17225" = "soxC (K17225)",
              "K17226" = "soxY (K17226)",
              "K17227" = "soxZ (K17227)",
              "K00958" = "sat (K00958)",
              "K00394" = "aprA (K00394)",
              "K00395" = "aprB (K00395)",
              "K11180" = "dsrA (K11180)",
              "K11181" = "dsrB (K11181)",
              "K15551" = "tauA (K15551)",
              "K15552" = "tauC (K15552)",
              "K10831" = "tauB (K10831)",
              "K03119" = "tauD (K03119)")

# Plot expression
my_cols <- colorRampPalette(c("#132B43", "#56B1F7"))(3) 

ggplot(df_list2$irc, aes(x = KO, y = Taxonomy, fill = Expression)) +
  geom_tile(colour = "black", size=0.10) +
  scale_fill_manual(values = my_cols, 
                    labels = c("Not found", "Gene found", "Gene expressed")) +
  scale_x_discrete(labels = x_labels) +
  scale_y_discrete(position = "right") +
  theme_dark() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.position = "left",
        legend.text = element_text(size = 12),
        legend.key.size = unit(2,"line"),
        legend.title.align = 1) + 
  facet_wrap(df_list2$irc$Process, nrow = 1, ncol = 5, scales = "free_x") +
  theme(strip.text = element_text(size = 12)) +
  guides(fill = guide_legend(label.position = "left"))

ggsave("car_sulfur_fam_f_20220923.svg", device = "svg", width = 37, height = 30, units = "cm", dpi = 300)
ggsave("irb_sulfur_fam_f_20220923.svg", device = "svg", width = 37, height = 31, units = "cm", dpi = 300)
ggsave("irc_sulfur_fam_f_20220923.svg", device = "svg", width = 37, height = 36, units = "cm", dpi = 300)


## Capacity for sulfur cycling in both sponges





