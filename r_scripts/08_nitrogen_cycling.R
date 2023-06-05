## 08 Nitrogen cycling

# Genome-centric metatranscriptomic analysis of nitrogen cycling

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

# Subset dataframes to nitrogen cycling KO's
nit <- lapply(ko_list, function(x) select(x, "mag_id", "Taxonomy", "K10944", "K10945",	
                                          "K10946",	"K10535",	"K00370", "K00371",	
                                          "K00374",	"K02567",	"K02568",	"K00368", 	
                                          "K15864",	"K04561",	"K02305",	"K00376",
                                          "K00362", "K00363", "K03385", "K15876"))

# then get KO's
ko <- get_ko(nit$Car)

# create vector for nitrogen processes
# using different processes of nitrogen cycling (i.e., fixation, nitrification etc)
# note: Denitrification/Reduction used in both processes
nit_cyc <- paste0(c(rep("Nitrification", 4), 
                    rep("Denitrification/Reduction", 5), 
                    rep("Denitrification", 5), 
                    rep("Nitrate_reduction", 4)))

# combine with KO list
nit_cyc_steps <- cbind(ko, nit_cyc) %>% as.data.frame()

# change dataframe to presence absence of KO (gene copy numbers no needed)
nit <- lapply(nit, pres_abs)

## Format expression data ----

# Subset expression data to Nitrogen KOs
exp_nit <- lapply(expression_data, function(x) exp_pathway(ko, x))

# Clean up annotation column
for(i in seq_along(exp_nit)) {
  exp_nit[[i]]$Annotation <- sub("\\,.*", "", exp_nit[[i]]$Annotation)
  exp_nit[[i]]$Annotation <- sub("&.*", "", exp_nit[[i]]$Annotation)
}

# Combine duplicate rows (after removing trailing annotations)
for(i in seq_along(exp_nit)) {
  exp_nit[[i]] <- ddply(exp_nit[[i]], 
                        c("MAG_id", "Annotation"), numcolwise(sum))
}

# Convert to wide format
for(i in seq_along(exp_nit)) {
  exp_nit[[i]] <- dcast(exp_nit[[i]], MAG_id ~ Annotation, 
                        value.var ="total")
  exp_nit[[i]][is.na(exp_nit[[i]])] <- 0
}

# Join the two dataframes
tmp <- list()
for(i in seq_along(exp_nit)) {
  names(exp_nit[[i]])[1] <- "mag_id"
  tmp[[i]] <- left_join(nit[[i]], exp_nit[[i]], by = "mag_id")
  tmp[[i]][is.na(tmp[[i]])] <- 0
}

# clean data frame
df_list <- lapply(tmp, expression_df)
names(df_list) <- c("car", "irb", "irc")

# Our final dataframe consists of presence/absence matrix, 
# where 0 is no gene found, 1 is gene found in metagenome but not expressed, 
# and 2 is gene in metagenome and expressed 


## Save data ----

saveRDS(df_list, "Formatted_data/nitrogen_exp_data_filtered_formatted_20220923")
df_list <- readRDS("Formatted_data/nitrogen_exp_data_filtered_formatted_20220923")

## Plot pathway heatmap ----

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

#irb
df_list$irc$mag_id <- factor(df_list$irc$mag_id, levels = irc_tree$tip.label)

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

# add nitrogen processes
for(i in seq_along(df_list)) {
  df_list[[i]]$Process <- nit_cyc_steps$nit_cyc[match(df_list[[i]]$KO, nit_cyc_steps$ko)]
}

# relevel steps
# process
for(i in seq_along(df_list)) {
  df_list[[i]]$Process <- factor(df_list[[i]]$Process, 
                                 levels = c("Denitrification/Reduction",
                                 "Denitrification",
                                 "Nitrate_reduction",
                                 "Nitrification"))
}

# Group by taxonomy
df_list2 <- list()
for(i in seq_along(df_list)) {
  df_list2[[i]] <- group_by(df_list[[i]], Process, Taxonomy, KO) %>% 
     dplyr::summarise(Expression = max(Expression))
  df_list2[[i]]$Expression <- df_list2[[i]]$Expression %>% as.factor
}

names(df_list2) <- c("car", "irb", "irc")

# edit x-axis labels
x_labels <- c("K10944" = "amoA (K10944)",
              "K10945" = "amoB (K10945)",
              "K10946" = "amoC (K10946)",
              "K10535" = "hao (K10535)",
              "K00370" = "narG (K00370)",
              "K00371" = "narH (K00371)",
              "K00374" = "narI (K00374)",
              "K02567" = "napA (K02567)",
              "K02568" = "napB (K02568)",
              "K00368" = "nirK (K00368)",
              "K15864" = "nirS (K15864)",
              "K04561" = "norB (K04561)",
              "K02305" = "norC (K02305)",
              "K00376" = "nosZ (K00376)",
              "K00362" = "nirB (K00362)",
              "K00363" = "nirD (K00363)",
              "K03385" = "nrfA (K03385)",
              "K15876" = "nrfH (K15876)")

# heatmap colours
my_cols <- colorRampPalette(c("#132B43", "#56B1F7"))(3) # blue palette (from ggplots continuous scale)

# plot
ggplot(df_list2$irc, aes(x = KO, y = Taxonomy, fill = Expression)) +
  geom_tile(colour = "black", size=0.10) +
  scale_fill_manual(values = my_cols, 
                    labels = c("Not found", "Gene found", "Gene expressed")) +
  scale_x_discrete(labels = x_labels) +
  scale_y_discrete(position = "right") +
  theme_dark() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.position = "left",
        legend.text = element_text(size = 12),
        legend.key.size = unit(2,"line"),
        legend.title.align = 1) + 
  facet_wrap(df_list2$irc$Process, nrow = 1, ncol = 5, scales = "free_x") +
  theme(strip.text = element_text(size = 12)) +
  guides(fill = guide_legend(label.position = "left"))


ggsave("car_nit_exp_fam_f_20220923.svg", device = "svg", width = 42, height = 30, units = "cm", dpi = 300)
ggsave("irb_nit_exp_fam_f_20220923.svg", device = "svg", width = 42, height = 30, units = "cm", dpi = 300)
ggsave("irc_nit_exp_fam_f_20220923.svg", device = "svg", width = 42, height = 30, units = "cm", dpi = 300)






