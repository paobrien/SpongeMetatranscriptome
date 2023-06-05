## 05 CAZY expression

# Genome-centric metatranscriptomic analysis of carbohydrate active enzymes

setwd("~/Documents/R/Metatranscriptomics/Metabolic_pathways")

library(tidyverse)
library(plyr)
library(reshape2)
library(viridis)
library(scales)
library(ggpubr)

## Load functions
source("~/Documents/R/Metatranscriptomics/Scripts/R_functions_metatranscriptomics.R")

## Import files ----

# CAZY
cazy_list <- readRDS(
  "~/Documents/R/Metagenomics/Metabolic_pathways/cazy_by_host_metagenome_list"
  )
cazy_list <- cazy_list[c("Car", "Irb", "Irc")]

# Expression data - HTseq
expression_data <- readRDS(
  "~/Documents/R/Metatranscriptomics/HTSeq/Output/sponge_combined_expression_data_filtered_20220923"
  )

# make CAZY datafram ----

# subset to GHs, CEs & PEs
cazy_list <- lapply(cazy_list, function(x) select(x, "mag_id", 
                                             "Taxonomy", 
                                             contains("GH"),
                                             contains("CE"),
                                             contains("PL")))
# get genes
genes <- get_ko(cazy_list$Car)


## format expression data ----

# subset expression data to cazys
exp_cazy <- lapply(expression_data, function(x) exp_pathway(genes, x))


# clean up annotation column (remove any KOs or PFs) 
for(i in seq_along(exp_cazy)) {
  exp_cazy[[i]]$Annotation <- sub("K.*\\,", "", exp_cazy[[i]]$Annotation)
  exp_cazy[[i]]$Annotation <- sub("PF.*\\,", "", exp_cazy[[i]]$Annotation)
}

# add full taxonomy string
for (i in seq_along(exp_cazy)) {
  # first match MAG columns
  colnames(cazy_list[[i]])[which(colnames(cazy_list[[i]]) == "mag_id")] <- "MAG_id" 
  # select and join columns from dfs
  exp_cazy[[i]] <- inner_join(select(cazy_list[[i]], Taxonomy, MAG_id), 
                              select(exp_cazy[[i]], Phylum, Class, Family, 
                                     MAG_id, Annotation, contains("rpkm")))
}

# clean tax strings in df
exp_cazy <- lapply(exp_cazy, clean_tax_strings2)

# format taxonomy
for(i in seq_along(exp_cazy)) {
  exp_cazy[[i]]$Taxonomy <- factor(exp_cazy[[i]]$Taxonomy, 
                                   levels = unique(exp_cazy[[i]]$Taxonomy))
}

## combine dataframes
# change col names
tmp <- exp_cazy
for (i in seq_along(tmp)) {
  colnames(tmp[[i]])[7:9] <- c("Reads_aligned_1",
                               "Reads_aligned_2", 
                               "Reads_aligned_3")
}

# add sponge
tmp$car$Sponge <- rep("P_foliascens", nrow(tmp$car))
tmp$irb$Sponge <- rep("I_microconulosa", nrow(tmp$irb))
tmp$irc$Sponge <- rep("I_ramosa", nrow(tmp$irc))

# combine
df.all <- bind_rows(tmp)

# melt data
# separate
for (i in seq_along(exp_cazy)) {
  exp_cazy[[i]] <- melt(exp_cazy[[i]], 
                        id.vars = c("Taxonomy", "MAG_id", "Annotation",
                                    "Phylum", "Class", "Family"), 
                        value.name = "Expression")
}

# all
df.all <- melt(df.all,
               id.vars = c("Taxonomy", "MAG_id", "Annotation", "Phylum", 
                           "Class", "Family",  "Sponge"), 
               value.name = "Expression")

# remove rows with 0 expression
for (i in seq_along(exp_cazy)) {
  exp_cazy[[i]] <- filter(exp_cazy[[i]], Expression > 0)
}

df.all <- filter(df.all, Expression > 0.1)


## plot expression ----

# transform data to bring out expression
for(i in seq_along(exp_cazy)) {
  exp_cazy[[i]]$Expression_log <- log(exp_cazy[[i]]$Expression+1)
}

df.all$Expression_log <- log(df.all$Expression+1)

# get range
# car
range(exp_cazy$car$Expression_log)
r_car <- c(0, 1.5, 4.41)

# irb
range(exp_cazy$irb$Expression_log)
r_irb <- c(0, 0.5, 3)

# irc
range(exp_cazy$irc$Expression_log)
r_irc <- c(0, 1, 2.3)

# all
range(df.all$Expression_log)
r_all <- c(0, 0.4, 4.41)

## plot relative abundance
# get colours
heat_col <- viridis(n = 3, 
                    alpha = 0.9, 
                    begin = 0, 
                    end = 0.95, 
                    option = "inferno", 
                    direction = 1)
show_col(heat_col)

# plot
#car
p1 <- ggplot(exp_cazy$car, aes(x = Annotation, y = Taxonomy, fill = Expression_log)) +
  geom_tile() +
  scale_fill_gradientn(colours = heat_col,
                       values = rescale(r_car),
                       name = "Expression (log)") +
  scale_y_discrete(position = "right") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "left",
        legend.key.size = unit(2,"line"),
        panel.background = element_rect(fill = "black"))
p1

#irb
p2 <- ggplot(exp_cazy$irb, aes(x = Annotation, y = Taxonomy, fill = Expression_log)) +
  geom_tile() +
  scale_fill_gradientn(colours = heat_col,
                       values = rescale(r_irb),
                       name = "Expression (log)") +
  scale_y_discrete(position = "right") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "left",
        legend.key.size = unit(2,"line"),
        panel.background = element_rect(fill = "black"))
p2

#irc
p3 <- ggplot(exp_cazy$irc, aes(x = Annotation, y = Taxonomy, fill = Expression_log)) +
  geom_tile() +
  scale_fill_gradientn(colours = heat_col,
                       values = rescale(r_irc),
                       name = "Expression (log)", ) +
  scale_y_discrete(position = "right") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "left",
        legend.key.size = unit(2,"line"),
        panel.background = element_rect(fill = "black"))
p3

#all
# get facet labels
facet_labels <- c("I. microconulosa", 
                  "I. ramosa", 
                  "P. foliascens")
names(facet_labels) <-c("I_microconulosa", 
                        "I_ramosa", 
                        "P_foliascens")

# get tax labels
tax_labels <- unite(df.all, "Tax_lab", 
                    c(Phylum, Class, Family), 
                    sep = ';') %>% select(Tax_lab)
tax_labels <- as.vector(tax_labels$Tax_lab)
tax_labels <- gsub("p__", "", tax_labels)
tax_labels <- gsub(";f__$", "", tax_labels)
names(tax_labels) <- df.all$Taxonomy

# remove low expressed cazymes
df.sub <- df.all[df.all$Expression >= 0.25, ]
df.sub$Taxonomy <- gsub("Latescibacterota;c__$", "Latescibacterota", df.sub$Taxonomy)

p4 <- ggplot(df.sub, aes(x = Annotation, y = Taxonomy, fill = Expression_log)) +
  geom_tile() +
  scale_fill_gradientn(colours = heat_col,
                       values = rescale(r_all),
                       name = "Expression \n(log RPKM)", ) +
  scale_y_discrete(labels = tax_labels, position = "right") +
  xlab("CAZy Annotation") +
  ylab("MAG Taxonomy") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16),
        legend.position = "left",
        legend.key.size = unit(2,"line"),
        strip.text = element_text(size = 16, face = "italic"),
        panel.background = element_rect(fill = "black"))

p4 + facet_wrap("Sponge", nrow = 3, ncol = 1, scales = "free_y",
                labeller = labeller(Sponge = facet_labels))

# save
ggsave("cazy_expression_all_f_short_tax.svg", device = "svg", width = 36, height = 45,
       units = "cm", dpi = 300)





## old pipline ----

# Combine duplicate rows (after removing trailing annotations)
for(i in seq_along(exp_cazy)) {
  exp_cazy[[i]] <- ddply(exp_cazy[[i]], 
                        c("MAG_id", "Annotation"), numcolwise(sum))
}

# Convert to wide format
for(i in seq_along(exp_cazy)) {
  exp_cazy[[i]] <- dcast(exp_cazy[[i]], MAG_id ~ Annotation, 
                        value.var ="total")
  exp_cazy[[i]][is.na(exp_cazy[[i]])] <- 0
}

# join with taxonomy
for(i in seq_along(exp_cazy)) {
  names(exp_cazy[[i]])[1] <- "mag_id"
  exp_cazy[[i]] <- left_join(select(cazy_list[[i]], mag_id, Taxonomy), exp_cazy[[i]], by = "mag_id")
  exp_cazy[[i]][is.na(exp_cazy[[i]])] <- 0
}

## Save data ----

# do this when happy with final df

## Plot pathway heatmap ----

# get mag phylogeny
car_tree <- keep.tip(full_tree, cazy_list$Car$mag_id)
irb_tree <- keep.tip(full_tree, cazy_list$Irb$mag_id)
irc_tree <- keep.tip(full_tree, cazy_list$Irc$mag_id)

# clean tax strings in df
for(i in seq_along(cazy_list)) {
    cazy_list[[i]]$Taxonomy <- gsub(";s__$", "", cazy_list[[i]]$Taxonomy)
    cazy_list[[i]]$Taxonomy <- gsub(";g__$", "", cazy_list[[i]]$Taxonomy)
    cazy_list[[i]]$Taxonomy <- gsub(";f__$", "", cazy_list[[i]]$Taxonomy)
    cazy_list[[i]]$Taxonomy <- gsub(";o__$", "", cazy_list[[i]]$Taxonomy)
    cazy_list[[i]]$Taxonomy <- gsub("d__Bacteria;p__", "", cazy_list[[i]]$Taxonomy)
}

# reduce tax strings to fam (if necessary to group plot by fam)
for(i in seq_along(cazy_list)) {
  cazy_list[[i]]$Taxonomy <- gsub(";g.*", "", cazy_list[[i]]$Taxonomy)
}

# Format dataframe
#car
cazy_list$Car <- cazy_list$Car[match(car_tree$tip.label, cazy_list$Car$mag_id),] # reorder by phylogeny

# irb
cazy_list$Irb <- cazy_list$Irb[match(irb_tree$tip.label, cazy_list$Irb$mag_id),] 

# irc
cazy_list$Irc <- cazy_list$Irc[match(irc_tree$tip.label, cazy_list$Irc$mag_id),] 

# relevel factors
# mag id
#car
cazy_list$Car$Car1 <- factor(cazy_list$Car$MAG_id, levels = car_tree$tip.label)

#irb
cazy_list$Irb$mag_id <- factor(cazy_list$Irb$mag_id, levels = irb_tree$tip.label)

#irb
cazy_list$Irc$mag_id <- factor(cazy_list$Irc$mag_id, levels = irc_tree$tip.label)

# taxonomy
for(i in seq_along(cazy_list)) {
  cazy_list[[i]]$Taxonomy <- factor(cazy_list[[i]]$Taxonomy, 
                                    levels = unique(cazy_list[[i]]$Taxonomy))
}

# melt data
for (i in seq_along(cazy_list)) {
  cazy_list[[i]] <- melt(cazy_list[[i]], 
                       id.vars = c("mag_id", "Taxonomy"), 
                       value.name = "Expression")
  colnames(cazy_list[[i]])[3] <- c("CAZY")
}

# Group by taxonomy
df_list2 <- list()
for(i in seq_along(cazy_list)) {
  df_list2[[i]] <- group_by(cazy_list[[i]], Taxonomy, CAZY) %>% 
    dplyr::summarise(Expression = max(Expression))
  df_list2[[i]]$Expression <- df_list2[[i]]$Expression %>% as.factor
}

names(df_list2) <- c("car", "irb", "irc")

# Plot expression
my_cols <- colorRampPalette(c("#132B43", "#56B1F7"))(3) # blue palette (from ggplots continuous scale)

ggplot(df_list2$car, aes(x = CAZY, y = Taxonomy, fill = Expression)) +
  geom_tile(colour = "black", size=0.10) +
  #scale_fill_manual(values = my_cols) +
  theme_dark() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.size = unit(2,"line"))# + 
 # facet_wrap(df_list2$car$Process, nrow = 1, ncol = 5, scales = "free_x") +
  #theme(strip.text = element_text(size = 12))





