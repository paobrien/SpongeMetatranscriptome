## 09 quantify nitrogen expression

# Quantify nitrogen expression to see if anaerobic and aerobic processes are
# occuring in different sample replicates or coupled within the same individual

setwd("~/Documents/R/Metatranscriptomics/Metabolic_pathways")

library(tidyverse)
library(ape)
library(plyr)
library(reshape2)
library(viridis)
library(scales)

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

## Make KO dataframe ----

# Subset dataframes to nitrogen cycling KO's
nit <- lapply(ko_list, function(x) select(x, "mag_id", "Taxonomy", "K10944", "K10945",	
                                          "K10946",	"K10535",	"K00370", "K00371",	
                                          "K00374",	"K02567",	"K02568",	"K00368", 	
                                          "K15864",	"K04561",	"K02305",	"K00376",
                                          "K00362", "K00363", "K03385", "K15876"))
# then get KO's
ko <- get_ko(nit$Car)

## Format expression data ----

# Subset expression data to Nitrogen KOs
exp_nit <- lapply(expression_data, function(x) exp_pathway(ko, x))

# Clean up annotation column
for(i in seq_along(exp_nit)) {
  exp_nit[[i]]$Annotation <- sub("\\,.*", "", exp_nit[[i]]$Annotation)
  exp_nit[[i]]$Annotation <- sub("&.*", "", exp_nit[[i]]$Annotation)
}

# Combine tax column (Phylum-MAG_id)
for(i in seq_along(exp_nit)) {
  exp_nit[[i]] <- unite(exp_nit[[i]], "Taxonomy", c(Phylum:Genus,MAG_id), sep = ';', remove = F)
  exp_nit[[i]]$Taxonomy <- gsub("p__", "", exp_nit[[i]]$Taxonomy)
}

# get new tax labels
tax_labels <- unite(exp_nit$irc, "Taxonomy", c(Phylum, Class, Family), sep = ';') %>%
  select(Taxonomy)

# Combine duplicate rows (after removing trailing annotations)
for(i in seq_along(exp_nit)) {
  exp_nit[[i]] <- ddply(exp_nit[[i]], 
                        c("MAG_id", "Annotation", "Taxonomy"), numcolwise(sum))
}

## Plot expression ----

# melt data
exp_nit.m <- list()
for (i in seq_along(exp_nit)) {
  exp_nit[[i]] <- select(exp_nit[[i]], !c(gene_length, total))
  exp_nit[[i]] <- select(exp_nit[[i]], !ends_with("raw"))
  exp_nit.m[[i]] <- melt(exp_nit[[i]], 
                        id.vars = c("MAG_id", "Annotation", "Taxonomy"), 
                        value.name = "Expression")
}
names(exp_nit.m) <- names(exp_nit)

# transform data to bring out expression
for(i in seq_along(exp_nit.m)) {
  exp_nit.m[[i]]$Expression_log <- log(exp_nit.m[[i]]$Expression+1)
}

# get range
# car
range(exp_nit.m$car$Expression_log)
r_car <- c(0, 0.4, 0.96)

# irb
range(exp_nit.m$irb$Expression_log)
r_irb <- c(0, 1, 3.6)

# irc
range(exp_nit.m$irc$Expression_log)
r_irc <- c(0, 0.4, 3.4)

# get colours
heat_col <- viridis(n = 3, 
                    alpha = 0.9, 
                    begin = 0, 
                    end = 0.95, 
                    option = "inferno", 
                    direction = 1)
show_col(heat_col)


# plot expression heatmaps
#car
p1 <- ggplot(exp_nit.m$car, aes(x = Annotation, y = Taxonomy, fill = Expression_log)) +
  geom_tile() +
  scale_fill_gradientn(colours = heat_col,
                       values = rescale(r_car),
                       name = "Expression \n(log RPKM)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.size = unit(2,"line"),
        panel.background = element_rect(fill = "black"))
p1 + facet_wrap(facets = "variable", nrow = 1, ncol = 3, scales = "fixed")

#irb
p2 <- ggplot(exp_nit.m$irb, aes(x = Annotation, y = Taxonomy, fill = Expression_log)) +
  geom_tile() +
  scale_fill_gradientn(colours = heat_col,
                       values = rescale(r_irb),
                       name = "Expression \n(log RPKM)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.key.size = unit(2,"line"),
        panel.background = element_rect(fill = "black"))
p2 + facet_wrap(facets = "variable", nrow = 1, ncol = 3, scales = "fixed")

#irc
# get facet labels
facet_labels <- c("Reads_aligned_irc1.rpkm" = "I. ramosa 1", 
                  "Reads_aligned_irc2.rpkm" = "I. ramosa 2", 
                  "Reads_aligned_irc3.rpkm" = "I. ramosa 3")

# edit x-axis labels
x_labels <- c("K10944" = "amoA (K10944)",
              "K10945" = "amoB (K10945)",
              "K10946" = "amoC (K10946)",
              "K10535" = "hao (K10535)",
              "K00370" = "narG (K00370)",
              "K00371" = "narH (K00371)",
              "K00374" = "narI (K00374)",
              "K00363" = "nirD (K00363)",
              "K00368" = "nirK (K00368)",
              "K04561" = "norB (K04561)")

#plot
p3 <- ggplot(exp_nit.m$irc, aes(x = Annotation, y = Taxonomy, fill = Expression_log)) +
  geom_tile() +
  scale_fill_gradientn(colours = heat_col,
                       values = rescale(r_irc),
                       name = "Expression \n(log RPKM)") +
  scale_x_discrete(labels = x_labels) +
  scale_y_discrete(labels = tax_labels, position = "right") +
  xlab("Gene Annotation") +
  ylab("MAG Taxonomy") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "left",
        legend.key.size = unit(2,"line"),
        panel.background = element_rect(fill = "black"),
        strip.background = element_rect(size = 0.5),
        strip.text = element_text(size = 14, face = "italic"),
        plot.margin = unit(c(1, 1, 2, 1), "lines"))
p3 + 
  coord_cartesian(ylim=c(1, 32), clip="off") +
  annotate("segment", x = 0.8, xend = 6.2, y = -6, yend = -6, color = "black", size = 0.8) +
  annotate("text", x = 3.6, y = -6.5, label = "anaerobic", size = 3.5) +
  annotate("segment", x = 6.8, xend = 10.2, y = -6, yend = -6, color = "black", size = 0.8) +
  annotate("text", x = 8.6, y = -6.5, label = "aerobic", size = 3.5) +
  facet_wrap(facets = "variable", 
                nrow = 1, 
                ncol = 3, 
                scales = "fixed",
                labeller = labeller(variable = facet_labels))

ggsave("nit_exp_quant_irc.svg", device = "svg", width = 30, height = 20, units = "cm", dpi = 300)

## Get aerobic vs anaerobic expression values ----

# add metabolism to df
exp_nit.m$irc <- exp_nit.m$irc %>%
  mutate(metabolism = case_when(
    Annotation %in% c("K10944", "K10945", "K10946", "K10535") ~ "aerobic",
    Annotation %in% c("K00370", "K00371", "K00374", "K00363", "K00368", "K04561") ~ "anaerobic",
    TRUE ~ NA_character_
  ))

# group data to calculate the mean and standard error of expression

df_summary <- exp_nit.m$irc %>%
  group_by(variable, metabolism) %>%
  dplyr::summarise(mean_expression = mean(Expression), 
                   se_expression = sd(Expression)/sqrt(n()),
                   total_expression = sum(Expression))
df_summary








