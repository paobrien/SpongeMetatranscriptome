## 02 MAG phylogeny

# Plot phylogeny of MAGs with relative abundnace heatmap

setwd("~/Documents/R/Metagenomics") 

library(tidyverse)
library(ape)
library(ggtree)
library(scales)
library(RColorBrewer)
#library(viridis)


## load and format data ----

# Load trees
bact <- ape::read.tree(
  "Data/GTDB/dastool_95.bac120.classify.tree"
  )
arch <- ape::read.tree(
  "Data/GTDB/dastool_95.ar122.classify.tree"
  )

full_tree <- bind.tree(bact, arch)

# load formatted data table
df_d <- read.table(
  "Bin_stats/95_ani/final_bins_stats_coverage_das_95ANI.tsv", 
  sep = '\t', 
  row.names = 1, 
  header = T, 
  strip.white = T)
df_d <- df_d[-416,] # remove unmapped data

# load coverage table
coverm <- read.table(
  "Data/Coverm/all_samples_genome_coverage_dastool_95ANI.txt", 
  sep = '\t', 
  row.names = 1, 
  header = T, 
  strip.white = T)
coverm <- coverm[-1,] # remove unmapped data

# reorder for heatmap
coverm <- coverm[,c("Car1", "Car2", "Car3", "Car4", "Car5", 
                    "IrB1", "IrB2", "IrB3", "IrB4", "IrB5", 
                    "Irc1", "Irc2", "Irc3", "Irc4", "Irc5", 
                    "CS69", "CS70", "CS70b", "CS71", "CS72", "CS73")]

## Subset tree to bins
full_tree <- keep.tip(full_tree, row.names(df_d))

# quick view
plot(full_tree, show.tip.label = F)

## plot with ggtree ----

# Create bin id (same as tip labels) column and move to col 1. 
# Necessary for ggtree %<+% operator to attach metadata

df_d$Bin_id <- rownames(df_d)
df_d <- select(df_d, Bin_id, everything())  # move bin_id to column 1

# get max and min values and summary
r_range <- range(coverm, na.rm = T)

# create values for heatmap
range_val <- c(r_range[1], 0.5, r_range[2])
heat_col <- c("black", "purple", "white")

# edit phylum labels
df_d$Phylum <- gsub("p__", "", df_d$Phylum)
df_d$Class <- gsub("c__", "", df_d$Class)
df_d$Phylum <- df_d$Phylum %>% as.factor()
df_d$Class <- df_d$Class %>% as.factor()

# get plot colours
tree_col <- colorRampPalette(brewer.pal(12, "Paired"))(22)
names(tree_col) <- levels(df_d$Phylum)

## relevel factors
# reorder rows to phylogeny
df_d <- df_d[match(full_tree$tip.label, rownames(df_d)),]
phylum <- df_d$Phylum %>% unique

# relevel phylum
df_d$Phylum <- factor(df_d$Phylum, levels = phylum)
levels(df_d$Phylum)

# plot tree
p <- ggtree(full_tree, layout = "fan", size = 0.4, 
            ladderize = F, open.angle = 15) %<+% df_d +
  geom_tiplab2(size=0.0, align=TRUE, linetype='dashed', 
               linesize=.5, offset = .2, 
               aes(label = Phylum, colour = Phylum))  +
  scale_color_manual(values = tree_col, limits = levels(df_d$Phylum)) +
  guides(color = guide_legend(ncol=1),
         override.aes = list(size = 5)) +
  theme(legend.position="right",
        #legend.key.size = unit(2, "lines"), 
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) 
p

# add heatmap
gheatmap(p, coverm, low = "black", high = "white", color = "black", 
         colnames = T, colnames_angle = 90, width = 0.3, 
         colnames_offset_y = 1, font.size = 0.001)  +
  scale_fill_gradientn(name = "Relative Abundance",
                       colours = heat_col,
                       values = rescale(range_val))

# save
ggsave("MAG_phylogeny_rel_abun5.pdf", 
       device = "pdf", width = 30, height = 30,
       units = "cm", dpi = 300)


