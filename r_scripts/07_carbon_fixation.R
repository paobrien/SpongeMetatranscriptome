## 07 Carbon fixation

# Genome-centric metatranscriptomic analysis of carbon fixation pathways

setwd("~/Documents/R/Metatranscriptomics/Metabolic_pathways")

library(tidyverse)
library(ape)
library(plyr)
library(reshape2)


## Load functions\
source("~/Documents/R/Metatranscriptomics/Scripts/R_functions_metatranscriptomics.R")

## Import files ----
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


## Make KO dataframes ----

## Autotrophic carbon fixation pathways 

# KOs and steps involved are based on KEGG

# subset data to pathway
# Reductive acetyl-CoA pathway (Wood-Ljungdahl pathway)
wl_path <- lapply(
  ko_list, function(x) select(x, "mag_id", "Taxonomy", "K00198", "K05299", 
"K15022", "K01938", "K01491","K00297","K15023", "K14138", "K00197", "K00194")
)
wl_ko <- get_ko(wl_path$Car)

wl_steps <- c("S_01", "S_02", "S_02", "S_03", "S_04",
           "S_05", "S_06", "S_07", "S_07", "S_07")

wl_ko_steps <- cbind(wl_ko, wl_steps) %>% as.data.frame

## Dicarboxylate/4-hydroxybutyrate cycle (DC-HB)
dc_hb_path <- lapply(
  ko_list, function(x) select(x, "mag_id", "Taxonomy", "K00169",	
                              "K00170",	"K00171",	"K00172",	"K01007", 
                              "K01595","K00024", "K01677",	"K01678",	
                              "K00239",	"K00240","K00241",	"K18860",	
                              "K01902",	"K01903", "K15038",	"K15017",	
                              "K14465",	"K14467", "K18861", "K14534",	
                              "K15016",	"K00626")
  )
dc_hb_ko <- get_ko(dc_hb_path$Car)

dc_hb_steps <-  c("S_01", "S_01", "S_01", "S_01", "S_02", "S_03",
                  "S_04", "S_05", "S_05", "S_06", "S_06", "S_06", 
                  "S_06", "S_07", "S_07", "S_08", "S_08", "S_09", 
                  "S_10", "S_10", "S_11", "S_12", "S_13")

dc_hb_ko_steps <- cbind(dc_hb_ko, dc_hb_steps) %>% as.data.frame


## 3-hydroxypropionate bicycle (3-HP)
hp_path <- lapply(
  ko_list, function(x) select(x,"mag_id", "Taxonomy", "K02160", 
                              "K01961", "K01962", "K01963", "K14468",
                              "K14469", "K15052", "K05606", "K01847", 
                              "K01848", "K01849", "K14471", "K14472",	
                              "K00239", "K00240", "K00241", "K01679",
                              "K08691", "K14449", "K14470", "K09709"))

hp_ko <- get_ko(hp_path$Car)

hp_steps <- c("S_01", "S_01", "S_01", "S_01", "S_02", "S_03", 
           "S_04", "S_05", "S_06", "S_06", "S_06", "S_07/10", 
           "S_07/10", "S_08", "S_08", "S_08", "S_09", "S_11", 
           "S_12", "S_13", "S_14")

hp_ko_steps <- cbind(hp_ko, hp_steps) %>% as.data.frame

## 3-hydroxypropionate/4-hydroxybutyrate cycle (HP-HB)
hp_hb_path <- lapply(
  ko_list, function(x) select(x, "mag_id", "Taxonomy","K01964",
                              "K15036",	"K15037",	"K15017",	
                              "K15039",	"K15018",	"K15019",
                              "K15020",	"K05606",	"K01848",	
                              "K01849",	"K15038",	"K14465",	
                              "K14466",	"K18861", "K14534",	
                              "K15016",	"K00626"))

hp_hb_steps <- c("S_01/07", "S_01/07", "S_01/07", 
                 "S_02/10", "S_03", "S_04", "S_05", 
                 "S_06", "S_08", "S_09", "S_09", 
                 "S_10", "S_11", "S_12", "S_12", 
                 "S_13", "S_14", "S_15")

hp_hb_ko <- get_ko(hp_hb_path$Car)

hp_hb_ko_steps <- cbind(hp_hb_ko, hp_hb_steps) %>% as.data.frame

## reductive citric acid cycle (rTCA)
rtca_path <- lapply(
  ko_list, function(x) select(x, "mag_id", "Taxonomy","K00169",
                              "K00170",	"K00171",	"K00172",	"K03737",	"K01007",
                              "K01006",	"K01595",	"K01959",	"K01960",	"K01958",
                              "K00024",	"K01676",	"K01679",	"K01677",	"K01678",	
                              "K00239",	"K00240",	"K00241",	"K00242",	"K00244",
                              "K00245",	"K00246",	"K00247",	"K01902",	"K01903",
                              "K00174",	"K00175",	"K00177",	"K00176",	"K00031",	
                              "K01681",	"K01682",	"K15230", "K15231", "K15232",	
                              "K15233",	"K15234"))

# note: removed "K18556",	"K18557","K18558",	"K18559",	"K18560"
# not in our data and not necessary


rtca_steps <- c("S_01", "S_01", "S_01", "S_01", "S_01", "S_02", "S_02", 
           "S_02", "S_02", "S_02", "S_02", "S_03", "S_04", "S_04", 
           "S_04", "S_04", "S_05", "S_05", "S_05", "S_05", "S_05",
           "S_05", "S_05", "S_05", "S_06", "S_06", "S_07", "S_07", 
           "S_07", "S_07", "S_08", "S_09", "S_09", "S_10", "S_10", 
           "S_10", "S_10", "S_10")

rtca_ko <- get_ko(rtca_path$Car)

rtca_ko_steps <- cbind(rtca_ko, rtca_steps) %>% as.data.frame


## Calvin-Benson-Bassham cycle (CBB)
cbb_path <- lapply(ko_list, function(x) select(x, "mag_id", "Taxonomy", "K00855",
                                          "K01601",	"K01602",	"K00927",	"K00150",
                                          "K00134",	"K01623",	"K01624",	"K03841",	
                                          "K02446", "K11532",	"K00615",	"K01623",	
                                          "K01624",	"K11532",	"K00615", "K01807",	
                                          "K01808"))

# note: removed "K01100", "K01086",	"K01086", "K05298",	
# not in data and function could be performed by other enzymes

cbb_ko <- get_ko(cbb_path$Car)

cbb_steps <- c("S_01", "S_02", "S_02", "S_03", 
               "S_04", "S_04", "S_05/08", 
               "S_05/08", "S_06/09", "S_06/09", 
               "S_06/09", "S_07/10", 
               "S_11", "S_11")

cbb_ko_steps <- cbind(cbb_ko, cbb_steps) %>% as.data.frame


# Combine pathways
path_names <- c("wl", "dc_hb", "hp", "hp_hb", "rTCA", "cbb")

cf_list <- list(wl_path, dc_hb_path, hp_path, hp_hb_path, rtca_path, cbb_path)
names(cf_list) <- path_names

# combine ko
ko <- list(wl_ko, dc_hb_ko, hp_ko, hp_hb_ko, rtca_ko, cbb_ko)
names(ko) <- path_names

# combine ko_steps 
ko_steps <- list(wl_ko_steps, dc_hb_ko_steps, hp_ko_steps, hp_hb_ko_steps,
                 rtca_ko_steps, cbb_ko_steps)

names(ko_steps) <- path_names

# change to presence / absence 
tmp_c <- list()
tmp_i <- list()
tmp_r <- list()
for(i in seq_along(cf_list)) {
# Carterio 
  tmp_c[[i]] <- ifelse(cf_list[[i]][["Car"]][,-c(1:2)] > 0, 1, 0)
  cf_list[[i]][["Car"]] <- cbind(cf_list[[i]][["Car"]][,1:2], tmp_c[[i]])
  print(head(cf_list[[i]][["Car"]]))
# Imicro
  tmp_i[[i]] <- ifelse(cf_list[[i]][["Irb"]][,-c(1:2)] > 0, 1, 0)
  cf_list[[i]][["Irb"]] <- cbind(cf_list[[i]][["Irb"]][,1:2], tmp_i[[i]])
  print(head(cf_list[[i]][["Irb"]]))
# Iramosa
  tmp_r[[i]] <- ifelse(cf_list[[i]][["Irc"]][,-c(1:2)] > 0, 1, 0)
  cf_list[[i]][["Irc"]] <- cbind(cf_list[[i]][["Irc"]][,1:2], tmp_r[[i]])
  print(head(cf_list[[i]][["Irb=c"]]))
}


## Format expression data ----

# Subset data and combine reads column
names(expression_data) <- c("Car", "Irb", "Irc")

# Get KO expression data
exp_cf <- list()
for(i in seq_along(expression_data)) {
  exp_cf[[i]] <- lapply(ko, function(x) exp_pathway(x, expression_data[[i]]))
}
names(exp_cf) <- c("car", "irb", "irc")

# Clean up annotation column
for(i in seq_along(exp_cf$car)) {
  exp_cf$car[[i]]$Annotation <- sub("\\,.*", "", exp_cf$car[[i]]$Annotation)
  exp_cf$car[[i]]$Annotation <- sub("&.*", "", exp_cf$car[[i]]$Annotation)
  exp_cf$irb[[i]]$Annotation <- sub("\\,.*", "", exp_cf$irb[[i]]$Annotation)
  exp_cf$irb[[i]]$Annotation <- sub("&.*", "", exp_cf$irb[[i]]$Annotation)
  exp_cf$irc[[i]]$Annotation <- sub("\\,.*", "", exp_cf$irc[[i]]$Annotation)
  exp_cf$irc[[i]]$Annotation <- sub("&.*", "", exp_cf$irc[[i]]$Annotation)
}


# Combine duplicate rows (after removing trailing annotations)
for (i in seq_along(exp_cf$car)) {
  exp_cf$car[[i]] <- ddply(exp_cf$car[[i]],c("MAG_id", "Annotation"),numcolwise(sum))
  exp_cf$irb[[i]] <- ddply(exp_cf$irb[[i]],c("MAG_id", "Annotation"),numcolwise(sum))
  exp_cf$irc[[i]] <- ddply(exp_cf$irc[[i]],c("MAG_id", "Annotation"),numcolwise(sum))
}

# Convert to wide format
for(i in seq_along(exp_cf$car)) {
  exp_cf$car[[i]] <- dcast(exp_cf$car[[i]], MAG_id ~ Annotation, value.var = "total")
  exp_cf$car[[i]][is.na(exp_cf$car[[i]])] <- 0  # change NAs to 0 (missing values meant there were no KOs for that MAG)
  exp_cf$irb[[i]] <- dcast(exp_cf$irb[[i]], MAG_id ~ Annotation, value.var = "total")
  exp_cf$irb[[i]][is.na(exp_cf$irb[[i]])] <- 0
  exp_cf$irc[[i]] <- dcast(exp_cf$irc[[i]], MAG_id ~ Annotation, value.var = "total")
  exp_cf$irc[[i]][is.na(exp_cf$irc[[i]])] <- 0
}

# Join the two dataframes
tmp_c <- list()
tmp_i <- list()
tmp_r <- list()
for(i in seq_along(exp_cf$car)) {
  #car
  names(exp_cf$car[[i]])[1] <- "mag_id"
  tmp_c[[i]] <- left_join(cf_list[[i]][["Car"]], exp_cf$car[[i]], by = "mag_id")
  tmp_c[[i]][is.na(tmp_c[[i]])] <- 0
  #irb
  names(exp_cf$irb[[i]])[1] <- "mag_id"
  tmp_i[[i]] <- left_join(cf_list[[i]][["Irb"]], exp_cf$irb[[i]], by = "mag_id")
  tmp_i[[i]][is.na(tmp_i[[i]])] <- 0
  #irc
  names(exp_cf$irc[[i]])[1] <- "mag_id"
  tmp_r[[i]] <- left_join(cf_list[[i]][["Irc"]], exp_cf$irc[[i]], by = "mag_id")
  tmp_r[[i]][is.na(tmp_r[[i]])] <- 0
}
tmp <- list(tmp_c, tmp_i, tmp_r)
names(tmp) <- c("car", "irb", "irc")

# clean up data frame
df_list_c <- list()
for(i in seq_along(tmp)) {
  df_list_c[[i]] <- lapply(tmp[[i]], expression_df)
}
names(df_list_c) <- c("car", "irb", "irc")

for(i in seq_along(df_list_c)) {
  names(df_list_c[[i]]) <- c("wl", "dc_hb", "hp", "hp_hb", "rTCA", "cbb")
}


# Our final dataframe consists of presence/absence matrix, 
# where 0 is no gene found, 1 is gene found in metagenome but not expressed, 
# and 2 is gene in metagenome and expressed 

## Save data ----

# carbon fixation data
saveRDS(df_list_c, "Formatted_data/carbon_fix_exp_data_filtered_20220923")
saveRDS(ko_steps, "Formatted_data/carbon_fix_steps_list_20220923")

df_list <- readRDS("Formatted_data/carbon_fix_exp_data_filtered_20220923")
ko_steps <-readRDS("Formatted_data/carbon_fix_steps_list_20220923") 

## Plot heatmap ----

# get mag phylogeny
car_tree <- keep.tip(full_tree, df_list$car$wl$mag_id)
irb_tree <- keep.tip(full_tree, df_list$irb$wl$mag_id)
irc_tree <- keep.tip(full_tree, df_list$irc$wl$mag_id)

# clean tax strings in df
for(i in seq_along(df_list)) {
  df_list[[i]] <- lapply(df_list[[i]], clean_tax_strings)
}

# reduce tax strings to fam
for(i in seq_along(df_list)) {
  df_list[[i]] <- lapply(df_list[[i]], reduce_tax_strings)
}

# Format dataframe
for(i in seq_along(df_list$car)) {
  df_list$car[[i]] <- df_list$car[[i]][match(car_tree$tip.label, 
                                             df_list$car$wl$mag_id),]
  df_list$irb[[i]] <- df_list$irb[[i]][match(irb_tree$tip.label, 
                                             df_list$irb$wl$mag_id),]
  df_list$irc[[i]] <- df_list$irc[[i]][match(irc_tree$tip.label,
                                             df_list$irc$wl$mag_id),]
}

## relevel factors
#mag id (to phylogeny)
for(i in seq_along(df_list$car)) {
  df_list$car[[i]]$mag_id <- factor(df_list$car[[i]]$mag_id, 
                                    levels = car_tree$tip.label)
  df_list$irb[[i]]$mag_id <- factor(df_list$irb[[i]]$mag_id, 
                                    levels = car_tree$tip.label)
  df_list$irc[[i]]$mag_id <- factor(df_list$irc[[i]]$mag_id, 
                                    levels = car_tree$tip.label)
}

#taxonomy
for (i in seq_along(df_list$car)) {
  df_list$car[[i]]$Tonomy <- factor(df_list$car[[i]]$Tonomy, 
                                    levels = unique(df_list$car[[i]]$Tonomy))
  df_list$irb[[i]]$Tonomy <- factor(df_list$irb[[i]]$Tonomy, 
                                    levels = unique(df_list$irb[[i]]$Tonomy))
  df_list$irc[[i]]$Tonomy <- factor(df_list$irc[[i]]$Tonomy, 
                                    levels = unique(df_list$irc[[i]]$Tonomy))
  names(df_list$car[[i]])[2] <- "Taxonomy"
  names(df_list$irb[[i]])[2] <- "Taxonomy"
  names(df_list$irc[[i]])[2] <- "Taxonomy"
}

# change to long format
for (i in seq_along(df_list$car)) {
  df_list$car[[i]] <- melt(df_list$car[[i]], 
                           id.vars = c("mag_id", "Taxonomy"), 
                           value.name = "Expression")
  colnames(df_list$car[[i]])[3] <- c("KO")
  df_list$irb[[i]] <- melt(df_list$irb[[i]], 
                           id.vars = c("mag_id", "Taxonomy"), 
                           value.name = "Expression")
  colnames(df_list$irb[[i]])[3] <- c("KO")
  df_list$irc[[i]] <- melt(df_list$irc[[i]], 
                           id.vars = c("mag_id", "Taxonomy"), 
                           value.name = "Expression")
  colnames(df_list$irc[[i]])[3] <- c("KO")
}

# Add pathway steps (steps from top of script)
for (i in seq_along(df_list$car)) {
  df_list$car[[i]]$Steps <- ko_steps[[i]][,2][match(df_list$car[[i]]$KO, 
                                                    ko_steps[[i]][,1])]
  df_list$irb[[i]]$Steps <- ko_steps[[i]][,2][match(df_list$irb[[i]]$KO, 
                                                    ko_steps[[i]][,1])]
  df_list$irc[[i]]$Steps <- ko_steps[[i]][,2][match(df_list$irc[[i]]$KO, 
                                                    ko_steps[[i]][,1])]
}

# Group by taxonomy (family level)
df_list2 <- list()
df_list2$car <- list()
df_list2$irb <- list()
df_list2$irc <- list()
for(i in seq_along(df_list$car)) {
  #car
  df_list2$car[[i]] <- group_by(df_list$car[[i]], Steps, Taxonomy, KO) %>% 
    dplyr::summarise(Expression = max(Expression)) 
  df_list2$car[[i]]$Expression <- df_list2$car[[i]]$Expression %>% as.factor
  #irb
  df_list2$irb[[i]] <- group_by(df_list$irb[[i]], Steps, Taxonomy, KO) %>% 
    dplyr::summarise(Expression = max(Expression)) 
  df_list2$irb[[i]]$Expression <- df_list2$irb[[i]]$Expression %>% as.factor
  #irc
  df_list2$irc[[i]] <- group_by(df_list$irc[[i]], Steps, Taxonomy, KO) %>% 
    dplyr::summarise(Expression = max(Expression)) 
  df_list2$irc[[i]]$Expression <- df_list2$irc[[i]]$Expression %>% as.factor
}

names(df_list2$car) <- path_names
names(df_list2$irb) <- path_names
names(df_list2$irc) <- path_names

# edit x-axis labels (for complete or near complete pathways)
hp_hb_labs <- c("K01964" = "K01964", 
                "K15036" = "K15036",
                "K15037" = "K15037",
                "K15017" = "K15017",
                "K15039" = "K15039",
                "K15018" = "K15018",
                "K15019" = "K15019",
                "K15020" = "K15020",
                "K05606" = "MCEE, epi (K05606)",
                "K01848" = "mcmA1 (K01848)",
                "K01849" = "mcmA2 (K01849)",
                "K15038" = "K15038",
                "K14465" = "K14465",
                "K14466" = "K14466",
                "K18861" = "K18861",
                "K14534" = "abfD (K14534)",
                "K15016" = "K15016",
                "K00626" = "atoB (K00626)")

rtca_labs <- c("K00169" = "porA (K00169)",
               "K00170" = "porB (K00170)",
               "K00171" = "porD (K00171)",
               "K00172" = "porC (K00172)",
               "K03737" = "nifJ (K03737)",
               "K01007" = "ppsA (K01007)",
               "K01006" = "ppdK (K01006)",
               "K01595" = "ppc (K01595)",
               "K01959" = "pycA (K01959)",
               "K01960" = "pycB (K01960)",
               "K01958" = "pyc (K01958)",
               "K00024" = "mdh (K00024)",
               "K01676" = "fumA, fumB (K01676)",
               "K01679" = "fumC, FH (K01679)",
               "K01677" = "fumA (K01677)",
               "K01678" = "fumB (K01678)",
               "K00239" = "sdhA, frdA (K00239)",
               "K00240" = "sdhB, frdB (K00240)",
               "K00241" = "sdhC, frdC (K00241)",
               "K00242" = "sdhD, frdD (K00242)",
               "K00244" = "frdA (K00244)",
               "K00245" = "frdB (K00245)",
               "K00246" = "frdC (K00246)",
               "K00247" = "frdD (K00247)",
               "K01902" = "sucD (K01902)",
               "K01903" = "sucC (K01903)",
               "K00174" = "korA, oorA, ofofA (K00174)",
               "K00175" = "korB, oorB, ofofB (K00175)",
               "K00177" = "korC, oorC (K00177)",
               "K00176" = "korD, oorD (K00176)",
               "K00031" = "IDH1, IDH2, icd (K00031)",
               "K01681" = "ACO, acnA (K01681)",
               "K01682" = "acnB (K01682)",
               "K15230" = "aclA (K15230)",
               "K15231" = "aclB (K15231)",
               "K15232" = "ccsA (K15232)",
               "K15233" = "ccsB (K15233)",
               "K15234" = "ccl (K15234)")

cbb_labs <- c("K00855" = "PRK, prkB (K00855)",
              "K01601" = "rbcl, cbbL (K01601)",
              "K01602" = "rbcS, cbbS (K01602)",
              "K00927" = "PGK, pgk (K00927)",
              "K00150" = "gap2, gapB (K00150)",
              "K00134" = "gapA (K00134)",
              "K01623" = "ALDO (K01623)",
              "K01624" = "FBA, fbaA (K01624)",
              "K03841" = "FBP, fbp (K03841)",
              "K02446" = "glpX (K02446)",
              "K11532" = "glpX-SEBP (K11532)",
              "K00615" = "fbp-SEBP (K00615)",
              "K01807" = "rpiA (K01807)",
              "K01808" = "rpiB (K01808)")

x_labels <- list(hp_hb_labs, rtca_labs, cbb_labs)
names(x_labels) <- c("hp_hb", "rtc", "cbb")

# heatmap colours
my_cols <- colorRampPalette(c("#132B43", "#56B1F7"))(3) 

# plot
ggplot(df_list2$irc$cbb, aes(x = KO, y = Taxonomy, fill = Expression)) +
  scale_x_discrete(labels = x_labels$cbb) +
  scale_y_discrete(position = "right") +
  geom_tile(colour = "black", size=0.10) +
  scale_fill_manual(values = my_cols, 
                    labels = c("Not found", "Gene found", "Gene expressed")) +
  theme_dark() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.position = "left",
        legend.text = element_text(size = 12),
        legend.key.size = unit(2,"line"),
        legend.title.align = 1) + 
  facet_wrap(df_list2$irc$cbb$Steps, nrow = 1, ncol = 15, scales = "free_x") +
  theme(strip.text = element_text(size = 10)) +
  guides(fill = guide_legend(label.position = "left"))

# hpbb
ggsave("hp_hb_car_exp_fam_20220924.svg", device = "svg", width = 45, height = 30, units = "cm", dpi = 300)
ggsave("hp_hb_irb_exp_fam_20220924.svg", device = "svg", width = 45, height = 30, units = "cm", dpi = 300)
ggsave("hp_hb_irc_exp_fam_20220924.svg", device = "svg", width = 45, height = 34, units = "cm", dpi = 300)

# rTCA
ggsave("rTCA_car_exp_fam_20220924.svg", device = "svg", width = 52, height = 30, units = "cm", dpi = 300)
ggsave("rTCA_irb_exp_fam_20220924.svg", device = "svg", width = 52, height = 30, units = "cm", dpi = 300)
ggsave("rTCA_irc_exp_fam_20220924.svg", device = "svg", width = 52, height = 34, units = "cm", dpi = 300)

# cbb
ggsave("cbb_car_exp_fam_20220924.svg", device = "svg", width = 38, height = 30, units = "cm", dpi = 300)
ggsave("cbb_irb_exp_fam_20220924.svg", device = "svg", width = 38, height = 30, units = "cm", dpi = 300)
ggsave("cbb_irc_exp_fam_20220924.svg", device = "svg", width = 38, height = 34, units = "cm", dpi = 300)







