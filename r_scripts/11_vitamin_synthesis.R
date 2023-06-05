## Genome-centric metatranscriptomic analysis of thiamine (Vitamin B1) metabolism

setwd("~/Documents/R/Metatranscriptomics/Metabolic_pathways")

library(tidyverse)
library(ape)
library(plyr)
library(reshape2)

## Load functions from 'functions_metatranscriptomics.R' script 

source("Scripts/R_functions_metatranscriptomics.R")

## Import files ----
ko_list <- readRDS("~/Documents/R/Metagenomics/Metabolic_pathways/ko_by_host_metagenome_list")
ko_list <- ko_list[c("Car", "Irb", "Irc")]

# Expression data - HTseq
expression_data <- readRDS("~/Documents/R/Metatranscriptomics/HTSeq/Output/sponge_combined_expression_data_filtered_20220923")

# Tree
bact <- ape::read.tree("~/Documents/R/Metagenomics/Data/GTDB/all_mags.bac120.classify.tree")
arch <- ape::read.tree("~/Documents/R/Metagenomics/Data/GTDB/all_mags.ar122.classify.tree")
full_tree <- bind.tree(bact, arch)


## Make KO dataframes ----

## Vitamin B1 biosynthesis (thiamine)
tm_path <- lapply(ko_list, function(x) select(x, "mag_id", "Taxonomy", "K02064",	
                                               "K02063", "K02062", "K03147", "K00941", 
                                               "K14153", "K03148", "K03151", "K03150", 
                                               "K03153", "K03149", "K10810", "K03146", 
                                               "K00788", "K00946", "K01077", "K00949"))

# note: "K14154" not in data and not key enzyme - removing

# get ko
tm_ko <- get_ko(tm_path$Car)

# Reactions - multiple reactions from different pathways combined to make thiamine
tm_rxn <- c(rep("Transport", 3), "R_01", 
            "R_02", "R_02/08", "R_03", 
            "R_03", "R_04", "R_05", "R_06",
            "R_06", "R_07", "R_08", "R_08",
            rep("R_09", 2))

#combine data
tm_ko_rxn <- cbind(tm_ko, tm_rxn) %>% as.data.frame



## Vitamin B2 biosynthesis (riboflavin)
rb_path <- lapply(ko_list, function(x) select(x, "mag_id", "Taxonomy", "K01497", 
                                              "K14652", "K11752", "K20862", "K00794", 
                                              "K00793", "K11753", "K02858"))

#note: removed K01498, K00082,K20860, K20861 as not in data and not essential enzyme

# get ko
rb_ko <- get_ko(rb_path$Car)

#get steps
rb_steps <- c("S_01", "S_01", "S_02", 
              "S_03", "S_04", "S_05", 
              "S_06", "S_01b")

#combine data
rb_ko_steps <- cbind(rb_ko, rb_steps) %>% as.data.frame


## Vitamin B5 (Pantothenate)
pt_path <- lapply(ko_list, function(x) select(x, "mag_id", "Taxonomy", "K00826",	"K00606",	
                                         "K00077", "K01579", "K18933", "K18966",
                                         "K17839", "K13367", "K00128", "K01918",
                                         "K13799"))
#get ko
pt_ko <- get_ko(pt_path$Car)

#get steps
pt_rxns <- c("R_01", "R_01", "R_01", 
             "R_02", "R_02", "R_02",
             "R_03", "R_03", "R_03",
             "R_04", "R_04")

#combine
pt_ko_rxns <- cbind(pt_ko, pt_rxns) %>% as.data.frame


## Vitamin B6 (Pyridoxine)
pd_path <- lapply(ko_list, function(x) select(x, "mag_id", "Taxonomy", "K03472",	
                                              "K03473", "K00831", "K00097", 
                                              "K03474", "K00275", "K06215", 
                                              "K08681", "K01733", "K00868"))

# get ko
pd_ko <- get_ko(pd_path$Car)

# rxn
pd_rxns <- c(rep("R_01", 6), "R_02", 
             "R_02", "R_03", "R_03")

#comine data
pd_ko_rxns <- cbind(pd_ko, pd_rxns) %>% as.data.frame


## Vitamin B7 (Biotin)
bt_path <- lapply(ko_list, function(x) select(x, "mag_id", "Taxonomy", "K03523",	
                                              "K16783", "K16784", "K16785", 
                                              "K16786", "K16787", "K16593", 
                                              "K01906", "K00652", "K00833", 
                                              "K01935", "K01012")) 

# note: removing "K19563", "K19562",

#get ko
bt_ko <- get_ko(bt_path$Car)

#get steps
bt_steps <- c(rep("Transporter_1", 3), 
              rep("Transporter_2", 3), 
              "S_01", "S_01", "S_02", 
              "S_03", "S_03", "S_04")

#combine
bt_ko_steps <- cbind(bt_ko, bt_steps) %>% as.data.frame


## Vitamin B12 (cobalamin)
cb_path <- lapply(ko_list, function(x) select(x, "mag_id", "Taxonomy", "K02303",	
                                              "K13542", "K02302", "K02304", "K03795",
                                              "K03394", "K05934",	"K13541", "K05936",
                                              "K02189", "K02188",	"K05895", "K00595",
                                              "K06042",	"K02224", "K03394", "K02229",
                                              "K05934", "K13541", "K05936",	"K02228",
                                              "K05895", "K00595", "K06042", "K02224",
                                              "K02230", "K09882", "K09883", "K00798", 
                                              "K19221", "K02232", "K02225", "K02227", 
                                              "K02231", "K00768", "K02226", "K02233"))


# get ko
cb_ko <- get_ko(cb_path$Car)

# get steps
cb_rxns <- c("R_01", "R_01", rep("R_02", 13),
             rep("R_03", 5), rep("R_04", 6), 
             "R_05", "R_05", "R_06")

#combine
cb_ko_rxns <- cbind(cb_ko, cb_rxns) %>% as.data.frame

## Make data lists ----

#get names
names.vb <- c("tm", "rb", "pt", "pd", "bt", "cb")

#combine kos
ko <- list(tm_ko, rb_ko, pt_ko, pd_ko, bt_ko, cb_ko)
names(ko) <- names.vb

#combine pathways
vb_path <- list(tm_path, rb_path, pt_path, pd_path, bt_path, cb_path)
names(vb_path) <- names.vb

#combine reactions/steps
vb_rxn <- list(tm_rxn, rb_steps, pt_rxns, pd_rxns, bt_steps, cb_rxns)
names(vb_rxn) <- names.vb

#combine ko + reactions/steps
vb_ko_rxn <- list(tm_ko_rxn, rb_ko_steps, pt_ko_rxns, pd_ko_rxns, bt_ko_steps, cb_ko_rxns)

# change to presence / absence
for (i in seq_along(vb_path)) {
  vb_path[[i]]$Car <- pres_abs(vb_path[[i]]$Car)
  vb_path[[i]]$Irb <- pres_abs(vb_path[[i]]$Irb)
  vb_path[[i]]$Irc <- pres_abs(vb_path[[i]]$Irc)
}


## Format expression data ----

names(expression_data) <- c("car", "irb", "irc")

# Get KO expression data
exp_vb <- list()
for(i in seq_along(expression_data)) {
  exp_vb[[i]] <- lapply(ko, function(x) exp_pathway(x, expression_data[[i]]))
}
names(exp_vb) <- c("car", "irb", "irc")


# Clean up annotation column
for(i in seq_along(exp_vb$car)) {
  exp_vb$car[[i]]$Annotation <- sub("\\,.*", "", exp_vb$car[[i]]$Annotation)
  exp_vb$irb[[i]]$Annotation <- sub("\\,.*", "", exp_vb$irb[[i]]$Annotation)
  exp_vb$irc[[i]]$Annotation <- sub("\\,.*", "", exp_vb$irc[[i]]$Annotation)
  exp_vb$car[[i]]$Annotation <- sub("\\&.*", "", exp_vb$car[[i]]$Annotation)
  exp_vb$irb[[i]]$Annotation <- sub("\\&.*", "", exp_vb$irb[[i]]$Annotation)
  exp_vb$irc[[i]]$Annotation <- sub("\\&.*", "", exp_vb$irc[[i]]$Annotation)
}

# Combine duplicate rows (after removing trailing annotations)
for (i in seq_along(exp_vb$car)) {
  exp_vb$car[[i]] <- ddply(exp_vb$car[[i]],c("MAG_id", "Annotation"),numcolwise(sum))
  exp_vb$irb[[i]] <- ddply(exp_vb$irb[[i]],c("MAG_id", "Annotation"),numcolwise(sum))
  exp_vb$irc[[i]] <- ddply(exp_vb$irc[[i]],c("MAG_id", "Annotation"),numcolwise(sum))
}

# Convert to wide format
for(i in seq_along(exp_vb$car)) {
  exp_vb$car[[i]] <- dcast(exp_vb$car[[i]], MAG_id ~ Annotation, value.var = "total")
  exp_vb$car[[i]][is.na(exp_vb$car[[i]])] <- 0  # change NAs to 0 (missing values meant there were no KOs for that MAG)
  exp_vb$irb[[i]] <- dcast(exp_vb$irb[[i]], MAG_id ~ Annotation, value.var = "total")
  exp_vb$irb[[i]][is.na(exp_vb$irb[[i]])] <- 0
  exp_vb$irc[[i]] <- dcast(exp_vb$irc[[i]], MAG_id ~ Annotation, value.var = "total")
  exp_vb$irc[[i]][is.na(exp_vb$irc[[i]])] <- 0
}

# Join the two dataframes
tmp_c <- list()
tmp_i <- list()
tmp_r <- list()
for(i in seq_along(exp_vb$car)) {
  # Car
  names(exp_vb$car[[i]])[1] <- "mag_id"
  tmp_c[[i]] <- left_join(vb_path[[i]][["Car"]], exp_vb$car[[i]], by = "mag_id")
  tmp_c[[i]][is.na(tmp_c[[i]])] <- 0
  #irb
  names(exp_vb$irb[[i]])[1] <- "mag_id"
  tmp_i[[i]] <- left_join(vb_path[[i]][["Irb"]], exp_vb$irb[[i]], by = "mag_id")
  tmp_i[[i]][is.na(tmp_i[[i]])] <- 0
  #irc
  names(exp_vb$irc[[i]])[1] <- "mag_id"
  tmp_r[[i]] <- left_join(vb_path[[i]][["Irc"]], exp_vb$irc[[i]], by = "mag_id")
  tmp_r[[i]][is.na(tmp_r[[i]])] <- 0
}
tmp <- list(tmp_c, tmp_i, tmp_r)
names(tmp) <- c("car", "irb", "irc")


# clean up data frame
df_list <- list()
for(i in seq_along(tmp)) {
  df_list[[i]] <- lapply(tmp[[i]], expression_df)
}
names(df_list) <- c("car", "irb", "irc")

for(i in seq_along(df_list)) {
  names(df_list[[i]]) <- names.vb
}

#Our final dataframe consists of presence/absence matrix, 
# where 0 is no gene found, 1 is gene found in metagenome but not expressed, 
# and 2 is gene in metagenome and expressed 

## Save data ----

# Carbon - note: can delete old car and irb dataframes
saveRDS(df_list, "Formatted_data/vb_exp_data_filtered_formatted_20220923")
df_list <- readRDS("Formatted_data/vb_exp_data_filtered_formatted_20220923")
df_list <- readRDS("Formatted_data/vb_exp_data_formatted") # debug


## Plot heatmap ----

# get mag phylogeny (edit to correct df)
car_tree <- keep.tip(full_tree, df_list$car$tm$mag_id)
irb_tree <- keep.tip(full_tree, df_list$irb$tm$mag_id)
irc_tree <- keep.tip(full_tree, df_list$irc$tm$mag_id)

# clean tax strings in df
for(i in seq_along(df_list)) {
  df_list[[i]] <- lapply(df_list[[i]], clean_tax_strings)
}

# reduce tax strings to fam (if necessary to group plot by fam)
for(i in seq_along(df_list)) {
  df_list[[i]] <- lapply(df_list[[i]], reduce_tax_strings)
}

# Format dataframe
df_car <- list()
df_irb <- list()
df_irc <- list()
for(i in seq_along(df_list$car)) {
  #car
  df_car[[i]] <- df_list$car[[i]][match(car_tree$tip.label, df_list$car[[i]]$mag_id),]
  #irb
  df_irb[[i]] <- df_list$irb[[i]][match(irb_tree$tip.label, df_list$irb[[i]]$mag_id),]
  #irc
  df_irc[[i]] <- df_list$irc[[i]][match(irc_tree$tip.label, df_list$irc[[i]]$mag_id),]
}


# relevel factors and format
for(i in seq_along(df_car)) {
  #car
  df_car[[i]]$mag_id <- factor(df_car[[i]]$mag_id, levels = car_tree$tip.label)
  df_car[[i]]$Tonomy <- factor(df_car[[i]]$Tonomy, levels = unique(df_car[[i]]$Tonomy))
  names(df_car[[i]])[2] <- "Taxonomy"
  #irb
  df_irb[[i]]$mag_id <- factor(df_irb[[i]]$mag_id, levels = irb_tree$tip.label)
  df_irb[[i]]$Tonomy <- factor(df_irb[[i]]$Tonomy, levels = unique(df_irb[[i]]$Tonomy))
  names(df_irb[[i]])[2] <- "Taxonomy"
  #irc
  df_irc[[i]]$mag_id <- factor(df_irc[[i]]$mag_id, levels = irc_tree$tip.label)
  df_irc[[i]]$Tonomy <- factor(df_irc[[i]]$Tonomy, levels = unique(df_irc[[i]]$Tonomy))
  names(df_irc[[i]])[2] <- "Taxonomy"
}

# melt dataframe
for(i in seq_along(df_car)) {
  #car
  df_car[[i]] <- melt(df_car[[i]], id.vars = c("mag_id", "Taxonomy"), value.name = "Expression")
  colnames(df_car[[i]])[3] <- c("KO")
  #irb
  df_irb[[i]] <- melt(df_irb[[i]], id.vars = c("mag_id", "Taxonomy"), value.name = "Expression")
  colnames(df_irb[[i]])[3] <- c("KO")
  #irc
  df_irc[[i]] <- melt(df_irc[[i]], id.vars = c("mag_id", "Taxonomy"), value.name = "Expression")
  colnames(df_irc[[i]])[3] <- c("KO")
}


# Add pathways 
for(i in seq_along(df_car)) {
  #car - match structure with ko id
  df_car[[i]]$Steps <- vb_ko_rxn[[i]][,2][match(df_car[[i]]$KO, vb_ko_rxn[[i]][,1])]
  #irb
  df_irb[[i]]$Steps <- vb_ko_rxn[[i]][,2][match(df_irb[[i]]$KO, vb_ko_rxn[[i]][,1])]
  #irc
  df_irc[[i]]$Steps <- vb_ko_rxn[[i]][,2][match(df_irc[[i]]$KO, vb_ko_rxn[[i]][,1])]
}
names(df_car) <- names.vb
names(df_irb) <- names.vb
names(df_irc) <- names.vb


# Group by taxonomy
df_car2 <- list()
df_irb2 <- list()
df_irc2 <- list()
for(i in seq_along(df_car)) {
  #car
  df_car2[[i]] <- group_by(df_car[[i]], Steps, Taxonomy, KO) %>% 
    dplyr::summarise(Expression = max(Expression)) %>% as.data.frame
  df_car2[[i]]$Expression <- df_car2[[i]]$Expression %>% as.factor
  df_car2[[i]]$Steps <- df_car2[[i]]$Steps %>% as.factor
  #irb
  df_irb2[[i]] <- group_by(df_irb[[i]], Steps, Taxonomy, KO) %>% 
    dplyr::summarise(Expression = max(Expression)) %>% as.data.frame
  df_irb2[[i]]$Expression <- df_irb2[[i]]$Expression %>% as.factor
  df_irb2[[i]]$Steps <- df_irb2[[i]]$Steps %>% as.factor
  #irc
  df_irc2[[i]] <- group_by(df_irc[[i]], Steps, Taxonomy, KO) %>% 
    dplyr::summarise(Expression = max(Expression)) %>% as.data.frame
  df_irc2[[i]]$Expression <- df_irc2[[i]]$Expression %>% as.factor
  df_irc2[[i]]$Steps <- df_irc2[[i]]$Steps %>% as.factor
}
names(df_car2) <- names.vb
names(df_irb2) <- names.vb
names(df_irc2) <- names.vb


# Plot expression

# edit x-axis labels

tm_labs <- c("K03147" = "thiC (K03147)",
             "K00941" = "thiD (K00941)",
             "K14153" = "thiDE (K14153)",
             "K03148" = "thiF (K03148)",
             "K03151" = "thiI (K03151)",
             "K03150" = "thiH (K03150)",
             "K03153" = "thiO (K03153)",
             "K03149" = "thiG (K03149)",
             "K10810" = "tenI (K10810)",
             "K03146" = "THI4 (K03146)",
             "K00788" = "thiE (K00788)",
             "K00946" = "thiL (K00946)",
             "K01077" = "phoA, phoB (K01077)",
             "K00949" = "thiN (K00949)",
             "K02064" = "thiB, tbpA (K02064)",
             "K02063" = "thiP (K02063)",
             "K02062" = "thiQ (K02062)")

rb_labs <- c("K01497" = "ribA (K01497)", 
             "K14652" = "ribBA (K14652)", 
             "K02858" = "ribB (K02858)", 
             "K11752" = "ribD (K11752)",
             "K20862" = "yigB (K20862)", 
             "K00794" = "ribH (K00794)", 
             "K00793" = "ribE (K00793)", 
             "K11753" = "ribF (K11753)")

pt_labs <- c("K00826" = "ilvE (K00826)",
             "K00606" = "panB (K00606)",
             "K00077" = "panE (K00077)",
             "K01579" = "panD (K01579)",
             "K18933" = "mfnA (K18933)",
             "K18966" = "GADL1, CSAD (K18966)",
             "K17839" = "PAO4, PAO3, PAO2 (K17839)",
             "K13367" = "FMS1 (K13367)",
             "K00128" = "ALDH (K00128)",
             "K01918" = "panC (K01918)",
             "K13799" = "panC-cmk (K13799)")

pd_labs <-c("K03472" = "epd (K03472)",
            "K03473" = "pdxB (K03473)",
            "K00831" = "serC (K00831)",
            "K00097" = "pdxA (K00097)",
            "K03474" = "pdxJ (K03474)",
            "K00275" = "pdxH (K00275)",
            "K06215" = "pdxS (K06215)",
            "K08681" = "pdxT (K08681)",
            "K01733" = "thrC (K01733)",
            "K00868" = "pdxK, pdxY (K00868)")

bt_labs <- c("K16593" = "bioI (K16593)",
             "K01906" = "bioW (K01906)",
             "K00652" = "bioF (K00652)",
             "K00833" = "bioA (K00833)",
             "K01935" = "bioD (K01935)",
             "K01012" = "bioB (K01012)",
             "K03523" = "bioY (K03523)",
             "K16783" = "bioN (K16783)",
             "K16784" = "bioM (K16784)",
             "K16785" = "ecfT (K16785)",
             "K16786" = "ecfA1 (K16786)",
             "K16787" = "ecfA2 (K16787)")

cb_labs <- c( "K02303" = "cobA (K02303)",
              "K13542" = "cobA-hemD (K13542)",
              "K02302" = "cysG (K02302)",
              "K02304" = "MET8 (K02304)",
              "K03795" = "cbiX (K03795)",
              "K03394" = "cobI-cbiL (K03394)",
              "K05934" = "cobJ (K05934)",
              "K13541" = "cbiGH-cobJ (K13541)",
              "K05936" = "cobM (K05936)",
              "K02189" = "cbiG (K02189)",
              "K02188" = "cbiD (K02188)",
              "K05895" = "cobK-cbiJ (K05895)",
              "K00595" = "cobL-cbiET (K00595)",
              "K06042" = "cobH-cbiC (K06042)",
              "K02224" = "cobB-cbiA (K02224)",
              "K02229" = "cobG (K02229)",
              "K02228" = "cobF (K02228)",
              "K02230" = "cobN (K02230)",
              "K09882" = "cobS (K09882)",
              "K09883" = "cobT (K09883)",
              "K00798" = "MMAB (K00798)",
              "K19221" = "cobA (K19221)",
              "K02232" = "cobQ (K02232)",
              "K02225" = "cobC1 (K02225)",
              "K02227" = "cbiB (K02227)",
              "K02231" = "cobP (K02231)",
              "K00768" = "cobU, cobT (K00768)",
              "K02226" = "cobC (K02226)",
              "K02233" = "cobS, cobV (K02233)")

x_labels <- list(tm_labs,
                 rb_labs,
                 pt_labs,
                 pd_labs,
                 bt_labs,
                 cb_labs)

names(x_labels) <- names.vb

# heatmap colours
my_cols <- colorRampPalette(c("#132B43", "#56B1F7"))(3)

ggplot(df_irc2$cb, aes(x = KO, y = Taxonomy, fill = Expression)) +
  geom_tile(colour = "black", size=0.10) +
  scale_x_discrete(labels = x_labels$cb) +
  scale_y_discrete(position = "right") +
  scale_fill_manual(values = my_cols, 
                    labels = c("Not found", "Gene found", "Gene expressed")) +
  theme_dark() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.position = "left",
        legend.text = element_text(size = 12),
        legend.key.size = unit(2,"line"),
        legend.title.align = 1) + 
  facet_wrap(df_irc2$cb$Steps, nrow = 1, ncol = 11, scales = "free_x") +
  theme(strip.text = element_text(size = 10)) +
  guides(fill = guide_legend(label.position = "left"))


#thiamine
ggsave("tm_car_exp_fam_f_20220924.svg", device = "svg", width = 44, height = 30, units = "cm", dpi = 300)
ggsave("tm_irb_exp_fam_f_20220924.svg", device = "svg", width = 44, height = 35, units = "cm", dpi = 300)
ggsave("tm_irc_exp_fam_f_20220924.svg", device = "svg", width = 44, height = 38, units = "cm", dpi = 300)

#riboflavin
ggsave("rb_car_exp_fam_f_20220924.svg", device = "svg", width = 35, height = 30, units = "cm", dpi = 300)
ggsave("rb_irb_exp_fam_f_20220924.svg", device = "svg", width = 35, height = 35, units = "cm", dpi = 300)
ggsave("rb_irc_exp_fam_f_20220924.svg", device = "svg", width = 35, height = 38, units = "cm", dpi = 300)

#pantothenate
ggsave("pt_car_exp_fam_f_20220924.svg", device = "svg", width = 35, height = 30, units = "cm", dpi = 300)
ggsave("pt_irb_exp_fam_f_20220924.svg", device = "svg", width = 35, height = 35, units = "cm", dpi = 300)
ggsave("pt_irc_exp_fam_f_20220924.svg", device = "svg", width = 35, height = 38, units = "cm", dpi = 300)

#pyridoxine 
ggsave("pd_car_exp_fam_f_20220924.svg", device = "svg", width = 35, height = 30, units = "cm", dpi = 300)
ggsave("pd_irb_exp_fam_f_20220924.svg", device = "svg", width = 35, height = 35, units = "cm", dpi = 300)
ggsave("pd_irc_exp_fam_f_20220924.svg", device = "svg", width = 35, height = 38, units = "cm", dpi = 300)

#biotin
ggsave("bt_car_exp_fam_f_20220924.svg", device = "svg", width = 39, height = 30, units = "cm", dpi = 300)
ggsave("bt_irb_exp_fam_f_20220924.svg", device = "svg", width = 39, height = 35, units = "cm", dpi = 300)
ggsave("bt_irc_exp_fam_f_20220924.svg", device = "svg", width = 39, height = 38, units = "cm", dpi = 300)

#cobalamin
ggsave("cb_car_exp_fam_f_20220924.svg", device = "svg", width = 55, height = 30, units = "cm", dpi = 300)
ggsave("cb_irb_exp_fam_f_20220924.svg", device = "svg", width = 55, height = 35, units = "cm", dpi = 300)
ggsave("cb_irc_exp_fam_f_20220924.svg", device = "svg", width = 55, height = 38, units = "cm", dpi = 300)



