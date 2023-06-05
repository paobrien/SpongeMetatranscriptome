## Functions used for genome-centric metatranscriptomic analysis

# load functions when necessary using souce("R_functions_metatranscriptomics)

# function to get ko for steps
get_ko <- function(df) {
  ko <- colnames(df)[3:length(df)]
  return(ko)
}

# change to presence-absence
pres_abs <- function(df) {
  tmp <- ifelse(df[, -c(1:2)] > 0, 1, 0)
  df <- cbind(df[,1:2], tmp)
  return(df)
}

# subset data and combine reads column
subset_combine <- function(df) {
  df$Reads_aligned_total <- rowSums(df[,grep("Reads_aligned", colnames(df))])
  df <- select(df, "MAG_id", "Annotation", "Reads_aligned_total")
}

# subset expression data to KOs of interest
exp_pathway <- function(ko, exp_data) {
  tmp.l <- list()
  for (i in seq_along(ko)) {
    print(ko[i])
    tmp <- exp_data[grep(ko[i], exp_data$Annotation),]
    tmp.l[[i]] <- tmp 
  }
  exp_combined <- do.call(rbind, tmp.l)
  return(exp_combined)
}


# function to clean up dataframe
# if KO is expressed in expression column, change metagenome column value to 2
# then remove old expression columns

expression_df <- function(df) {
  i3 <- vector()
  for(i in grep(pattern = "\\.x", colnames(df), value = T)) {
    print(i)
    i2 <- gsub(".x", ".y", i)
    print(i2)
    df[,i] <- ifelse(df[,i2] > 0, 2, df[,i])
    i3 <- c(i3, i2)
  }
  df_final <- df[, !colnames(df) %in% i3]
  colnames(df_final) <- gsub(".x", "", colnames(df_final))
  return(df_final)
}


# clean taxonomy strings in dataframe
# $ matches end of string
clean_tax_strings <- function(df) {
  df$Tonomy <- gsub(";s__$", "", df$Tonomy)
  df$Tonomy <- gsub(";g__$", "", df$Tonomy)
  df$Tonomy <- gsub(";f__$", "", df$Tonomy)
  df$Tonomy <- gsub(";o__$", "", df$Tonomy)
  df$Tonomy <- gsub("d__Bacteria;p__", "", df$Tonomy)
  df$Tonomy <- gsub("d__Archaea", "Archaea", df$Tonomy)
  return(df)
} 

clean_tax_strings2 <- function(df) {
  df$Taxonomy <- gsub(";s__$", "", df$Taxonomy)
  df$Taxonomy <- gsub(";g__$", "", df$Taxonomy)
  df$Taxonomy <- gsub(";f__$", "", df$Taxonomy)
  df$Taxonomy <- gsub(";o__$", "", df$Taxonomy)
  df$Taxonomy <- gsub("d__Bacteria;p__", "", df$Taxonomy)
  df$Taxonomy <- gsub("p__", "", df$Taxonomy)
  df$Taxonomy <- gsub("d__Archaea", "Archaea", df$Tonomy)
  return(df)
} 


# To group by family, remove all strings after ;g (gives same row name from fam)
reduce_tax_strings <- function(df) {
  df$Tonomy <- gsub(";g.*", "", df$Tonomy)
  return(df)
}

reduce_tax_strings_cazy <- function(df) {
  df$Taxonomy <- gsub(";g.*", "", df$Taxonomy)
  return(df)
}









