library(dplyr)
library(openxlsx)
library(readxl)

convert_columns <- function(df) {
  if ("chr" %in% names(df)) df$chr <- as.character(df$chr)
  if ("pos" %in% names(df)) df$pos <- as.numeric(df$pos)
  if ("ref" %in% names(df)) df$ref <- as.character(df$ref)
  if ("alt" %in% names(df)) df$alt <- as.character(df$alt)
  if ("rsid" %in% names(df)) df$rsid <- as.character(df$rsid)
  if ("pvalue" %in% names(df)) df$pvalue <- as.numeric(df$pvalue)
  if ("fdr" %in% names(df)) df$fdr <- as.numeric(df$fdr)
  if ("disease" %in% names(df)) df$disease <- as.character(df$disease)
  if ("celltype" %in% names(df)) df$celltype <- as.character(df$celltype)
  if ("genome" %in% names(df)) df$genome <- as.character(df$genome)
  if ("GWASorQTL" %in% names(df)) df$GWASorQTL <- as.character(df$GWASorQTL)
  if ("gene_start" %in% names(df)) df$gene_start <- as.numeric(df$gene_start)
  if ("gene_end" %in% names(df)) df$gene_end <- as.numeric(df$gene_end)
  if ("gene_strand" %in% names(df)) df$gene_strand <- as.character(df$gene_strand)
  if ("insideGene" %in% names(df)) df$insideGene <- as.character(df$insideGene)
  if ("distancetoGene" %in% names(df)) df$distancetoGene <- as.numeric(df$distancetoGene)
  if ("gene_name" %in% names(df)) df$gene_name <- as.character(df$gene_name)
  if ("log2FC" %in% names(df)) df$log2FC <- as.numeric(df$log2FC)
  
  return(df)
}

datas=list.files('genetics_processed_clean/')
datas=datas[grep(".xlsx$",datas)]

dataset_names = character(0)  # Initialize an empty character vector to store the dataset-Names

for (name in datas) {
  df <- read_xlsx(paste0("genetics_processed_clean/", name))
  df <- convert_columns(df)
  
  cleaned_name <- gsub(".xlsx$", "", name)
  cleaned_name <- gsub("[^[:alnum:]]", "_", cleaned_name)
  
  dataset_names <- c(dataset_names, cleaned_name)
  
  assign(cleaned_name, df)
}

combined_dataset <- bind_rows(
  mget(dataset_names),
  .id = "source_dataset"
)

combined_dataset$fdr[is.na(combined_dataset$fdr)]=1
combined_dataset$log2FC[is.na(combined_dataset$log2FC)]=0

combined_dataset$chr <- gsub("^chr", "", combined_dataset$chr)

write.csv(combined_dataset, "combined_dataset.csv", row.names = FALSE)

combined_dataset <- combined_dataset %>%
  relocate(source_dataset, .after = last_col())


save(combined_dataset, file = "/Users/javlon_nizomov/Documents/Research/PHHP Honors Thesis/ML_Project/combined_dataset.RData")
