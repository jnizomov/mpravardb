# library(dplyr)
# library(openxlsx)
# library(readxl)
# 
# # --------------- #
# 
# # Convert columns to appropriate data types
# 
# convert_columns <- function(df) {
#   if ("chr" %in% names(df)) df$chr <- as.character(df$chr)
#   if ("pos" %in% names(df)) df$pos <- as.numeric(df$pos)
#   if ("ref" %in% names(df)) df$ref <- as.character(df$ref)
#   if ("alt" %in% names(df)) df$alt <- as.character(df$alt)
#   if ("rsid" %in% names(df)) df$rsid <- as.character(df$rsid)
#   if ("pvalue" %in% names(df)) df$pvalue <- as.numeric(df$pvalue)
#   if ("fdr" %in% names(df)) df$fdr <- as.numeric(df$fdr)
#   if ("disease" %in% names(df)) df$disease <- as.character(df$disease)
#   if ("celltype" %in% names(df)) df$celltype <- as.character(df$celltype)
#   if ("genome" %in% names(df)) df$genome <- as.character(df$genome)
#   if ("GWASorQTL" %in% names(df)) df$GWASorQTL <- as.character(df$GWASorQTL)
#   if ("gene_start" %in% names(df)) df$gene_start <- as.numeric(df$gene_start)
#   if ("gene_end" %in% names(df)) df$gene_end <- as.numeric(df$gene_end)
#   if ("gene_strand" %in% names(df)) df$gene_strand <- as.character(df$gene_strand)
#   if ("insideGene" %in% names(df)) df$insideGene <- as.character(df$insideGene)
#   if ("distancetoGene" %in% names(df)) df$distancetoGene <- as.numeric(df$distancetoGene)
#   if ("gene_name" %in% names(df)) df$gene_name <- as.character(df$gene_name)
#   if ("log2FC" %in% names(df)) df$log2FC <- as.numeric(df$log2FC)
# 
#   return (df)
# }
# 
# # Process .xlsx files in genetics_processed_clean
# 
# datas <- list.files('genetics_processed_clean/')
# datas <- datas[grep(".xlsx$",datas)]
# 
# dataset_names <- character(0)
# 
# for (name in datas) {
#   df <- read_xlsx(paste0("genetics_processed_clean/", name))
#   df <- convert_columns(df)
# 
#   cleaned_name <- gsub(".xlsx$", "", name)
#   cleaned_name <- gsub("[^[:alnum:]]", "_", cleaned_name)
# 
#   dataset_names <- c(dataset_names, cleaned_name)
# 
#   assign(cleaned_name, df)
# }
# 
# rm(df)
# 
# # Combine all datasets into one
# 
# combined_dataset <- bind_rows(
#   mget(dataset_names),
#   .id = "source_dataset"
# )
# 
# # Empty cells initialized to 1 and 0 for FDR and log2FC respectively
# 
# combined_dataset$fdr[is.na(combined_dataset$fdr)]=1
# combined_dataset$log2FC[is.na(combined_dataset$log2FC)]=0
# 
# # Clean chromosome column (i.e. convert "chr1" to just "1")
# 
# combined_dataset$chr <- gsub("^chr", "", combined_dataset$chr)
# 
# combined_dataset <- combined_dataset %>%
#   relocate(source_dataset, .after = last_col())
# 
# # ---------- #
# #   Export   #
# # ---------- #
# 
# save(combined_dataset, file = "combined_dataset.RData")
# 
# # Create analysis dataset for ML analysis section
# 
# analysis_dataset <- combined_dataset %>%
#   dplyr::mutate(keyword = case_when(
#     source_dataset == "Functional_regulatory_variants_implicate_distinct_transcriptional_networks_in_dementia" ~ "dementia",
#     source_dataset == "Genome_wide_functional_screen_of_30_UTR_variants_uncovers_causal_variants_for_human_disease_and_evolution" ~ "evolution3UTR",
#     source_dataset == "Transcriptional_regulatory_convergence_across_functional_MDD_risk_variants_identified_by_massively_parallel_reporter_assays" ~ "MDD",
#     source_dataset == "Saturation_mutagenesis_of_twenty_disease_associated_regulatory_elements_at_single_base_pair_resolution_GRCh37_ALL" ~ "mutagenesis",
#     source_dataset == "Prioritization_of_autoimmune_disease_associated_genetic_variants_that_perturb_regulatory_element_activity_in_T_cells" ~ "autoimmune",
#     source_dataset == "Massively_parallel_reporter_assays_and_variant_scoring_identified_functional_variants_and_target_genes_for_melanoma_loci_and_highlighted_cell_type_specificity" ~ "melanoma",
#     TRUE ~ NA
#   ))
# 
# analysis_dataset <- analysis_dataset %>%
#   dplyr::filter(!is.na(keyword))