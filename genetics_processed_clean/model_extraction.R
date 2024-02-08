# extracts models in models folder and splits by dashes to creates df of keyword, disease/celltype, and model type 

all_models <- list.files(path = "models/", pattern = "\\.rds$", full.names = TRUE)

extract_model_details <- function(filename) {
  clean_name <- sub("\\.rds$", "", filename)
  parts <- strsplit(clean_name, "-")[[1]]
  
  keyword <- parts[1]
  
  model_type <- parts[length(parts)]
  disease_celltype <- paste(parts[-c(1, length(parts))], collapse = "-")
  
  return (list(keyword = keyword, disease_celltype = disease_celltype, model_type = model_type))
}

model_details_list <- lapply(all_models, function(model) {
  filename <- basename(model) 
  extract_model_details(filename)
})

models_df <- do.call(rbind, lapply(model_details_list, function(x) as.data.frame(t(x), stringsAsFactors = FALSE)))
models_df <- as.data.frame(lapply(models_df, trimws), stringsAsFactors = FALSE)
colnames(models_df) <- c("keyword", "disease_celltype", "model_type")
row.names(models_df) <- NULL

#save(models_df, file = "models_df.RData")