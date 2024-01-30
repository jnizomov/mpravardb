library(shiny)
library(dplyr)
library(randomForest)
library(data.table)

load("combined_dataset.RData")

source("themes.R", local = TRUE)
#source("predictions.R", local = TRUE)

source('genetics_processed_clean/prediction.lib.R')
source('genetics_processed_clean/genetics_lib.R')

#MPRAVarDB.shinny.io

getFilteredDataset = function(input_chr, input_start, input_end, input_disease, input_cell) {
  filtered_data <- combined_dataset
  
  if (input_cell != "" && !is.na(input_cell)) {
    filtered_data <- filtered_data %>% 
      dplyr::filter(celltype == input_cell)
  }
  
  if (input_disease != "" && !is.na(input_disease)) {
    if (input_disease == "NA") {
      filtered_data <- filtered_data %>%
        dplyr::filter(is.na(disease))
    } else {
      filtered_data <- filtered_data %>% 
        dplyr::filter(disease == input_disease)
    }
  }
  
  if (!is.na(input_chr)) {
    filtered_data <- filtered_data %>% 
      dplyr::filter(chr == input_chr)
  }
  
  if (!is.na(input_start) && !is.na(input_end)) {
    filtered_data <- filtered_data %>% 
      dplyr::filter(pos >= input_start & pos <= input_end)
  }
  
  return (filtered_data)
}

server <- function(input, output, session) {
 # bs_themer() # preview different themes in real time, remember to change the 'theme =' in the ui.R file to right one 
  
  # Update browse buttons to right colors
  
  runjs("$('#file1').parent().removeClass('btn btn-default').addClass('btn btn-dark');")
  runjs("$('#file2').parent().removeClass('btn btn-default').addClass('btn btn-dark');")
  
  filtered_data <- combined_dataset
  fileSubmitted <- NA
  isFileAvailable <- FALSE
  
  observeEvent(input$paper, {
    dt <- combined_dataset[combined_dataset$source_dataset == input$paper, ]
    
    # Update the choices in the 'dc' (disease/cell type) input
    updateSelectInput(session, "dc",
                      choices = c("Select a disease/cell type:" = "", 
                                  unique(dt$celltype)))
  })
  
  output$exampleDownload <- downloadHandler(
    filename = function(){
      paste("file_input_examples/RShinyInputTest","txt",sep=".")
    },
    content = function(file){
      file.copy("file_input_examples/RShinyInputTest.txt", file)
    }
  )
  
  observeEvent(input$loadFile, {
    isFileAvailable <- !is.null(input$file1$datapath)
    
    if (!isFileAvailable) {
      shinyjs::runjs('$("#checkmark").attr("src", "error.svg");')
      shinyjs::runjs('$("#checkmark").css("visibility", "visible");')
      return();
    } else {
      shinyjs::runjs('$("#checkmark").attr("src", "Checkmark.svg");')
      shinyjs::runjs('$("#checkmark").css("visibility", "visible");')
    }
    
    fileSubmitted <<- fread(input$file1$datapath, fill=FALSE)
    
    output$filePreview <- renderTable({
      return (head(fileSubmitted))
    })
    
    output$previewUI <- renderUI({
      if (grepl("\\.txt$", input$file1$name)) {
        fluidPage(
          h1(strong("Preview your file"), style = "font-size:20px; margin-top: 5px"),
          div(class = "scrollable-table", tableOutput("filePreview"))
        )
      }
    })
  })
  
  observeEvent(input$queryData, {
    isFileAvailable <- !is.null(input$file1$datapath)
    
    if (isFileAvailable) {
      aggregated_data <- data.frame()
      
      for(i in 1:nrow(fileSubmitted)) {
        row <- fileSubmitted[i,]
        
        disease_column = NA # optional 
        cell_line_column = NA # optional 
        
        if (ncol(fileSubmitted) >= 4) {
          disease_column <- row[[4]]
        }
        
        if (ncol(fileSubmitted) >= 5) {
          cell_line_column <- row[[5]]
        }
        
        filter_result <- getFilteredDataset(row[[1]], row[[2]], row[[3]], disease_column, cell_line_column)
        
        aggregated_data <- rbind(aggregated_data, filter_result)
      }
      
      filtered_data <- aggregated_data
    } else {
      filtered_data <- getFilteredDataset(input$chr, input$start_position, input$end_position, input$disease, input$cell_line)
    }
    
    output$table <- renderDT({filtered_data}, options=list(scrollX=T))
  })
  
  observeEvent(input$reset, {
    updateNumericInput(session, "cell_line", value = "")
    updateNumericInput(session, "disease", value = "")
    updateNumericInput(session, "chr", value = "")
    updateNumericInput(session, "start_position", value = "")
    updateNumericInput(session, "end_position", value = "")
    updateNumericInput(session, "file1", value = "")
    
    shinyjs::runjs('$("#checkmark").css("visibility", "hidden");')
  })
  
  ### ANALYSIS PAGE SERVER CODE
  
  analysisFileLoaded <- FALSE
  
  uploadedAnalysisFile <- reactive({
    if (is.null(input$file2)) {
      return (NULL)
    }
    
    return (fread(input$file2$datapath, header = TRUE, sep = "\t"))
  })
  
  observeEvent(input$loadAnalysisFile, {
    y_hats <- vector("numeric", length = nrow(analysisFile))
    
    keyword = ""
    
    if (input$paper == "Functional_regulatory_variants_implicate_distinct_transcriptional_networks_in_dementia") {
      keyword = "dementia"
    } else if (input$paper == "Genome_wide_functional_screen_of_30_UTR_variants_uncovers_causal_variants_for_human_disease_and_evolution") {
      keyword = "evolution3UTR"
    } else if (input$paper == "Transcriptional_regulatory_convergence_across_functional_MDD_risk_variants_identified_by_massively_parallel_reporter_assays") {
      keyword = "MDD"  
    } else if (input$paper == "Saturation_mutagenesis_of_twenty_disease_associated_regulatory_elements_at_single_base_pair_resolution_GRCh37_ALL") {
      keyword == "mutagenesis"
    } else if (input$paper == "Prioritization_of_autoimmune_disease_associated_genetic_variants_that_perturb_regulatory_element_activity_in_T_cells") {
      keyword == "autoimmune"
    } else if (input$paper == "Massively_parallel_reporter_assays_and_variant_scoring_identified_functional_variants_and_target_genes_for_melanoma_loci_and_highlighted_cell_type_specificity") {
      keyword == "melanoma"
    }
    
    if (input$analysisFileType == "BED") {
      data <- fread(input$file2$datapath, fill=FALSE)
      
      probabilities <- getProbability(data, paste0(keyword, '-', input$dc, '-', input$model), input$model, input$analysisFileType)
    } else if (input$analysisFileType == 'FASTA') {
      data <- readDNAStringSet(input$file2$datapath)
      
      print(data)
      fac
      probabilities <- getProbability(data, paste0(keyword, '-', input$dc, '-', input$model), input$model, input$analysisFileType)
    }
    
    output$probabilitiesOutput <- renderDT({
      data.frame(Probability = probabilities)
    })
  })
}

