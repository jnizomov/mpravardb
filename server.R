library(shiny)
library(dplyr)
library(data.table)
library(shinybusy)

source('genetics_processed_clean/prediction.lib.R')
source('genetics_processed_clean/genetics_lib.R')

load("combined_dataset.RData")

# ----------------------------- #

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

spinList <- c("circle", "fading-circle", "half-circle", "double-bounce") # apt spin themes

server <- function(input, output, session) {
  # Update the browse buttons from gray to black
  
  runjs("$('#file1').parent().removeClass('btn btn-default').addClass('btn btn-dark');")
  runjs("$('#file2').parent().removeClass('btn btn-default').addClass('btn btn-dark');")
  
  filtered_data <- combined_dataset
  fileSubmitted <- NA
  isFileAvailable <- FALSE
  
  observeEvent(input$paper, {
    dt <- combined_dataset[combined_dataset$source_dataset == input$paper, ]
    
    # Update the choices in the 'dc' (disease/cell type) input
    
    updateSelectInput(
      session, 
      "dc",
      choices = c("Select a disease/cell type:" = "", unique(dt$celltype))
    )
  })
  
  output$exampleDownload <- downloadHandler(
    filename = "Example.txt",
    
    content = function(file){
      file.copy("file_input_examples/Example.txt", file)
    }
  )
  
  observeEvent(input$loadFile, {
    isFileAvailable <- !is.null(input$file1$datapath)
    
    if (!isFileAvailable) {
      shinyjs::runjs('$("#checkmark").attr("src", "error.svg");')
      shinyjs::runjs('$("#checkmark").css("visibility", "visible");')
      
      return ();
    } else {
      shinyjs::runjs('$("#checkmark").attr("src", "checkmark.svg");')
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
    if (!is.null(input$file1$datapath)) { # Checks if the input file exists
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
  
  ### -------------- ###
  ###    ANALYSIS    ###
  ### -------------- ###
  
  analysisFileLoaded <- FALSE
  
  uploadedAnalysisFile <- reactive({
    return (if (is.null(input$file2)) NULL else (fread(input$file2$datapath, header = TRUE, sep = "\t")))
  })
  
  observeEvent(input$loadAnalysisFile, {
    #fileUploaded <- fread(input$file2$datapath, header = TRUE, sep = "\t")
    #y_hats <- vector("numeric", length = nrow(fileUploaded))
    
    show_modal_spinner(
      spin = "half-circle",
      color = "#2372CA",
    )
    
    selected_keyword <- analysis_dataset$keyword[which(analysis_dataset$source_dataset == input$paper)[1]]
    
    if (length(selected_keyword) == 0 || is.na(selected_keyword)) {
      selected_keyword <- NA_character_ 
    }
    
    model.filename <- paste0(selected_keyword, '-', input$dc, '-', input$model)
    
    data <- if (input$analysisFileType == "BED") 
                fread(input$file2$datapath, fill = FALSE) 
            else 
                readDNAStringSet(input$file2$datapath)
    
    data$probabilities <- getProbability(data, model.filename, input$model, input$analysisFileType)
    
    remove_modal_spinner()
    
    output$probabilitiesOutput <- renderDT({data})
  })
}

