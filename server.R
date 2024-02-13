library(shiny)
library(dplyr)
library(data.table)
library(shinybusy)
library(Biostrings)

load("combined_dataset.RData")
load("analysis_dataset.RData")

source('genetics_processed_clean/prediction.lib.R')
source('genetics_processed_clean/model_extraction.R')
#source('genetics_processed_clean/genetics_lib.R')

# ----------------------------- #

getPaperByKeyword = function(keyword) {
  return (analysis_dataset$source_dataset[which(analysis_dataset$keyword == keyword)[1]]);
}

getKeywordByPaper = function(paper) {
  return (analysis_dataset$keyword[which(analysis_dataset$source_dataset == paper)[1]]);
}

getFilteredDataset = function(input_chr, input_start, input_end, input_disease, input_cell) {
  # Additional server logic as needed
  
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

# -------------------------- #
#       Server Logic 
# -------------------------- #

getFileNameFromInputs <- function(disease, cell_line, chr, start_position, end_position) {
  parts <- c()

  if (!is.na(disease) && !is.null(disease) && disease != "") parts <- c(parts, disease)
  if (!is.na(cell_line) && !is.null(cell_line) && cell_line != "") parts <- c(parts, cell_line)
  if (!is.na(chr) && !is.null(chr) && chr != "") parts <- c(parts, chr)
  if (!is.na(start_position) && !is.null(start_position) && start_position != "") parts <- c(parts, start_position)
  if (!is.na(end_position) && !is.null(end_position) && end_position != "") parts <- c(parts, end_position)
  
  filename <- ifelse(length(parts) > 0, paste(parts, collapse = "_"), "mpra_data")
  
  return (paste0(filename, ".csv"));
}

server <- function(input, output, session) {
  # Update the browse buttons from gray to black to match theme
  
  runjs("$('#customFile').parent().removeClass('btn btn-default').addClass('btn btn-dark');")
  runjs("$('#analysisFile').parent().removeClass('btn btn-default').addClass('btn btn-dark');")
  
  browse_filtered_data <- reactiveVal(combined_dataset)
  custom_filtered_data <- reactiveVal(data.frame())
  
  fileSubmitted <- reactiveVal(NA)
  customFileExists <- reactiveVal(FALSE)
  
  browseQueryRun <- reactiveVal(FALSE)
  customQueryRun <- reactiveVal(FALSE)
  probabilitiesQueryRun <- reactiveVal(FALSE)
  
  ### ------------------- ###
  ###    DATABASE TAB     ###
  ### ------------------- ###
  
  # Browse Database Section
  
  output$browseDownloadTable <- downloadHandler(
    filename = function() {
      getFileNameFromInputs(input$disease, input$cell_line, input$chr, input$start_position, input$end_position)
    },
    
    content = function(file) {
      write.csv(browse_filtered_data(), file, row.names = FALSE)
    }
  )
  
  output$browseDownloadButton <- renderUI({
    if(browseQueryRun()) {
      downloadButton(
        outputId = "browseDownloadTable",
        label = "Download",
        class = "btn btn-dark",
        style = "width:100%;"
      )
    }
  })
  
  output$browseQueryTable <- renderUI({
    if(browseQueryRun()) {
      renderDT(
        browse_filtered_data(), 
        options = list(scrollX = TRUE)
      )
    }
  })
  
  observeEvent(input$browseQuery, {
    browse_filtered_data(getFilteredDataset(input$chr, input$start_position, input$end_position, input$disease, input$cell_line))
    
    browseQueryRun(TRUE)
  })
  
  observeEvent(input$cell_line, {
    if (input$cell_line == "") {
      return();
    }
    
    disease_choices = unique(combined_dataset$disease[combined_dataset$celltype == input$cell_line])
    
    if (input$disease %in% disease_choices) {
      updateSelectInput(
        session,
        "disease",
        choices = c("Select a disease" = "", disease_choices),
        selected = input$disease
      )
    } else {
      updateSelectInput(
        session,
        "disease",
        choices = c("Select a disease" = "", disease_choices),
        selected = NULL
      )
    }
  })
  
  observeEvent(input$disease, {
    if (input$disease == "") {
  
      return();
    }
    
    cell_choices = unique(combined_dataset$celltype[combined_dataset$disease == input$disease])
    
    if (input$cell_line %in% cell_choices) {
      updateSelectInput(
        session,
        "cell_line",
        choices = c("Select a cell line" = "", cell_choices),
        selected = input$cell_line
      )
    } else {
      updateSelectInput(
        session,
        "cell_line",
        choices = c("Select a cell line" = "", cell_choices),
        selected = NULL
      )
    }
  })
  
  observeEvent(input$reset, {
    updateNumericInput(session, "cell_line", value = "")
    updateNumericInput(session, "disease", value = "")
    updateNumericInput(session, "chr", value = "")
    updateNumericInput(session, "start_position", value = "")
    updateNumericInput(session, "end_position", value = "")
    updateNumericInput(session, "customFile", value = "")
    
    browseQueryRun(FALSE)
    
    updateSelectInput(
      session,
      "disease",
      choices = c("Select a disease" = "", unique(combined_dataset$disease[!combined_dataset$disease %in% c("None", "NA")])),
    )
    
    updateSelectInput(
      session,
      "cell_line",
      choices = c("Select a cell line" = "", unique(combined_dataset$celltype)),
    )
    
    shinyjs::runjs('$("#checkmark").css("visibility", "hidden");')
  })

  
  # Upload a Custom File Section
  
  output$exampleDownload <- downloadHandler(
    filename = "Example.txt",
    
    content = function(file){
      file.copy("file_input_examples/Example.txt", file)
    }
  )
  
  shinyjs::runjs('
    $("#customFile").on("click", function() {
      $("#checkmark").css("visibility", "hidden");
    });
  ')
  
  observeEvent(input$loadFile, {
    customFileExists(!is.null(input$customFile$datapath))
    customQueryRun(FALSE)
    
    if (!customFileExists()) {
      shinyjs::runjs('$("#checkmark").attr("src", "error.svg");')
      shinyjs::runjs('$("#checkmark").css("visibility", "visible");')
      
      return ();
    } else {
      shinyjs::runjs('$("#checkmark").attr("src", "checkmark.svg");')
      shinyjs::runjs('$("#checkmark").css("visibility", "visible");')
    }
    
    fileSubmitted(fread(input$customFile$datapath, fill=FALSE))
    
    output$filePreview <- renderTable({
      return (head(fileSubmitted()))
    })
    
    output$previewUI <- renderUI({
      if (grepl("\\.txt$", input$customFile$name)) {
        fluidPage(
          h1(strong("Preview your file"), style = "font-size:20px; margin-top: 5px"),
          div(class = "scrollable-table", tableOutput("filePreview"))
        )
      }
    })
  })
  
  output$customDownloadTable <- downloadHandler(
    filename = paste0(ifelse (customFileExists(), input$customFile$name, "custom_query"), ".csv"),
    
    content = function(file) {
      write.csv(custom_filtered_data(), file, row.names = FALSE)
    }
  )
  
  output$customDownloadButton <- renderUI({
    if(customQueryRun()) {
      downloadButton(
        outputId = "customDownloadTable",
        label = "Download",
        class = "btn btn-dark",
        style = "width:100%;"
      )
    }
  })
  
  output$customQueryTable <- renderUI({
    if(customQueryRun()) {
      renderDT(
        custom_filtered_data(), 
        options = list(scrollX = TRUE)
      )
    }
  })
  
  observeEvent(input$customQuery, {
    if (!is.null(input$customFile$datapath)) { 
      aggregated_data <- data.frame()
      
      for(i in 1:nrow(fileSubmitted())) {
        row <- fileSubmitted()[i,]
        
        disease_column = NA # optional 
        cell_line_column = NA # optional 
        
        if (ncol(fileSubmitted()) >= 4) {
          disease_column <- row[[4]]
        }
        
        if (ncol(fileSubmitted()) >= 5) {
          cell_line_column <- row[[5]]
        }
        
        filter_result <- getFilteredDataset(row[[1]], row[[2]], row[[3]], disease_column, cell_line_column)
        
        aggregated_data <- rbind(aggregated_data, filter_result)
      }
      
      custom_filtered_data(aggregated_data)
    }
    
    customQueryRun(TRUE)
  })
  
  ### ------------------- ###
  ###    ANALYSIS TAB     ###
  ### ------------------- ###
  
  model.filename <- reactiveVal("")
  probabilitiesOutput <- reactiveVal(data.frame())
  
  observeEvent(input$paper, {
    dt <- analysis_dataset[analysis_dataset$source_dataset == input$paper, ]
    
    # Update the choices in the 'disease_celltype' input based on selected paper
    
    updateSelectInput(
      session, 
      "disease_celltype",
      choices = c("Select a disease/cell type:" = "", unique(models_df$disease_celltype[models_df$keyword == getKeywordByPaper(input$paper)]))
    )
  })
  
  observeEvent(input$disease_celltype, {
    # Update the choices in the 'model' input based on selected disease_celltype (validation)
    
    updateSelectInput(
      session, 
      "model",
      choices = c("Select a model type:" = "", 
                  unique(models_df$model_type[models_df$keyword == getKeywordByPaper(input$paper) & models_df$disease_celltype == input$disease_celltype])
      )
    )
  }) 
  
  observeEvent(input$loadAnalysisFile, {
    probabilitiesQueryRun(FALSE)
    
    show_modal_spinner(
      spin = "half-circle",
      color = "#2372CA",
    )
    
    selected_keyword <- getKeywordByPaper(input$paper)
    
    if (length(selected_keyword) == 0 || is.na(selected_keyword)) {
      selected_keyword <- NA_character_ 
    }
    
    model.filename(paste0(selected_keyword, '-', input$disease_celltype, '-', input$model))
    
    if (input$analysisFileType == "BED") {
      data <- fread(input$analysisFile$datapath, header = TRUE, sep = "\t", fill = FALSE) 
      probabilities <- getProbability(data, model.filename(), input$model, input$analysisFileType, input$genome)
      
      data$probabilities <- probabilities
    } else if (input$analysisFileType == "FASTA") {
      data <- readDNAStringSet(input$analysisFile$datapath)
      probabilities <- getProbability(data, model.filename(), input$model, input$analysisFileType, input$genome)
      
      data <- data.frame(sequence = data, probabilities)
    }
    
    probabilitiesOutput(data) # sets reactive probabilityOutput to data for downloading purposes
    
    remove_modal_spinner()
    
    probabilitiesQueryRun(TRUE)
  })
  
  output$probabilitiesDownloadTable <- downloadHandler(
    filename = function() {
      paste0(model.filename(), ".csv")
    },
    
    content = function(file) {
      write.csv(probabilitiesOutput(), file, row.names = FALSE)
    }
  )
  
  output$probabilitiesDownloadButton <- renderUI({
    if(probabilitiesQueryRun()) {
      downloadButton(
        outputId = "probabilitiesDownloadTable",
        label = "Download",
        class = "btn btn-dark",
        style = "width:100%;"
      )
    }
  })
  
  output$probabilitiesOutput <- renderUI({
    if(probabilitiesQueryRun()) {
      renderDT(
        probabilitiesOutput(), 
        options = list(scrollX = TRUE)
      )
    }
  })
}

