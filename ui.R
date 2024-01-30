library(DT)
library(shinyjs)
library(shinythemes)
library(bslib)
library(dplyr)
library(shinyWidgets)

load("combined_dataset.RData")

source("themes.R", local = TRUE)

project_title <- "MRPAVarDB"

database_page <- tabPanel(
  title = div(icon("database", style = "margin-right: 5px;"), "Database"),
  
  radioGroupButtons(
    inputId = "Id01",
    label = "",
    choices = c("Browse Database", "Upload a Custom File"),
    justified = TRUE,
    #status = "btn btn-secondary"
  ),
  
  tags$div(class = "rounded-grey-square",
    conditionalPanel(
      condition = "input.Id01 == 'Browse Database'",
      tabPanel(
        title = "Custom Inputs",
        fluidRow(
          column(width = 6,
             selectInput("disease", "Disease:", choices = c("Select a disease" = "", unique(combined_dataset$disease)), width = '100%')
          ),
          
          column(width = 6, 
             #tags$div(class = "header-spacing"),
             selectInput("cell_line", "Cell Line:", choices = c("Select a cell line" = "", unique(combined_dataset$celltype)), width = '100%'),
          )
        ),
        fluidRow(
          column(width = 4,
              numericInput("chr", "Chromosome:", value = NULL, width = '100%'),
          ),
          
          column(width = 4,
              numericInput("start_position", "Start Position:", value = NULL, width = '100%'),
          ),
          
          column(width = 4,
              numericInput("end_position", "End Position:", value = NULL, width = '100%'),
          ),
        ),
        
        actionButton("reset", "Reset", class = "btn btn-dark")
      ),
    ),
    conditionalPanel(
      condition = "input.Id01 == 'Upload a Custom File'",
      tabPanel(
        title = "Upload file",
        fluidRow(
          column(width = 7, 
             tags$style(".shiny-file-input-progress {display: none}"),
             tags$b("File Instructions"),
             tags$div(class = "text-style",
               HTML("There are three required columns and two optional columns."),
               p(""),
               tags$div(class = "center-div",
                  tags$table(
                    tags$thead(
                      # Additional row for "Required" and "Optional" headers
                      tags$tr(
                        tags$th("Required", colspan = 3),
                        tags$th(tags$span(style = "color: gray;", "Optional"), colspan = 2)
                      ),
                      tags$tr(
                        tags$th("chr"),
                        tags$th("start"),
                        tags$th("end"),
                        tags$th(tags$span(style = "color: gray;", "disease")),
                        tags$th(tags$span(style = "color: gray;", "celltype"))
                      )
                    ),
                    
                    tags$tbody(
                      tags$tr(
                        tags$td("chr1"),
                        tags$td(1),
                        tags$td(3),
                        tags$td("Disease1"),
                        tags$td("Type1")
                      ),
                      tags$tr(
                        tags$td("chr2"),
                        tags$td(4),
                        tags$td(6),
                        tags$td("Disease2"),
                        tags$td("Type2")
                      ),
                      tags$tr(
                        tags$td("chr3"),
                        tags$td(7),
                        tags$td(9),
                        tags$td("Disease3"),
                        tags$td("Type3")
                      )
                    )
                  ),
                p(""),
               )
             ),
             tags$div(
               class = "text-style",
               HTML("Need to see an example file?<br>"),
             ),
             tags$div(style = "margin-bottom: 12px;"),
             downloadButton("exampleDownload", "Download example .txt", class = "btn btn-dark"),
             tags$div(style = "margin-bottom: 12px;"),
             tags$hr(),
             fileInput("file1", "Upload a .txt file", accept = c(".txt")),
             includeCSS("www/app.css"),
             actionButton("loadFile", "Load File", class = "btn btn-dark"),
             img(id = "checkmark", src = "", height = 20, width = 20, margin = 10, style = "visibility: hidden; margin-left: 5px; margin-bottom: 5px"),
             p(""),
          ),
          
          column(width = 5, 
             tags$div(class = "preview-back",
                tags$b("Preview file"),
                tags$div(class = "scrollable-table", tableOutput("filePreview"))
             )
          ),
        ),
      )
    ),
  ),

  p(""),
  actionButton("queryData", "Run Query", width = '100%', class = "btn btn-primary", icon = icon("check")),
  tags$hr(),
  
  fluidRow(
    column(width = 12,  
      tags$div(class = "scrollable-table", DTOutput("table", fill = TRUE))
    )
  )
)

analysis_page <- tabPanel(
  title = div(icon("chart-bar", style = "margin-right: 5px;"), "Analysis"),
  
  p(""),
  p("Upload a file below containing chromosome (chr) and starting (start) and ending (end) positions to run an analysis using a pre-trained Random Forest model and get probabilities output."),
  p(""),
  
  tags$div(
    class = "rounded-grey-square",
    selectInput("genome", tags$h6(strong("Choose a genome:")), choices = c("Select a genome" = "", unique(combined_dataset$genome)), width = '100%'),
  ),
  
  p(""),
  
  tags$div(
    class = "rounded-grey-square",
    selectInput("paper", tags$h6(strong("Choose a paper:")), choices = c("Select a paper" = "", unique(combined_dataset$source_dataset)), width = '100%'),
    selectInput("dc", tags$h6(strong("Choose a disease/cell type:")), choices = c("Select a disease/cell type:" = "", c("Please select a paper first"), width = '100%')),
    selectInput("model", tags$h6(strong("Choose a model:")), choices = c("Select a model type:" = "", c("motif" = "motif", "3mer" = "3mer")), width = '100%'),
  ),
  
  p(""),
  
  tags$div(class = "rounded-grey-square",
    p(""),
    tags$h6(strong("Choose your file type (BED/FASTA):")),
    selectInput("analysisFileType", label = NULL, choices = c("Select a genome" = "", c("BED" = "BED", "FASTA" = "FASTA")), width = '100%'),
    tags$hr(),
    fileInput('file2', tags$h6(strong('Upload a .txt file')), accept = c(".txt")),
  ),
  
  p(""),
  actionButton("loadAnalysisFile", "Generate Predictions", width = '100%', icon = icon("check"), class = "btn btn-primary"),
  tags$hr(),
  p(""),
  
  DTOutput("probabilitiesOutput"),
)

about_page <- tabPanel(
  title = div(icon("info-circle", style = "margin-right: 5px;"), "About"),
  
  tags$h5(strong("Background")),
  tags$hr(),
  p(""),
  p(paste0(project_title, " is a massively-parallel reporter assay project that allows you to dynamically filter and explore an integrated genomic dataset. You can explore the dataset and filter by chromosome number, chromosome positions, disease, cell type, or even input your own file for filtering.")),
  p("There is also an analysis tool for predicting how certain DNA sequences or genetic variants affect gene expression based on our aggregated MRPA data. We trained several machine-learning models using a Random Forest approach for providing prediction scores."),
  tags$img(src = "test.png", height = "520px", width = "900px"),
  
  tags$h5(strong("Instructions")),
  tags$hr(),
  p("Some text here."),
  p("Text here.")
)


ui <- fluidPage(
  #theme = bs_theme(),
  
  useShinyjs(),

  tags$head(
    tags$style(
      HTML("
          .custom-header {
            display: inline-block;
            border: 1px solid rgb(225, 225, 225);
            border-radius: 5px;
            padding: 12px;
            align-items: center;
            margin-bottom: 15px;
          }
          
          .custom-header h1 {
            line-height: 0.8;
            margin: 0;
          }
          
          .custom-header h2 {
            fill: solid rgb(255, 0, 0)
          }
          
          .header-spacing {
            margin-bottom: 52.5px;
          }
          
          .scrollable-table {
            max-height: 700px;
            overflow-y: auto;
            overflow-x: auto;
          }
    
          .text-style {
            font-size: 15px;
            margin-top: 5px;
            margin-bottom: 5px;
          }
          
          .btn-dark {
            background-color: #373A3C !important;
            border-color: #373A3C !important;
            color: #yourTextColor !important;
          }
          
          .rounded-grey-square {
            background-color: #FAFAFA;
            border-radius: 5px;
            padding: 10px;
            border: 1.25px solid #DCDCDC;  /* border */
          }
        
          .rounded-grey-square2 {
            background-color: #F5F5F5;
            border-radius: 5px;
            padding: 10px;
            border: 1.25px solid #DCDCDC;  /* border */
          }
        
          .center-div {
            width: 450px; /* or whatever width you want */
            text-align: center; /* optional, for text within the div */
          }
        
          table {
            width: 100%;
            border-collapse: separate;
            border-spacing: 15px 0;
          }
        
          .preview-back {
            background-color: #FFFFFF;
            border-radius: 5px;
            padding: 10px;
            border: 1px solid #DCDCDC;  /* border */
            width: 100%;
            height: 100%;
          }
        
          .navbar-brand {
            font-weight: 500; /* Make the text bold */
            font-size: larger; /* Increase the font size */
            /* Add other styles as needed */
          }
    
          .navbar.navbar-default {
            background-color: #373A3C !important;
          }
      ),
    ")),
  ),
  
  fluidRow(
    navbarPage(
      theme = custom, #theme = bs_theme(bootswatch = "united"),
      inverse = TRUE,
      id = "myNavbar", 
      title = project_title,
      about_page,
      database_page,
      analysis_page
    ),
  )
)


