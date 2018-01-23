x <- c("shiny", "DT", "shinyjs")
lapply(x, require, character.only=TRUE)

instructions <- "instructions.html"
pastAnalyses <- gsub(".Rdata", "", list.files(paste0(getwd(), "/archive/")))
if (length(pastAnalyses) == 0) pastAnalyses <- NULL
if (is.null(pastAnalyses)){
  pastSamples_placeholder <- "(No archived analyses yet.)"
} else {
  pastSamples_placeholder <- "(Select a prior analysis)"
}

shinyUI(
  fluidPage(
    useShinyjs(),
    tags$head(
      tags$style(HTML("
                      .nav-tabs>li {width: 50%; text-align: center;}
                      .nav-pills>li {position: relative; display: inline; width: 49%; text-align: center;}
                      .shiny-download-link {text-align: center; width: 100%;}
                      .btn {width: 100%;}
                      body {text-align: justify;}
                      .modal-body>btn {width: 49%;}
                      "))
    ),
    
    # Application title
    titlePanel("Isotope Dilution Assistant v0.9"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
      sidebarPanel(
        h4("1. Apply processing values to:"),
        tabsetPanel(id = "apply", type="pills",
                    tabPanel("All",
                             br(),
                             sliderInput("buffer_all",
                                         "Buffer size: (default 10 points)",
                                         min = 0,
                                         max = 50,
                                         value = 10),
                             sliderInput("tolerance_all",
                                         "Tolerance: (default 10% change)",
                                         min = 0,
                                         max = 0.5,
                                         value = 0.1),
                             sliderInput("expansion_all",
                                         "Expansion size: (default 20 points)",
                                         min = 10,
                                         max = 50,
                                         value = 20)
                    ), # Close tabPanel "All"
                    tabPanel("Single Sample",
                             br(),
                             sliderInput("buffer_single",
                                         "Buffer size: (default 10 points)",
                                         min = 0,
                                         max = 50,
                                         value = 10),
                             sliderInput("tolerance_single",
                                         "Tolerance: (default 10% change)",
                                         min = 0,
                                         max = 0.5,
                                         value = 0.1),
                             sliderInput("expansion_single",
                                         "Expansion size: (default 20 points)",
                                         min = 10,
                                         max = 50,
                                         value = 20)
                    ) # Close tabPanel "Single Sample"
        ), # Close tabesetPanel "apply"
        hr(),
        h4("2. Select an instrument file:"),
        fileInput("fn",
                  NULL,
                  multiple = FALSE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv"),
                  buttonLabel = "Browse...",
                  placeholder = "No file selected"),
        downloadButton("example", "Get an example data set", width="100%"),
        hr(),
        h4("3. Download results:"),
        downloadButton("saveCSVSummary", "Summary Only (.csv)", width="100%"),
        downloadButton("saveCSVProcVals", "Processing Values (.csv)", width="100%"),
        #downloadButton('saveMSXL', "All Results (.xlsx)", width="100%")#,
        #downloadButton('saveGraphs', "Graphs Only (.jpg)", icon('download'), width="100%")
        hr(),
        h4("4. Archive:"),
        selectizeInput("pastSamples", 
                       label = NULL, 
                       choices = pastAnalyses, 
                       multiple = TRUE, 
                       selected = NULL, 
                       options = list(placeholder = pastSamples_placeholder,
                                      maxItems = 1)),
        actionButton("archive", 
                     "Click to archive this analysis", 
                     icon=icon("toggle-off"), 
                     width="100%")
      ), # Close sidebarPanel
      
      # Show a plot of the generated distribution
      mainPanel(
        tabsetPanel(id = "main",
                    tabPanel("Instructions",
                             includeHTML(instructions),
                             img(src="A product of NIST DTD.jpg", align="right")
                             ), # Close tabPanel "Instructions"
                    tabPanel("Results",
                             tabsetPanel(id = "results", type="pills",
                                         tabPanel("Summary Results",
                                                  br(),
                                                  DT::dataTableOutput("IDAsummary")),
                                         tabPanel("Details",
                                                  h4("Quality Overview"),
                                                  p("Selected sample is highlighted with a green diamond."),
                                                  plotOutput("IDAquality",
                                                             click = "click_quality"),
                                                  hr(),
                                                  h4("Sample Investigation"),
                                                  p("Select a sample or click a point in the quality overview:"),
                                                  selectInput("sample",
                                                              NULL,
                                                              choices = c("Please load a file first."),
                                                              width = "100%",
                                                              selectize = FALSE),
                                                  tabsetPanel(id = "details", type="pills",
                                                              tabPanel("Graphs",
                                                                       br(),
                                                                       plotOutput("rawPlot", 
                                                                                  height="200px",
                                                                                  brush=brushOpts("brush_raw", 
                                                                                                  direction = "x",
                                                                                                  resetOnNew = TRUE)),
                                                                       hr(),
                                                                       plotOutput("ratioPlot", 
                                                                                  height="200px",
                                                                                  brush=brushOpts("brush_ratio", 
                                                                                                  direction = "x",
                                                                                                  resetOnNew = TRUE)),
                                                                       #textOutput("brushText"),
                                                                       actionButton("saveNewTime", 
                                                                                    "Draw on the ratio plot to select a new stability region.",
                                                                                    icon=icon("floppy-o"))),
                                                              tabPanel("Values",
                                                                       br(),
                                                                       DT::dataTableOutput("IDAvalues"))
                                                  ) # Close tabsetPanel 'details'
                                         ) # Close tabPanel "Quality Overview"
                             ) # Close tabsetPanel "results"
                    ) # Close tabPanel "Results"
        ) # Close tabsetPanel "main"
      ) # Close mainPanel
    ) # Close sidebarLayout
  ) # Close fluidPage
) # Close shinyUI
