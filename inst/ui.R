# packs <- c("shiny", "DT", "shinyjs", "ggrepel")
# packs_false <- packs[-which(packs %in% installed.packages())]
# if (length(packs_false) > 0) {
#   install.packages(pkgs = packs_false, dependencies = TRUE)
# }
# lapply(packs, library, character.only=TRUE)

instructions <- "instructions.html"
disclaimer   <- "disclaimer.html"
pastAnalyses <- gsub(".RDS", "",
                     list.files(
                       path = file.path(getwd(), "archive"),
                       pattern = ".RDS$"
                     )
)
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
                      #quality_controls > .btn {width: 100px; border: 0; padding: 0; color: #337ab7; float: right;}
                      #quality_controls > .form-group {width: 150px; display: inline-block;}
                      #quality_controls > .form-group > .checkbox {margin: 0;}
                      .hidden {position: absolute; z-index: 1;}
                      .overlay {position: fixed; z-index: 3; opacity: 0.85; top: 0; left: 0; width: 100%; height: 100%; background-color: White; color: Black;}
                      .overlay>img {position: fixed; top: calc(50vh - 150px); left: calc(50vw - 100px); width: 200px; height: 200px; opacity: 1;}
                      .overlay>h3 {position: fixed; top: calc(50vh + 20px); left: calc(50vw - 100px); width: 200px; color: Black; opacity: 1; text-align: center}
                      #dtd_logo {float: left; text-align: left;}
                      #dtd_logo > p {margin-bottom: 0; font-size: xx-small;}
                      #dtd_logo > img {height: 36px; margin-bottom: 5px;}
                      #dtd_logo p:first-child {color: #337ab7; font-size: small;}
                      #dtd_contact_disclaimer {float: right; text-align: right; max-width: 200px;}
                      #dtd_contact_disclaimer > * {font-size: small;}
                      "))
    ),
    div(id="download_mask", class="hidden", img(src="processing.gif"), h3("Building .xlsx file...")),
    div(id="processing_mask", class="hidden", img(src="processing.gif")),
    
    # Application title
    titlePanel("Isotope Dilution Assistant"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
      sidebarPanel(
        tagList(
          if(TESTING) {
            actionButton("inspect", "Inspect App")
          } else {
            NULL
          }),
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
        downloadButton("saveMSXL", "All Results (.xlsx)", width="100%"),
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
                             div(id = "dtd_logo",
                                 p("A data tool product from"),
                                 img(src="nist_logo.png"),
                                 # p("U.S. Department of Commerce"),
                                 # p("National Institute of Standards and Technology"),
                                 p("Material Measurement Laboratory"),
                                 p("Chemical Sciences Division"),
                                 p("Chemical Informatics Group")
                             ),
                             div(id = "dtd_contact_disclaimer",
                                 p(tags$strong("IDA v1.0.1, last updated on 2023-10-18")),
                                 p(id = "contact",
                                   "Please direct any questions to ",
                                   a(href = "mailto::jared.ragland@nist.gov?subject=[IDA] Question",
                                     "jared.ragland@nist.gov")),
                                 p(id = "disclaimer",
                                   "This software published under the ",
                                   a(href = "https://www.nist.gov/disclaimer",
                                     "NIST Software Disclaimer"))
                             )
                             # img(src="A product of NIST DTD.png", align="right")
                    ), # Close tabPanel "Instructions"
                    tabPanel("Results",
                             tabsetPanel(id = "results", type="pills",
                                         tabPanel("Summary Results",
                                                  br(),
                                                  p("Display of Mean, StDev, and RSD values has been truncated to four decimal places. Full precision values will be available in the download."),
                                                  DT::dataTableOutput("IDAsummary")),
                                         tabPanel("Details",
                                                  h4("Quality Overview"),
                                                  p("Selected sample is highlighted with a green diamond."),
                                                  div(id = "quality_controls",
                                                      checkboxInput("repel_labels", "Space out labels", value = FALSE),
                                                      actionButton("zoom_in", "Zoom In", icon = icon("search-plus")),
                                                      actionButton("zoom_out", "Zoom Out", icon = icon("search-minus"))
                                                  ),
                                                  plotOutput("IDAquality",
                                                             click = "click_quality",
                                                             brush = "brush_quality"),
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
                                                                                    icon=icon("pen-square"))),
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
