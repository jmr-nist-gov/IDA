# Server function ----
shinyServer(function(session, input, output) {
  
  # Inputs ----
  file_selected     <- reactive(input$fn)
  apply_to          <- reactive(input$apply)
  buffer_all        <- reactive(input$buffer_all)
  tolerance_all     <- reactive(input$tolerance_all)
  expansion_all     <- reactive(input$expansion_all)
  buffer_single     <- reactive(input$buffer_single)
  tolerance_single  <- reactive(input$tolerance_single)
  expansion_single  <- reactive(input$expansion_single)
  current_sample    <- reactive({
    if (is.null(IDA_result)) {
      return(NULL)
    } else {
      which(IDA_result$Eval$Sample %in% input$sample)
    }
  })
  redraw <- reactiveVal(1)
  IDA_result <- reactiveValues(
    Eval      = NULL,
    Info      = NULL,
    Processed = NULL,
    Raw       = NULL,
    Graphs    = NULL,
    Quality   = NULL,
    isotopes  = NULL
  )
  slide_trigger <- reactiveVal(TRUE)
  change_archive_list <- reactiveVal(TRUE)
  disable("saveCSVSummary")
  disable("saveCSVProcVals")
  disable("saveMSXL")
  
  # Main function ----
  # Generates IDA object on file load
  observeEvent({
    input$fn
    input$buffer_all
    input$tolerance_all
    input$expansion_all
  }, {
    if (!is.null(input$fn)) {
      tmp <- isolate({
        IDA(
          input$fn$datapath,
          input$buffer_all,
          input$tolerance_all,
          input$expansion_all,
          draw_stable_bounds = FALSE
        )
      })
      lapply(names(tmp),
             function(x) {
               IDA_result[[x]] <- tmp[[x]]
             })
      sample_list <- as.vector(IDA_result$Eval$Sample)
      updateSelectInput(session,
                        "sample",
                        choices = sample_list,
                        selected = head(sample_list, 1))
      updateTabsetPanel(session,
                        'main',
                        selected = "Results")
      updateTabsetPanel(session,
                        'results',
                        selected = "Details")
      slide_trigger(!slide_trigger())
      updateSliderInput(session,
                        'buffer_single',
                        value = input$buffer_all)
      updateSliderInput(session,
                        'tolerance_single',
                        value = input$tolerance_all)
      updateSliderInput(session,
                        'expansion_single',
                        value = input$expansion_all)
      slide_trigger(!slide_trigger())
      redraw(redraw() + 1)
      enable("saveCSVSummary")
      enable("saveCSVProcVals")
      enable("saveMSXL")
    }
  })
  
  # Single sample parameter update ----
  observeEvent({
    input$buffer_single
    input$tolerance_single
    input$expansion_single
  }, {
    req(IDA_result$Processed)
    if (slide_trigger() & !is.null(IDA_result$Processed)) {
      x <- which(IDA_result$Eval$Sample == input$sample)
      # Replace processed, including
      IDA_temp <- IDA_calc(
        IDA_result$Raw[[x]],
        input$buffer_single,
        input$tolerance_single,
        input$expansion_single
      )
      # Replace processed attributes
      attributes(IDA_result$Processed[[x]]) <- attributes(IDA_temp)
      # Replace info start and end times
      IDA_result$Info[x, 6] <- attr(IDA_temp, "stable.start")
      IDA_result$Info[x, 7] <- attr(IDA_temp, "stable.end")
      IDA_result$Info[x, 8] <- attr(IDA_temp, "proc.buffer")
      IDA_result$Info[x, 9] <- attr(IDA_temp, "proc.tolerance")
      IDA_result$Info[x, 10] <- attr(IDA_temp, "proc.expansion")
      # Replace eval values
      IDA_result$Eval[x, 4:6] <- IDA_summary(IDA_temp,
                                              attr(IDA_temp, "stable.start"),
                                              attr(IDA_temp, "stable.end"))
      # Replace graphs
      IDA_result$Graphs[[x]] <- IDA_graphs(
        IDA_result$Processed[[x]],
        IDA_temp,
        names(IDA_result$Processed)[x],
        IDA_result$isotopes,
        draw_stable_bounds = FALSE
      )
      IDA_result$Quality <- IDA_quality(IDA_result$Eval)
      redraw(redraw() + 1)
      temp_sample <- isolate(input$sample)
      updateSelectInput(session,
                        "sample",
                        selected = temp_sample)
    }
  })
  
  # Update outputs ----
  # Observe changes to redraw() and recreate outputs.  This can now be fired from anywhere by changing redraw().
  observeEvent(redraw(), {
    if (!is.null(IDA_result$Processed)) {
      output$IDAsummary <- DT::renderDataTable({
        if (is.null(IDA_result$Eval)) {
          NULL
        } else {
          full_join(IDA_result$Eval, IDA_result$Info)
        }
      })
      if(input$sample %in% IDA_result$Eval$Sample){
        i <- which(IDA_result$Eval$Sample == input$sample)
      } else {
        i <- 1
      }
      ii <- dim(IDA_result$Eval)[1] - i + 1
      output$IDAquality <- renderPlot({
        IDA_result$Quality +
          annotate(
            "point",
            x = IDA_result$Eval$RSD[i],
            y = ii,
            size = 4,
            pch = 5,
            stroke = 2,
            colour = 'forestgreen'
          )
        
      })
      output$IDAvalues <- DT::renderDataTable({
        if (is.null(IDA_result$Eval)) {
          NULL
        } else {
          raw <- as.data.frame(IDA_result$Raw[[current_sample()]])[, -1]
          names(raw) <- gsub("X", "", names(raw))
          proc <- as.data.frame(IDA_result$Processed[[current_sample()]])
          names(proc) <- gsub("X", "", names(proc))
          names(proc)[2:3] <- paste(names(proc)[2:3], " corrected")
          x <- data.frame(cbind(proc$Time,
                                proc$Ratio,
                                proc[, 2:3],
                                raw))
          nam_list <- c("Time",
                        "Ratio",
                        gsub("X", "", names(proc)[2:3]),
                        gsub("X", "", names(raw)))
          names(x) <- nam_list
          x
        }
      })
      if (input$sample == "Please load a file first.") {
        i <- 1
      } else {
        i <- current_sample()
      }
      stable_start <-attr(IDA_result$Processed[[i]], "stable.start")
      stable_end <- attr(IDA_result$Processed[[i]], "stable.end")
      stability_caption <- paste0(
        "\nSelected stability region is marked by red lines (",
        stable_start,
        "-",
        stable_end,
        " s)."
      )
      output$rawPlot <- renderPlot({
        IDA_result$Graphs[[i]]$Signal +
          geom_vline(xintercept = stable_start, colour = 'red') +
          geom_vline(xintercept = stable_end, colour = 'red')
      })
      output$ratioPlot <- renderPlot({
        IDA_result$Graphs[[i]]$Ratio +
          geom_vline(xintercept = stable_start, colour = 'red') +
          geom_vline(xintercept = stable_end, colour = 'red') +
          labs(caption = stability_caption)
      })
    }
  })
  
  # File required ----
  # Make certain user doesn't select "Single Sample" prior to a file being loaded.
  observeEvent(input$apply, {
    if (input$apply == "Single Sample" & is.null(IDA_result$Processed)) {
      shinyjs::info("Please load a sample set first.")
      updateTabsetPanel(session, 'apply', selected = "All")
    }
  })
  
  # Sample selection ----
  # Observe sample selection change to highlight the given sample in the quality plot
  observeEvent(input$sample, {
    if (!is.null(IDA_result$Processed)) {
      i <- which(IDA_result$Eval$Sample == input$sample)
      ii <- dim(IDA_result$Eval)[1] - i + 1
      output$IDAquality <- renderPlot({
        IDA_result$Quality +
          annotate(
            "point",
            x = IDA_result$Eval$RSD[i],
            y = ii,
            size = 4,
            pch = 5,
            stroke = 2,
            colour = 'forestgreen'
          )
      })
      # updateTabsetPanel(session, 'apply', selected = "Single Sample")
      slide_trigger(!slide_trigger())
      updateSliderInput(session,
                        'buffer_single',
                        value = attr(IDA_result$Processed[[current_sample()]], "proc.buffer"))
      updateSliderInput(
        session,
        'tolerance_single',
        value = attr(IDA_result$Processed[[current_sample()]], "proc.tolerance")
      )
      updateSliderInput(
        session,
        'expansion_single',
        value = attr(IDA_result$Processed[[current_sample()]], "proc.expansion")
      )
      slide_trigger(!slide_trigger())
      stable_start <-
        attr(IDA_result$Processed[[i]], "stable.start")
      stable_end <- attr(IDA_result$Processed[[i]], "stable.end")
      stability_caption <- paste0(
        "\nSelected stability region is marked by red lines (",
        stable_start,
        "-",
        stable_end,
        " s)."
      )
      output$rawPlot <- renderPlot({
        IDA_result$Graphs[[i]]$Signal +
          geom_vline(xintercept = stable_start, colour = 'red') +
          geom_vline(xintercept = stable_end, colour = 'red')
      })
      output$ratioPlot <- renderPlot({
        IDA_result$Graphs[[i]]$Ratio +
          geom_vline(xintercept = stable_start, colour = 'red') +
          geom_vline(xintercept = stable_end, colour = 'red') +
          labs(caption = stability_caption)
      })
    }
  })
  
  # Qual plot click ----
  # Observe quality plot click event and change selected sample to nearest point.
  observe({
    if (!is.null(input$click_quality)) {
      i <- round(input$click_quality$y)
      i <- dim(IDA_result$Eval)[1] - i + 1
      selected_sample <- as.character(IDA_result$Eval$Sample[i])
      updateSelectInput(session, 'sample', selected = selected_sample)
    }
  })
  
  # Ratio plot brush ----
  # Observe ratio plot brushing event and provide the user real-time feedback on performance.
  observe({
    if (!is.null(input$brush_ratio)) {
      x_min <- round(input$brush_ratio$xmin, 3)
      x_max <- round(input$brush_ratio$xmax, 3)
      output$rawPlot <- renderPlot(
        IDA_result$Graphs[[current_sample()]]$Signal +
          annotate(
            "rect",
            xmin = x_min,
            xmax = x_max,
            ymin = -Inf,
            ymax = Inf,
            fill = 'lightskyblue',
            alpha = 0.25
          )
      )
      btn_text <-
        paste0("Save ", x_min, "-", x_max, " s as the stability region.")
      updateActionButton(session, 'saveNewTime', label = btn_text, icon=icon("floppy-o"))
    }
  })
  
  # Save the new stability time region ----
  observeEvent(input$saveNewTime, {
    if (is.null(input$brush_ratio)) {
      info("Please select a region of the ratio plot.")
    } else {
      x_min <- round(input$brush_ratio$xmin, 3)
      x_max <- round(input$brush_ratio$xmax, 3)
      i <- which(IDA_result$Eval$Sample %in% input$sample)
      IDA_result$Eval[i, 4:6] <- IDA_summary(IDA_result$Processed[[i]], x_min, x_max)
      attr(IDA_result$Processed[[i]], "stable.start") <- x_min
      IDA_result$Info$`Stable Time Start`[i] <- x_min
      attr(IDA_result$Processed[[i]], "stable.end") <- x_max
      IDA_result$Info$`Stable Time End`[i] <- x_max
      IDA_result$Quality <- IDA_quality(IDA_result$Eval)
      redrawit <- redraw() + 1
      redraw(redrawit)
      updateActionButton(session, 'saveNewTime', label = "Draw on the ratio plot to select a new stability region.")
    }
  })
  
  # Load an archived analysis ----
  observeEvent(input$pastSamples, {
    if (change_archive_list()) {
      if (!is.null(input$pastSamples)) {
        if (!is.null(IDA_result)) {
          showModal(
            modalDialog(
              title = "This will remove the current analysis and load one from the archive.",
              actionButton('loadArchive', "Proceed", icon = icon("check")),
              footer = modalButton("Cancel", icon = icon("close")),
              fade = FALSE
            )
          )
        } else {
          load(file.path(getwd(), "archive", input$pastSamples[1], ".Rdata"))
          IDA_result <- IDA_result
          sample_list <- as.vector(IDA_result$Eval$Sample)
          updateSelectInput(
            session,
            "sample",
            choices = sample_list,
            selected = head(sample_list, 1)
          )
          updateTabsetPanel(session,
                            'main',
                            selected = "Results")
          updateTabsetPanel(session,
                            'results',
                            selected = "Details")
          updateActionButton(session,
                             'archive',
                             label = "Archived, click to update",
                             icon = icon("toggle-on"))
          redrawit <- redraw() + 1
          redraw(redrawit)
          reset("fn")
          enable("saveCSVSummary")
          enable("saveCSVProcVals")
          enable("saveMSXL")
        }
      }
    }
  })
  
  # Load from archive ----
  # Load an archived analysis from the modal dialog (same code snippet as above but triggers from modal dialogue)
  observeEvent(input$loadArchive, {
    load(paste0(getwd(), "/archive/", input$pastSamples[1], ".Rdata"))
    IDA_result <- IDA_result
    sample_list <- as.vector(IDA_result$Eval$Sample)
    updateSelectInput(session,
                      "sample",
                      choices = sample_list,
                      selected = head(sample_list, 1))
    updateTabsetPanel(session,
                      'main',
                      selected = "Results")
    updateTabsetPanel(session,
                      'results',
                      selected = "Details")
    updateActionButton(session,
                       'archive',
                       label = "Archived, click to update",
                       icon = icon("toggle-on"))
    reset("fn")
    removeModal()
    redrawit <- redraw() + 1
    redraw(redrawit)
  })
  
  # Archive the current analysis ----
  observeEvent(input$archive, {
    if (is.null(IDA_result)) {
      alert("Please load a file for analysis.")
    } else {
      fname <-
        paste0(Sys.Date(),
               " - ",
               gsub(".csv", "", file_selected()$name),
               ".Rdata")
      destination <- file.path(getwd(), "archive", fname)
      updateActionButton(session,
                         'archive',
                         label = "Archived, click to update",
                         icon = icon("toggle-on"))
      out <- reactiveValuesToList(IDA_result)
      save(out, file = destination)
      alert(paste0("The current analysis has been archived as ", fname, "."))
      change_archive_list(!change_archive_list())
      pastAnalyses <- gsub(".Rdata", "", list.files(file.path(getwd(), "archive")))
      selected_archive <- isolate(input$pastSamples)
      updateSelectizeInput(session,
                           'pastSamples',
                           choices = pastAnalyses,
                           selected = selected_archive, 
                           options = list(placeholder = "(Select a prior analysis)",
                                          maxItems = 1))
      change_archive_list(!change_archive_list())
    }
  })
  
  # Download Example ----
  # Provide an example file from the app directory to provider users with the expected input format and to have a file to explore the tool.
  output$example <- downloadHandler(
    filename = "IDA Example - SRM 2778.csv",
    content = function(file) {
      # Example file
      file_example <-
        read_csv(paste0(getwd(), "/SRM 2778.csv"), col_types = "ccc")
      names(file_example) <- gsub("X", "", names(file_example))
      write.csv(
        file_example,
        quote = FALSE,
        row.names = FALSE,
        na = "",
        file
      )
      shinyjs::info("Drag this file to the Browse... input to get started.")
    }
  )
  
  ## Download summary ----
  output$saveCSVSummary <- downloadHandler(
    filename = function() {
      paste0(
        gsub(".csv", "", file_selected()$name),
        " - processed ",
        Sys.Date(),
        " Summary.csv")
    },
    content = function(file) {
      write.csv(
        IDA_result$Eval,
        quote = TRUE,
        row.names = FALSE,
        na = "",
        file
      )
    }
  )
  
  # Download processing values ----
  output$saveCSVProcVals <- downloadHandler(
    filename = function() {
      paste0(
        gsub(".csv", "", file_selected()$name),
        " - processed ",
        Sys.Date(),
        " Processing Values.csv")
    },
    content = function(file) {
      write.csv(
        IDA_result$Info,
        quote = TRUE,
        row.names = FALSE,
        na = "",
        file
      )
    }
  )
  
  # Download as XL ----
  output$saveMSXL <- downloadHandler(
    filename = function() {
      if (is.null(file_selected())){
        paste0(input$pastSamples, " - processed ", Sys.Date(), ".xlsx")
      } else {
        paste0(file_selected()$name, " - processed ", Sys.Date(), ".xlsx")
      }
    },
    content = function(file) {
      addClass("mask", "overlay")
      removeClass("mask", "hidden")
      if (is.null(file_selected())){
        fname <- paste0(input$pastSamples, " - processed ", Sys.Date(), ".xlsx")
      } else {
        fname <- paste0(file_selected()$name, " - processed ", Sys.Date(), ".xlsx")
      }
      saveWorkbook(pack_as_excel(IDA_result, draw_stable_bounds=TRUE), fname, overwrite = TRUE)
      file.copy(fname, file)
      file.remove(fname)
      addClass("mask", "hidden")
      removeClass("mask", "overlay")
    }
  )
}) 
