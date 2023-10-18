# Server function ----
shinyServer(function(session, input, output) {
  
  # Inputs ----
  file_selected     <- reactive(input$fn)
  analysis_name     <- reactiveVal(NULL)
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
      tmp <- which(IDA_result$Eval$Sample %in% input$sample)
      if (length(tmp) == 1) {
        return(tmp)
      } else {
        return(1)
      }
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
  load_archive <- reactiveVal(0)
  is_archived <- reactiveVal(FALSE)
  is_zoomed <- reactiveVal(FALSE)
  disable("saveCSVSummary")
  disable("saveCSVProcVals")
  disable("saveMSXL")
  
  # Testing ----
  observeEvent(input$inspect, { browser() })
  
  # Main function ----
  # Generates IDA object on file or archive load
  observeEvent(input$fn, {
    is_archived(FALSE)
    analysis_name(tools::file_path_sans_ext(input$fn$name))
    addClass("processing_mask", "overlay")
    removeClass("processing_mask", "hidden")
    tmp <- isolate({
      IDA(
        raw_data = input$fn$datapath,
        buffer = input$buffer_all,
        tolerance = input$tolerance_all,
        expansion = input$expansion_all,
        draw_stable_bounds = FALSE,
        is_file = TRUE
      )
    })
    tmp$Quality <- list(
      tmp$Quality,
      IDA_quality(tmp$Eval, repel = TRUE)
    )
    lapply(names(tmp),
           function(x) {
             IDA_result[[x]] <- tmp[[x]]
           })
    sample_list <- as.vector(IDA_result$Eval$Sample)
    enable("saveCSVSummary")
    enable("saveCSVProcVals")
    enable("saveMSXL")
    redraw(redraw() + 1)
    updateTabsetPanel(session,
                      'main',
                      selected = "Results")
    updateTabsetPanel(session,
                      'results',
                      selected = "Details")
    updateSelectInput(session,
                      "sample",
                      choices = sample_list,
                      selected = head(sample_list, 1))
    addClass("processing_mask", "hidden")
    removeClass("processing_mask", "overlay")
  })
  
  observeEvent({
    input$buffer_all
    input$tolerance_all
    input$expansion_all
  }, {
    req(slide_trigger(), any(!is.null(input$fn), is_archived()))
    addClass("processing_mask", "overlay")
    removeClass("processing_mask", "hidden")
    use_direct <- !is.null(IDA_result$Raw) || is_archived()
    if (is.null(IDA_result$Raw)) {
      raw_data <- input$fn$datapath
    } else {
      raw_data <- IDA_result$Raw
    }
    tmp <- isolate({
      IDA(
        raw_data = raw_data,
        buffer = input$buffer_all,
        tolerance = input$tolerance_all,
        expansion = input$expansion_all,
        draw_stable_bounds = FALSE,
        is_file = !use_direct
      )
    })
    tmp$Quality <- list(
      tmp$Quality,
      IDA_quality(tmp$Eval, repel = TRUE)
    )
    lapply(names(tmp),
           function(x) {
             IDA_result[[x]] <- tmp[[x]]
           })
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
    enable("saveCSVSummary")
    enable("saveCSVProcVals")
    enable("saveMSXL")
    redraw(redraw() + 1)
    addClass("processing_mask", "hidden")
    removeClass("processing_mask", "overlay")
  })
  
  # Single sample parameter update ----
  observeEvent({
    input$buffer_single
    input$tolerance_single
    input$expansion_single
  }, {
    req(IDA_result$Processed, slide_trigger())
    addClass("processing_mask", "overlay")
    removeClass("processing_mask", "hidden")
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
    IDA_result$Eval[x, 2:4] <- IDA_summary(IDA_temp,
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
    IDA_result$Quality <- list(
      IDA_quality(IDA_result$Eval),
      IDA_quality(IDA_result$Eval, repel = TRUE)
    )
    redraw(redraw() + 1)
    addClass("processing_mask", "hidden")
    removeClass("processing_mask", "overlay")
  })
  
  # Update outputs ----
  # Observe changes to redraw() and recreate outputs.  This can now be fired from anywhere by changing redraw().
  observeEvent(redraw(), {
    if (!is.null(IDA_result$Processed)) {
      output$IDAsummary <- DT::renderDataTable({
        if (is.null(IDA_result$Eval)) {
          NULL
        } else {
          # full_join(IDA_result$Eval, IDA_result$Info)
          DT::datatable(
            data = full_join(IDA_result$Eval, IDA_result$Info) |>
              mutate(across(c(Mean, StDev, RSD), ~ round(.x, digits = 4))),
            rownames = FALSE,
            extensions = c("Responsive")
          )
        }
      })
      if(input$sample %in% IDA_result$Eval$Sample){
        i <- which(IDA_result$Eval$Sample == input$sample)
      } else {
        i <- 1
      }
      ii <- dim(IDA_result$Eval)[1] - i + 1
      out <- IDA_result$Quality[[1 + isolate(input$repel_labels)]] +
        annotate(
          "point",
          x = IDA_result$Eval$RSD[i],
          y = ii,
          size = 4,
          pch = 5,
          stroke = 2,
          colour = 'forestgreen'
        )
      if (is_zoomed()) {
        xlims <- isolate(
          c(input$brush_quality$xmin, input$brush_quality$xmax)
        )
        ylims <- isolate(
          c(input$brush_quality$ymin, input$brush_quality$ymax)
        )
        out <- out +
          coord_cartesian(xlim = xlims, ylim = ylims)
      }
      output$IDAquality <- renderPlot(out)
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
      if (input$sample == "Please load a file first." || length(current_sample()) == 0) {
        i <- 1
      } else {
        i <- current_sample()
      }
      stable_start <- attr(IDA_result$Processed[[i]], "stable.start")
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
  observeEvent(input$repel_labels, {
    redraw(redraw() + 1)
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
      redraw(redraw() + 1)
    }
  })
  
  # Qual plot click ----
  # Observe quality plot click event and change selected sample to nearest point.
  observeEvent(input$click_quality, {
    req(input$click_quality)
    selected_sample <- nearPoints(
      df = IDA_result$Quality[[1 + isolate(input$repel_labels)]]$data,
      coordinfo = input$click_quality,
      maxpoints = 1,
      threshold = 10
    )
    if (nrow(selected_sample) == 1) {
      updateSelectInput(session, 'sample', selected = selected_sample$Sample)
    }
  })
  
  # Qual plot brush ----
  # Observe the quality plot zoom button events and zoom to the brushed area
  observeEvent(input$zoom_out, {
    is_zoomed(FALSE)
    i <- which(IDA_result$Eval$Sample == input$sample)
    ii <- dim(IDA_result$Eval)[1] - i + 1
    output$IDAquality <- renderPlot({
      IDA_result$Quality[[1 + isolate(input$repel_labels)]]+
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
  })
  observeEvent(input$zoom_in, {
    is_zoomed(TRUE)
    if (is.null(input$brush_quality)) {
      alert("Select a region in the quality plot to zoom in and click this button again.")
    }
    req(input$brush_quality)
    i <- which(IDA_result$Eval$Sample == input$sample)
    ii <- dim(IDA_result$Eval)[1] - i + 1
    xlims <- isolate(
      c(input$brush_quality$xmin, input$brush_quality$xmax)
    )
    ylims <- isolate(
      c(input$brush_quality$ymin, input$brush_quality$ymax)
    )
    output$IDAquality <- renderPlot({
      IDA_result$Quality[[1 + isolate(input$repel_labels)]] +
        coord_cartesian(xlim = xlims,
                        ylim = ylims) +
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
      updateActionButton(session, 'saveNewTime', label = btn_text, icon=icon("save"))
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
      IDA_result$Eval[i, 2:4] <- IDA_summary(IDA_result$Processed[[i]], x_min, x_max)
      attr(IDA_result$Processed[[i]], "stable.start") <- x_min
      IDA_result$Info$`Stable Time Start`[i] <- x_min
      attr(IDA_result$Processed[[i]], "stable.end") <- x_max
      IDA_result$Info$`Stable Time End`[i] <- x_max
      IDA_result$Quality <- list(
        IDA_quality(IDA_result$Eval),
        IDA_quality(IDA_result$Eval, repel = TRUE)
      )
      redrawit <- redraw() + 1
      redraw(redrawit)
      updateActionButton(session, 'saveNewTime', label = "Draw on the ratio plot to select a new stability region.")
    }
  })
  
  # Archived analysis modal ----
  observeEvent(input$pastSamples, {
    if (change_archive_list()) {
      if (!is.null(input$pastSamples)) {
        if (!is.null(IDA_result$Raw)) {
          showModal(
            modalDialog(
              title = "This will remove the current analysis and load one from the archive.",
              actionButton('loadArchive', "Proceed", icon = icon("check")),
              footer = modalButton("Cancel", icon = icon("times")),
              fade = FALSE
            )
          )
        } else {
          load_archive(load_archive() + 1)
          reset("fn")
          enable("saveCSVSummary")
          enable("saveCSVProcVals")
          enable("saveMSXL")
        }
      }
    }
  })
  
  # Confirm archive load ----
  observeEvent(input$loadArchive, {
    load_archive(load_archive() + 1)
  })
  
  # Load from archive ----
  observeEvent(load_archive(), {
    req(load_archive() > 0)
    addClass("processing_mask", "overlay")
    removeClass("processing_mask", "hidden")
    fname <- paste0(input$pastSamples[1], ".RDS")
    analysis_name(tools::file_path_sans_ext(fname))
    archive <- readRDS(file.path(getwd(), "archive", fname))
    validate(
      need(all(names(IDA_result %in% names(archive))),
           message = "This archive does not contain all necessary information.")
    )
    if (length(archive$Quality) > 2) {
      archive$Quality <- list(
        archive$Quality,
        IDA_quality(archive$Eval, repel = TRUE)
      )
    }
    lapply(names(IDA_result),
           function(x) {
             IDA_result[[x]] <- archive[[x]]
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
    updateActionButton(session,
                       'archive',
                       label = "Archived, click to update",
                       icon = icon("toggle-on"))
    reset("fn")
    removeModal()
    redrawit <- redraw() + 1
    redraw(redrawit)
    is_archived(TRUE)
    addClass("processing_mask", "hidden")
    removeClass("processing_mask", "overlay")
  })
  
  # Archive the current analysis ----
  observeEvent(input$archive, {
    if (is.null(IDA_result$Eval)) {
      alert("Please load a file for analysis.")
    } else {
      req(analysis_name())
      if (is_archived()) {
        fname <- analysis_name()
      } else {
        fname <- paste0(Sys.Date(), " - ", analysis_name())
      }
      destination <- file.path(getwd(), "archive",
                               paste0(fname,
                                      ".RDS"))
      updateActionButton(session,
                         'archive',
                         label = "Archived, click to update",
                         icon = icon("toggle-on"))
      out <- reactiveValuesToList(IDA_result)
      saveRDS(out, file = destination)
      if (is_archived()) {
        alert(paste0("This analysis has been updated as ", fname, "."))
      } else {
        alert(paste0("The current analysis has been archived as ", fname, "."))
      }
      change_archive_list(!change_archive_list())
      pastAnalyses <- gsub(".RDS", "", list.files(file.path(getwd(), "archive"), pattern = ".RDS"))
      selected_archive <- isolate(input$pastSamples)
      if (!is_archived()) {
        updateSelectizeInput(session,
                             'pastSamples',
                             choices = pastAnalyses,
                             selected = selected_archive, 
                             options = list(placeholder = "(Select a prior analysis)",
                                            maxItems = 1))
      }
      change_archive_list(!change_archive_list())
    }
  })
  
  # ______________________________----
  # Downloads ----
  ## Example file ----
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
  
  ## Current summary ----
  output$saveCSVSummary <- downloadHandler(
    filename = function() {
      out <- paste0(
        tools::file_path_sans_ext(analysis_name()),
        " - processed ",
        Sys.Date(),
        " Summary.csv")
      return(out)
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
  
  ## All processing values ----
  output$saveCSVProcVals <- downloadHandler(
    filename = function() {
      paste0(
        analysis_name(),
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
  
  ## All as MSXL ----
  output$saveMSXL <- downloadHandler(
    filename = function() {
      paste0(analysis_name(), " - processed ", Sys.Date(), ".xlsx")
    },
    content = function(file) {
      addClass("download_mask", "overlay")
      removeClass("download_mask", "hidden")
      IDA_result %>%
        reactiveValuesToList() %>%
        pack_as_excel(draw_stable_bounds=TRUE) %>%
        saveWorkbook("tmp", overwrite = TRUE)
      file.copy("tmp", file)
      file.remove("tmp")
      addClass("download_mask", "hidden")
      removeClass("download_mask", "overlay")
    }
  )
}) 
