


# Setup things ----

library(shiny)
library(shinyDirectoryInput)
library(shinythemes)
library(plotly)
library(data.table)
library(xml2)
library(base64enc)
source("RaMS_custom.R")
R
#Debugging things
# input <- list(mz=118.0865, ppm=5, directory=r"(G:\My Drive\FalkorFactor\mzMLs\pos\MSMS)")
# files_to_load <- function()c("G:\\My Drive\\FalkorFactor\\mzMLs\\pos\\MSMS/180205_Poo_TruePooPos_dda1.mzML",
#                              "G:\\My Drive\\FalkorFactor\\mzMLs\\pos\\MSMS/190715_Poo_TruePooFK180310_DDApos50.mzML")
# files_to_load <- function()c("Z:\1_QEdata\Will\MSMS_mzMLs\200929_Poo_TruePooMortality_DDApos50/mzML")
# mzml_data <- grabMzmlData(files_to_load()[[1]])
# clean_filenames <- gsub("\\.mzML$", "", basename(files_to_load()[[1]]))
# sutom_MS1_data <- cbind(mzml_data[[1]], filename=clean_filenames)
# sutom_MS2_data <- cbind(mzml_data[[2]], filename=clean_filenames)
# mzml_data <- grabMzmlData(files_to_load()[[2]])
# clean_filenames <- gsub("\\.mzML$", "", basename(files_to_load())[[2]])
# sutom_MS1_data <- rbind(sutom_MS1_data, cbind(mzml_data[[1]], filename=clean_filenames))
# sutom_MS2_data <- rbind(sutom_MS2_data, cbind(mzml_data[[2]], filename=clean_filenames))
# current_MS1_data <- function()sutom_MS1_data
# current_MS2_data <- function()sutom_MS2_data


# UI ----

ui <- fluidPage(
  sidebarLayout(
    sidebarPanel(
      numericInput(inputId = "mz", label = "Enter a mass of interest:", value = 118.0865),
      numericInput(inputId = "ppm", label = "Enter instrument ppm:", value = 5),
      directoryInput('directory', label = "Choose a directory", 
                     value = r"(G:\My Drive\FalkorFactor\mzMLs\pos\MSMS)"),
      br(),
      strong("Files to load:"),
      tableOutput("files_to_load"),
      actionButton(inputId = "loadem", label = "Load selected files"),
      br(),
      br(),
      strong("Files found:"),
      tableOutput("found_files"),
      style="margin-top: 20px;"
    ),
    mainPanel(
      # textOutput("files_loaded"),
      plotlyOutput("MS1_chrom", height="300px"),
      plotlyOutput("MS2_chrom", height="300px"),
      style="margin-top: 20px;"
    )
  ),
  includeScript("detect_click.js"),
  theme = "sandstone.mod.css"
)






server <- function(input, output, session){
  files_to_load <- reactiveVal()
  files_loaded <- reactiveVal()
  current_MS1_data <- reactiveVal()
  current_MS2_data <- reactiveVal()
  
  output$found_files <- renderTable({
    dir <- readDirectoryInput(session, 'directory')
    mzml_paths <- list.files(dir, pattern = ".mzML", full.names = TRUE)
    if(!length(files_to_load())){
      mzml_files <- basename(mzml_paths)
    } else {
      mzml_files <- setdiff(basename(mzml_paths), basename(files_to_load()))
    }
    if(length(files_loaded())){
      mzml_files <- setdiff(mzml_files, basename(files_loaded()))
    }
    
    if(!length(mzml_files)){
      return(data.frame())
    }
    linked_names <- paste0("<a href='#' onclick='detect_click(this)'>",
                         mzml_files, "</a>")
    `colnames<-`(data.frame(linked_names), NULL)
  }, sanitize.text.function = function(x) x)
  
  output$files_to_load <- renderTable({
    if(!length(files_to_load())){
      `colnames<-`(data.frame(), NULL)
    } else {
      `colnames<-`(data.frame(basename(files_to_load())), NULL)
    }
  })
  
  output$files_loaded <- renderText({
    req(files_loaded())
    paste("Files loaded:", paste(basename(files_loaded()), collapse = ", "))
  })
  
  output$MS1_chrom <- renderPlotly({
    if(req(current_MS1_data())){
      EIC <- current_MS1_data()[mz%between%pmppm(input$mz, input$ppm)]
      
      MS2_scans <- current_MS2_data()[premz%between%pmppm(input$mz, input$ppm)]
      MS2_scans <- unique(MS2_scans[,c("rt", "voltage", "filename")])
      
      
      nearest_MS1 <- lapply(split(MS2_scans, seq_len(nrow(MS2_scans))), function(MS2_scan){
        EIC_file <- EIC[filename==MS2_scan$filename]
        ms1_scan_data <- EIC_file[which.min(abs(EIC_file$rt-MS2_scan$rt)), c("rt", "int")]
        names(ms1_scan_data) <- c("rt_to_plot", "int_to_plot")
        cbind(ms1_scan_data, MS2_scan[,c("rt", "voltage", "filename")])
      })
      nearest_MS1 <- do.call(what = rbind, nearest_MS1)
      nearest_MS1$label <- paste0(
        'Retention time: ', round(nearest_MS1$rt_to_plot, digits = 3), ' min<br>',
        'Intensity: ', nearest_MS1$int_to_plot, '<br>',
        'Voltage: ', nearest_MS1$voltage, '<br>',
        'Filename: ', nearest_MS1$filename
      )
      
      plot_ly(source="MS1") %>%
        add_trace(data = EIC, type="scatter", 
                  x=~rt, y=~int, color=~filename,
                  mode="lines", hoverinfo='skip') %>%
        add_trace(data = nearest_MS1, type="scatter", 
                  x=~rt_to_plot, y=~int_to_plot, 
                  color=~filename, text=~label,
                  customdata=~rt,
                  hovertemplate=paste0('%{text}<extra></extra>'),
                  mode="markers", showlegend = FALSE) %>%
        layout(legend = list(orientation="h", x=0.5, y=100), showlegend=TRUE)
    } else {
      plot_ly(x=1, y=1, mode="text", type="scatter", text="No data found :/")
    }
  })
  
  output$MS2_chrom <- renderPlotly({
    selected_scan <- event_data(event = "plotly_click", source = "MS1")$x
    req(selected_scan)
    
    mass_MS2_data <- current_MS2_data()[premz%between%pmppm(input$mz, input$ppm)]
    MS2_scan_rts <- unique(mass_MS2_data$rt)
    nearest_MS2_rt <- MS2_scan_rts[which.min(abs(MS2_scan_rts-selected_scan))]
    MS2_data <- mass_MS2_data[rt==nearest_MS2_rt]
    MS2_data$int <- MS2_data$int/max(MS2_data$int)*100
    split_MS2 <- split(MS2_data, seq_len(nrow(MS2_data)))
    split_MS2 <- lapply(split_MS2, function(row){
      rbind(data.frame(rt=row$rt, premz=row$premz, fragmz=row$fragmz, 
                       int=0, voltage=row$voltage, filename=row$filename), 
            row, 
            data.frame(rt=row$rt, premz=row$premz, fragmz=row$fragmz, 
                       int=0, voltage=row$voltage, filename=row$filename))
    })
    comb_MS2 <- do.call(what = rbind, split_MS2)
    
    plot_ly(source="MS2") %>%
      add_trace(data=comb_MS2, x=~fragmz, y=~int,
                type = "scatter", mode="lines",
                line=list(color="black"), hoverinfo="skip") %>%
      add_trace(data=MS2_data, x=~fragmz, y=~int,
                type="scatter", mode="markers",
                marker=list(color="black")) %>%
      layout(xaxis = list(range=c(0, input$mz*1.1)),
             showlegend=FALSE)
  })
  
  #Observe clicked text
  observeEvent(input$clicked_text, {
    files_to_load(c(files_to_load(), paste(readDirectoryInput(session, 'directory'), 
                                           input$clicked_text, sep = "/")))
  })
  
  #Observe file load button
  observeEvent(input$loadem, {
    req(files_to_load())
    new_MS1_data <- list()
    new_MS2_data <- list()
    n_files <- length(files_to_load())
    withProgress({
      for(i in seq_along(files_to_load())){
        mzml_data <- grabMzmlData(files_to_load()[[i]])
        new_MS1_data[[i]] <- mzml_data[[1]]
        new_MS2_data[[i]] <- mzml_data[[2]]
      }
    }, message = "Loading your data...", value = 0)

    clean_filenames <- gsub("\\.mzML$", "", basename(files_to_load()))
    
    new_MS1_data <- mapply(cbind, new_MS1_data, clean_filenames, SIMPLIFY = FALSE)
    new_MS1_data <- do.call(what = rbind, new_MS1_data)
    new_MS1_data <- `colnames<-`(new_MS1_data, c("rt", "mz", "int", "filename"))
    
    new_MS2_data <- mapply(cbind, new_MS2_data, clean_filenames, SIMPLIFY = FALSE)
    new_MS2_data <- do.call(what = rbind, new_MS2_data)
    new_MS2_data <- `colnames<-`(new_MS2_data, c("rt", "premz", "fragmz", "int", 
                                                 "voltage", "filename"))
    
    current_MS1_data(rbind(current_MS1_data(), new_MS1_data))
    current_MS2_data(rbind(current_MS2_data(), new_MS2_data))
    
    files_loaded(c(files_loaded(), files_to_load()))
    files_to_load(NULL)
  })
  
  #Observe triple button click
  observeEvent(
    ignoreNULL = TRUE,
    eventExpr = {
      input$directory
    },
    handlerExpr = {
      if (input$directory > 0) {
        # launch the directory selection dialog with initial path read from the widget
        path = choose.dir(default = readDirectoryInput(session, 'directory'))
        
        # update the widget value
        updateDirectoryInput(session, 'directory', value = path)
      }
    }
  )
}

shinyApp(ui = ui, server = server)
