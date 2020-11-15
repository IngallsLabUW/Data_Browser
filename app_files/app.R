


# Setup things ----

library(shiny)
library(shinyDirectoryInput)
library(plotly)
library(data.table)
library(xml2)
library(base64enc)
source("app_files/RaMS_custom.R")

# Debugging things ----
# input <- list(mz=118.0865, ppm=5, directory=r"(G:\My Drive\FalkorFactor\mzMLs\pos\MSMS)")
# files_to_load <- function()c("G:\\My Drive\\FalkorFactor\\mzMLs\\pos\\MSMS/180205_Poo_TruePooPos_dda1.mzML",
#                              "G:\\My Drive\\FalkorFactor\\mzMLs\\pos\\MSMS/190715_Poo_TruePooFK180310_DDApos50.mzML")
# files_to_load <- function()c("Z:\1_QEdata\Will\MSMS_mzMLs\200929_Poo_TruePooMortality_DDApos50/mzML")
# files_to_load <- function()"G:\\My Drive\\FalkorFactor\\mzMLs\\pos\\190715_Poo_TruePooFK180310_Half1.mzML"
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

# files_to_load <- function()r"(Z:\1_QEdata\LTC\DATA\HILIC\190718_DepthProfiles_FK180310\MSMS\190715_Poo_TruePooFK180310_DDAneg50.mzXML)"
# mzxml_data <- grabMzxmlData(files_to_load())
# current_MS1_data <- function()mzxml_data[[1]]
# current_MS2_data <- function()mzxml_data[[2]]

# UI ----

ui <- fluidPage(
  # Include custom CSS so it looks pretty
  tags$head(includeCSS("app_files/sandstone.mod.css")),
  
  # Default sidebar layout
  sidebarLayout(
    sidebarPanel(
      numericInput(inputId = "mz", label = "Enter a mass of interest:", value = 118.0865),
      numericInput(inputId = "ppm", label = "Enter instrument ppm:", value = 5),
      directoryInput('directory', label = "Choose a directory", 
                     value = "~"),
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
      plotlyOutput("MS1_chrom", height="300px"),
      plotlyOutput("MS2_chrom", height="300px"),
      style="margin-top: 20px;"
    )
  ),
  # Include the Javascript necessary to detect text clicks
  includeScript("app_files/detect_click.js")
)




# Server ----
# General scheme: 
#  Load in files
#  Render MS1 chromatogram
#  When a point is clicked on the MS1 chrom, render MS2 chromatogram
# Watch for new files to be loaded, new masses added, etc.

server <- function(input, output, session){
  # Define our reactive values
  # files_to_load and files_loaded produce character vectors containing
  # paths to the mzML files
  # files_to_load contains paths to the "staged" files
  # files_loaded tracks all the files that have already been loaded and
  #   removes them from the list presented to the user
  files_to_load <- reactiveVal()
  files_loaded <- reactiveVal()
  
  # current_MS1_data and current_MS2_data produce data.table objects
  # current_MS1_data has four columns: rt, mz, int, and filename
  # current_MS2_data has six columns: rt, premz, fragmz, int, filename, and voltage
  current_MS1_data <- reactiveVal()
  current_MS2_data <- reactiveVal()
  
  # Render the table of mzML files that are found in the given directory
  output$found_files <- renderTable({
    # This function is from shinyDirectoryInput and returns a file path
    dir <- readDirectoryInput(session, 'directory')
    
    # Find all the .mzML files and collect their full paths
    # I've chosen to retain the full paths in the files_to_load and files_loaded
    # variables to avoid sloppy paste()s but it does mean I use a lot of 
    # basename() to get just the important part.
    mzml_paths <- list.files(dir, pattern = ".mzML|.mzXML", full.names = TRUE)
    if(!length(files_to_load())){
      mzml_files <- basename(mzml_paths)
    } else {
      # If files have already been staged, remove them from the full list
      mzml_files <- setdiff(basename(mzml_paths), basename(files_to_load()))
    }
    if(length(files_loaded())){
      # If files have already been loaded, remove them from the full list
      mzml_files <- setdiff(mzml_files, basename(files_loaded()))
    }
    
    if(!length(mzml_files)){
      #Handle cases when no mzML files are found in the directory
      return(data.frame())
    }
    
    # Combine file names with HTML hyperlink text that activates the
    # detect_click javascript
    linked_names <- paste0("<a href='#' onclick='detect_click(this)'>",
                         mzml_files, "</a>")
    
    # Remove column names so the output looks pretty
    `colnames<-`(data.frame(linked_names), NULL)
    
    # Fake sanitize the text so the raw HTML <a> tags don't get escaped
  }, sanitize.text.function = function(x) x)
  
  # Render the "staged" files area
  output$files_to_load <- renderTable({
    if(!length(files_to_load())){
      f2load <- data.frame()
    } else {
      f2load <- data.frame(basename(files_to_load()))
    }
    colnames(f2load) <- NULL
    f2load
  })
  
  # Use plotly to render the MS1 chromatogram (top chart)
  output$MS1_chrom <- renderPlotly({
    # Check whether an MS1 data exists - if not, don't bother plotting
    req(current_MS1_data())
    
    # Grab Extracted Ion Chromatogram - i.e. points within a given
    # mass window. pmppm() takes a mass and a ppm accuracy and returns
    # an upper and lower mass bound, and data.table's between()
    # does a binary search to locate the relevant data points.
    EIC <- current_MS1_data()[mz%between%pmppm(input$mz, input$ppm)]
    if(!nrow(EIC)){
      return(plot_ly(x=1, y=1, text="No data found :/", type = "scatter", mode="text"))
    }
    
    # Create the plotly, with x=retention time, y=intensity, and color by file
    p <- plot_ly(source="MS1") %>%
      add_trace(data = EIC, type="scatter", 
                x=~rt, y=~int, color=~filename,
                mode="lines", hoverinfo='skip') %>%
      layout(legend = list(orientation="h", x=0.5, y=100), showlegend=TRUE)
    # Add the legend at top middle, in horizontal layout
    
    # If the loaded files contain MSMS data, add points on the chromatogram
    # that can be clicked on
    if(!is.null(current_MS2_data())){
      # Collect MSMS data that corresponds to the parent mass (premz)
      MS2_scans <- current_MS2_data()[premz%between%pmppm(input$mz, input$ppm)]
      
      # Fun jujitsu code
      # Since the MS2 scans follow the MS1 scan that got plotted, we need to
      # find the nearest MS1 scan and add the MS2 point to the MS1 coordinates
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
        'Voltage: ', nearest_MS1$voltage, 'V<br>',
        'Filename: ', nearest_MS1$filename
      )
      
      # Add MS2 data points to the plot, if they exist
      p <- p %>%
        add_trace(data = nearest_MS1, type="scatter", 
                  x=~rt_to_plot, y=~int_to_plot, 
                  color=~filename, text=~label,
                  customdata=~rt,
                  hovertemplate=paste0('%{text}<extra></extra>'),
                  mode="markers", showlegend = FALSE)
    }
    # Return the final plotly object
    p
  })
  
  # Annoying code to avoid warning about plotly not finding the
  # MS1 plot before it's even been rendered
  # Basically just checks whether we have the MS1 plot and returns information
  # about which MS2 point was clicked on
  relayout_data <- reactive({
    req(input$loadem)
    event_data(event = "plotly_click", source = "MS1")
  })
  
  # Use plotly to render the MS2 chromatogram (bottom chart)
  # after a point has been clicked on in the top chart
  output$MS2_chrom <- renderPlotly({
    # Grab the point that's been chosen - this could probably be handled
    # better with the customdata argument but I haven't figured out how yet
    selected_scan <- relayout_data()$x
    req(selected_scan)
    
    # Grab MS2 data that's close to the precursor mass (premz)
    mass_MS2_data <- current_MS2_data()[premz%between%pmppm(input$mz, input$ppm)]
    
    # Fun jiujitsu code to again find that weird "closest" MS2 scan to the 
    # chosen MS1
    MS2_scan_rts <- unique(mass_MS2_data$rt)
    nearest_MS2_rt <- MS2_scan_rts[which.min(abs(MS2_scan_rts-selected_scan))]
    MS2_data <- mass_MS2_data[rt==nearest_MS2_rt]
    # Normalize out of 100, per convention
    MS2_data$int <- MS2_data$int/max(MS2_data$int)*100
    
    # Gross code to make the nice histogram plot
    # Basically sets an intensity of zero immediately before and after each
    # data point to make the spike vertical
    split_MS2 <- split(MS2_data, seq_len(nrow(MS2_data)))
    split_MS2 <- lapply(split_MS2, function(row){
      rbind(data.frame(rt=row$rt, premz=row$premz, fragmz=row$fragmz, 
                       int=0, voltage=row$voltage, filename=row$filename), 
            row, 
            data.frame(rt=row$rt, premz=row$premz, fragmz=row$fragmz, 
                       int=0, voltage=row$voltage, filename=row$filename))
    })
    comb_MS2 <- do.call(what = rbind, split_MS2)
    
    # Create the plotly object and return it.
    # Use lines to create "lollipop" sticks and markers to create 
    # mouseover-sensitive "lollipops"
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
  
  # Observe clicked text
  # When text is clicked, move it into the "files to load" object and paste
  # the current directory name to it
  observeEvent(input$clicked_text, {
    files_to_load(c(files_to_load(), paste(readDirectoryInput(session, 'directory'), 
                                           input$clicked_text, sep = "/")))
  })
  
  # Observe file load button
  # When clicked, load the files which are currently in files_to_load
  # Use a Shiny-specific withProgress bar
  # grabMzmlData is the magical function here - check out the RaMS_custom.R
  # script for details or ask Will
  observeEvent(input$loadem, {
    req(files_to_load())
    new_MS1_data <- list()
    new_MS2_data <- list()
    n_files <- length(files_to_load())
    withProgress({
      for(i in seq_along(files_to_load())){
        setProgress(value = (i-1)/n_files, detail = basename(files_to_load()[[i]]))
        if(grepl("mzML", files_to_load()[[i]])){
          mzml_data <- grabMzmlData(files_to_load()[[i]])
        } else if(grepl("mzXML", files_to_load()[[i]])) {
          mzml_data <- grabMzxmlData(files_to_load()[[i]])
        } else {
          print(files_to_load()[[i]])
          print("How on earth did you sneak a non-mz(X)ML file in here?")
          req(FALSE)
        }
        new_MS1_data[[i]] <- mzml_data[[1]]
        new_MS2_data[[i]] <- mzml_data[[2]]
      }
      setProgress(value = n_files, detail = "Done!")
    }, message = "Loading your data...", value = 0)

    clean_filenames <- gsub("\\.mzML|\\.mzXML$", "", basename(files_to_load()))
    
    # Add filenames to new data and clean up a little
    new_MS1_data <- mapply(cbind, new_MS1_data, clean_filenames, SIMPLIFY = FALSE)
    new_MS1_data <- do.call(what = rbind, new_MS1_data)
    new_MS1_data <- `colnames<-`(new_MS1_data, c("rt", "mz", "int", "filename"))
    if(!is.null(new_MS2_data)){
      new_MS2_data <- mapply(cbind, new_MS2_data, clean_filenames, SIMPLIFY = FALSE)
      new_MS2_data <- do.call(what = rbind, new_MS2_data)
      new_MS2_data <- `colnames<-`(new_MS2_data, c("rt", "premz", "fragmz", "int", 
                                                   "voltage", "filename"))
    }
    
    # Combine any existing MS1 and MS2 data with the new data
    current_MS1_data(rbind(current_MS1_data(), new_MS1_data))
    current_MS2_data(rbind(current_MS2_data(), new_MS2_data))
    
    # Add files_to_load to files_loaded and clear out files_to_load
    files_loaded(c(files_loaded(), files_to_load()))
    files_to_load(NULL)
  })
  
  #Observe triple button click
  # Magic code, no idea how this works - directly from shinyDirectoryInput vignette
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


# Run app ----

# Launch the default web browser
# URL used is localhost (127.0.0.1, port 1313)
browseURL("http://127.0.0.1:1313/")

# Set the port and launch the app
shinyApp(ui = ui, server = server, options = list(port=1313))
