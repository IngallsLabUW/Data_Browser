


# Setup things ----

reqd_packages <- c("shiny", "plotly", "devtools", "data.table", "RaMS")
missing_packages <- !reqd_packages%in%installed.packages()[,"Package"]
if(any(missing_packages)){
  message("Whoops! Looks like not all required packages are installed.")
  message(paste("Packages missing:", paste(
    reqd_packages[missing_packages], collapse = ", ")
  ))
  message("I can try to install them for you, or you can do it yourself.")
  message("1: Install them for you\n2: Do it yourself")
  choice <- "neither"
  while(!choice%in%c("1", "2")){
    message("Choose 1 or 2 by entering that number here: ")
    if(interactive()){
      choice <- readline()
    } else {
      choice <- readLines(con = "stdin", n = 1)
    }
  }
  if(choice=="1"){
    message("Okay, I'll try to install them.")
    message(paste("Installing packages:", paste(
      reqd_packages[missing_packages], collapse = ", ")
    ))
    cran_to_install <- setdiff(reqd_packages[missing_packages], c("RaMS"))
    if(length(cran_to_install)>0)install.packages(cran_to_install, 
                                                  repos='https://cran.rstudio.com/')
    if(missing_packages[5]){
      devtools::install_github('wkumler/RaMS')
    }
  } else if (choice=="2"){
    message("Okay, I'll let you do it yourself. Press Enter to exit.")
    if(interactive()){
      readline()
    } else {
      readLines(con = "stdin", n = 1)
    }
  } else {
    message("How on earth did you do that?")
    readLines(con = "stdin", n = 1)
  }
}

library(shiny)
library(plotly)
library(data.table)
library(RaMS)


metadata <- fread("metadata.csv")
names(metadata) <- tolower(gsub(" ", "_", names(metadata)))
metadata <- metadata[project_name!=""&file_name!=""&file_path!=""]


stan_data <- fread("Ingalls_Lab_Standards_NEW.csv")
stan_data <- stan_data[order(Compound.Name_figure),
                       c("Compound.Name_figure", "m.z", "z", "Fraction1")]
stan_data <- stan_data[!is.na(m.z)]
stan_data$polarity <- ifelse(stan_data$z==1, "pos", "neg")
stan_data$split <- paste0(stan_data$Fraction1, " (", stan_data$polarity, ")")

stan_choices <- `names<-`(stan_data$m.z, stan_data$Compound.Name_figure)
stan_choices <- split(stan_choices, stan_data$split)
stan_choices <- stan_choices[c("HILICPos (pos)", "HILICNeg (neg)", "CyanoDCM (pos)",
                               "CyanoAq (pos)", "CyanoAq (neg)")]

files <- data.frame(file_name=gsub(".mzML.*", "", list.files("MSMS files")))
# files <- files[!files$file_name%in%c(files_to_load(), files_loaded()),,drop=FALSE]
mzml_files <- merge(files, metadata, by = "file_name")[
  c("file_name", setdiff(names(metadata), c("file_path", "file_name")))
]


# Debugging things ----
# input <- list(mz=118.0865, ppm=5, directory=r"(G:\My Drive\FalkorFactor\mzMLs\pos\MSMS)")
# files_to_load <- function()c("G:\\My Drive\\FalkorFactor\\mzMLs\\pos\\MSMS/180205_Poo_TruePooPos_dda1.mzML",
#                              "G:\\My Drive\\FalkorFactor\\mzMLs\\pos\\MSMS/190715_Poo_TruePooFK180310_DDApos50.mzML")

# UI ----

ui <- fluidPage(
  # Include custom CSS so it looks pretty
  tags$head(includeCSS("sandstone.mod.css")),
  
  fluidRow(
    sidebarLayout(
      sidebarPanel(
        strong("Choose files to load from table at bottom:"),
        tableOutput("files_to_load"),
        actionButton(inputId = "loadem", label = "Load selected files"),
        br(),
        br(),
        selectInput(inputId = "stan_mz", label = "Pick a standard:",
                    choices = stan_choices, multiple = FALSE, selectize = TRUE),
        numericInput(inputId = "mz", label = "or enter a mass directly:", value = 118.0865),
        numericInput(inputId = "ppm", label = "Enter instrument ppm:", value = 5),
        checkboxGroupInput(inputId = "polarity", label = "Particular polarity?",
                           choices = c(Positive="pos", Negative="neg"), 
                           selected = c(Positive="pos", Negative="neg"),
                           inline = TRUE),
        checkboxGroupInput(inputId = "which_cruises", label = "Which cruises?",
                           choices = sort(unique(metadata$project_name)),
                           selected = sort(unique(metadata$project_name))),
        actionButton(inputId = "clearchecks", label = "Uncheck all"),
        style="margin: 10px;",
        width = 3
      ),
      mainPanel(
        plotlyOutput("MS1_chrom", height="300px"),
        plotlyOutput("MS2_chrom", height="300px"),
        style="margin-top: 20px;"
      )
    )
  ),
  fluidRow(
    DT::dataTableOutput("found_files"),
    style="margin: 20px;"
  )
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
  
  current_mz <- reactiveVal()
  
  observeEvent(input$stan_mz, {
    current_mz(as.numeric(input$stan_mz))
    updateNumericInput(session=session, inputId = "mz", 
                       value = as.numeric(input$stan_mz))
  })
  observeEvent(input$mz, {
    current_mz(input$mz)
  })
  observeEvent(input$clearchecks, {
    print("Clearing checks")
    updateCheckboxGroupInput(inputId = "which_cruises", selected = character(0))
  })
  
  output$found_files <- DT::renderDataTable({
    mzml_files <- mzml_files[mzml_files$project_name%in%input$which_cruises,]
    mzml_files <- mzml_files[!mzml_files$file_name%in%files_loaded(),]
    mzml_files <- mzml_files[tolower(mzml_files$polarity)%in%input$polarity,]
    mzml_files
  }, escape = FALSE, options = list(
    scrollX = TRUE, 
    pageLength = 5,
    lengthMenu=list(c(5, 10, -1), c('5', '10', 'All')))
  )
  
  # Render the "staged" files area
  output$files_to_load <- renderTable({
    if(!length(files_to_load())){
      f2load <- data.frame("No files selected")
    } else {
      f2load <- data.frame(basename(files_to_load()))
    }
    colnames(f2load) <- NULL
    f2load
  })
  
  # Observe clicked text
  # When text is clicked, move it into the "files to load" object and paste
  # the current directory name to it
  observeEvent(input$found_files_rows_selected, {
    mzml_files <- mzml_files[mzml_files$project_name%in%input$which_cruises,]
    mzml_files <- mzml_files[!mzml_files$file_name%in%files_loaded(),]
    
    files_to_load(mzml_files$file_name[input$found_files_rows_selected])
  })
  
  # Observe file load button
  # When clicked, load the files which are currently in files_to_load
  # Use a Shiny-specific withProgress bar
  observeEvent(input$loadem, {
    req(files_to_load())
    filepaths <- paste0("MSMS files/", files_to_load(), ".mzML.gz")
    new_MS1_data <- list()
    new_MS2_data <- list()
    n_files <- length(filepaths)
    withProgress({
      for(i in seq_along(files_to_load())){
        setProgress(value = (i-1)/n_files, detail = basename(files_to_load()[[i]]))
        mzml_data <- grabMSdata(filepaths[[i]], grab_what = c("MS1", "MS2"))
        current_MS1_data(rbind(current_MS1_data(), mzml_data$MS1))
        current_MS2_data(rbind(current_MS2_data(), mzml_data$MS2))
      }
      setProgress(value = n_files, detail = "Done!")
    }, message = "Loading your data...", value = 0)

    clean_filenames <- gsub("\\.mzML|\\.mzXML$", "", basename(files_to_load()))
    
    # Add files_to_load to files_loaded and clear out files_to_load
    files_loaded(c(files_loaded(), files_to_load()))
    files_to_load(NULL)
  })
  
  # Use plotly to render the MS1 chromatogram (top chart)
  output$MS1_chrom <- renderPlotly({
    # Check whether an MS1 data exists - if not, don't bother plotting
    req(current_MS1_data())
    req(current_mz())
    
    # Grab Extracted Ion Chromatogram - i.e. points within a given
    # mass window. pmppm() takes a mass and a ppm accuracy and returns
    # an upper and lower mass bound, and data.table's between()
    # does a binary search to locate the relevant data points.
    EIC <- current_MS1_data()[mz%between%pmppm(current_mz(), input$ppm)]
    if(!nrow(EIC)){
      return(plot_ly(x=1, y=1, text="No data found :/", type = "scatter", 
                     mode="text") %>% layout(yaxis=list(title=""), xaxis=list(title="")))
    }

    # Create the plotly, with x=retention time, y=intensity, and color by file
    p <- plot_ly(source="MS1") %>%
      add_trace(data = EIC, type="scatter", 
                x=~rt, y=~int, color=~filename,
                mode="lines", hoverinfo='skip') %>%
      layout(legend = list(orientation="h", x=0.5, y=100), showlegend=TRUE,
             xaxis=list(title = "Retention time"), yaxis=list(title = "Intensity")
             )
    # Add the legend at top middle, in horizontal layout
    
    # If the loaded files contain MSMS data, add points on the chromatogram
    # that can be clicked on
    if(!is.null(current_MS2_data())){
      # Collect MSMS data that corresponds to the parent mass (premz)
      MS2_scans <- current_MS2_data()[premz%between%pmppm(current_mz(), input$ppm)]
      if(nrow(MS2_scans)>0){
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
    mass_MS2_data <- current_MS2_data()[premz%between%pmppm(current_mz(), input$ppm)]
    if(!nrow(mass_MS2_data)){
      return(plot_ly(x=1, y=1, text="No data found :/", type = "scatter", 
                     mode="text") %>% layout(yaxis=list(title=""), xaxis=list(title="")))
    }
    
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
      layout(xaxis = list(range=c(0, current_mz()*1.1)),
             showlegend=FALSE)
  })
  
}


# Run app ----

# Launch the default web browser
# URL used is localhost (127.0.0.1, port 1313)
browseURL("http://127.0.0.1:1313/")

# Set the port and launch the app
shinyApp(ui = ui, server = server, options = list(port=1313))
