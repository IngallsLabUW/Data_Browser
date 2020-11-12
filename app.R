
# Setup things ----

library(shiny)
library(shinyDirectoryInput)
library(RaMS)

ui <- fluidPage(
  titlePanel("TBD"),
  sidebarLayout(
    sidebarPanel(
      numericInput(inputId = "mz", label = "Choose a mass of interest:", value = 118.0865),
      directoryInput('directory', label = "Choose a directory", 
                     value = r"(G:\My Drive\FalkorFactor\mzMLs\pos\MSMS)"),
      h4("Files to load:"),
      tableOutput("files_to_load"),
      actionButton(inputId = "loadem", label = "Load files"),
      h4("Files found:"),
      tableOutput("found_files")
    ),
    mainPanel(
      
    )
  ),
  includeScript("detect_click.js")
)

server <- function(input, output, session){
  files_to_load <- reactiveVal()
  
  output$found_files <- renderTable({
    dir <- readDirectoryInput(session, 'directory')
    mzml_files <- list.files(dir, pattern = ".mzML")
    mzml_files <- setdiff(mzml_files, files_to_load())
    
    if(!length(mzml_files)){
      return(data.frame())
    }
    linked_names <- paste0("<a href='#' onclick='detect_click(this)'>",
                         mzml_files, "</a>")
    `colnames<-`(data.frame(linked_names), NULL)
  }, sanitize.text.function = function(x) x)
  
  output$files_to_load <- renderTable({
    `colnames<-`(data.frame(files_to_load()), NULL)
  })
  
  
  #Observe clicked text
  observeEvent(input$clicked_text, {
    files_to_load(c(files_to_load(), input$clicked_text))
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
