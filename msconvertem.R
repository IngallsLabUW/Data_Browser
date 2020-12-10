library(tidyverse)
library(pbapply)

outfile <- "G:/Shared drives/Ingalls Lab/Collaborative_Projects/Database_work/Data_Browser/MSMS files"
already_converted <- list.files(path = outfile) %>% str_replace(".mzML.gz", "")

file_data <- read.csv("metadata.csv") %>%
  filter(File.path!=""&File.name!="") %>%
  select(File.path, File.name) %>%
  filter(!File.name%in%already_converted)
file_names <- paste0(file_data$File.path, "/", file_data$File.name, ".mzXML")
file_exists <- file.exists(file_names)


success_vec <- pbsapply(file_names[file_exists], function(filename){
  msconvert_command <- paste0(
    'msconvert --mzML -gzv "', filename, '" -o "', outfile, '"'
  )
  system(msconvert_command)
})

if(any(success_vec)){
  message(paste("Unable to convert", names(success_vec[success_vec])))
}
