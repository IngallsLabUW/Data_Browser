# Data_Browser
A small R app designed to visualize mzML data from past Ingalls Lab experiments.

## Philosophy

It shouldn't be hard to open up a mass-spec file and browse around the data inside. But it is. 
Often the libraries required are difficult to install, the functions are slow and memory-inefficient,
and the process produces arcane R objects from which it's difficult to extract *m/z*, retention time,
or intensity data. 

This R Shiny-based app is meant to resolve some of these problems by providing a simple user-interface to our 
mass-spec data, allowing the user to browse intuitively through files on their computer or connected drives, load those
into R automagically, and search in those files by mass. The app is optimized for MSMS data, with
scans that collected fragmentation information at a given mass denoted by points on the chromatogram,
and clicking on those will pull up the fragments and their relative intensities in the bottom plot.

Loading the files into R is the slowest part of this process, but at a few seconds per file it's still
much faster than the typical `readMSdata` used by `xcms` or `MSnbase`. This is done via a custom
XML parser written around Hadley's `xml2` library that also has the memory advantage of using external
pointers until the data itself is needed. Switching between masses should be super fast thanks to `data.table`'s
subsetting and binary search. While there's currently no support for mzXML files, this is a top
priority for future development.

## Requirements
Several R libraries are required, but none of which should be difficult to install.
  - shiny
  - shinyDirectoryInput (handles intuitive file selection)
  - plotly (for making interactive plots)
  - data.table (for super-fast indexing)
  - xml2 (for parsing the mzML file)
  - base64enc (for parsing *m/z* and intensity data encoded as binary arrays)

Additionally, this app sources the custom script file "RaMS_custom.R" which
borrows several functions from the upcoming `RaMS` package and optimizes
them for MSMS data extraction.

Also, `shinyDirectoryInput` requires Powershell 3.0 or newer.

## Usage
The app can be launched on Windows simply by double-clicking the "runDataBrowser.cmd"
command. This should initialize both an R script in the background and a browser
window to interact with the app. Support for Mac and Linux upcoming shortly.

First, locate the mzML files of interest by clicking on the "..." button of the
directory input window. Files ending in .mzML are automatically detected and rendered
below. Clicking on a file name will "stage" it for loading, and multiple files can
be loaded simultaneously this way.

Once the files render, click around the upper plot or enter a new mass in the sidebar.
If a given file includes MSMS data, points on the upper plot will appear on the 
chromatogram and clicking on those will render the lower plot of fragments collected
near that scan.

Additional files can be added at any time in the same way, and Plotly should
assign them new colors. The legends in Plotly charts are also interactive - click
on a single trace to hide it, and double click to select only that one. Zooming
is also supported.

## Bugs
Absolutely full of 'em. Let me know if you find some, and what you were doing to cause it!