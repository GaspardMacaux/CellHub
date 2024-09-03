############################## Server ############################## 

############################## Library ############################## 

#Interface
library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinyjs)
library(shinyFiles)
library(fontawesome)
library(jsonlite)
library(tools)
library(R.utils)

#Graphics and tables
library(plotly)
library(stringr)
library(dplyr)
library(patchwork)
library(Matrix)
library(ggplot2)
library(ggrepel)
library(readr)
library(DT)
library(cowplot)
library(data.table)
library(httr)
library(rhdf5)
library(tidyverse)
library(reticulate)
library(colourpicker)
library(viridis)
library(uwot)

#Single cell
library(Seurat)
library(Signac)
library(SeuratWrappers)
library(SeuratDisk)
library(SeuratData)
library(monocle3)
library(anndata)
library(clustree)

#Genomes         
library(EnsDb.Hsapiens.v86)
library(biovizBase)

############################## Source Files ############################## 
# Importer les données
options(shiny.maxRequestSize = 20000*1024^2)  # 2000 MB

# Source server files:
source("single_dataset_server.R")
source("multiple_datasets_server.R")
source("trajectory_server.R")
source("main_ui.R")

server <- function(input, output, session) {
  single_dataset_server(input, output, session)
  multiple_datasets_server(input, output, session)
  trajectory_server(input, output, session)
}

# Exécuter l'application
options(shiny.launch.browser = TRUE)

shinyApp(ui = ui, server = server)
