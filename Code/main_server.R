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
library(uwot)


#Single cell
library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
library(SeuratData)
library(monocle3)
library(SingleCellExperiment)


library(anndata)
library(clustree)
library(circlize)
library(viridis)
library(callr)
library(EnsDb.Hsapiens.v86)
library(biovizBase)
library(biomaRt)
library(DoubletFinder)
library(nichenetr)
library(multinichenetr)

############################## Source Files ##############################
# Importer les données
options(shiny.maxRequestSize = 20000*1024^2)  # 2000 MB

# Source server files:
source("single_dataset_server.R")
source("multiple_datasets_server.R")
source("trajectory_server.R")
source("multinichenetr_server.R")
source("main_ui.R")

server <- function(input, output, session) {
  single_dataset_server(input, output, session)
  multiple_datasets_server(input, output, session)
  trajectory_server(input, output, session)
  multinichenetr_server(input, output, session)

}

# Exécuter l'application
options(shiny.launch.browser = TRUE)

shinyApp(ui = ui, server = server)
