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
library(CellChat)
library(usethis)
library(circlize)
library(callr)
library(viridis)
library(VennDiagram)
library(grid)


#Single cell
library(Seurat)
library(SeuratWrappers)
library(SeuratDisk)
library(SeuratData)
library(monocle3)
library(SingleCellExperiment)
library(harmony)
library(anndata)
library(clustree)
library(EnsDb.Hsapiens.v86)
library(biovizBase)
library(biomaRt)
library(DoubletFinder)
library(CellChat)
library(VennDiagram)
library(grid)
library(openxlsx)
library(arrow)
library(glmGamPoi)

library(RBioFormats)
library(magick)
library(tiff)

library(png)
library(jpeg)
library(EBImage)
library(rJava)
library(readODS)

############################## Source Files ##############################
# Importer les données
# Augmenter la limite de mémoire vectorielle
options(shiny.maxRequestSize = 300000*1024^2)  # Augmente à 30GB par fichier
options(future.globals.maxSize = 6000000 * 1024^2)  # 100GB
Sys.setenv('R_MAX_VSIZE'=20000000000000)



# Source server files:
source("single_dataset_server_newarch.R")
source("multiple_datasets_server_newarch.R")
source("cellchat_server_newarch.R")
source("spatial_transcriptomic_server_newarch.R")
source("trajectory_server_newarch.R")
source("multinichenet_server_newarch.R")
#Source function:
source("functions/data_loading_functions_newarch.R")
source("functions/workspace_management_newarch.R")
source("functions/download_functions_newarch.R")
source("functions/venn_diagram_functions_newarch.R")
source("functions/ui_update_functions_newarch.R")
source("functions/cells_genes_expressions_newarch.R")
source("functions/integration_functions_newarch.R")

#Source UI:
source("main_ui_newarch.R")

server <- function(input, output, session) {
  single_dataset_server(input, output, session)
  multiple_datasets_server(input, output, session)
  trajectory_server(input, output, session)
  multinichenet_server(input, output, session)
  cellchat_server(input, output, session)
  spatial_transcriptomic_server (input, output, session)
    
}

# Exécuter l'application
options(shiny.launch.browser = TRUE)

shinyApp(ui = ui, server = server)
