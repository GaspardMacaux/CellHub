# download_functions_newarch.R

generateFileName <- function(object_name = "Dataset", plot_name = NULL, file_type = "pdf") {
  # Generate standardized file name for downloads
  # Args:
  #   object_name: Name of the dataset/object
  #   plot_name: Name of the plot (optional)
  #   file_type: File extension (pdf, png, tiff, csv, rds, etc.)
  # Returns:
  #   Formatted file name with timestamp
  
  # Clean object name (remove special characters)
  clean_object_name <- gsub("[^A-Za-z0-9_-]", "_", object_name)
  
  # Create base name
  if (!is.null(plot_name)) {
    clean_plot_name <- gsub("[^A-Za-z0-9_-]", "_", plot_name)
    base_name <- paste(clean_object_name, clean_plot_name, sep = "_")
  } else {
    base_name <- clean_object_name
  }
  
  # Add timestamp
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  file_name <- paste(base_name, timestamp, sep = "_")
  
  # Add extension
  final_name <- paste0(file_name, ".", file_type)
  
  return(final_name)
}

downloadPlot <- function(plot_object, object_name, plot_name, file_type = "pdf", 
                         width = 10, height = 8, dpi = 300) {
  # Download plot with specified parameters
  # Args:
  #   plot_object: ggplot or other plot object
  #   object_name: Name of the dataset
  #   plot_name: Name of the plot
  #   file_type: "pdf", "png", or "tiff"
  #   width: Plot width in inches
  #   height: Plot height in inches
  #   dpi: DPI for raster files (PNG, TIFF)
  # Returns:
  #   List with filename and content function
  
  req(plot_object)
  
  tryCatch({
    # Generate file name
    filename <- generateFileName(object_name, plot_name, file_type)
    
    # Create content function based on file type
    content <- function(file) {
      ggsave(file, plot = plot_object, device = file_type, width = width, height = height, dpi = dpi)
    }
    
    content_type <- switch(file_type,
                           "pdf" = "application/pdf",
                           "png" = "image/png",
                           "tiff" = "image/tiff",
                           "jpeg" = "image/jpeg"
    )
    
    return(list(
      filename = filename,
      content = content,
      contentType = content_type
    ))
    
  }, error = function(e) {
    message(paste("Error creating plot download:", e$message))
    return(NULL)
  })
}

downloadCSV <- function(data_object, object_name, data_name = "data") {
  # Download data as CSV file
  # Args:
  #   data_object: Data frame or matrix to download
  #   object_name: Name of the dataset
  #   data_name: Name of the data (e.g., "metadata", "markers", etc.)
  # Returns:
  #   List with filename and content function
  
  req(data_object)
  
  tryCatch({
    # Convert to data frame if needed
    if (is.matrix(data_object)) {
      data_df <- as.data.frame(data_object)
    } else if (is.data.frame(data_object)) {
      data_df <- data_object
    } else {
      stop("Data object must be a data frame or matrix")
    }
    
    # Generate file name
    filename <- generateFileName(object_name, data_name, "csv")
    
    # Create content function
    content <- function(file) {
      write.csv(data_df, file, row.names = TRUE)
    }
    
    return(list(
      filename = filename,
      content = content,
      contentType = "text/csv"
    ))
    
  }, error = function(e) {
    message(paste("Error creating CSV download:", e$message))
    return(NULL)
  })
}


downloadOpenDocumentSpreadsheet <- function(data_lists, object_name, data_name = "gene_lists") {
  # Download multiple data sets as OpenDocument Spreadsheet with multiple sheets
  # Args:
  #   data_lists: Named list of data frames or vectors
  #   object_name: Name of the dataset
  #   data_name: Name of the data type
  # Returns:
  #   List with filename and content function
  
  req(data_lists)
  
  tryCatch({
    # Generate file name with .ods extension
    filename <- generateFileName(object_name, data_name, "ods")
    
    # Create content function
    content <- function(file) {
      # Load required library
      if (!requireNamespace("readODS", quietly = TRUE)) {
        stop("Package 'readODS' is required but not installed. Please install it with: install.packages('readODS')")
      }
      library(readODS)
      
      # Prepare data for ODS format
      sheet_data <- list()
      
      # Create data frame for each gene list
      for (list_name in names(data_lists)) {
        data_content <- data_lists[[list_name]]
        
        # Create clean sheet name (ODS sheet names have restrictions)
        sheet_name <- gsub("[\\/:*?\"<>|\\[\\]]", "_", list_name)
        sheet_name <- gsub("&", "and", sheet_name)  # Replace & with 'and'
        sheet_name <- substr(sheet_name, 1, 30)  # Limit length
        
        # Ensure unique sheet names
        original_name <- sheet_name
        counter <- 1
        while (sheet_name %in% names(sheet_data)) {
          sheet_name <- paste0(substr(original_name, 1, 27), "_", counter)
          counter <- counter + 1
        }
        
        # Handle different data types
        if (is.vector(data_content)) {
          if (length(data_content) > 0) {
            df_to_write <- data.frame(
              Gene = as.character(data_content), 
              Gene_Set = rep(list_name, length(data_content)),
              stringsAsFactors = FALSE
            )
          } else {
            df_to_write <- data.frame(
              Gene = "No genes found",
              Gene_Set = list_name,
              stringsAsFactors = FALSE
            )
          }
        } else if (is.data.frame(data_content)) {
          df_to_write <- data_content
        } else {
          df_to_write <- data.frame(
            Gene = "Invalid data type",
            Gene_Set = list_name,
            stringsAsFactors = FALSE
          )
        }
        
        sheet_data[[sheet_name]] <- df_to_write
      }
      
      # Add index sheet
      index_data <- data.frame(
        Sheet_Name = names(sheet_data),
        Original_Name = names(data_lists),
        Gene_Count = sapply(data_lists, function(x) {
          if (is.vector(x)) length(x) else if (is.data.frame(x)) nrow(x) else 0
        }),
        stringsAsFactors = FALSE
      )
      sheet_data[["00_Index"]] <- index_data
      
      # Write ODS file - add debug messages
      message("About to write ODS file with ", length(sheet_data), " sheets")
      message("File path: ", file)
      
      readODS::write_ods(sheet_data, path = file)
      
      message("ODS file written successfully")
      
      # Verify file was created
      if (file.exists(file)) {
        message("File verified: ", file.size(file), " bytes")
      } else {
        stop("ODS file was not created!")
      }
    }
    
    return(list(
      filename = filename,
      content = content,
      contentType = "application/octet-stream"
    ))
    
  }, error = function(e) {
    message("ERROR in downloadOpenDocumentSpreadsheet: ", e$message)
    return(NULL)  # <-- C'est ça qui cause le problème !
  })
}



downloadSeuratObject <- function(seurat_object, object_name, show_modal = TRUE) {
  # Download Seurat object as RDS file
  # Args:
  #   seurat_object: Seurat object to download
  #   object_name: Name of the dataset
  #   show_modal: Whether to show modal dialog during download
  # Returns:
  #   List with filename and content function
  
  req(seurat_object)
  
  tryCatch({
    # Validate Seurat object
    if (!inherits(seurat_object, "Seurat")) {
      stop("Object is not a valid Seurat object")
    }
    
    # Generate file name
    filename <- generateFileName(object_name, "seurat_object", "rds")
    
    # Create content function with modal management
    content <- function(file) {
      if (show_modal) {
        showModal(modalDialog(
          title = "Please Wait",
          "Preparing the Seurat object for download...",
          footer = NULL,
          easyClose = FALSE
        ))
      }
      
      tryCatch({
        saveRDS(seurat_object, file)
        if (show_modal) {
          removeModal()
        }
      }, error = function(e) {
        if (show_modal) {
          removeModal()
        }
        stop(e$message)
      })
    }
    
    return(list(
      filename = filename,
      content = content,
      contentType = "application/octet-stream"
    ))
    
  }, error = function(e) {
    message(paste("Error creating Seurat object download:", e$message))
    return(NULL)
  })
}

createDownloadHandler <- function(reactive_data, object_name_reactive, data_name, 
                                  download_type, plot_params = NULL, show_modal = FALSE) {
  # Create a downloadHandler for Shiny
  # Args:
  #   reactive_data: Reactive expression containing the data/plot
  #   object_name_reactive: Reactive expression containing object name
  #   data_name: Name of the data/plot
  #   download_type: "csv", "plot", "seurat", "ods"
  #   plot_params: List with plot parameters (file_type, width, height, dpi)
  #   show_modal: Whether to show modal for Seurat downloads
  # Returns:
  #   downloadHandler function
  
  downloadHandler(
    filename = function() {
      message("=== DEBUG: Starting filename generation ===")
      message("download_type: ", download_type)
      
      tryCatch({
        # Get object name without causing reactivity issues
        obj_name <- isolate({
          if (is.reactive(object_name_reactive)) {
            object_name_reactive()
          } else {
            object_name_reactive
          }
        })
        message("obj_name: ", obj_name)
        
        # Get data name
        name <- isolate({
          if (is.reactive(data_name)) data_name() else data_name
        })
        message("name: ", name)
        
        if (download_type == "plot") {
          file_type <- isolate({
            if (!is.null(plot_params$file_type)) {
              if (is.reactive(plot_params$file_type)) plot_params$file_type() else plot_params$file_type
            } else "pdf"
          })
          result <- generateFileName(obj_name, name, file_type)
          message("Generated filename (plot): ", result)
          return(result)
        } else if (download_type == "csv") {
          result <- generateFileName(obj_name, name, "csv")
          message("Generated filename (csv): ", result)
          return(result)
        } else if (download_type == "ods") {
          result <- generateFileName(obj_name, name, "ods")
          message("Generated filename (ods): ", result)
          return(result)
        } else if (download_type == "seurat") {
          result <- generateFileName(obj_name, "seurat_object", "rds")
          message("Generated filename (seurat): ", result)
          return(result)
        }
        
        message("No matching download_type found!")
        return("no_match.txt")
        
      }, error = function(e) {
        message("ERROR in filename generation: ", e$message)
        return("download_error.txt")
      })
    },
    
    content = function(file) {
      tryCatch({
        # Get data
        data <- if (is.reactive(reactive_data)) {
          reactive_data()
        } else {
          reactive_data
        }
        req(data)
        
        # Get object name
        obj_name <- if (is.reactive(object_name_reactive)) {
          object_name_reactive()
        } else {
          object_name_reactive
        }
        
        # Get data name
        name <- if (is.reactive(data_name)) data_name() else data_name
        
        if (download_type == "plot") {
          # Get plot parameters
          file_type <- if (!is.null(plot_params$file_type)) {
            if (is.reactive(plot_params$file_type)) plot_params$file_type() else plot_params$file_type
          } else "pdf"
          
          width <- if (!is.null(plot_params$width)) {
            if (is.reactive(plot_params$width)) plot_params$width() else plot_params$width
          } else 10
          
          height <- if (!is.null(plot_params$height)) {
            if (is.reactive(plot_params$height)) plot_params$height() else plot_params$height
          } else 8
          
          dpi <- if (!is.null(plot_params$dpi)) {
            if (is.reactive(plot_params$dpi)) plot_params$dpi() else plot_params$dpi
          } else 300
          
          # Check if it's a Venn diagram (gList) or regular plot
          if (inherits(data, "gList")) {
            # For Venn diagrams, use grid plotting
            if (file_type == "pdf") {
              pdf(file, width = width, height = height)
            } else if (file_type == "png") {
              png(file, width = width, height = height, units = "in", res = dpi)
            } else if (file_type == "tiff") {
              tiff(file, width = width, height = height, units = "in", res = dpi, compression = "lzw")
            } else if (file_type == "jpeg") {
              jpeg(file, width = width, height = height, units = "in", res = dpi, quality = 100)
            }
            grid.draw(data)
            dev.off()
          } else {
            # For regular ggplots, use ggsave
            ggsave(file, plot = data, device = file_type, width = width, height = height, dpi = dpi)
          }
          
        } else if (download_type == "csv") {
          # Convert to data frame if needed
          if (is.matrix(data)) {
            data_df <- as.data.frame(data)
          } else if (is.data.frame(data)) {
            data_df <- data
          } else {
            stop("Data object must be a data frame or matrix")
          }
          write.csv(data_df, file, row.names = TRUE)
          
        } else if (download_type == "ods") {
          # Download as OpenDocument Spreadsheet
          req(data)
          download_info <- downloadOpenDocumentSpreadsheet(data, obj_name, name)
          download_info$content(file)
          
        } else if (download_type == "seurat") {
          # Validate Seurat object
          if (!inherits(data, "Seurat")) {
            stop("Object is not a valid Seurat object")
          }
          
          if (show_modal) {
            showModal(modalDialog(
              title = "Please Wait",
              "Preparing the Seurat object for download...",
              footer = NULL,
              easyClose = FALSE
            ))
          }
          
          saveRDS(data, file)
          
          if (show_modal) {
            removeModal()
          }
        }
        
      }, error = function(e) {
        if (download_type == "seurat" && show_modal) {
          removeModal()
        }
        message(paste("Download error:", e$message))
        showNotification(paste("Download failed:", e$message), type = "error")
      })
    },

  )
}



getObjectNameForDownload <- function(seurat_object = NULL, input_name = NULL, default_name = "Dataset") {
  # Get object name for download file naming
  # Args:
  #   seurat_object: Seurat object (to get project name)
  #   input_name: Input field with dataset name
  #   default_name: Default name if others fail
  # Returns:
  #   Object name for file naming
  
  tryCatch({
    # Try to get from Seurat object project name
    if (!is.null(seurat_object) && !is.null(seurat_object@project.name)) {
      return(seurat_object@project.name)
    }
    
    # Try to get from input field
    if (!is.null(input_name) && input_name != "") {
      return(input_name)
    }
    
    # Return default
    return(default_name)
    
  }, error = function(e) {
    return(default_name)
  })
}