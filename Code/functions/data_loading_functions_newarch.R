# data_loading_functions_newarch.R

############################## Utility Functions ##############################

getDatasetFileName <- function(file_inputs = NULL, default_name = "Dataset") {
  # Extract dataset file name from various input sources
  # Args:
  #   file_inputs: List of file inputs to check (e.g., list(input$file, input$load_seurat_file))
  #   default_name: Default name if no file found
  
  tryCatch({
    file_name <- ""
    
    # If specific inputs provided, check them
    if (!is.null(file_inputs)) {
      for (file_input in file_inputs) {
        if (!is.null(file_input) && 
            !is.null(file_input$name) && 
            nchar(file_input$name) > 0) {
          file_name <- tools::file_path_sans_ext(basename(file_input$name))
          break  # Use first valid file found
        }
      }
    }
    
    # Return file name or default
    if (nchar(file_name) == 0) {
      return(default_name)
    } else {
      return(file_name)
    }
    
  }, error = function(e) {
    message(paste("Error extracting file name:", e$message))
    return(default_name)
  })
}

getClusters <- function(seurat_object) {
  # Get clusters from Seurat object
  # Args:
  #   seurat_object: Seurat object
  # Returns:
  #   Vector of cluster levels
  
  if (is.null(seurat_object)) {
    return(NULL)
  }
  
  tryCatch({
    return(levels(Idents(seurat_object)))
  }, error = function(e) {
    message(paste("Error getting clusters:", e$message))
    return(NULL)
  })
}

find10XFolder <- function(base_dir) {
  # Find folder containing 10X files
  # Args:
  #   base_dir: Base directory to search in
  # Returns:
  #   Path to folder with 10X files or NULL
  
  message(paste("Searching for 10X files in:", base_dir))
  
  # List all folders recursively
  folders <- list.dirs(base_dir, recursive = TRUE, full.names = TRUE)
  
  # Search each folder for 10X files
  for (folder in folders) {
    # Check for compressed files (.gz)
    compressed_files <- list.files(
      folder,
      pattern = "barcodes\\.tsv\\.gz$|features\\.tsv\\.gz$|matrix\\.mtx\\.gz$|genes\\.tsv\\.gz$",
      full.names = TRUE
    )
    
    # Check for uncompressed files
    uncompressed_files <- list.files(
      folder,
      pattern = "barcodes\\.tsv$|features\\.tsv$|matrix\\.mtx$|genes\\.tsv$",
      full.names = TRUE
    )
    
    # Check if there are enough files (compressed or uncompressed)
    if (length(compressed_files) >= 3 || length(uncompressed_files) >= 3) {
      message(paste("Found 10X files in:", folder))
      
      # If files are not compressed, compress them
      if (length(uncompressed_files) >= 3 && length(compressed_files) < 3) {
        message("Compressing uncompressed 10X files...")
        for (file in uncompressed_files) {
          if (!file.exists(paste0(file, ".gz"))) {
            message(paste("Compressing:", file))
            R.utils::gzip(file, destname = paste0(file, ".gz"), overwrite = TRUE, remove = FALSE)
          }
        }
      }
      
      return(folder)
    }
  }
  return(NULL)
}

compressFiles <- function(files) {
  # Compress uncompressed files
  # Args:
  #   files: Vector of file paths to compress
  # Returns:
  #   Vector of compressed file paths
  
  gz_files <- sapply(files, function(file) {
    gz_file <- paste0(file, ".gz")
    if (!file.exists(gz_file)) {
      message(paste("Compressing file:", file))
      R.utils::gzip(file, destname = gz_file, overwrite = TRUE)
      message(paste("Compressed", file, "to", gz_file))
    } else {
      message(gz_file, "already exists.")
    }
    return(gz_file)
  })
  return(gz_files)
}

extractAndProcessZip <- function(zip_path, dataset_name = NULL) {
  # Extract and process ZIP file with handling for nested ZIP files
  # Args:
  #   zip_path: Path to ZIP file
  #   dataset_name: Name for dataset (optional)
  # Returns:
  #   Path to folder containing 10X files or NULL
  
  # Create unique identifier for this dataset
  if (is.null(dataset_name)) {
    dataset_name <- format(Sys.time(), "%Y%m%d_%H%M%S")
  }
  
  # Create specific temporary folder for this dataset
  temp_base_dir <- file.path(tempdir(), "single_dataset")
  dir.create(temp_base_dir, showWarnings = FALSE, recursive = TRUE)
  
  dataset_dir <- file.path(temp_base_dir, paste0("dataset_", dataset_name))
  dir.create(dataset_dir, showWarnings = FALSE, recursive = TRUE)
  
  message(paste("Extracting zip file to unique directory:", dataset_dir))
  
  # Extract ZIP file
  tryCatch({
    unzip(zip_path, exdir = dataset_dir)
    
    # Search for nested ZIP files
    nested_zips <- list.files(dataset_dir, pattern = "\\.zip$", recursive = TRUE, full.names = TRUE)
    
    while (length(nested_zips) > 0) {
      for (nested_zip in nested_zips) {
        message(paste("Extracting nested zip file:", nested_zip))
        nested_dir <- file.path(dirname(nested_zip), tools::file_path_sans_ext(basename(nested_zip)))
        dir.create(nested_dir, showWarnings = FALSE)
        unzip(nested_zip, exdir = nested_dir)
        file.remove(nested_zip)  # Remove ZIP after extraction
      }
      # Check if there are remaining nested ZIP files
      nested_zips <- list.files(dataset_dir, pattern = "\\.zip$", recursive = TRUE, full.names = TRUE)
    }
    
    # Search for folder containing 10X files
    valid_folder <- find10XFolder(dataset_dir)
    
    if (is.null(valid_folder)) {
      stop("No complete set of 10X files found in the uploaded ZIP")
    }
    
    return(valid_folder)
    
  }, error = function(e) {
    message(paste("Error extracting zip:", e$message))
    return(NULL)
  })
}

############################## Main Loading Functions ##############################

loadRawData <- function(file_path, dataset_type, species) {
  # Load raw data from either ZIP (10X files) or H5 format
  # Args:
  #   file_path: Path to ZIP or H5 file
  #   dataset_type: Type of dataset ("snRNA" or "multiome") 
  #   species: Species for mitochondrial pattern
  # Returns:
  #   Seurat object
  
  # Get file extension
  file_extension <- tools::file_ext(file_path)
  
  # Get mitochondrial pattern
  mt_pattern <- switch(species,
                       "human" = "^MT-",
                       "mouse" = "^mt-|^Mt-", 
                       "rat" = "^Mt-",
                       "^MT-"
  )
  
  seuratObj <- NULL
  
  if (file_extension == "zip") {
    # Handle ZIP file with 10X format
    seuratObj <- loadFromZip(file_path, dataset_type, mt_pattern)
    
  } else if (file_extension %in% c("h5", "hdf5")) {
    # Handle H5 file
    seuratObj <- loadFromH5(file_path, dataset_type, mt_pattern)
    
  } else {
    stop("Unsupported file format. Please provide ZIP or H5 file.")
  }
  
  if (is.null(seuratObj)) {
    stop("Failed to create Seurat object from raw data")
  }
  
  # Add metadata
  seuratObj@meta.data$dataset_type <- dataset_type
  seuratObj@meta.data$species <- species
  seuratObj@meta.data$source_format <- file_extension
  
  return(seuratObj)
}

loadFromZip <- function(zip_path, dataset_type, mt_pattern) {
  # Load data from ZIP file containing 10X files
  
  # Extract ZIP and find 10X files
  valid_folder <- extractAndProcessZip(zip_path)
  
  if (is.null(valid_folder)) {
    stop("Could not extract and process ZIP file")
  }
  
  # Load 10X data
  message("Reading 10X data...")
  data_10x <- Read10X(valid_folder)
  
  # Create Seurat object based on dataset type
  if (dataset_type == "snRNA") {
    seuratObj <- CreateSeuratObject(
      counts = data_10x,
      project = "RawData"
    )
    
  } else if (dataset_type == "multiome") {
    if (is.list(data_10x) && !is.matrix(data_10x)) {
      # Multimodal data
      rna_assay <- NULL
      if ("Gene Expression" %in% names(data_10x)) {
        rna_assay <- "Gene Expression"
      } else {
        rna_candidates <- grep("RNA|Gene|Expression", names(data_10x), ignore.case = TRUE, value = TRUE)
        if (length(rna_candidates) > 0) {
          rna_assay <- rna_candidates[1]
        }
      }
      
      if (is.null(rna_assay)) {
        stop("Could not identify RNA assay in multimodal data")
      }
      
      seuratObj <- CreateSeuratObject(
        counts = data_10x[[rna_assay]],
        project = "RawData"
      )
      
      # Add other assays
      for (assay_name in names(data_10x)) {
        if (assay_name != rna_assay) {
          seuratObj[[assay_name]] <- CreateAssayObject(counts = data_10x[[assay_name]])
        }
      }
      
    } else {
      # Single assay
      seuratObj <- CreateSeuratObject(
        counts = data_10x,
        project = "RawData"
      )
    }
  }
  
  # Add mitochondrial percentage
  seuratObj[["percent.mt"]] <- PercentageFeatureSet(seuratObj, pattern = mt_pattern)
  
  return(seuratObj)
}

loadFromH5 <- function(h5_path, dataset_type, mt_pattern) {
  # Load data from H5 file
  
  message("Reading H5 file...")
  h5_data <- Read10X_h5(h5_path)
  
  # Create Seurat object based on dataset type
  if (dataset_type == "snRNA") {
    # ✅ SIMPLE: No filtering conditions
    seuratObj <- CreateSeuratObject(
      counts = h5_data,
      project = "RawData"
    )
    
  } else if (dataset_type == "multiome") {
    if (is.list(h5_data) && !is.matrix(h5_data)) {
      # Multimodal H5
      rna_assay <- NULL
      if ("Gene Expression" %in% names(h5_data)) {
        rna_assay <- "Gene Expression"
      } else {
        rna_candidates <- grep("RNA|Gene|Expression", names(h5_data), ignore.case = TRUE, value = TRUE)
        if (length(rna_candidates) > 0) {
          rna_assay <- rna_candidates[1]
        }
      }
      
      if (is.null(rna_assay)) {
        stop("Could not identify RNA assay in H5 file")
      }
      
      rna_data <- h5_data[[rna_assay]]
      
      # ✅ SIMPLE: No filtering conditions
      seuratObj <- CreateSeuratObject(
        counts = rna_data,
        project = "RawData"
      )
      
      # Add other assays
      for (assay_name in names(h5_data)) {
        if (assay_name != rna_assay) {
          seuratObj[[assay_name]] <- CreateAssayObject(counts = h5_data[[assay_name]])
        }
      }
      
    } else {
      # Single assay H5
      # ✅ SIMPLE: No filtering conditions
      seuratObj <- CreateSeuratObject(
        counts = h5_data,
        project = "RawData"
      )
    }
  }
  
  # Add mitochondrial percentage
  seuratObj[["percent.mt"]] <- PercentageFeatureSet(seuratObj, pattern = mt_pattern)
  
  return(seuratObj)
}
  


loadSeuratObject <- function(rds_path, add_dataset_column = FALSE, dataset_name = NULL, clean_before = FALSE, module_type = "single") {
  # Load processed Seurat object from RDS file with enhanced functionality
  # Args:
  #   rds_path: Path to RDS file
  #   add_dataset_column: Whether to add/verify dataset column in metadata
  #   dataset_name: Name for dataset column (if NULL, extracted from filename)
  #   clean_before: Whether to clean workspace before loading
  #   module_type: Module type for cleaning ("single", "multiple", etc.)
  # Returns:
  #   Seurat object
  
  if (clean_before) {
    cleanWorkspace(module_type)
  }
  
  message("Loading Seurat object from RDS...")
  
  # Load RDS file
  seuratObj <- readRDS(rds_path)
  
  # Validate Seurat object
  if (!inherits(seuratObj, "Seurat")) {
    stop("File is not a valid Seurat object")
  }
  
  # ✅ KEEP YOUR COMPREHENSIVE SEURAT V5 HANDLING
  if (packageVersion("Seurat") >= "5.0.0") {
    message("Detected Seurat V5+, checking for layers to join...")
    tryCatch({
      # Handle RNA assay first
      if ("RNA" %in% names(seuratObj@assays)) {
        rna_layers <- names(seuratObj[["RNA"]]@layers)
        if (length(rna_layers) > 1) {
          message(paste("Found", length(rna_layers), "RNA layers, joining..."))
          seuratObj[["RNA"]] <- JoinLayers(seuratObj[["RNA"]])
          message("RNA layers joined successfully")
        } else {
          message("Only one RNA layer found, no joining needed")
        }
      }
      
      # ✅ YOUR ADDITION: Handle other assays too
      for (assay_name in names(seuratObj@assays)) {
        if (assay_name != "RNA") {
          assay_layers <- names(seuratObj[[assay_name]]@layers)
          if (length(assay_layers) > 1) {
            message(paste("Found", length(assay_layers), "layers in", assay_name, "assay, joining..."))
            seuratObj[[assay_name]] <- JoinLayers(seuratObj[[assay_name]])
            message(paste(assay_name, "layers joined successfully"))
          }
        }
      }
    }, error = function(e) {
      message(paste("Warning: Could not join layers:", e$message))
    })
  } else {
    message("Detected older Seurat version, no layer joining needed")
  }
  
  # ✅ ADD: Handle dataset column for multiple datasets workflow
  if (add_dataset_column) {
    if (!"dataset" %in% colnames(seuratObj@meta.data)) {
      message("Adding 'dataset' column to metadata")
      
      if (is.null(dataset_name)) {
        dataset_name <- tools::file_path_sans_ext(basename(rds_path))
        if (is.null(dataset_name) || dataset_name == "") {
          dataset_name <- "Dataset"
        }
      }
      
      seuratObj$dataset <- dataset_name
      message(paste("Set dataset name to:", dataset_name))
    } else {
      datasets <- unique(seuratObj$dataset)
      message(paste("Found existing datasets:", paste(datasets, collapse = ", ")))
    }
  }
  
  message(paste("Loaded Seurat object:", ncol(seuratObj), "cells,", nrow(seuratObj), "genes"))
  
  return(seuratObj)
}


generateInitialPlot <- function(seurat_obj, remove_axes = FALSE, remove_legend = FALSE) {
  # Generate initial clustering plot with color restoration
  # Args:
  #   seurat_obj: Seurat object
  #   remove_axes: Remove axes from plot
  #   remove_legend: Remove legend from plot
  # Returns:
  #   List with plot and updated seurat object
  
  if (is.null(seurat_obj) || !"umap" %in% names(seurat_obj@reductions)) {
    return(NULL)
  }
  
  # Handle cluster colors
  if (!is.null(seurat_obj@misc$cluster_colors)) {
    message("Restoring cluster colors from Seurat object.")
    print(seurat_obj@misc$cluster_colors)
    colors <- seurat_obj@misc$cluster_colors
  } else {
    message("No cluster colors found. Generating default colors.")
    n_clusters <- length(levels(Idents(seurat_obj)))
    colors <- rainbow(n_clusters)
    names(colors) <- levels(Idents(seurat_obj))
    seurat_obj@misc$cluster_colors <- colors
  }
  
  # Create plot
  plot <- DimPlot(seurat_obj, group.by = "ident", label = FALSE) + 
    ggtitle("") +
    scale_color_manual(values = colors)
  
  if (remove_axes) { plot <- plot + NoAxes() }
  if (remove_legend) { plot <- plot + NoLegend() }
  
  return(list(plot = plot, seurat_obj = seurat_obj))
}



