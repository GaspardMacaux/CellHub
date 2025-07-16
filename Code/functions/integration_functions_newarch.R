# integration_functions_newarch.R


############################## Seurat v5 Layer Management ##############################

prepareSeuratV5ForIntegration <- function(seurat_object) {
  # Prepare Seurat v5 object for integration by managing layers properly
  # Args:
  #   seurat_object: Seurat object to prepare
  # Returns:
  #   Prepared Seurat object with proper layer structure
  
  if (packageVersion("Seurat") >= "5.0.0") {
    message("Preparing Seurat v5 object for integration...")
    
    tryCatch({
      # Join layers in RNA assay if multiple layers exist
      if ("RNA" %in% names(seurat_object@assays)) {
        rna_layers <- names(seurat_object[["RNA"]]@layers)
        if (length(rna_layers) > 1) {
          message(paste("Found", length(rna_layers), "RNA layers, joining..."))
          seurat_object[["RNA"]] <- JoinLayers(seurat_object[["RNA"]])
        }
      }
      
      # Ensure proper layer structure
      if ("RNA" %in% names(seurat_object@assays)) {
        # Check if counts layer exists
        if (!"counts" %in% names(seurat_object[["RNA"]]@layers)) {
          message("Adding counts layer...")
          seurat_object[["RNA"]]@layers$counts <- seurat_object[["RNA"]]@layers[[1]]
        }
        
        # Check if data layer exists
        if (!"data" %in% names(seurat_object[["RNA"]]@layers)) {
          message("Adding data layer...")
          seurat_object[["RNA"]]@layers$data <- seurat_object[["RNA"]]@layers[[1]]
        }
      }
      
      message("Seurat v5 object prepared successfully")
      
    }, error = function(e) {
      message(paste("Warning: Could not prepare Seurat v5 object:", e$message))
    })
  }
  
  return(seurat_object)
}

cleanIntegratedSeuratV5 <- function(seurat_object) {
  # Clean integrated Seurat v5 object to resolve layer issues
  # Args:
  #   seurat_object: Integrated Seurat object
  # Returns:
  #   Cleaned Seurat object
  
  if (packageVersion("Seurat") >= "5.0.0") {
    message("Cleaning integrated Seurat v5 object...")
    
    tryCatch({
      # Join layers in all assays
      for (assay_name in names(seurat_object@assays)) {
        assay_layers <- names(seurat_object[[assay_name]]@layers)
        if (length(assay_layers) > 1) {
          message(paste("Joining layers in", assay_name, "assay"))
          seurat_object[[assay_name]] <- JoinLayers(seurat_object[[assay_name]])
        }
      }
      
      message("Seurat v5 object cleaned successfully")
      
    }, error = function(e) {
      message(paste("Warning: Could not clean Seurat v5 object:", e$message))
    })
  }
  
  return(seurat_object)
}



############################## Data Preprocessing Functions ##############################

preprocessRawDataset <- function(file_path, dataset_type, species, dataset_name, 
                                 qc_params) {
  # Preprocess a single raw dataset for integration
  # Args:
  #   file_path: Path to the uploaded file
  #   dataset_type: "snRNA_merge", "multiome_merge", or "seurat_object_merge"
  #   species: Species for mitochondrial pattern
  #   dataset_name: Name for the dataset
  #   qc_params: List with min_features_merge, max_features_merge, max_mt_percent_merge
  # Returns:
  #   Processed Seurat object ready for integration
  
  message(paste("Starting preprocessing for dataset:", dataset_name))
  
  # Define mitochondrial pattern based on species
  mt_pattern <- switch(species,
                       "mouse" = "^mt-",
                       "human" = "^MT-",
                       "rat" = "^Mt-",
                       "^mt-")
  
  seurat_object <- NULL
  
  tryCatch({
    if (dataset_type == "seurat_object_merge") {
      # Load pre-processed Seurat object
      message("Loading pre-processed Seurat object")
      seurat_object <- loadSeuratObject(file_path, add_dataset_column = FALSE)
      
      if (!inherits(seurat_object, "Seurat")) {
        stop("Loaded file is not a valid Seurat object")
      }
      DefaultAssay(seurat_object) <- "RNA"
      
    } else {
      # Process raw data
      message("Processing raw 10X data")
      
      # Create temporary directory for extraction
      temp_dir <- createTempDirectory(paste0("dataset_", dataset_name))
      
      # Extract data
      unzip(file_path, exdir = temp_dir)
      data_10x <- Read10X(temp_dir)
      
      # Create Seurat object based on data type
      if (dataset_type == "snRNA_merge") {
        seurat_object <- CreateSeuratObject(
          counts = data_10x,
          project = dataset_name,
          min.cells = 3,
          min.features = qc_params$min_features_merge
        )
      } else if (dataset_type == "multiome_merge") {
        # Handle multiome data - extract Gene Expression
        if (is.list(data_10x) && "Gene Expression" %in% names(data_10x)) {
          seurat_object <- CreateSeuratObject(
            counts = data_10x$`Gene Expression`,
            project = dataset_name,
            min.cells = 3,
            min.features = qc_params$min_features_merge
          )
        } else {
          stop("Could not find Gene Expression data in multiome file")
        }
      }
      
      # Clean up temporary directory
      unlink(temp_dir, recursive = TRUE)
    }
    
    # Apply QC and preprocessing steps
    seurat_object <- applyQCAndPreprocessing(
      seurat_object = seurat_object,
      mt_pattern = mt_pattern,
      qc_params = qc_params,
      dataset_name = dataset_name,
      is_preprocessed = (dataset_type == "seurat_object_merge")
    )
    
    message(paste("Dataset", dataset_name, "preprocessing completed successfully"))
    return(seurat_object)
    
  }, error = function(e) {
    # Clean up on error
    if (exists("temp_dir") && dir.exists(temp_dir)) {
      unlink(temp_dir, recursive = TRUE)
    }
    stop(paste("Error preprocessing dataset", dataset_name, ":", e$message))
  })
}

applyQCAndPreprocessing <- function(seurat_object, mt_pattern, qc_params, 
                                    dataset_name, is_preprocessed = FALSE) {
  # Apply QC filtering and standard preprocessing steps with v5 compatibility
  # Args:
  #   seurat_object: Seurat object to process
  #   mt_pattern: Mitochondrial gene pattern
  #   qc_params: List with min_features_merge, max_features_merge, max_mt_percent_merge
  #   dataset_name: Dataset identifier
  #   is_preprocessed: Whether object is already processed
  # Returns:
  #   QC-filtered and preprocessed Seurat object
  
  message("Applying QC and preprocessing steps")
  
  # Prepare for Seurat v5 if needed
  seurat_object <- prepareSeuratV5ForIntegration(seurat_object)
  
  # Add mitochondrial percentage if not present
  if (!"percent.mt" %in% colnames(seurat_object@meta.data)) {
    seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = mt_pattern)
  }
  
  # Log pre-filter statistics
  pre_filter_cells <- ncol(seurat_object)
  message(paste("Pre-filter cell count:", pre_filter_cells))
  
  # Apply QC filters using UI parameters
  seurat_object <- subset(seurat_object,
                          subset = nFeature_RNA > qc_params$min_features_merge &
                            nFeature_RNA < qc_params$max_features_merge &
                            percent.mt < qc_params$max_mt_percent_merge)
  
  post_filter_cells <- ncol(seurat_object)
  message(paste("Post-filter cell count:", post_filter_cells))
  message(paste("Retained", round(post_filter_cells/pre_filter_cells * 100, 2), "% of cells"))
  
  # Validate that cells remain after filtering
  if (post_filter_cells == 0) {
    stop(paste("No cells remain after QC filtering for dataset", dataset_name, 
               "- consider relaxing QC parameters"))
  }
  
  # Apply standard preprocessing if not already done
  if (!is_preprocessed) {
    message("Applying normalization, scaling, and PCA")
    
    # Use layer-aware functions for Seurat v5
    seurat_object <- NormalizeData(seurat_object, verbose = FALSE)
    seurat_object <- FindVariableFeatures(seurat_object,
                                          selection.method = "vst",
                                          nfeatures = 4000,
                                          verbose = FALSE)
    seurat_object <- ScaleData(seurat_object,
                               features = rownames(seurat_object),
                               verbose = FALSE)
    seurat_object <- RunPCA(seurat_object,
                            features = VariableFeatures(seurat_object),
                            verbose = FALSE)
  }
  
  # Add dataset metadata
  seurat_object$dataset <- dataset_name
  seurat_object$orig.ident <- dataset_name
  
  return(seurat_object)
}

############################## Batch Processing Function ##############################

processDatasetsForIntegration <- function(file_inputs, dataset_types, dataset_names, 
                                          species, qc_params) {
  # Process multiple datasets for integration
  # Args:
  #   file_inputs: List of file paths
  #   dataset_types: Vector of dataset types
  #   dataset_names: Vector of dataset names
  #   species: Species identifier
  #   qc_params: QC parameters from UI
  # Returns:
  #   List of processed Seurat objects
  
  if (length(file_inputs) != length(dataset_types) || 
      length(file_inputs) != length(dataset_names)) {
    stop("Mismatch in number of files, types, and names")
  }
  
  message(paste("Processing", length(file_inputs), "datasets for integration"))
  
  seurat_list <- list()
  
  for (i in 1:length(file_inputs)) {
    message(paste("Processing dataset", i, "of", length(file_inputs)))
    
    processed_object <- preprocessRawDataset(
      file_path = file_inputs[[i]],
      dataset_type = dataset_types[i],
      species = species,
      dataset_name = dataset_names[i],
      qc_params = qc_params
    )
    
    seurat_list[[i]] <- processed_object
  }
  
  # Validate all objects
  validateSeuratList(seurat_list)
  
  message("All datasets processed successfully")
  return(seurat_list)
}

############################## Integration Functions ##############################

performDataIntegration <- function(seurat_list, integration_method = "standard") {
  # Perform data integration using specified method
  # Args:
  #   seurat_list: List of preprocessed Seurat objects
  #   integration_method: "standard" or "simple"
  # Returns:
  #   Integrated Seurat object
  
  if (length(seurat_list) < 2) {
    stop("At least two datasets are required for integration")
  }
  
  message(paste("Starting", integration_method, "integration with", length(seurat_list), "datasets"))
  
  # Validate input objects
  validateSeuratList(seurat_list)
  
  if (integration_method == "simple") {
    return(performSimpleMerge(seurat_list))
  } else {
    return(performStandardIntegration(seurat_list))
  }
}

performStandardIntegration <- function(seurat_list) {
  # Perform standard Seurat integration workflow with v5 compatibility
  # Args:
  #   seurat_list: List of preprocessed Seurat objects
  # Returns:
  #   Integrated Seurat object with "integrated" assay
  
  message("Starting standard Seurat integration")
  
  # Prepare objects for Seurat v5 integration
  seurat_list <- lapply(seurat_list, prepareSeuratV5ForIntegration)
  
  # Calculate dataset statistics for parameter adjustment
  dataset_sizes <- sapply(seurat_list, function(obj) ncol(obj))
  min_cells <- min(dataset_sizes)
  total_cells <- sum(dataset_sizes)
  
  message(paste("Dataset sizes:", paste(dataset_sizes, collapse = ", ")))
  message(paste("Minimum cells:", min_cells, "Total cells:", total_cells))
  
  # Adjust integration parameters based on dataset sizes
  integration_params <- calculateIntegrationParameters(min_cells)
  
  message(paste("Using integration parameters:",
                "k.filter =", integration_params$k_filter,
                "k.score =", integration_params$k_score,
                "k.anchor =", integration_params$k_anchor,
                "k.weight =", integration_params$k_weight))
  
  # Select integration features
  message("Selecting integration features")
  features <- SelectIntegrationFeatures(
    object.list = seurat_list,
    nfeatures = 2000
  )
  
  # Find integration anchors
  message("Finding integration anchors")
  anchors <- FindIntegrationAnchors(
    object.list = seurat_list,
    dims = 1:30,
    k.filter = integration_params$k_filter,
    k.score = integration_params$k_score,
    k.anchor = integration_params$k_anchor,
    anchor.features = features
  )
  
  # Integrate data
  message("Integrating datasets")
  integrated_object <- IntegrateData(
    anchorset = anchors,
    dims = 1:30,
    features.to.integrate = features,
    k.weight = integration_params$k_weight
  )
  
  # Clean integrated object for v5 compatibility
  integrated_object <- cleanIntegratedSeuratV5(integrated_object)
  
  # Set default assay and clean up metadata
  DefaultAssay(integrated_object) <- "integrated"
  integrated_object <- cleanIntegratedMetadata(integrated_object)
  
  message("Standard integration completed successfully")
  return(integrated_object)
}


performSimpleMerge <- function(seurat_list) {
  # Perform simple merge without integration, preserving original clusters
  # Args:
  #   seurat_list: List of preprocessed Seurat objects
  # Returns:
  #   Merged Seurat object with preserved cluster information
  
  message("Starting simple merge without integration")
  
  # Prepare objects for merging - preserve cluster information
  for (i in 1:length(seurat_list)) {
    dataset_name <- paste0("Dataset_", i)
    seurat_list[[i]]$dataset_origin <- dataset_name
    
    # Preserve existing cluster information if available
    if ("seurat_clusters" %in% colnames(seurat_list[[i]]@meta.data)) {
      seurat_list[[i]]$original_clusters <- paste0(dataset_name, "_C", seurat_list[[i]]$seurat_clusters)
      message(paste("Preserved clusters for", dataset_name))
    }
  }
  
  # Perform merge
  merged_object <- merge(
    seurat_list[[1]], 
    y = seurat_list[-1],
    add.cell.ids = paste0("Dataset_", 1:length(seurat_list)),
    project = "Simple_Merged_Datasets"
  )
  
  # Set appropriate identities
  if ("original_clusters" %in% colnames(merged_object@meta.data)) {
    Idents(merged_object) <- "original_clusters"
    message("Set identities to original clusters")
  } else {
    Idents(merged_object) <- "dataset_origin"
    message("Set identities to dataset origin")
  }
  
  # Mark merge method
  merged_object$merge_method <- "simple_merge"
  
  message(paste("Simple merge completed:", ncol(merged_object), "total cells"))
  return(merged_object)
}

############################## Helper Functions ##############################

calculateIntegrationParameters <- function(min_cells) {
  # Calculate appropriate integration parameters based on dataset size
  # Args:
  #   min_cells: Minimum number of cells across datasets
  # Returns:
  #   List of integration parameters
  
  # Conservative parameter calculation to avoid errors
  k_filter <- min(30, max(5, floor(min_cells / 15)))
  k_score <- min(20, max(5, floor(min_cells / 20)))
  k_anchor <- min(5, max(2, floor(min_cells / 50)))
  k_weight <- min(30, max(5, floor(min_cells / 15)))
  
  return(list(
    k_filter = k_filter,
    k_score = k_score,
    k_anchor = k_anchor,
    k_weight = k_weight
  ))
}

cleanIntegratedMetadata <- function(seurat_object) {
  # Clean and standardize metadata for integrated object
  # Args:
  #   seurat_object: Integrated Seurat object
  # Returns:
  #   Seurat object with cleaned metadata
  
  # Standardize dataset column name
  metadata <- seurat_object@meta.data
  if ("orig.ident" %in% colnames(metadata) && !"dataset" %in% colnames(metadata)) {
    colnames(metadata)[colnames(metadata) == "orig.ident"] <- "dataset"
    seurat_object <- AddMetaData(seurat_object, metadata = metadata)
  }
  
  return(seurat_object)
}

validateSeuratList <- function(seurat_list) {
  # Validate that all objects in list are proper Seurat objects
  # Args:
  #   seurat_list: List of objects to validate
  # Returns:
  #   TRUE if all valid, throws error otherwise
  
  if (length(seurat_list) == 0) {
    stop("Empty Seurat object list provided")
  }
  
  for (i in 1:length(seurat_list)) {
    if (is.null(seurat_list[[i]])) {
      stop(paste("Seurat object", i, "is NULL"))
    }
    
    if (!inherits(seurat_list[[i]], "Seurat")) {
      stop(paste("Object", i, "is not a valid Seurat object"))
    }
    
    if (ncol(seurat_list[[i]]) == 0) {
      stop(paste("Seurat object", i, "has no cells"))
    }
  }
  
  message("All Seurat objects validated successfully")
  return(TRUE)
}

extractQCParams <- function(input) {
  # Extract QC parameters from Shiny input
  # Args:
  #   input: Shiny input object
  # Returns:
  #   List of QC parameters
  
  return(list(
    min_features_merge = input$min_features_merge,
    max_features_merge = input$max_features_merge,
    max_mt_percent_merge = input$max_mt_percent_merge
  ))
}

############################## Metadata Management Functions ##############################

addDatasetMetadata <- function(seurat_object, metadata_fields, metadata_values, dataset_name) {
  # Add custom metadata to a specific dataset
  # Args:
  #   seurat_object: Integrated Seurat object
  #   metadata_fields: Vector of metadata field names
  #   metadata_values: Vector of metadata values for this dataset
  #   dataset_name: Name of the dataset to update
  # Returns:
  #   Seurat object with updated metadata
  
  if (length(metadata_fields) != length(metadata_values)) {
    stop("Metadata fields and values must have the same length")
  }
  
  # Find cells belonging to this dataset
  dataset_cells <- which(seurat_object@meta.data$dataset == dataset_name)
  
  if (length(dataset_cells) == 0) {
    warning(paste("No cells found for dataset:", dataset_name))
    return(seurat_object)
  }
  
  # Add metadata for each field
  for (i in 1:length(metadata_fields)) {
    field_name <- metadata_fields[i]
    field_value <- metadata_values[i]
    
    # Initialize field if it doesn't exist
    if (!field_name %in% colnames(seurat_object@meta.data)) {
      seurat_object@meta.data[[field_name]] <- NA
    }
    
    # Set values for this dataset
    seurat_object@meta.data[dataset_cells, field_name] <- field_value
  }
  
  message(paste("Added metadata for dataset:", dataset_name))
  return(seurat_object)
}

processMetadataFromUI <- function(seurat_object, input, num_fields) {
  # Process metadata inputs from UI and add to Seurat object
  # Args:
  #   seurat_object: Integrated Seurat object
  #   input: Shiny input object
  #   num_fields: Number of metadata fields to process
  # Returns:
  #   Seurat object with updated metadata
  
  datasets <- unique(seurat_object@meta.data$dataset)
  
  for (j in 1:num_fields) {
    metadata_field_name <- input[[paste0("metadata_name_", j)]]
    
    if (is.null(metadata_field_name) || metadata_field_name == "") {
      next  # Skip empty fields
    }
    
    # Initialize field in metadata
    if (!metadata_field_name %in% colnames(seurat_object@meta.data)) {
      seurat_object@meta.data[[metadata_field_name]] <- NA
    }
    
    # Process each dataset
    for (dataset in datasets) {
      metadata_field_value <- input[[paste0("metadata_value_", dataset, "_", j)]]
      
      if (!is.null(metadata_field_value) && metadata_field_value != "") {
        dataset_cells <- which(seurat_object@meta.data$dataset == dataset)
        seurat_object@meta.data[dataset_cells, metadata_field_name] <- metadata_field_value
      }
    }
  }
  
  message("Processed metadata from UI successfully")
  return(seurat_object)
}

############################## Memory Management Functions ##############################

cleanupIntegrationMemory <- function() {
  # Clean up memory after integration process
  # Returns:
  #   Nothing (performs garbage collection)
  
  # Force garbage collection
  gc()
  
  # Clean temporary directories related to integration
  temp_dirs <- list.dirs(tempdir(), full.names = TRUE, recursive = FALSE)
  integration_dirs <- temp_dirs[grepl("dataset_|integration_", basename(temp_dirs))]
  
  for (dir in integration_dirs) {
    if (dir.exists(dir)) {
      unlink(dir, recursive = TRUE)
      message(paste("Cleaned temporary directory:", basename(dir)))
    }
  }
  
  message("Memory cleanup completed")
}

monitorIntegrationProgress <- function(current_step, total_steps, step_name) {
  # Monitor and log integration progress
  # Args:
  #   current_step: Current step number
  #   total_steps: Total number of steps
  #   step_name: Name of current step
  # Returns:
  #   Progress percentage
  
  progress_pct <- round((current_step / total_steps) * 100, 1)
  message(paste("Integration Progress:", progress_pct, "% -", step_name))
  
  return(progress_pct)
}