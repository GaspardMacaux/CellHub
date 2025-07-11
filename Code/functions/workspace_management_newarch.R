# workspace_management_newarch.R

cleanReactiveVariables <- function(variable_patterns = NULL, exclude_patterns = NULL) {
  # Clean reactive variables based on patterns
  # Args:
  #   variable_patterns: Vector of patterns to match variable names
  #   exclude_patterns: Vector of patterns to exclude from cleaning
  
  # Get all objects in current environment
  all_objects <- ls(envir = .GlobalEnv)
  
  # If no patterns specified, clean common reactive variable patterns
  if (is.null(variable_patterns)) {
    variable_patterns <- c("_plot$", "_data$", "_seurat$", "_object$", "gene_list", 
                           "subset_", "normalized_", "scaled_", "feature_", 
                           "neighbors_", "clustering_", "spatial_", "multiome_")
  }
  
  # Find objects matching patterns
  objects_to_clean <- character(0)
  for (pattern in variable_patterns) {
    matching_objects <- grep(pattern, all_objects, value = TRUE, ignore.case = TRUE)
    objects_to_clean <- c(objects_to_clean, matching_objects)
  }
  
  # Remove excluded patterns
  if (!is.null(exclude_patterns)) {
    for (exclude_pattern in exclude_patterns) {
      objects_to_clean <- objects_to_clean[!grepl(exclude_pattern, objects_to_clean, ignore.case = TRUE)]
    }
  }
  
  # Remove duplicates
  objects_to_clean <- unique(objects_to_clean)
  
  # Clean the objects
  cleaned_count <- 0
  for (obj_name in objects_to_clean) {
    if (exists(obj_name, envir = .GlobalEnv)) {
      tryCatch({
        obj <- get(obj_name, envir = .GlobalEnv)
        
        # ✅ FIX: Add proper type checking for reactive objects
        if (is.function(obj)) {
          # Check if it's a reactive function by testing for reactive class
          if (inherits(obj, "reactive") || inherits(obj, "reactiveVal")) {
            obj(NULL)
            cleaned_count <- cleaned_count + 1
          } else {
            # Skip regular functions like server functions
            message(paste("Skipping function:", obj_name, "- not a reactive object"))
          }
        } else if (is.list(obj) && "features" %in% names(obj)) {
          # For reactiveValues like gene_list
          obj$features <- NULL
          cleaned_count <- cleaned_count + 1
        }
      }, error = function(e) {
        message(paste("Could not clean object:", obj_name, "-", e$message))
      })
    }
  }
  message(paste("Cleaned", cleaned_count, "reactive variables"))
}

cleanTempDirectories <- function(prefix_patterns = c("single_dataset", "spatial", "multiome")) {
  # Clean temporary directories created by the app
  
  temp_dirs <- list.dirs(tempdir(), full.names = TRUE, recursive = FALSE)
  
  # Filter directories matching our patterns
  app_temp_dirs <- character(0)
  for (pattern in prefix_patterns) {
    matching_dirs <- temp_dirs[grepl(pattern, basename(temp_dirs))]
    app_temp_dirs <- c(app_temp_dirs, matching_dirs)
  }
  
  # Remove directories
  removed_count <- 0
  for (dir in app_temp_dirs) {
    if (dir.exists(dir)) {
      unlink(dir, recursive = TRUE)
      message(paste("Removed temp directory:", basename(dir)))
      removed_count <- removed_count + 1
    }
  }
  
  message(paste("Removed", removed_count, "temporary directories"))
}

cleanWorkspace <- function(module = "all", custom_patterns = NULL) {
  # General workspace cleaning function
  # Args:
  #   module: "single", "spatial", "multiome", or "all"
  #   custom_patterns: Additional patterns to clean
  
  message(paste("Cleaning workspace for module:", module))
  
  patterns <- switch(module,
                     "single" = c("single_dataset", "_plot$", "_data$", "gene_list", "subset_", "normalized_", "scaled_"),
                     "spatial" = c("spatial_", "tissue_", "spot_"),
                     "multiome" = c("multiome_", "atac_", "rna_", "peak_"),
                     # ✅ ADD: Support for multiple datasets
                     "multiple" = c("multiple_datasets_", "_merge$", "_plot_merge$", "clustering_plot_merge", 
                                    "feature_plot_merge", "vln_plot_merge", "dot_plot_merge", "ridge_plot_merge", 
                                    "heatmap_plot_multidataset", "seurat_objects", "data_loaded", "merged_gene_tables"),
                     "all" = c("_plot$", "_data$", "_seurat$", "_object$", "gene_list", "subset_", "normalized_", 
                               "scaled_", "feature_", "neighbors_", "clustering_", "spatial_", "multiome_", 
                               "multiple_datasets_", "_merge$")
  )
  
  if (!is.null(custom_patterns)) {
    patterns <- c(patterns, custom_patterns)
  }
  
  cleanReactiveVariables(variable_patterns = patterns)
  
  # Clean temp directories
  temp_patterns <- switch(module,
                          "multiple" = "multiple_datasets",
                          "single" = "single_dataset", 
                          "spatial" = "spatial",
                          "multiome" = "multiome",
                          c("single_dataset", "multiple_datasets", "spatial", "multiome")
  )
  
  cleanTempDirectories(prefix_patterns = temp_patterns)
  
  gc()
  message(paste("Workspace cleaning completed for module:", module))
}

createTempDirectory <- function(prefix = "analysis") {
  # Create a timestamped temporary directory
  
  temp_dir <- file.path(tempdir(), paste0(prefix, "_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  message(paste("Created temp directory:", temp_dir))
  
  return(temp_dir)
}