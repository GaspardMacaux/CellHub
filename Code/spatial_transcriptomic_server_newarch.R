# spatial_server.R - Server functions for Spatial Transcriptomics

spatial_transcriptomic_server <- function(input, output, session) {
  
  ############################## Reactive Variables ##############################
  
  # Main spatial object storage
  spatial_obj <- reactiveVal(NULL)
  sketching_applied <- reactiveVal(FALSE)
  available_assays <- reactiveVal(c("Spatial"))
  # Simple assay selection storage
  selected_assay <- reactiveVal("Spatial")
  
  
  # Processing state flags
  processing_states <- reactiveValues(
    qc_applied = FALSE,
    normalized = FALSE,
    pca_computed = FALSE,
    clustered = FALSE
  )
  
  # Individual state flags (for compatibility)
  spatial_normalized <- reactiveVal(FALSE)
  spatial_qc_applied <- reactiveVal(FALSE)
  spatial_pca_computed <- reactiveVal(FALSE)
  spatial_clustered <- reactiveVal(FALSE)
  
  # Plots for download
  spatial_plots <- reactiveValues()
  
  ############################## Helper Functions ##############################
  
  # Update available assays when object changes
  observe({
    req(spatial_obj())
    obj <- spatial_obj()
    assays <- names(obj@assays)
    available_assays(assays)
    
    # Update assay selection menus
    updateSelectInput(session, "spatial_assay_select", choices = assays, selected = assays[1])
    updateSelectInput(session, "spatial_viz_assay", choices = assays, selected = if("RNA" %in% assays) "RNA" else assays[1])
  })
  
  
  # Update object and state in one function
  update_spatial_object <- function(new_obj, state_updates = list()) {
    spatial_obj(new_obj)
    for (state_name in names(state_updates)) {
      processing_states[[state_name]] <- state_updates[[state_name]]
      # Also update individual flags
      if (state_name == "normalized") spatial_normalized(state_updates[[state_name]])
      if (state_name == "qc_applied") spatial_qc_applied(state_updates[[state_name]])
      if (state_name == "pca_computed") spatial_pca_computed(state_updates[[state_name]])
      if (state_name == "clustered") spatial_clustered(state_updates[[state_name]])
    }
  }
  
  # Get current processing status
  get_processing_status <- function() {
    obj <- spatial_obj()
    if (is.null(obj)) return("No data loaded")
    
    status_parts <- c()
    if (processing_states$qc_applied) status_parts <- c(status_parts, "QC")
    if (processing_states$normalized) status_parts <- c(status_parts, "Normalized")
    if (sketching_applied()) status_parts <- c(status_parts, "Sketched")
    if (processing_states$pca_computed) status_parts <- c(status_parts, "PCA")
    if (processing_states$clustered) status_parts <- c(status_parts, "Clustered")
    
    if (length(status_parts) == 0) return("Raw data")
    return(paste(status_parts, collapse = " → "))
  }
  
  # Add debug information display
  output$spatial_debug_info <- renderText({
    req(spatial_obj())
    
    obj <- spatial_obj()
    current_assay <- input$spatial_viz_assay %||% DefaultAssay(obj)
    
    debug_info <- paste(
      "Current assay:", current_assay,
      "\nAvailable assays:", paste(names(obj@assays), collapse = ", "),
      "\nDefault assay:", DefaultAssay(obj),
      "\nGenes in current assay:", nrow(obj[[current_assay]]),
      "\nHas spatial images:", length(obj@images) > 0
    )
    
    return(debug_info)
  })
  # Get spatial dataset name
  get_spatial_dataset_name <- function() {
    if (!is.null(spatial_obj())) {
      if ("orig.ident" %in% colnames(spatial_obj()@meta.data)) {
        return(unique(spatial_obj()@meta.data$orig.ident)[1])
      }
    }
    return("SpatialDataset")
  }
  
  # Validate spatial data structure
  validate_spatial_data <- function(seurat_obj) {
    required_assays <- c("Spatial")
    has_spatial <- any(required_assays %in% names(seurat_obj@assays))
    has_images <- length(seurat_obj@images) > 0
    
    has_coordinates <- TRUE
    if (has_images) {
      tryCatch({
        if (class(seurat_obj@images[[1]])[1] == "VisiumV2") {
          coords <- GetTissueCoordinates(seurat_obj@images[[1]])
          has_coordinates <- !is.null(coords) && nrow(coords) > 0
        } else {
          has_coordinates <- !is.null(seurat_obj@images[[1]]@coordinates)
        }
      }, error = function(e) {
        message("Could not validate coordinates, assuming they exist")
        has_coordinates <- TRUE
      })
    }
    
    return(list(
      valid = has_spatial,
      has_images = has_images,
      has_coordinates = has_coordinates,
      message = if (!has_spatial) "Missing Spatial assay" else "Valid spatial data"
    ))
  }
  
  # Safely extract ZIP file
  safe_extract_zip <- function(zip_path, extract_dir) {
    if (!file.exists(zip_path)) {
      stop("ZIP file does not exist at the specified path")
    }
    
    file_size <- file.info(zip_path)$size
    if (is.na(file_size) || file_size == 0) {
      stop("ZIP file appears to be empty or corrupted")
    }
    
    tryCatch({
      zip_contents <- unzip(zip_path, list = TRUE)
      if (nrow(zip_contents) == 0) {
        stop("ZIP file contains no files")
      }
      
      extracted_files <- unzip(zip_path, exdir = extract_dir)
      return(list(success = TRUE, contents = zip_contents, extracted_files = extracted_files))
      
    }, error = function(e) {
      stop(paste("Failed to process ZIP file:", e$message))
    })
  }
  
  # Organize spatial files into expected structure
  organize_spatial_files <- function(data_dir) {
    all_files <- list.files(data_dir, full.names = TRUE, recursive = FALSE)
    file_names <- basename(all_files)
    
    spatial_dir <- file.path(data_dir, "spatial")
    if (!dir.exists(spatial_dir)) {
      dir.create(spatial_dir, recursive = TRUE)
    }
    
    # Define file patterns for different spatial files
    file_patterns <- list(
      expression_h5 = list(
        patterns = c("filtered_feature_bc_matrix\\.h5", ".*feature.*matrix.*\\.h5", ".*\\.h5"),
        destination = "filtered_feature_bc_matrix.h5",
        target_dir = data_dir
      ),
      scalefactors = list(
        patterns = c("scalefactors_json.*\\.json", ".*scalefactors.*\\.json", ".*\\.json"),
        destination = "scalefactors_json.json",
        target_dir = spatial_dir
      ),
      tissue_positions_csv = list(
        patterns = c("tissue_positions.*\\.csv", ".*positions.*\\.csv", ".*tissue.*\\.csv"),
        destination = "tissue_positions_list.csv",
        target_dir = spatial_dir
      ),
      tissue_positions_parquet = list(
        patterns = c("tissue_positions.*\\.parquet", ".*positions.*\\.parquet"),
        destination = "tissue_positions.parquet",
        target_dir = spatial_dir
      ),
      hires_image = list(
        patterns = c("tissue_hires_image.*\\.png", ".*hires.*\\.png", ".*high.*\\.png"),
        destination = "tissue_hires_image.png",
        target_dir = spatial_dir
      ),
      lowres_image = list(
        patterns = c("tissue_lowres_image.*\\.png", ".*lowres.*\\.png", ".*low.*\\.png"),
        destination = "tissue_lowres_image.png",
        target_dir = spatial_dir
      )
    )
    
    # Match and move files
    moved_files <- list()
    for (file_type in names(file_patterns)) {
      pattern_info <- file_patterns[[file_type]]
      
      for (pattern in pattern_info$patterns) {
        matching_files <- grep(pattern, file_names, ignore.case = TRUE, value = TRUE)
        
        if (length(matching_files) > 0) {
          source_file <- file.path(data_dir, matching_files[1])
          dest_file <- file.path(pattern_info$target_dir, pattern_info$destination)
          
          if (file.copy(source_file, dest_file, overwrite = TRUE)) {
            moved_files[[file_type]] <- dest_file
            message(paste("Organized:", matching_files[1], "->", pattern_info$destination))
            break
          }
        }
      }
    }
    
    return(moved_files)
  }
  
  # Detect and organize spatial data structure
  detect_spatial_data_structure_flexible <- function(data_dir) {
    # First, try to organize files
    organized_files <- organize_spatial_files(data_dir)
    
    # Check for spatial directory
    spatial_dir <- file.path(data_dir, "spatial")
    if (!dir.exists(spatial_dir)) {
      return(list(valid = FALSE, message = "Could not create or find spatial directory"))
    }
    
    # Check required files
    required_files <- c(
      expression = "filtered_feature_bc_matrix.h5",
      scalefactors = file.path("spatial", "scalefactors_json.json")
    )
    
    missing_files <- c()
    found_files <- c()
    
    for (file_name in required_files) {
      file_path <- file.path(data_dir, file_name)
      if (file.exists(file_path)) {
        found_files <- c(found_files, file_name)
      } else {
        missing_files <- c(missing_files, file_name)
      }
    }
    
    # Check for tissue positions (either format)
    tissue_pos_parquet <- file.path(spatial_dir, "tissue_positions.parquet")
    tissue_pos_csv <- file.path(spatial_dir, "tissue_positions_list.csv")
    
    tissue_positions_found <- FALSE
    tissue_positions_format <- "none"
    
    if (file.exists(tissue_pos_parquet)) {
      tissue_positions_found <- TRUE
      tissue_positions_format <- "parquet"
    } else if (file.exists(tissue_pos_csv)) {
      tissue_positions_found <- TRUE
      tissue_positions_format <- "csv"
    }
    
    if (!tissue_positions_found) {
      missing_files <- c(missing_files, "tissue_positions file")
    }
    
    # Check for images (optional)
    has_images <- file.exists(file.path(spatial_dir, "tissue_hires_image.png")) ||
      file.exists(file.path(spatial_dir, "tissue_lowres_image.png"))
    
    # Determine if structure is valid
    valid <- length(missing_files) == 0
    
    if (valid) {
      message <- paste("Valid spatial data structure created. Found:", paste(found_files, collapse = ", "))
    } else {
      message <- paste("Missing required files:", paste(missing_files, collapse = ", "))
    }
    
    return(list(
      valid = valid,
      format = "H5",
      tissue_positions_format = tissue_positions_format,
      has_images = has_images,
      organized_files = organized_files,
      message = message
    ))
  }
  
  ############################## Load Spatial Dataset ##############################
  
  # Enhanced spatial data loading function
  observeEvent(input$spatial_file, {
    req(input$spatial_file)
    
    message("=== SPATIAL FILE UPLOAD STARTED ===")
    
    tryCatch({
      showModal(modalDialog(
        title = "Loading Spatial Data",
        "Processing your spatial transcriptomics data...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      file_info <- input$spatial_file
      temp_dir <- file.path(tempdir(), paste0("spatial_", format(Sys.time(), "%Y%m%d_%H%M%S")))
      dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
      
      # Extract ZIP
      message("Extracting ZIP file...")
      extracted_files <- unzip(file_info$datapath, exdir = temp_dir)
      message(paste("Extracted files:", paste(basename(extracted_files), collapse = ", ")))
      
      # Get all files recursively to handle nested structures
      all_files <- list.files(temp_dir, full.names = TRUE, recursive = TRUE)
      all_basenames <- basename(all_files)
      
      message(paste("All extracted files:", paste(all_basenames, collapse = ", ")))
      
      # Enhanced file detection with multiple patterns
      h5_files <- grep("(filtered_feature_bc_matrix|matrix|expression).*\\.h5$", all_files, value = TRUE, ignore.case = TRUE)
      json_files <- grep("scalefactors.*\\.json$", all_files, value = TRUE, ignore.case = TRUE)
      csv_files <- grep("tissue_positions.*\\.csv$", all_files, value = TRUE, ignore.case = TRUE)
      parquet_files <- grep("tissue_positions.*\\.parquet$", all_files, value = TRUE, ignore.case = TRUE)
      hires_images <- grep("tissue_hires_image.*\\.png$", all_files, value = TRUE, ignore.case = TRUE)
      lowres_images <- grep("tissue_lowres_image.*\\.png$", all_files, value = TRUE, ignore.case = TRUE)
      
      # Alternative patterns for different resolutions
      if (length(h5_files) == 0) {
        h5_files <- grep("\\.h5$", all_files, value = TRUE)
        message("Using broader H5 file detection")
      }
      
      if (length(json_files) == 0) {
        json_files <- grep("\\.json$", all_files, value = TRUE)
        message("Using broader JSON file detection")
      }
      
      if (length(csv_files) == 0) {
        csv_files <- grep("positions.*\\.csv$", all_files, value = TRUE, ignore.case = TRUE)
        message("Using broader CSV file detection")
      }
      
      if (length(parquet_files) == 0) {
        parquet_files <- grep("positions.*\\.parquet$", all_files, value = TRUE, ignore.case = TRUE)
        message("Using broader Parquet file detection")
      }
      
      message(paste("Enhanced detection - H5:", length(h5_files), "JSON:", length(json_files), 
                    "CSV:", length(csv_files), "Parquet:", length(parquet_files)))
      
      if (length(h5_files) == 0) {
        message("Available files:", paste(all_basenames, collapse = ", "))
        stop("No H5 matrix files found in uploaded data. Please check file structure.")
      }
      
      # Create spatial directory structure
      spatial_dir <- file.path(temp_dir, "spatial")
      if (!dir.exists(spatial_dir)) {
        dir.create(spatial_dir, recursive = TRUE)
      }
      
      # Copy H5 file to root with standard name
      h5_target <- file.path(temp_dir, "filtered_feature_bc_matrix.h5")
      if (!file.exists(h5_target)) {
        file.copy(h5_files[1], h5_target, overwrite = TRUE)
        message(paste("Copied H5 file:", basename(h5_files[1]), "to root"))
      }
      
      # Handle JSON file
      if (length(json_files) > 0) {
        json_target <- file.path(spatial_dir, "scalefactors_json.json")
        if (!file.exists(json_target)) {
          file.copy(json_files[1], json_target, overwrite = TRUE)
          message(paste("Copied JSON file:", basename(json_files[1])))
        }
      } else {
        # Create minimal scalefactors if missing (common issue with some datasets)
        json_target <- file.path(spatial_dir, "scalefactors_json.json")
        minimal_scalefactors <- list(
          tissue_hires_scalef = 1.0,
          tissue_lowres_scalef = 1.0,
          fiducial_diameter_fullres = 1.0,
          spot_diameter_fullres = 1.0
        )
        write(jsonlite::toJSON(minimal_scalefactors, auto_unbox = TRUE), file = json_target)
        message("Created minimal scalefactors file")
      }
      
      # Handle tissue positions - prefer Parquet over CSV
      positions_copied <- FALSE
      if (length(parquet_files) > 0) {
        parquet_target <- file.path(spatial_dir, "tissue_positions.parquet")
        if (!file.exists(parquet_target)) {
          file.copy(parquet_files[1], parquet_target, overwrite = TRUE)
          message(paste("Copied Parquet positions file:", basename(parquet_files[1])))
          positions_copied <- TRUE
        }
      }
      
      if (!positions_copied && length(csv_files) > 0) {
        csv_target <- file.path(spatial_dir, "tissue_positions_list.csv")
        if (!file.exists(csv_target)) {
          file.copy(csv_files[1], csv_target, overwrite = TRUE)
          message(paste("Copied CSV positions file:", basename(csv_files[1])))
          positions_copied <- TRUE
        }
      }
      
      if (!positions_copied) {
        message("Warning: No tissue positions file found. This may cause issues with spatial visualization.")
      }
      
      # Handle image files
      if (length(hires_images) > 0) {
        hires_target <- file.path(spatial_dir, "tissue_hires_image.png")
        if (!file.exists(hires_target)) {
          file.copy(hires_images[1], hires_target, overwrite = TRUE)
          message("Copied hires image")
        }
      }
      
      if (length(lowres_images) > 0) {
        lowres_target <- file.path(spatial_dir, "tissue_lowres_image.png")
        if (!file.exists(lowres_target)) {
          file.copy(lowres_images[1], lowres_target, overwrite = TRUE)
          message("Copied lowres image")
        }
      }
      
      # Enhanced loading with error handling for different formats
      message("Loading spatial data...")
      
      # Try different loading approaches
      spatial_data <- NULL
      loading_success <- FALSE
      
      # Approach 1: Standard Load10X_Spatial
      if (!loading_success) {
        tryCatch({
          message("Attempting standard Load10X_Spatial...")
          spatial_data <- Load10X_Spatial(
            data.dir = temp_dir,
            filename = "filtered_feature_bc_matrix.h5",
            assay = "Spatial",
            slice = "slice1",
            filter.matrix = TRUE,
            to.upper = FALSE
          )
          loading_success <- TRUE
          message("Standard loading successful")
        }, error = function(e) {
          message(paste("Standard loading failed:", e$message))
        })
      }
      
      # Approach 2: Try without filtering
      if (!loading_success) {
        tryCatch({
          message("Attempting Load10X_Spatial without filtering...")
          spatial_data <- Load10X_Spatial(
            data.dir = temp_dir,
            filename = "filtered_feature_bc_matrix.h5",
            assay = "Spatial",
            slice = "slice1",
            filter.matrix = FALSE,
            to.upper = FALSE
          )
          loading_success <- TRUE
          message("Loading without filtering successful")
        }, error = function(e) {
          message(paste("Loading without filtering failed:", e$message))
        })
      }
      
      # Approach 3: Manual loading for problematic datasets
      if (!loading_success) {
        tryCatch({
          message("Attempting manual H5 loading...")
          # Load expression matrix
          h5_data <- Read10X_h5(h5_target)
          
          # Create Seurat object
          spatial_data <- CreateSeuratObject(
            counts = h5_data,
            assay = "Spatial",
            project = "SpatialData"
          )
          
          # Try to add spatial information if available
          if (positions_copied) {
            # This is a simplified approach - may need adjustment for your specific data
            message("Adding spatial coordinates...")
            # Additional spatial coordinate loading would go here
          }
          
          loading_success <- TRUE
          message("Manual loading successful")
        }, error = function(e) {
          message(paste("Manual loading failed:", e$message))
        })
      }
      
      if (!loading_success || is.null(spatial_data)) {
        stop("All loading methods failed. Please check your data format and structure.")
      }
      
      message(paste("Spatial data loaded successfully:", ncol(spatial_data), "spots,", nrow(spatial_data), "genes"))
      
      # Check if we have actual data
      if (ncol(spatial_data) == 0) {
        stop("No cells/spots found in the loaded data. This may indicate a format issue or empty dataset.")
      }
      
      # Add mitochondrial genes based on species
      if (input$spatial_species_choice == "mouse") {
        spatial_data[["percent.mt"]] <- PercentageFeatureSet(spatial_data, pattern = "^mt-|^Mt-")
      } else if (input$spatial_species_choice == "human") {
        spatial_data[["percent.mt"]] <- PercentageFeatureSet(spatial_data, pattern = "^MT-")
      } else if (input$spatial_species_choice == "rat") {
        spatial_data[["percent.mt"]] <- PercentageFeatureSet(spatial_data, pattern = "^Mt-")
      }
      
      # Update object
      update_spatial_object(spatial_data, list(
        qc_applied = FALSE, 
        normalized = FALSE, 
        pca_computed = FALSE, 
        clustered = FALSE
      ))
      
      # Success message
      n_spots <- ncol(spatial_data)
      n_genes <- nrow(spatial_data)
      
      showNotification(
        paste0("Spatial data loaded successfully! Spots: ", format(n_spots, big.mark = ","), 
               ", Genes: ", format(n_genes, big.mark = ",")), 
        type = "message", 
        duration = 8
      )
      removeModal()
      
      message("=== SPATIAL LOADING COMPLETED ===")
      
    }, error = function(e) {
      message(paste("=== SPATIAL LOADING ERROR ===", e$message))
      removeModal()
      showNotification(
        paste("Error loading spatial data:", e$message, 
              "\nPlease check your file structure and ensure all required files are included."), 
        type = "error",
        duration = 15
      )
    })
  })
  
  ############################## Dataset Information Display ##############################
  
  # Info boxes for dataset overview
  output$spatial_spots_count <- renderInfoBox({
    count <- if (!is.null(spatial_obj())) ncol(spatial_obj()) else 0
    infoBox(
      title = "Spots/Cells",
      value = format(count, big.mark = ","),
      icon = icon("circle"),
      color = "blue"
    )
  })
  
  output$spatial_genes_count <- renderInfoBox({
    count <- if (!is.null(spatial_obj())) nrow(spatial_obj()) else 0
    infoBox(
      title = "Genes",
      value = format(count, big.mark = ","),
      icon = icon("dna"),
      color = "green"
    )
  })
  
  output$spatial_samples_count <- renderInfoBox({
    count <- if (!is.null(spatial_obj())) {
      length(unique(spatial_obj()@meta.data$orig.ident))
    } else 0
    infoBox(
      title = "Samples",
      value = count,
      icon = icon("flask"),
      color = "orange"
    )
  })
  
  # Dataset summary table
  output$spatial_dataset_summary <- renderDT({
    req(spatial_obj())
    
    obj <- spatial_obj()
    
    summary_data <- data.frame(
      Metric = c(
        "Total Spots/Cells",
        "Total Genes", 
        "Available Assays",
        "Sketched Data",
        "Median Genes per Spot",
        "Median UMI per Spot", 
        "Mean Mitochondrial %",
        "Has Tissue Image",
        "Processing Status"
      ),
      Value = c(
        format(ncol(obj), big.mark = ","),
        format(nrow(obj), big.mark = ","),
        paste(names(obj@assays), collapse = ", "),
        if (sketching_applied()) "Yes" else "No",
        round(median(obj$nFeature_Spatial), 0),
        format(round(median(obj$nCount_Spatial), 0), big.mark = ","),
        paste0(round(mean(obj$percent.mt, na.rm = TRUE), 2), "%"),
        if (length(obj@images) > 0) "Yes" else "No",
        get_processing_status()
      ),
      stringsAsFactors = FALSE
    )
    
    datatable(summary_data, options = list(pageLength = 10, searching = FALSE, lengthChange = FALSE, info = FALSE), rownames = FALSE)
  })
  
  
  # Load processed Seurat object from RDS file
  observeEvent(input$load_spatial_seurat_file, {
    req(input$load_spatial_seurat_file)
    
    message("=== LOADING SPATIAL SEURAT OBJECT FROM RDS ===")
    
    tryCatch({
      showModal(modalDialog(
        title = "Loading Spatial Object",
        "Reading your saved spatial Seurat object...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      # Read the RDS file
      file_path <- input$load_spatial_seurat_file$datapath
      message(paste("Loading file from:", file_path))
      
      # Load the object
      loaded_obj <- readRDS(file_path)
      message("RDS file loaded successfully")
      
      # Validate that it's a Seurat object
      if (!inherits(loaded_obj, "Seurat")) {
        stop("The loaded file is not a valid Seurat object")
      }
      
      # Check if it has spatial data
      has_spatial <- FALSE
      if ("Spatial" %in% names(loaded_obj@assays)) {
        has_spatial <- TRUE
        message("Found Spatial assay")
      }
      
      # Check for images
      has_images <- length(loaded_obj@images) > 0
      if (has_images) {
        message(paste("Found", length(loaded_obj@images), "spatial image(s)"))
      } else {
        message("No spatial images found")
      }
      
      # Get basic info
      n_spots <- ncol(loaded_obj)
      n_genes <- nrow(loaded_obj)
      
      message(paste("Object contains:", n_spots, "spots and", n_genes, "genes"))
      
      # Check what processing has been done
      has_pca <- "pca" %in% names(loaded_obj@reductions)
      has_umap <- "umap" %in% names(loaded_obj@reductions)
      has_clusters <- "seurat_clusters" %in% colnames(loaded_obj@meta.data)
      has_normalized <- FALSE
      
      # Check if data is normalized by looking at the assay data
      if ("Spatial" %in% names(loaded_obj@assays)) {
        assay_data <- GetAssayData(loaded_obj, assay = "Spatial", slot = "data")
        # If max value is much larger than raw counts would be, it's likely normalized
        max_val <- max(assay_data@x[1:min(1000, length(assay_data@x))])
        has_normalized <- max_val < 50  # Normalized data typically has smaller values
      }
      
      # Add mitochondrial percentage if not present
      if (!"percent.mt" %in% colnames(loaded_obj@meta.data)) {
        message("Adding mitochondrial percentage...")
        if (input$spatial_species_choice == "mouse") {
          loaded_obj[["percent.mt"]] <- PercentageFeatureSet(loaded_obj, pattern = "^mt-|^Mt-")
        } else if (input$spatial_species_choice == "human") {
          loaded_obj[["percent.mt"]] <- PercentageFeatureSet(loaded_obj, pattern = "^MT-")
        } else if (input$spatial_species_choice == "rat") {
          loaded_obj[["percent.mt"]] <- PercentageFeatureSet(loaded_obj, pattern = "^Mt-")
        }
      }
      
      # Update the spatial object and processing states
      update_spatial_object(loaded_obj, list(
        qc_applied = has_normalized,  # Assume QC was done if normalized
        normalized = has_normalized,
        pca_computed = has_pca,
        clustered = has_clusters
      ))
      
      # Update sketching status if sketch assay exists
      if ("sketch" %in% names(loaded_obj@assays)) {
        sketching_applied(TRUE)
        message("Found sketch assay")
      }
      
      # Update available assays
      assay_names <- names(loaded_obj@assays)
      available_assays(assay_names)
      updateSelectInput(session, "spatial_assay_select", choices = assay_names, selected = DefaultAssay(loaded_obj))
      updateSelectInput(session, "spatial_viz_assay", choices = assay_names, selected = DefaultAssay(loaded_obj))
      
      # Success message
      info_parts <- c()
      if (has_spatial) info_parts <- c(info_parts, "Spatial data")
      if (has_images) info_parts <- c(info_parts, "Images")
      if (has_normalized) info_parts <- c(info_parts, "Normalized")
      if (has_pca) info_parts <- c(info_parts, "PCA")
      if (has_umap) info_parts <- c(info_parts, "UMAP")
      if (has_clusters) info_parts <- c(info_parts, "Clusters")
      
      success_msg <- paste0(
        "Spatial object loaded successfully!\n",
        "Spots: ", format(n_spots, big.mark = ","), "\n",
        "Genes: ", format(n_genes, big.mark = ","), "\n",
        if (length(info_parts) > 0) paste("Contains:", paste(info_parts, collapse = ", ")) else ""
      )
      
      removeModal()
      
      message("=== SPATIAL OBJECT LOADING COMPLETED ===")
      
    }, error = function(e) {
      message(paste("=== SPATIAL OBJECT LOADING ERROR ===", e$message))
      removeModal()
      showNotification(
        paste("Error loading spatial object:", e$message), 
        type = "error",
        duration = 10
      )
    })
  })
  
  
  
  ############################## Spatial QC & Normalization ##############################
  
  # Info box showing current spot count
  output$spatial_spots_info <- renderInfoBox({
    count <- if (!is.null(spatial_obj())) ncol(spatial_obj()) else 0
    infoBox(
      title = "Current Spots",
      value = format(count, big.mark = ","),
      icon = icon("circle"),
      color = "blue",
      width = 12
    )
  })
  
  # Generate QC plots
  # Generate QC plots - FIXED VERSION with proper layer handling
  observeEvent(input$spatial_qc_plots, {
    req(spatial_obj())
    
    tryCatch({
      obj <- spatial_obj()
      options(warn = -1)  # Suppress warnings for plotting
      
      # Check available layers in the Spatial assay
      spatial_assay <- obj[["Spatial"]]
      available_layers <- names(spatial_assay@layers)
      message(paste("Available layers in Spatial assay:", paste(available_layers, collapse = ", ")))
      
      # Use counts layer explicitly for QC metrics
      if ("counts" %in% available_layers) {
        message("Using 'counts' layer for QC calculations")
      }
      
      # Spatial tissue plot with proper layer specification
      output$spatial_tissue_plot <- renderPlot({
        if (length(obj@images) > 0) {
          tryCatch({
            # Force raster=FALSE to avoid the rasterization warning for QC
            SpatialFeaturePlot(obj, features = "nCount_Spatial", 
                               pt.size.factor = 1.2, 
                               alpha = c(0.1, 1),
                               raster = FALSE) +  # Disable rasterization for QC
              theme(legend.position = "bottom") +
              ggtitle("UMI Count Distribution")
          }, error = function(e) {
            message(paste("Spatial plot error:", e$message))
            # Fallback plot if spatial plot fails
            ggplot(obj@meta.data, aes(x = nCount_Spatial, y = nFeature_Spatial)) +
              geom_point(alpha = 0.6, size = 0.5) +
              labs(title = "Spatial QC Metrics", x = "UMI Count", y = "Gene Count") +
              theme_minimal()
          })
        } else {
          ggplot(obj@meta.data, aes(x = nCount_Spatial, y = nFeature_Spatial)) +
            geom_point(alpha = 0.6, size = 0.5) +
            labs(title = "Spatial QC Metrics", x = "UMI Count", y = "Gene Count") +
            theme_minimal()
        }
      })
      
      # Standard QC violin plots with raster control
      output$spatial_vlnplot <- renderPlot({
        VlnPlot(obj, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"),
                ncol = 3, pt.size = 0, raster = FALSE) +  # Disable raster for QC
          theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
      })
      
      output$spatial_qc_violin <- renderPlot({
        VlnPlot(obj, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"),
                ncol = 1, pt.size = 0, raster = FALSE) &  # Disable raster for QC
          theme(axis.text.x = element_text(angle = 45, hjust = 1))
      })
      
      # Feature scatter plots with raster control for large datasets
      output$spatial_feature_scatter <- renderPlot({
        n_spots <- ncol(obj)
        use_raster <- n_spots > 50000  # Only rasterize for very large datasets
        
        plot1 <- FeatureScatter(obj, feature1 = "nCount_Spatial", feature2 = "percent.mt", 
                                pt.size = 0.5, raster = use_raster) + NoLegend()
        plot2 <- FeatureScatter(obj, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial", 
                                pt.size = 0.5, raster = use_raster) + NoLegend()
        plot1 + plot2
      })
      
      output$spatial_scatter1 <- renderPlot({
        n_spots <- ncol(obj)
        use_raster <- n_spots > 50000
        FeatureScatter(obj, feature1 = "nCount_Spatial", feature2 = "percent.mt", 
                       pt.size = 0.5, raster = use_raster)
      })
      
      output$spatial_scatter2 <- renderPlot({
        n_spots <- ncol(obj)
        use_raster <- n_spots > 50000
        FeatureScatter(obj, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial", 
                       pt.size = 0.5, raster = use_raster)
      })
      
      options(warn = 0)  # Restore warnings
      showNotification("QC plots generated successfully!", type = "message")
      
    }, error = function(e) {
      options(warn = 0)
      showNotification(paste("Error generating QC plots:", e$message), type = "error")
    })
  })
  
  # Analyze data distribution to suggest optimal filters
  analyze_data_distribution <- function(obj) {
    features <- obj$nFeature_Spatial
    counts <- obj$nCount_Spatial
    mt <- obj$percent.mt
    
    # Clean NA values
    features[is.na(features)] <- 0
    counts[is.na(counts)] <- 0
    mt[is.na(mt)] <- 0
    
    # Calculate suggested filter thresholds
    suggested_filters <- list(
      min_features = max(1, round(quantile(features[features > 0], 0.25, na.rm = TRUE) * 0.5)),
      max_features = min(round(quantile(features[features > 0], 0.95, na.rm = TRUE)), 1000),
      min_counts = max(1, round(quantile(counts[counts > 0], 0.25, na.rm = TRUE) * 0.5)),
      max_counts = min(round(quantile(counts[counts > 0], 0.95, na.rm = TRUE)), 5000),
      max_mt = min(round(quantile(mt, 0.95, na.rm = TRUE)), 50)
    )
    
    # Estimate how many spots would pass these filters
    would_pass <- features >= suggested_filters$min_features & 
      features <= suggested_filters$max_features &
      counts >= suggested_filters$min_counts &
      counts <= suggested_filters$max_counts &
      mt <= suggested_filters$max_mt
    
    spots_passing <- sum(would_pass, na.rm = TRUE)
    
    return(list(suggested_filters = suggested_filters, spots_passing = spots_passing))
  }
  
  # Apply QC filters to spatial data
  observeEvent(input$apply_spatial_qc, {
    req(spatial_obj())
    
    tryCatch({
      showModal(modalDialog(
        title = "Applying QC Filters",
        "Analyzing and filtering spatial data...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      obj <- spatial_obj()
      original_spots <- ncol(obj)
      
      # Analyze data distribution for adaptive filtering
      distribution_analysis <- analyze_data_distribution(obj)
      
      # Remove completely empty spots first
      empty_spots <- obj$nFeature_Spatial == 0 & obj$nCount_Spatial == 0
      empty_spots[is.na(empty_spots)] <- TRUE
      
      if (sum(!empty_spots) > 0) {
        obj <- obj[, !empty_spots]
      }
      
      # Determine which filters to use
      user_filters_would_pass <- sum(
        obj$nFeature_Spatial >= input$spatial_nFeature_range[1] & 
          obj$nFeature_Spatial <= input$spatial_nFeature_range[2] &
          obj$nCount_Spatial >= input$spatial_nCount_range[1] &
          obj$nCount_Spatial <= input$spatial_nCount_range[2] &
          obj$percent.mt <= input$spatial_mt_max,
        na.rm = TRUE
      )
      
      # Use adaptive filters if user filters would remove too many spots
      if (user_filters_would_pass < 1000 && distribution_analysis$spots_passing > user_filters_would_pass) {
        filters <- distribution_analysis$suggested_filters
        filter_type <- "adaptive"
      } else {
        filters <- list(
          min_features = input$spatial_nFeature_range[1],
          max_features = input$spatial_nFeature_range[2],
          min_counts = input$spatial_nCount_range[1],
          max_counts = input$spatial_nCount_range[2],
          max_mt = input$spatial_mt_max
        )
        filter_type <- "user-defined"
      }
      
      # Apply filters
      keep_cells <- obj$nFeature_Spatial >= filters$min_features & 
        obj$nFeature_Spatial <= filters$max_features &
        obj$nCount_Spatial >= filters$min_counts &
        obj$nCount_Spatial <= filters$max_counts &
        obj$percent.mt <= filters$max_mt
      
      keep_cells[is.na(keep_cells)] <- FALSE
      
      # Fallback to minimal filters if nothing passes
      if (sum(keep_cells) == 0) {
        keep_cells <- obj$nFeature_Spatial > 0 & obj$nCount_Spatial > 0 & obj$percent.mt < 90
        keep_cells[is.na(keep_cells)] <- FALSE
        filter_type <- "minimal"
      }
      
      obj <- obj[, keep_cells]
      
      update_spatial_object(obj, list(qc_applied = TRUE))
      
      filtered_spots <- ncol(obj)
      removal_rate <- round((1 - filtered_spots/original_spots) * 100, 1)
      
      success_message <- paste0(
        "QC completed! Filter type: ", filter_type, "\n",
        "Remaining: ", format(filtered_spots, big.mark = ","), " spots (", 100-removal_rate, "%)"
      )
      
      showNotification(success_message, type = "message", duration = 10)
      removeModal()
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("QC failed:", e$message), type = "error")
    })
  })
  
  # Normalize spatial data WITH DEBUG LOGS
  observeEvent(input$normalize_spatial_data, {
    message("=== NORMALIZE BUTTON CLICKED ===")
    
    if (is.null(spatial_obj())) {
      message("ERROR: spatial_obj() is NULL")
      showNotification("No spatial data loaded!", type = "error")
      return()
    }
    
    message("spatial_obj() exists, proceeding...")
    
    tryCatch({
      obj <- spatial_obj()
      message("Successfully retrieved spatial object")
      
      n_spots <- ncol(obj)
      message(paste("Dataset size:", n_spots, "spots"))
      
      showModal(modalDialog(
        title = "Normalizing Data",
        div(
          paste("Normalizing", format(n_spots, big.mark = ","), "spots..."),
          br(),
          tags$small("This may take several minutes")
        ),
        easyClose = FALSE,
        footer = NULL
      ))
      
      message("Modal dialog shown")
      message("=== STARTING NORMALIZATION ===")
      
      # Force garbage collection
      gc()
      
      # Suppress warnings
      options(warn = -1)
      
      normalization_success <- FALSE
      method_used <- "none"
      
      # Choose normalization strategy based on dataset size
      if (n_spots > 10000) {
        message("Large dataset - using memory-optimized normalization")
        
        tryCatch({
          message("Starting NormalizeData...")
          obj <- NormalizeData(obj, assay = "Spatial", 
                               normalization.method = "LogNormalize",
                               scale.factor = 10000, 
                               verbose = FALSE)
          message("NormalizeData completed")
          
          message("Starting FindVariableFeatures...")
          obj <- FindVariableFeatures(obj, assay = "Spatial",
                                      selection.method = "vst",
                                      nfeatures = min(1000, input$spatial_var_features),
                                      verbose = FALSE)
          message("FindVariableFeatures completed")
          
          normalization_success <- TRUE
          method_used <- "Standard normalization"
          message("Standard normalization SUCCESS")
        })
      }
      
      options(warn = 0)
      
      if (!normalization_success) {
        message("ERROR: All normalization methods failed!")
        stop("All normalization methods failed")
      }
      
      message("Generating variable features plot...")
      # Generate variable features plot
      output$spatial_variable_features <- renderPlot({
        tryCatch({
          var_features <- VariableFeatures(obj)
          if (length(var_features) > 0) {
            assay_to_use <- if ("SCT" %in% names(obj@assays)) "SCT" else "Spatial"
            VariableFeaturePlot(obj, assay = assay_to_use) + 
              ggtitle(paste("Variable Features -", method_used))
          } else {
            ggplot() + geom_text(aes(x=0.5, y=0.5, label="No variable features"), size=6) + theme_void()
          }
        }, error = function(e) {
          ggplot() + geom_text(aes(x=0.5, y=0.5, label="Normalization completed"), size=6) + theme_void()
        })
      })
      
      message("Updating spatial object...")
      # Update object and state
      update_spatial_object(obj, list(normalized = TRUE))
      message("Spatial object updated")
      
      # Force garbage collection after normalization
      gc()
      
      message(paste("=== NORMALIZATION COMPLETED ===", method_used))
      showNotification(paste("Normalization completed using:", method_used), 
                       type = "message", duration = 8)
      removeModal()
      
    }, error = function(e) {
      message(paste("=== NORMALIZATION ERROR ===", e$message))
      options(warn = 0)
      removeModal()
      gc()
      
      showNotification(paste("Normalization failed:", e$message, 
                             "\nTry applying stricter QC filters to reduce dataset size"), 
                       type = "error")
    })
  })
  
  
  # Auto-adjust sketch size based on dataset size
  observe({
    req(spatial_obj())
    n_spots <- ncol(spatial_obj())
    
    # Suggest reasonable sketch sizes based on dataset size
    suggested_sketch <- if (n_spots > 500000) {
      100000  # Large datasets: sketch to 100k
    } else if (n_spots > 200000) {
      50000   # Medium-large: sketch to 50k
    } else if (n_spots > 100000) {
      30000   # Medium: sketch to 30k
    } else {
      min(20000, floor(n_spots * 0.5))  # Small: sketch to 50% or 20k max
    }
    
    # Update the input only if it hasn't been manually changed
    current_value <- input$sketch_ncells %||% 50000
    if (current_value == 50000) {  # Only auto-adjust if still at default
      updateNumericInput(session, "sketch_ncells", value = suggested_sketch)
    }
  })
  
  
  
  
  
  
  
  ########################Sketching########################
  # Apply sketching with improved error handling and method selection
  observeEvent(input$apply_sketching, {
    req(spatial_obj())
    
    if (!processing_states$normalized) {
      showNotification("Please normalize the data first!", type = "warning")
      return()
    }
    
    tryCatch({
      obj <- spatial_obj()
      n_cells_current <- ncol(obj)
      n_cells_sketch <- input$sketch_ncells
      
      if (n_cells_sketch >= n_cells_current) {
        showNotification("Sketch size should be smaller than current dataset size", type = "warning")
        return()
      }
      
      showModal(modalDialog(
        title = "Applying Sketching",
        div(
          paste("Sketching from", format(n_cells_current, big.mark = ","), 
                "to", format(n_cells_sketch, big.mark = ","), "cells..."),
          br(),
          tags$small("This will create a representative subset and normalize it")
        ),
        easyClose = FALSE,
        footer = NULL
      ))
      
      message("=== STARTING SKETCHING ===")
      
      # Try sketching with timeout and fallback method
      sketch_method <- input$sketch_method
      sketching_success <- FALSE
      
      # First, try the selected method with timeout protection
      if (sketch_method == "LeverageScore") {
        message("Attempting LeverageScore sketching...")
        
        # For very large datasets, LeverageScore can be too slow
        if (n_cells_current > 200000) {
          message("Dataset too large for LeverageScore, switching to Uniform")
          sketch_method <- "Uniform"
        } else {
          tryCatch({
            # Set a reasonable timeout for LeverageScore
            obj_sketched <- R.utils::withTimeout({
              SketchData(object = obj, ncells = n_cells_sketch, 
                         method = "LeverageScore", sketched.assay = "sketch")
            }, timeout = 300)  # 5 minute timeout
            
            obj <- obj_sketched
            sketching_success <- TRUE
            message("LeverageScore sketching completed successfully")
            
          }, error = function(e) {
            message(paste("LeverageScore failed:", e$message))
            if (grepl("timeout|slow", e$message, ignore.case = TRUE)) {
              message("LeverageScore too slow, falling back to Uniform method")
              sketch_method <- "Uniform"
            } else {
              stop(e$message)
            }
          })
        }
      }
      
      # If LeverageScore failed or Uniform was selected, use Uniform method
      if (!sketching_success) {
        message("Using Uniform sketching method...")
        tryCatch({
          obj <- SketchData(object = obj, ncells = n_cells_sketch, 
                            method = "Uniform", sketched.assay = "sketch")
          sketching_success <- TRUE
          sketch_method <- "Uniform"  # Update method name for reporting
          message("Uniform sketching completed successfully")
          
        }, error = function(e) {
          stop(paste("Both sketching methods failed:", e$message))
        })
      }
      
      if (!sketching_success) {
        stop("Sketching failed with all available methods")
      }
      
      message("Sketching completed, now normalizing sketch assay...")
      
      # Normalize and scale the sketch assay automatically
      obj <- NormalizeData(obj, assay = "sketch", normalization.method = "LogNormalize", 
                           scale.factor = 10000, verbose = FALSE)
      message("Sketch normalization completed")
      
      # Find variable features for sketch assay
      obj <- FindVariableFeatures(obj, assay = "sketch", selection.method = "vst", 
                                  nfeatures = min(2000, input$spatial_var_features), verbose = FALSE)
      message("Variable features found for sketch")
      
      # Scale the sketch assay data
      var_features_sketch <- VariableFeatures(obj, assay = "sketch")
      if (length(var_features_sketch) > 0) {
        obj <- ScaleData(obj, assay = "sketch", features = var_features_sketch, verbose = FALSE)
        message("Sketch scaling completed")
      }
      
      # IMPORTANT: Set sketch as the default assay
      DefaultAssay(obj) <- "sketch"
      
      # Update object and state
      update_spatial_object(obj, list(sketched = TRUE))
      sketching_applied(TRUE)
      selected_assay("sketch")  # Update our reactive value
      
      # Update available assays and set sketch as selected in ALL menus
      new_assays <- names(obj@assays)
      available_assays(new_assays)
      updateSelectInput(session, "spatial_analysis_assay", choices = new_assays, selected = "sketch")
      updateSelectInput(session, "spatial_viz_assay", choices = new_assays, selected = "sketch")
      updateSelectInput(session, "spatial_default_assay", choices = new_assays, selected = "sketch")
      updateSelectInput(session, "spatial_assay_select", choices = new_assays, selected = "sketch")
      
      success_message <- paste0(
        "Sketching completed successfully!\n",
        "Original: ", format(n_cells_current, big.mark = ","), " cells\n",
        "Sketched: ", format(n_cells_sketch, big.mark = ","), " cells\n",
        "Method used: ", sketch_method, "\n",
        "Default assay set to: sketch"
      )
      
      showNotification(success_message, type = "message", duration = 10)
      removeModal()
      
    }, error = function(e) {
      removeModal()
      message(paste("Sketching error:", e$message))
      showNotification(paste("Sketching failed:", e$message, 
                             "\nTry using Uniform method or reducing sketch size"), 
                       type = "error")
    })
  })
  # Apply assay selection
  observeEvent(input$apply_assay_selection, {
    req(spatial_obj())
    
    tryCatch({
      obj <- spatial_obj()
      chosen_assay <- input$spatial_analysis_assay
      
      if (chosen_assay %in% names(obj@assays)) {
        selected_assay(chosen_assay)
        showNotification(paste("Analysis assay set to:", chosen_assay), type = "message", duration = 5)
      } else {
        showNotification("Selected assay not found", type = "error")
      }
      
    }, error = function(e) {
      showNotification(paste("Error setting assay:", e$message), type = "error")
    })
  })
  
  # Update available assays when object changes
  observe({
    req(spatial_obj())
    obj <- spatial_obj()
    assays <- names(obj@assays)
    
    # Update assay selection menu
    updateSelectInput(session, "spatial_analysis_assay", choices = assays, selected = selected_assay())
  })
  
  
  # Display current assay information - FIXED VERSION
  output$spatial_assay_info <- renderText({
    req(spatial_obj())
    obj <- spatial_obj()
    
    assay_info <- paste("Available assays:", paste(names(obj@assays), collapse = ", "))
    
    if (sketching_applied()) {
      # Get the actual number of cells in the sketch assay
      sketch_ncells <- ncol(obj[["sketch"]])
      sketch_info <- paste("Sketched data available with", format(sketch_ncells, big.mark = ","), "cells")
      assay_info <- paste(assay_info, "\n", sketch_info)
    }
    
    selected_assay_name <- input$spatial_viz_assay
    if (!is.null(selected_assay_name) && selected_assay_name %in% names(obj@assays)) {
      # Get number of features and cells for the selected assay
      n_features <- nrow(obj[[selected_assay_name]])
      n_cells <- ncol(obj[[selected_assay_name]])
      current_info <- paste("Current assay:", selected_assay_name, "with", 
                            format(n_features, big.mark = ","), "features and",
                            format(n_cells, big.mark = ","), "cells")
      assay_info <- paste(assay_info, "\n", current_info)
    }
    
    return(assay_info)
  })
  
  
  
  
  
  # Keep these functions for dynamic recommendations
  output$sketching_recommendations <- renderText({
    req(spatial_obj())
    
    if (!processing_states$normalized) {
      "⚠️ Please normalize the data first before sketching"
    } else {
      n_spots <- ncol(spatial_obj())
      
      if (n_spots > 200000) {
        "⚠️ Very large dataset detected. Sketching highly recommended. Use 'Uniform' method for faster processing."
      } else if (n_spots > 100000) {
        "💡 Large dataset detected. Sketching recommended. 'LeverageScore' provides better representation but may be slower."
      } else if (n_spots > 50000) {
        "💡 Medium dataset. Sketching optional but can speed up downstream analysis."
      } else {
        "✅ Small dataset. Sketching not necessary but available if desired."
      }
    }
  })
  
  output$sketch_size_recommendation <- renderText({
    req(spatial_obj())
    n_spots <- ncol(spatial_obj())
    current_sketch_size <- input$sketch_ncells %||% 50000
    
    if (!processing_states$normalized) {
      "Sketching requires normalized data"
    } else if (current_sketch_size >= n_spots) {
      "⚠️ Sketch size should be smaller than total dataset size"
    } else {
      percentage <- round((current_sketch_size / n_spots) * 100, 1)
      paste("✓ Will sketch", percentage, "% of the data")
    }
  })
  
  # Set Default Assay function
  observeEvent(input$set_default_assay, {
    req(spatial_obj())
    
    tryCatch({
      obj <- spatial_obj()
      selected_assay <- input$spatial_default_assay
      
      if (selected_assay %in% names(obj@assays)) {
        # Set as default assay in the Seurat object
        DefaultAssay(obj) <- selected_assay
        
        # Update our reactive values
        spatial_obj(obj)
        selected_assay(selected_assay)
        
        showNotification(paste("Default assay set to:", selected_assay), type = "message", duration = 5)
        
        # Update all assay selection menus to reflect the change
        updateSelectInput(session, "spatial_viz_assay", selected = selected_assay)
        updateSelectInput(session, "spatial_assay_select", selected = selected_assay)
        
      } else {
        showNotification("Selected assay not found", type = "error")
      }
      
    }, error = function(e) {
      showNotification(paste("Error setting assay:", e$message), type = "error")
    })
  })
  
  # Update available assays when object changes
  observe({
    req(spatial_obj())
    obj <- spatial_obj()
    assays <- names(obj@assays)
    current_assay <- DefaultAssay(obj)
    
    # Update assay selection menu
    updateSelectInput(session, "spatial_default_assay", choices = assays, selected = current_assay)
    
    # Update current assay info
    output$current_assay_info <- renderText({
      paste("Current:", current_assay)
    })
  })
  
  output$spatial_current_assay_box <- renderInfoBox({
    req(spatial_obj())
    obj <- spatial_obj()
    
    current_assay <- selected_assay()
    if (current_assay %in% names(obj@assays)) {
      n_features <- nrow(obj[[current_assay]])
      color <- if (current_assay == "sketch") "purple" else if (current_assay == "SCT") "green" else "blue"
    } else {
      n_features <- 0
      color <- "red"
    }
    
    infoBox(
      title = "Selected Assay",
      value = current_assay,
      subtitle = paste(format(n_features, big.mark = ","), "features"),
      icon = icon("layer-group"),
      color = color,
      width = 12
    )
  })
  
  
  
  
  
  ############################## PCA & Dimensionality Reduction ##############################
  
  # Run PCA on spatial data
  observeEvent(input$spatial_scale_pca, {
    req(spatial_obj())
    
    if (!processing_states$normalized) {
      showNotification("Please normalize the data first!", type = "warning")
      return()
    }
    
    tryCatch({
      showModal(modalDialog(title = "Running PCA", "Computing principal components...", easyClose = FALSE, footer = NULL))
      
      obj <- spatial_obj()
      assay_for_pca <- selected_assay()
      
      # Determine scaling strategy
      if ("sketch" %in% names(obj@assays) && assay_for_pca == "sketch") {
        message("Using pre-scaled sketch assay for PCA")
      } else if ("SCT" %in% names(obj@assays) && assay_for_pca == "SCT") {
        message("Using SCT assay for PCA")
      } else {
        # Scale data for other assays
        var_features <- VariableFeatures(obj, assay = assay_for_pca)
        if (length(var_features) > 0) {
          message(paste("Scaling", assay_for_pca, "assay with", length(var_features), "variable features"))
          obj <- ScaleData(obj, assay = assay_for_pca, features = var_features, verbose = FALSE)
        } else {
          stop(paste("No variable features found for", assay_for_pca, "assay"))
        }
      }
      
      # Calculate optimal number of PCs  
      n_spots <- ncol(obj)
      n_features <- nrow(obj[[assay_for_pca]])
      max_pcs <- if (n_spots > 100000) 30 else if (n_spots > 50000) 40 else 50
      max_pcs <- min(max_pcs, n_spots - 1, n_features - 1)
      
      message(paste("Running PCA with", max_pcs, "components on", assay_for_pca, "assay"))
      obj <- RunPCA(obj, assay = assay_for_pca, npcs = max_pcs, verbose = FALSE)
      
      output$spatial_elbow_plot <- renderPlot({
        ElbowPlot(obj, ndims = max_pcs) + 
          ggtitle(paste("PCA Elbow Plot -", assay_for_pca, "assay")) +
          labs(subtitle = paste(n_spots, "spots,", max_pcs, "components"))
      })
      
      update_spatial_object(obj, list(pca_computed = TRUE))
      showNotification(paste("PCA completed using", assay_for_pca, "assay"), type = "message")
      removeModal()
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("PCA failed:", e$message), type = "error")
    })
  })
  
  
  ############################## Clustering ##############################
  
  # Combined Neighbors + UMAP function
  observeEvent(input$spatial_run_neighbors_umap, {
    message("=== NEIGHBORS + UMAP PIPELINE STARTED ===")
    
    if (is.null(spatial_obj())) {
      showNotification("No spatial data loaded!", type = "error")
      return()
    }
    
    if (!processing_states$pca_computed) {
      showNotification("Please run PCA first!", type = "warning")
      return()
    }
    
    tryCatch({
      obj <- spatial_obj()
      
      # Check PCA availability
      if (!"pca" %in% names(obj@reductions)) {
        stop("PCA reduction not found in object")
      }
      
      # Get parameters from UI inputs
      cluster_dims <- input$spatial_cluster_dims %||% 15
      umap_dims <- input$spatial_umap_dims %||% 30
      
      available_pcs <- ncol(obj@reductions$pca@cell.embeddings)
      max_dims_neighbors <- min(cluster_dims, available_pcs)
      max_dims_umap <- min(umap_dims, available_pcs)
      
      showModal(modalDialog(
        title = "Running Neighbors + UMAP",
        div(
          "Step 1: Computing neighbors...", br(),
          "Step 2: Computing UMAP...", br(),
          tags$small(paste("Neighbors:", max_dims_neighbors, "dims | UMAP:", max_dims_umap, "dims"))
        ),
        easyClose = FALSE,
        footer = NULL
      ))
      
      # Step 1: Find neighbors
      message(paste("Step 1: Finding neighbors with", max_dims_neighbors, "dimensions..."))
      obj <- FindNeighbors(obj, dims = 1:max_dims_neighbors, verbose = FALSE)
      message("FindNeighbors computation completed")
      
      # Step 2: Run UMAP
      message(paste("Step 2: Computing UMAP with", max_dims_umap, "dimensions..."))
      obj <- RunUMAP(obj, dims = 1:max_dims_umap, verbose = FALSE)
      message("RunUMAP computation completed")
      
      # Update object
      spatial_obj(obj)
      message("Spatial object updated with neighbors and UMAP")
      
      # Generate plots
      output$spatial_umap_clusters <- renderPlot({
        DimPlot(obj, reduction = "umap", raster = FALSE) + 
          ggtitle(paste("UMAP (", max_dims_umap, "dims) - Ready for clustering")) +
          theme(plot.title = element_text(hjust = 0.5))
      })
      
      output$spatial_tissue_clusters <- renderPlot({
        if (length(obj@images) > 0) {
          tryCatch({
            SpatialDimPlot(obj) + 
              ggtitle("Spatial View - Ready for clustering") +
              theme(plot.title = element_text(hjust = 0.5))
          }, error = function(e) {
            ggplot() + 
              geom_text(aes(x = 0.5, y = 0.5, label = "Spatial coordinates available"), size = 6) +
              theme_void()
          })
        } else {
          ggplot() + 
            geom_text(aes(x = 0.5, y = 0.5, label = "No spatial coordinates"), size = 6) +
            theme_void()
        }
      })
      
      success_message <- paste0(
        "Neighbors + UMAP completed!\n",
        "Neighbors: ✓ (", max_dims_neighbors, " dims)\n",
        "UMAP: ✓ (", max_dims_umap, " dims)\n",
        "Ready for clustering"
      )
      
      showNotification(success_message, type = "message", duration = 8)
      removeModal()
      
      message("=== NEIGHBORS + UMAP PIPELINE COMPLETED ===")
      
    }, error = function(e) {
      message(paste("=== NEIGHBORS + UMAP ERROR ===", e$message))
      removeModal()
      showNotification(paste("Neighbors + UMAP failed:", e$message), type = "error", duration = 10)
    })
  })
  
  # Find clusters function
  observeEvent(input$spatial_find_clusters, {
    message("=== FIND CLUSTERS BUTTON CLICKED ===")
    
    if (is.null(spatial_obj())) {
      message("ERROR: spatial_obj() is NULL")
      showNotification("No spatial data loaded!", type = "error")
      return()
    }
    
    obj <- spatial_obj()
    
    # Check if neighbors have been computed
    if (is.null(obj@graphs) || length(obj@graphs) == 0) {
      message("ERROR: No neighbor graph found")
      showNotification("Please run 'Run Neighbors + UMAP' first!", type = "warning")
      return()
    }
    
    message("spatial_obj() exists and neighbors computed, proceeding...")
    
    tryCatch({
      # Get and validate parameters
      resolution <- as.numeric(input$spatial_resolution)
      algorithm <- as.integer(input$spatial_cluster_algorithm)
      
      # Validate inputs
      if (is.na(resolution) || resolution <= 0) {
        resolution <- 0.5
        message("Invalid resolution, using default 0.5")
      }
      if (is.na(algorithm) || !algorithm %in% 1:4) {
        algorithm <- 1
        message("Invalid algorithm, using default Louvain")
      }
      
      message(paste("Using resolution:", resolution, "and algorithm:", algorithm))
      
      showModal(modalDialog(
        title = "Finding Clusters",
        div(
          paste("Clustering spots with resolution", resolution, "..."),
          br(),
          tags$small("This may take a few minutes for large datasets")
        ),
        easyClose = FALSE,
        footer = NULL
      ))
      
      message("=== STARTING FIND CLUSTERS ===")
      
      # Get the correct graph name
      graph_name <- names(obj@graphs)[1]
      message(paste("Using graph:", graph_name))
      
      # Find clusters with explicit parameters
      message("Starting FindClusters computation...")
      obj <- FindClusters(obj, 
                          resolution = resolution,
                          algorithm = algorithm,
                          graph.name = graph_name,
                          verbose = FALSE)
      message("FindClusters computation completed")
      
      # Ensure cluster identities are properly formatted
      if ("seurat_clusters" %in% colnames(obj@meta.data)) {
        obj@meta.data$seurat_clusters <- as.factor(obj@meta.data$seurat_clusters)
        message("Cluster identities converted to factors")
      }
      
      # Update our reactive object and state
      update_spatial_object(obj, list(clustered = TRUE))
      message("Spatial object updated with clusters")
      
      n_clusters <- length(unique(obj$seurat_clusters))
      message(paste("Found", n_clusters, "clusters"))
      
      success_message <- paste0(
        "Clustering completed successfully!\n",
        "Found ", n_clusters, " clusters\n",
        "Resolution: ", resolution, "\n",
        "Algorithm: ", c("Louvain", "Louvain (multilevel)", "SLM")[algorithm]
      )
      
      showNotification(success_message, type = "message", duration = 10)
      removeModal()
      
      message("=== FIND CLUSTERS COMPLETED ===")
      
    }, error = function(e) {
      message(paste("=== FIND CLUSTERS ERROR ===", e$message))
      removeModal()
      
      showNotification(paste("Clustering failed:", e$message), 
                       type = "error", duration = 12)
    })
  })
  
  # Input validation observer
  observe({
    # Validate spatial resolution input
    if (!is.null(input$spatial_resolution)) {
      if (is.na(input$spatial_resolution) || input$spatial_resolution <= 0) {
        updateNumericInput(session, "spatial_resolution", value = 0.5)
        showNotification("Invalid resolution value, reset to 0.5", type = "warning", duration = 3)
      }
    }
    
    # Validate clustering algorithm input
    if (!is.null(input$spatial_cluster_algorithm)) {
      if (is.na(input$spatial_cluster_algorithm) || !input$spatial_cluster_algorithm %in% 1:4) {
        updateSelectInput(session, "spatial_cluster_algorithm", selected = 1)
        showNotification("Invalid algorithm selection, reset to Louvain", type = "warning", duration = 3)
      }
    }
  })
  
  # Generate clustering plots after successful clustering
  observe({
    req(spatial_obj())
    
    if (processing_states$clustered) {
      obj <- spatial_obj()
      
      if ("seurat_clusters" %in% colnames(obj@meta.data)) {
        n_clusters <- length(unique(obj$seurat_clusters))
        
        # Generate UMAP clustering plot
        output$spatial_umap_clusters <- renderPlot({
          if ("umap" %in% names(obj@reductions)) {
            p <- DimPlot(obj, reduction = "umap", group.by = "seurat_clusters", 
                         label = input$spatial_plot_labels, pt.size = 0.8, 
                         label.size = 4, raster = FALSE)
            if (!input$spatial_plot_labels) p <- p + NoLegend()
            p + ggtitle(paste("UMAP Clusters (n =", n_clusters, ")")) +
              theme(plot.title = element_text(hjust = 0.5, size = 14)) +
              guides(color = guide_legend(override.aes = list(size = 3)))
          } else {
            ggplot() + 
              geom_text(aes(x = 0.5, y = 0.5, label = "UMAP not available\nRun UMAP first"), 
                        size = 6, color = "red") +
              theme_void()
          }
        })
        
        # Generate spatial clustering plot
        output$spatial_tissue_clusters <- renderPlot({
          if (length(obj@images) > 0) {
            tryCatch({
              SpatialDimPlot(obj, group.by = "seurat_clusters", 
                             label = input$spatial_plot_labels,
                             pt.size.factor = 1.2, alpha = c(0.3, 1),
                             label.size = 3) +
                ggtitle(paste("Spatial Clusters (n =", n_clusters, ")")) +
                theme(plot.title = element_text(hjust = 0.5, size = 14),
                      legend.text = element_text(size = 10)) +
                guides(fill = guide_legend(override.aes = list(size = 3)))
            }, error = function(e) {
              message(paste("Spatial plot error:", e$message))
              ggplot() + 
                geom_text(aes(x = 0.5, y = 0.5, 
                              label = paste("Spatial plot unavailable\n", n_clusters, "clusters found")), 
                          size = 6) +
                theme_void()
            })
          } else {
            ggplot() + 
              geom_text(aes(x = 0.5, y = 0.5, 
                            label = paste("No spatial coordinates\n", n_clusters, "clusters found")), 
                        size = 6) +
              theme_void()
          }
        })
        
        # Generate cluster statistics table
        output$spatial_cluster_stats <- renderDT({
          tryCatch({
            cluster_stats <- obj@meta.data %>%
              group_by(seurat_clusters) %>%
              summarise(
                n_spots = n(),
                mean_genes = round(mean(nFeature_Spatial, na.rm = TRUE), 1),
                mean_umis = round(mean(nCount_Spatial, na.rm = TRUE), 1),
                mean_mt = round(mean(percent.mt, na.rm = TRUE), 2),
                .groups = 'drop'
              )
            
            datatable(cluster_stats, 
                      options = list(pageLength = 15, scrollX = TRUE, dom = 'tip'),
                      rownames = FALSE) %>%
              formatStyle(columns = 1:ncol(cluster_stats), fontSize = '12px')
          }, error = function(e) {
            message(paste("Cluster stats error:", e$message))
            datatable(data.frame(Message = "Statistics unavailable"))
          })
        })
      }
    } else {
      # Show placeholder plots when not clustered
      output$spatial_umap_clusters <- renderPlot({
        if (!is.null(spatial_obj())) {
          obj <- spatial_obj()
          if ("umap" %in% names(obj@reductions)) {
            DimPlot(obj, reduction = "umap", raster = FALSE) + 
              ggtitle("UMAP (Ready for clustering)") +
              theme(plot.title = element_text(hjust = 0.5))
          } else {
            ggplot() + 
              geom_text(aes(x = 0.5, y = 0.5, label = "Run neighbors & UMAP first"), size = 6) +
              theme_void()
          }
        } else {
          ggplot() + 
            geom_text(aes(x = 0.5, y = 0.5, label = "No data loaded"), size = 6) +
            theme_void()
        }
      })
      
      output$spatial_tissue_clusters <- renderPlot({
        if (!is.null(spatial_obj())) {
          obj <- spatial_obj()
          if (length(obj@images) > 0) {
            ggplot() + 
              geom_text(aes(x = 0.5, y = 0.5, label = "Ready for clustering"), size = 6) +
              theme_void()
          } else {
            ggplot() + 
              geom_text(aes(x = 0.5, y = 0.5, label = "No spatial coordinates"), size = 6) +
              theme_void()
          }
        } else {
          ggplot() + 
            geom_text(aes(x = 0.5, y = 0.5, label = "No data loaded"), size = 6) +
            theme_void()
        }
      })
    }
  })
  
  # Status displays for processing state
  output$spatial_processing_status <- renderInfoBox({
    status <- get_processing_status()
    
    if (grepl("Clustered", status)) {
      color <- "green"
      icon_name <- "check-circle"
    } else if (grepl("PCA", status)) {
      color <- "blue"
      icon_name <- "chart-line"
    } else if (grepl("Normalized", status)) {
      color <- "yellow"
      icon_name <- "balance-scale"
    } else {
      color <- "red"
      icon_name <- "exclamation-triangle"
    }
    
    infoBox(
      title = "Processing Status",
      value = status,
      icon = icon(icon_name),
      color = color,
      width = 12
    )
  })
  
  # Individual status outputs
  output$spatial_neighbors_status <- renderText({
    if (!is.null(spatial_obj())) {
      obj <- spatial_obj()
      if (!is.null(obj@graphs) && length(obj@graphs) > 0) {
        "✓ Neighbors computed"
      } else {
        "⚪ Neighbors not computed"
      }
    } else {
      "⚪ No data loaded"
    }
  })
  
  output$spatial_clusters_status <- renderText({
    if (!is.null(spatial_obj())) {
      obj <- spatial_obj()
      if ("seurat_clusters" %in% colnames(obj@meta.data)) {
        n_clusters <- length(unique(obj$seurat_clusters))
        paste("✓ Clusters found:", n_clusters)
      } else {
        "⚪ Clusters not computed"
      }
    } else {
      "⚪ No data loaded"
    }
  })
  
  output$spatial_pca_status <- renderText({
    if (processing_states$pca_computed) {
      "✓ PCA completed"
    } else {
      "⚪ PCA not computed"
    }
  })
  
  output$spatial_umap_status <- renderText({
    if (!is.null(spatial_obj())) {
      obj <- spatial_obj()
      if ("umap" %in% names(obj@reductions)) {
        "✓ UMAP completed"
      } else {
        "⚪ UMAP not computed"
      }
    } else {
      "⚪ No data loaded"
    }
  })
  
  
  # Add this helper function to validate spatial data integrity
  validate_spatial_clustering_data <- function(obj) {
    issues <- list()
    
    # Check cluster assignments
    if ("seurat_clusters" %in% colnames(obj@meta.data)) {
      cluster_data <- obj@meta.data$seurat_clusters
      na_count <- sum(is.na(cluster_data))
      inf_count <- sum(is.infinite(cluster_data))
      
      if (na_count > 0) issues <- append(issues, paste(na_count, "NA cluster assignments"))
      if (inf_count > 0) issues <- append(issues, paste(inf_count, "infinite cluster assignments"))
    }
    
    # Check spatial coordinates
    if (length(obj@images) > 0) {
      tryCatch({
        coords <- GetTissueCoordinates(obj@images[[1]])
        if (!is.null(coords)) {
          coord_na <- sum(is.na(coords))
          coord_inf <- sum(is.infinite(as.matrix(coords)))
          
          if (coord_na > 0) issues <- append(issues, paste(coord_na, "NA coordinates"))
          if (coord_inf > 0) issues <- append(issues, paste(coord_inf, "infinite coordinates"))
        }
      }, error = function(e) {
        issues <- append(issues, "Coordinate extraction failed")
      })
    }
    
    return(list(valid = length(issues) == 0, issues = issues))
  }
  
  # Helper function to check if spatial coordinates are available
  check_spatial_coordinates <- function(obj) {
    if (length(obj@images) == 0) {
      return(list(available = FALSE, message = "No spatial images found"))
    }
    
    tryCatch({
      coords <- GetTissueCoordinates(obj@images[[1]])
      if (is.null(coords) || nrow(coords) == 0) {
        return(list(available = FALSE, message = "No tissue coordinates found"))
      }
      
      # Check if coordinates match the cells
      cell_names <- colnames(obj)
      coord_names <- rownames(coords)
      
      overlap <- intersect(cell_names, coord_names)
      if (length(overlap) == 0) {
        return(list(available = FALSE, message = "No matching coordinates for cells"))
      }
      
      return(list(available = TRUE, 
                  message = paste("Coordinates available for", length(overlap), "spots"),
                  coords = coords))
      
    }, error = function(e) {
      return(list(available = FALSE, message = paste("Error accessing coordinates:", e$message)))
    })
  }
  
  # Update the spatial clustering plot to handle missing coordinates better
  observe({
    req(spatial_obj())
    
    if (processing_states$clustered) {
      obj <- spatial_obj()
      
      if ("seurat_clusters" %in% colnames(obj@meta.data)) {
        n_clusters <- length(unique(obj$seurat_clusters))
        
        # Generate spatial clustering plot with better error handling
        output$spatial_tissue_clusters <- renderPlot({
          coord_check <- check_spatial_coordinates(obj)
          
          if (!coord_check$available) {
            # Show informative message
            ggplot() + 
              geom_text(aes(x = 0.5, y = 0.5, 
                            label = paste("Spatial visualization unavailable\n", 
                                          coord_check$message, "\n",
                                          n_clusters, "clusters found")), 
                        size = 5, color = "orange") +
              theme_void() +
              labs(title = "Use UMAP view for cluster visualization")
          } else {
            tryCatch({
              SpatialDimPlot(obj, group.by = "seurat_clusters", 
                             label = input$spatial_plot_labels,
                             pt.size.factor = 1.2, alpha = c(0.3, 1),
                             label.size = 3) +
                ggtitle(paste("Spatial Clusters (n =", n_clusters, ")")) +
                theme(plot.title = element_text(hjust = 0.5, size = 14),
                      legend.text = element_text(size = 10)) +
                guides(fill = guide_legend(override.aes = list(size = 3)))
            }, error = function(e) {
              message(paste("Spatial plot error:", e$message))
              ggplot() + 
                geom_text(aes(x = 0.5, y = 0.5, 
                              label = paste("Spatial plot error\n", 
                                            "Try using UMAP view\n",
                                            n_clusters, "clusters found")), 
                          size = 5, color = "red") +
                theme_void()
            })
          }
        })
      }
    }
  })
  
  ##########################Interactive Visualization##############################
  
  # Interactive spatial plot - FIXED VERSION pour sketch assay
  # Interactive spatial plot - MANUAL REFRESH ONLY
  output$interactive_spatial_plot <- renderPlotly({
    # Only trigger on refresh button click
    input$refresh_interactive_spatial
    
    # Don't auto-trigger on other inputs
    # input$spatial_interactive_feature  # REMOVED
    # input$spatial_interactive_pt_size  # REMOVED
    # etc.
    
    message("=== MANUAL REFRESH SPATIAL PLOT ===")
    
    req(spatial_obj())
    
    tryCatch({
      obj <- spatial_obj()
      
      if (length(obj@images) == 0 || !"seurat_clusters" %in% colnames(obj@meta.data)) {
        p <- ggplot() + 
          geom_text(aes(x = 0.5, y = 0.5, label = "Run clustering first"), size = 8) + 
          theme_void()
        return(plotly::ggplotly(p))
      }
      
      # Get current parameters (only when refresh is clicked)
      feature_to_show <- isolate(input$spatial_interactive_feature) %||% "seurat_clusters"
      max_points <- isolate(input$spatial_interactive_max_points) %||% 10000
      pt_size <- isolate(input$spatial_interactive_pt_size) %||% 0.8
      alpha_val <- isolate(input$spatial_interactive_alpha) %||% 0.7
      
      message(paste("Feature:", feature_to_show, "Max points:", max_points))
      
      # Get coordinates and prepare data
      coords <- GetTissueCoordinates(obj@images[[1]])
      current_assay <- DefaultAssay(obj)
      
      # Handle sketch vs spatial assay
      if (current_assay == "sketch") {
        sketch_cells <- colnames(obj[["sketch"]])
        coords_available <- rownames(coords)
        common_cells <- intersect(sketch_cells, coords_available)
        
        if (length(common_cells) == 0) {
          p <- ggplot() + 
            geom_text(aes(x = 0.5, y = 0.5, 
                          label = "Switch to 'Spatial' assay\nfor spatial visualization"), 
                      size = 6, color = "orange") + 
            theme_void()
          return(plotly::ggplotly(p))
        }
        
        coords <- coords[common_cells, ]
        meta_subset <- obj@meta.data[common_cells, ]
      } else {
        coords_available <- rownames(coords)
        common_cells <- intersect(colnames(obj), coords_available)
        coords <- coords[common_cells, ]
        meta_subset <- obj@meta.data[common_cells, ]
      }
      
      # Check if feature exists
      if (!feature_to_show %in% colnames(meta_subset)) {
        p <- ggplot() + 
          geom_text(aes(x = 0.5, y = 0.5, label = paste("Feature", feature_to_show, "not available")), 
                    size = 6, color = "red") + 
          theme_void()
        return(plotly::ggplotly(p))
      }
      
      # Create plot data
      plot_data <- data.frame(
        x = coords[,1],
        y = coords[,2],
        feature_value = meta_subset[[feature_to_show]],
        umi = meta_subset$nCount_Spatial,
        genes = meta_subset$nFeature_Spatial,
        mt = meta_subset$percent.mt,
        spot_id = rownames(coords),
        stringsAsFactors = FALSE
      )
      
      # Add cluster info if available
      if ("seurat_clusters" %in% colnames(meta_subset)) {
        plot_data$cluster <- as.factor(meta_subset$seurat_clusters)
      } else {
        plot_data$cluster <- "Unknown"
      }
      
      # Remove incomplete cases
      plot_data <- plot_data[complete.cases(plot_data[, c("x", "y", "feature_value")]), ]
      
      # Smart subsampling
      if (nrow(plot_data) > max_points) {
        set.seed(123)
        if (feature_to_show == "seurat_clusters" && length(unique(plot_data$cluster)) > 1) {
          # Stratified sampling by cluster
          plot_data <- plot_data %>%
            group_by(cluster) %>%
            sample_n(min(n(), ceiling(max_points / length(unique(plot_data$cluster))))) %>%
            ungroup() %>%
            as.data.frame()
        } else {
          # Random sampling
          plot_data <- plot_data[sample(nrow(plot_data), max_points), ]
        }
      }
      
      # Create hover text
      plot_data$text <- paste(
        "Spot:", plot_data$spot_id,
        "<br>", feature_to_show, ":", plot_data$feature_value,
        "<br>Cluster:", plot_data$cluster,
        "<br>UMI:", format(plot_data$umi, big.mark = ","),
        "<br>Genes:", plot_data$genes,
        "<br>MT%:", round(plot_data$mt, 2)
      )
      
      # Create plot based on feature type
      if (feature_to_show == "seurat_clusters") {
        # Categorical plot for clusters
        p <- ggplot(plot_data, aes(x = x, y = y, color = as.factor(feature_value), text = text)) +
          geom_point(size = pt_size, alpha = alpha_val) +
          scale_y_reverse() +
          theme_minimal() +
          ggtitle(paste("Spatial Clusters (", nrow(plot_data), "spots)")) +
          theme(plot.title = element_text(hjust = 0.5, size = 12)) +
          labs(x = "Spatial X", y = "Spatial Y", color = "Cluster")
        
      } else {
        # Continuous plot for other features
        p <- ggplot(plot_data, aes(x = x, y = y, color = as.numeric(feature_value), text = text)) +
          geom_point(size = pt_size, alpha = alpha_val) +
          scale_color_viridis_c(name = feature_to_show) +
          scale_y_reverse() +
          theme_minimal() +
          ggtitle(paste("Spatial", feature_to_show, "(", nrow(plot_data), "spots)")) +
          theme(plot.title = element_text(hjust = 0.5, size = 12)) +
          labs(x = "Spatial X", y = "Spatial Y")
      }
      
      # Convert to plotly
      plotly_obj <- plotly::ggplotly(p, tooltip = "text") %>%
        plotly::layout(
          dragmode = "pan", 
          hovermode = "closest"
        ) %>%
        plotly::config(
          displayModeBar = TRUE,
          scrollZoom = TRUE,
          modeBarButtonsToRemove = c("select2d", "lasso2d")
        )
      
      return(plotly_obj)
      
    }, error = function(e) {
      message(paste("Plot error:", e$message))
      p <- ggplot() + 
        geom_text(aes(x = 0.5, y = 0.5, label = "Plot error - check data"), 
                  size = 6, color = "red") + 
        theme_void()
      plotly::ggplotly(p)
    })
  })
  
  # Info display
  output$spatial_plot_info <- renderText({
    req(spatial_obj())
    
    obj <- spatial_obj()
    current_assay <- DefaultAssay(obj)
    n_cells <- ncol(obj)
    max_points <- input$spatial_interactive_max_points %||% 10000
    
    info_text <- paste(
      "Current assay:", current_assay,
      "\nTotal cells:", format(n_cells, big.mark = ","),
      "\nMax display points:", format(max_points, big.mark = ",")
    )
    
    if (n_cells > max_points) {
      info_text <- paste(info_text, "\n⚡ Will be subsampled for performance")
    }
    
    return(info_text)
  })
  
  
  
  # Update gene choices when spatial data is loaded
  observe({
    req(spatial_obj())
    
    if (processing_states$normalized || !is.null(spatial_obj())) {
      obj <- spatial_obj()
      gene_choices <- rownames(obj)
      
      # Update gene selection picker
      updatePickerInput(session, "spatial_genes_select",
                        choices = gene_choices,
                        selected = NULL)
      
      # Update cluster choices if clusters exist
      if ("seurat_clusters" %in% colnames(obj@meta.data)) {
        cluster_choices <- unique(obj@meta.data$seurat_clusters)
        cluster_choices <- sort(as.character(cluster_choices))
        
        updateSelectInput(session, "spatial_cluster_order",
                          choices = cluster_choices,
                          selected = cluster_choices)
      }
    }
  })
  
  observe({
    req(spatial_obj())
    
    obj <- spatial_obj()
    
    # Update gene choices for all assays
    all_genes <- character(0)
    current_assay <- input$spatial_viz_assay
    
    if (!is.null(current_assay) && current_assay %in% names(obj@assays)) {
      all_genes <- rownames(obj[[current_assay]])
    } else {
      # Use default assay genes
      all_genes <- rownames(obj[[DefaultAssay(obj)]])
    }
    
    message(paste("Updating gene choices. Available genes:", length(all_genes)))
    
    # Update gene selection picker
    updatePickerInput(session, "spatial_genes_select",
                      choices = all_genes,
                      selected = NULL)
    
    # Update assay choices to include all available assays
    available_assays_list <- names(obj@assays)
    current_default <- DefaultAssay(obj)
    
    updateSelectInput(session, "spatial_viz_assay", 
                      choices = available_assays_list, 
                      selected = current_default)
    
    message(paste("Available assays:", paste(available_assays_list, collapse = ", ")))
    message(paste("Default assay:", current_default))
  })
  
  
  
  # Sync between text inputs and picker selection - SIMPLIFIED AND SAFE
  observeEvent(input$spatial_genes_select, {
    # Only update if the picker has actual selections
    if (!is.null(input$spatial_genes_select) && length(input$spatial_genes_select) > 0) {
      
      # Simple replacement - less complex, more stable
      gene_text <- paste(input$spatial_genes_select, collapse = ", ")
      
      # Update text inputs safely
      updateTextInput(session, "spatial_genes_text", value = gene_text)
      updateTextInput(session, "spatial_violin_genes_text", value = gene_text)
      updateTextInput(session, "spatial_dotplot_genes_text", value = gene_text)
    }
  }, ignoreInit = TRUE)
  
  # Reactive values to store generated plots
  spatial_feature_plot_reactive <- reactiveVal(NULL)
  spatial_violin_plot_reactive <- reactiveVal(NULL)
  spatial_dot_plot_reactive <- reactiveVal(NULL)
  
  # Generate spatial feature plots - AVEC JUSTE LUMINOSITE AJOUTEE
  observeEvent(input$show_spatial_features, {
    genes_from_text <- character(0)
    if (!is.null(input$spatial_genes_text) && nchar(trimws(input$spatial_genes_text)) > 0) {
      genes_from_text <- trimws(strsplit(input$spatial_genes_text, ",")[[1]])
      genes_from_text <- genes_from_text[genes_from_text != ""]
    }
    
    if (length(genes_from_text) == 0) {
      showNotification("Please enter at least one gene", type = "warning")
      return()
    }
    
    req(spatial_obj())
    
    tryCatch({
      obj <- spatial_obj()
      selected_genes <- genes_from_text
      current_assay <- input$spatial_viz_assay %||% "Spatial"
      
      message(paste("=== SPATIAL FEATURE PLOT ==="))
      message(paste("Using assay:", current_assay))
      message(paste("Selected genes:", paste(selected_genes, collapse = ", ")))
      
      # Set the assay temporarily
      original_assay <- DefaultAssay(obj)
      DefaultAssay(obj) <- current_assay
      
      # Validate genes exist in current assay
      available_genes <- rownames(obj[[current_assay]])
      valid_genes <- selected_genes[selected_genes %in% available_genes]
      missing_genes <- selected_genes[!selected_genes %in% available_genes]
      
      if (length(missing_genes) > 0) {
        showNotification(paste("Missing genes:", paste(missing_genes, collapse = ", ")), 
                         type = "warning", duration = 8)
      }
      
      if (length(valid_genes) == 0) {
        showNotification("None of the selected genes found in current assay", type = "error")
        DefaultAssay(obj) <- original_assay
        return()
      }
      
      showModal(modalDialog(
        title = "Generating Spatial Feature Plots",
        paste("Creating plots for", length(valid_genes), "genes using", current_assay, "assay..."),
        easyClose = FALSE, footer = NULL
      ))
      
      # Generate the spatial feature plot with dynamic sizing
      output$spatial_feature_plots <- renderPlot({
        if (length(obj@images) == 0) {
          ggplot() + 
            geom_text(aes(x = 0.5, y = 0.5, label = "No spatial coordinates available"), size = 8) + 
            theme_void()
        } else {
          tryCatch({
            n_genes <- length(valid_genes)
            
            # Get image size factor from slider
            image_size_factor <- input$spatial_image_size %||% 1
            
            # Get brightness factor
            brightness_factor <- input$spatial_brightness %||% 1
            
            # Determine number of columns
            if (input$spatial_combine_plots && n_genes > 1) {
              ncol_val <- input$spatial_plot_ncol %||% min(3, n_genes)
            } else {
              ncol_val <- 1
            }
            
            message(paste("Creating spatial plot with", n_genes, "genes in", ncol_val, "columns"))
            message(paste("Image size factor:", image_size_factor))
            message(paste("Brightness factor:", brightness_factor))
            
            # Create plots based on number of genes
            if (n_genes == 1) {
              # Single gene plot
              p <- SpatialFeaturePlot(
                object = obj, 
                features = valid_genes[1],
                crop = input$spatial_crop %||% TRUE,
                slot = "data",
                pt.size.factor = (input$spatial_pt_size %||% 1.6) * image_size_factor,
                alpha = c(0.1, (input$spatial_alpha %||% 1) * brightness_factor),
                image.alpha = 1,
                stroke = 0.25 * image_size_factor
              ) +
                ggtitle(paste("Spatial Expression:", valid_genes[1], "(", current_assay, ")")) +
                theme(plot.title = element_text(hjust = 0.5, size = 14 * image_size_factor))
              
            } else if (input$spatial_combine_plots) {
              # Multi-gene combined plot
              p <- SpatialFeaturePlot(
                object = obj, 
                features = valid_genes,
                crop = input$spatial_crop %||% TRUE,
                slot = "data",
                ncol = ncol_val,
                combine = TRUE,
                pt.size.factor = (input$spatial_pt_size %||% 1.6) * image_size_factor,
                alpha = c(0.1, (input$spatial_alpha %||% 1) * brightness_factor),
                image.alpha = 0.8,
                stroke = 0.25 * image_size_factor
              )
              
            } else {
              # Show first gene only if not combining
              p <- SpatialFeaturePlot(
                object = obj, 
                features = valid_genes[1],
                crop = input$spatial_crop %||% TRUE,
                slot = "data",
                pt.size.factor = (input$spatial_pt_size %||% 1.6) * image_size_factor,
                alpha = c(0.1, (input$spatial_alpha %||% 1) * brightness_factor),
                image.alpha = 1,
                stroke = 0.25 * image_size_factor
              ) +
                ggtitle(paste("Spatial Expression:", valid_genes[1], "(", current_assay, ")")) +
                theme(plot.title = element_text(hjust = 0.5, size = 14 * image_size_factor))
            }
            
            return(p)
            
          }, error = function(e) {
            message(paste("Spatial plot error:", e$message))
            ggplot() + 
              geom_text(aes(x = 0.5, y = 0.5, label = paste("Error:", substr(e$message, 1, 50))), 
                        size = 4, color = "red") +
              theme_void()
          })
        }
      }, 
      height = function() {
        size_factor <- input$spatial_image_size %||% 1
        n_genes <- length(valid_genes)
        
        if (n_genes == 1) {
          return(600 * size_factor)
        } else {
          ncol_val <- input$spatial_plot_ncol %||% min(3, n_genes)
          nrow_val <- ceiling(n_genes / ncol_val)
          return(min(1200, 400 * nrow_val * size_factor))
        }
      })
      
      # Store plot parameters for download
      spatial_plots$feature_plots <- list(
        genes = valid_genes, 
        assay = current_assay,
        parameters = list(
          pt_size = (input$spatial_pt_size %||% 1.6) * (input$spatial_image_size %||% 1), 
          alpha = (input$spatial_alpha %||% 1) * (input$spatial_brightness %||% 1),
          crop = input$spatial_crop %||% TRUE, 
          combine = input$spatial_combine_plots,
          image_size = input$spatial_image_size %||% 1
        )
      )
      
      # Restore original assay
      DefaultAssay(obj) <- original_assay
      
      removeModal()
      showNotification(paste("Spatial plots generated for", length(valid_genes), "genes"), type = "message")
      spatial_feature_plot_reactive(p)
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error generating spatial plots:", e$message), type = "error")
    })
  })
  
  # Generate spatial dot plot by clusters
  # Update spatial_dotplot to use annotated clusters
  observeEvent(input$show_spatial_dotplot, {
    genes_from_text <- character(0)
    if (!is.null(input$spatial_dotplot_genes_text) && nchar(trimws(input$spatial_dotplot_genes_text)) > 0) {
      genes_from_text <- trimws(strsplit(input$spatial_dotplot_genes_text, ",")[[1]])
      genes_from_text <- genes_from_text[genes_from_text != ""]
    }
    
    if (length(genes_from_text) == 0) {
      showNotification("Please enter at least one gene for dot plot", type = "warning")
      return()
    }
    
    req(spatial_obj())
    
    tryCatch({
      obj <- spatial_obj()
      
      # Use annotated clusters if available
      cluster_col <- if ("annotated_clusters" %in% colnames(obj@meta.data)) {
        "annotated_clusters"
      } else if ("seurat_clusters" %in% colnames(obj@meta.data)) {
        "seurat_clusters"
      } else {
        showNotification("Please run clustering first", type = "warning")
        return()
      }
      
      selected_genes <- genes_from_text
      selected_assay <- input$spatial_viz_assay %||% "Spatial"
      
      # Validate genes exist in assay
      available_genes <- rownames(obj[[selected_assay]])
      valid_genes <- selected_genes[selected_genes %in% available_genes]
      
      if (length(valid_genes) == 0) {
        showNotification("Selected genes not found in current assay", type = "error")
        return()
      }
      
      # Set assay temporarily
      original_assay <- DefaultAssay(obj)
      DefaultAssay(obj) <- selected_assay
      
      # Generate plot
      p <- DotPlot(obj, 
                   features = valid_genes, 
                   group.by = cluster_col,  # Use the appropriate cluster column
                   scale = input$spatial_dot_scale %||% TRUE, 
                   scale.by = input$spatial_dot_scale_by %||% "radius",
                   dot.min = (input$spatial_dot_min %||% 0)/100, 
                   dot.scale = input$spatial_dot_max %||% 6,
                   cluster.idents = input$spatial_dot_cluster_idents %||% FALSE, 
                   cols = input$spatial_dot_colors %||% "RdYlBu") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), 
              plot.title = element_text(hjust = 0.5)) +
        labs(x = "Genes", y = "Clusters") + 
        coord_flip()
      
      # Store and display
      spatial_dot_plot_reactive(p)
      output$spatial_dotplot <- renderPlot({ p })
      
      # Restore original assay
      DefaultAssay(obj) <- original_assay
      
      showNotification("Dot plot generated", type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error generating dot plot:", e$message), type = "error")
    })
  })
  
  # Generate spatial violin plot by clusters - PROPERLY FIXED
  observeEvent(input$show_spatial_violin, {
    genes_from_text <- character(0)
    if (!is.null(input$spatial_violin_genes_text) && nchar(trimws(input$spatial_violin_genes_text)) > 0) {
      genes_from_text <- trimws(strsplit(input$spatial_violin_genes_text, ",")[[1]])
      genes_from_text <- genes_from_text[genes_from_text != ""]
    }
    
    if (length(genes_from_text) == 0) {
      showNotification("Please enter at least one gene for violin plot", type = "warning")
      return()
    }
    
    req(spatial_obj())
    
    tryCatch({
      obj <- spatial_obj()
      
      # Use annotated clusters if available
      cluster_col <- if ("annotated_clusters" %in% colnames(obj@meta.data)) {
        "annotated_clusters"
      } else if ("seurat_clusters" %in% colnames(obj@meta.data)) {
        "seurat_clusters"
      } else {
        showNotification("Please run clustering first", type = "warning")
        return()
      }
      
      selected_genes <- genes_from_text
      selected_assay <- input$spatial_viz_assay %||% "Spatial"
      
      # Validate genes exist in assay
      available_genes <- rownames(obj[[selected_assay]])
      valid_genes <- selected_genes[selected_genes %in% available_genes]
      missing_genes <- selected_genes[!selected_genes %in% available_genes]
      
      if (length(missing_genes) > 0) {
        showNotification(paste("Missing genes:", paste(missing_genes, collapse = ", ")), 
                         type = "warning", duration = 8)
      }
      
      if (length(valid_genes) == 0) {
        showNotification("Selected genes not found in current assay", type = "error")
        return()
      }
      
      # Set assay temporarily
      original_assay <- DefaultAssay(obj)
      DefaultAssay(obj) <- selected_assay
      
      message(paste("Creating violin plot with genes:", paste(valid_genes, collapse = ", ")))
      message(paste("Using cluster column:", cluster_col))
      message(paste("Using assay:", selected_assay))
      
      # Create violin plot with CORRECT parameters
      p <- VlnPlot(
        object = obj,
        features = valid_genes,
        assay = selected_assay,          # Specify assay explicitly
        group.by = cluster_col,          # Group by cluster column
        pt.size = if (input$spatial_violin_points) input$spatial_violin_pt_size else 0,
        log = input$spatial_violin_log,
        ncol = if (length(valid_genes) > 1) min(3, length(valid_genes)) else NULL,
        combine = TRUE,
        fill.by = "ident",               # Color by identity (clusters)
        layer = "data",                  # Use processed data layer
        flip = FALSE,                    # Keep standard orientation
        same.y.lims = FALSE,            # Allow different y-axis for each gene
        raster = FALSE                   # Don't rasterize for better quality
      )
      
      # Apply custom theming and labels
      if (length(valid_genes) == 1) {
        # Single gene - enhance appearance
        p <- p + 
          ggtitle(paste("Expression of", valid_genes[1], "across", cluster_col)) +
          labs(
            x = paste("Clusters (", cluster_col, ")"),
            y = paste0("Expression Level", if(input$spatial_violin_log) " (log)" else ""),
            subtitle = paste("Assay:", selected_assay)
          ) +
          theme(
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 12),
            axis.title.x = element_text(size = 14, face = "bold"),
            axis.title.y = element_text(size = 14, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
            axis.text.y = element_text(size = 12),
            legend.title = element_text(size = 14, face = "bold"),
            legend.text = element_text(size = 12),
            legend.position = "right"
          )
        
      } else {
        # Multiple genes - use patchwork for better control
        p <- p & theme(
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size = 12, face = "bold"),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10),
          strip.text = element_text(size = 12, face = "bold")  # Gene names
        )
        
        # Add overall annotation
        p <- p + plot_annotation(
          title = paste("Gene Expression across", cluster_col),
          subtitle = paste("Assay:", selected_assay, "| Genes:", paste(valid_genes, collapse = ", ")),
          theme = theme(
            plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 12)
          )
        )
      }
      
      # Apply cluster order if specified
      cluster_order <- input$spatial_cluster_order
      if (!is.null(cluster_order) && length(cluster_order) > 0) {
        # Get current clusters and apply order
        current_clusters <- levels(as.factor(obj@meta.data[[cluster_col]]))
        valid_order <- cluster_order[cluster_order %in% current_clusters]
        
        if (length(valid_order) > 0) {
          # Add missing clusters at the end
          missing_clusters <- setdiff(current_clusters, valid_order)
          final_order <- c(valid_order, missing_clusters)
          
          # Apply the order
          p <- p + scale_x_discrete(limits = final_order)
        }
      }
      
      # Store in reactive value for reuse
      spatial_violin_plot_reactive(p)
      
      # Display the plot
      output$spatial_violin_plot <- renderPlot({
        spatial_violin_plot_reactive()
      })
      
      # Restore original assay
      DefaultAssay(obj) <- original_assay
      
      
    }, error = function(e) {
      message(paste("Violin plot error:", e$message))
      showNotification(paste("Error generating violin plot:", e$message), type = "error")
    })
  })
  
  # Download handlers using stored plots
  output$download_spatial_features <- downloadHandler(
    filename = function() {
      paste0("spatial_features_", Sys.Date(), ".", input$spatial_plot_format %||% "png")
    },
    content = function(file) {
      p <- spatial_feature_plot_reactive()
      
      if (!is.null(p)) {
        plot_width <- input$spatial_plot_width %||% 10
        plot_height <- plot_width * 0.8
        
        ggsave(file, plot = p, 
               width = plot_width, 
               height = plot_height,
               dpi = input$spatial_export_dpi %||% 300, 
               device = input$spatial_plot_format %||% "png")
      }
    }
  )
  
  output$download_spatial_violin <- downloadHandler(
    filename = function() {
      paste0("spatial_violin_", Sys.Date(), ".", input$spatial_plot_format %||% "png")
    },
    content = function(file) {
      p <- spatial_violin_plot_reactive()
      
      if (!is.null(p)) {
        plot_width <- input$spatial_plot_width %||% 10
        plot_height <- plot_width * 0.6
        
        ggsave(file, plot = p, 
               width = plot_width, 
               height = plot_height,
               dpi = input$spatial_export_dpi %||% 300,
               device = input$spatial_plot_format %||% "png")
      }
    }
  )
  
  output$download_spatial_dotplot <- downloadHandler(
    filename = function() {
      paste0("spatial_dotplot_", Sys.Date(), ".", input$spatial_plot_format %||% "png")
    },
    content = function(file) {
      p <- spatial_dot_plot_reactive()
      
      if (!is.null(p)) {
        plot_width <- input$spatial_plot_width %||% 10
        plot_height <- plot_width * 0.8
        
        ggsave(file, plot = p, 
               width = plot_width, 
               height = plot_height,
               dpi = input$spatial_export_dpi %||% 300,
               device = input$spatial_plot_format %||% "png")
      }
    }
  )
  
  # Save spatial object
  output$save_spatial_object <- downloadHandler(
    filename = function() {
      paste0("spatial_object_", Sys.Date(), ".rds")
    },
    content = function(file) {
      showModal(modalDialog(
        title = "Saving Spatial Object",
        "Saving your spatial data object...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      req(spatial_obj())
      saveRDS(spatial_obj(), file)
      
      removeModal()
      showNotification("Spatial object saved successfully!", type = "message")
    }
  )
  
  
  
  ############################## Cluster Annotation ##############################
  
  # Reactive value to store cluster names
  spatial_cluster_names <- reactiveVal(NULL)
  
  # Initialize cluster names when clusters are available
  observe({
    req(spatial_obj())
    
    if (processing_states$clustered && "seurat_clusters" %in% colnames(spatial_obj()@meta.data)) {
      obj <- spatial_obj()
      unique_clusters <- sort(unique(as.character(obj$seurat_clusters)))
      
      # Initialize with default names if not already set
      if (is.null(spatial_cluster_names())) {
        default_names <- setNames(paste0("Cluster_", unique_clusters), unique_clusters)
        spatial_cluster_names(default_names)
      }
    }
  })
  
  # Generate dynamic input fields for cluster renaming
  output$spatial_cluster_rename_inputs <- renderUI({
    req(spatial_obj())
    req(processing_states$clustered)
    
    obj <- spatial_obj()
    if (!"seurat_clusters" %in% colnames(obj@meta.data)) {
      return(p("No clusters found. Please run clustering first.", style = "color: #666;"))
    }
    
    unique_clusters <- sort(unique(as.character(obj$seurat_clusters)))
    current_names <- spatial_cluster_names()
    
    # Create input fields for each cluster
    input_list <- lapply(unique_clusters, function(cluster) {
      current_name <- if (!is.null(current_names) && cluster %in% names(current_names)) {
        current_names[[cluster]]
      } else {
        paste0("Cluster_", cluster)
      }
      
      div(
        class = "form-group",
        style = "margin-bottom: 10px;",
        tags$label(paste0("Cluster ", cluster, ":"), style = "font-weight: bold; color: #333;"),
        textInput(
          inputId = paste0("spatial_cluster_name_", cluster),
          label = NULL,
          value = current_name,
          placeholder = paste0("Name for cluster ", cluster)
        )
      )
    })
    
    do.call(tagList, input_list)
  })
  
  # Apply all cluster names - FIXED VERSION
  observeEvent(input$apply_all_spatial_names, {
    message("=== APPLY CLUSTER NAMES STARTED ===")
    
    req(spatial_obj())
    
    tryCatch({
      obj <- spatial_obj()
      
      # Check if clusters exist
      if (!"seurat_clusters" %in% colnames(obj@meta.data)) {
        showNotification("No clusters found in the object!", type = "error")
        return()
      }
      
      unique_clusters <- sort(unique(as.character(obj$seurat_clusters)))
      message(paste("Found clusters:", paste(unique_clusters, collapse = ", ")))
      
      # Collect all new names
      new_names <- list()
      for (cluster in unique_clusters) {
        input_id <- paste0("spatial_cluster_name_", cluster)
        new_name <- input[[input_id]]
        
        if (!is.null(new_name) && nchar(trimws(new_name)) > 0) {
          new_names[[cluster]] <- trimws(new_name)
          message(paste("Cluster", cluster, "->", new_names[[cluster]]))
        } else {
          new_names[[cluster]] <- paste0("Cluster_", cluster)
          message(paste("Cluster", cluster, "-> (default)", new_names[[cluster]]))
        }
      }
      
      # Update reactive value
      spatial_cluster_names(new_names)
      message("Reactive values updated")
      
      # Create new annotated clusters column
      obj$annotated_clusters <- as.factor(obj$seurat_clusters)
      
      # Get the current levels
      current_levels <- levels(obj$annotated_clusters)
      message(paste("Current levels:", paste(current_levels, collapse = ", ")))
      
      # Create new level names in the same order
      new_level_names <- sapply(current_levels, function(x) {
        if (x %in% names(new_names)) {
          return(new_names[[x]])
        } else {
          return(paste0("Cluster_", x))
        }
      })
      
      # Apply new names
      levels(obj$annotated_clusters) <- new_level_names
      message("New cluster names applied to factor levels")
      
      # Update the spatial object
      spatial_obj(obj)
      message("Spatial object updated")
      
      # FIXED: Remove the problematic showNotification call or fix it
      
      message("=== APPLY CLUSTER NAMES COMPLETED ===")
      
    }, error = function(e) {
      message(paste("=== APPLY CLUSTER NAMES ERROR ===", e$message))
      showNotification(paste("Error applying cluster names:", e$message), type = "error", duration = 10)
    })
  })
  
  # Reset cluster names
  observeEvent(input$reset_spatial_names, {
    req(spatial_obj())
    
    obj <- spatial_obj()
    unique_clusters <- sort(unique(as.character(obj$seurat_clusters)))
    
    # Reset to default names
    default_names <- setNames(paste0("Cluster_", unique_clusters), unique_clusters)
    spatial_cluster_names(default_names)
    
    # Remove annotated_clusters column if it exists
    if ("annotated_clusters" %in% colnames(obj@meta.data)) {
      obj$annotated_clusters <- NULL
      spatial_obj(obj)
    }
    
  })
  
  # Live update plots when cluster names change
  observe({
    req(spatial_obj())
    req(spatial_cluster_names())
    
    obj <- spatial_obj()
    
    # Determine which cluster column to use
    cluster_col <- if ("annotated_clusters" %in% colnames(obj@meta.data)) {
      "annotated_clusters"
    } else {
      "seurat_clusters"
    }
    
    # Dans la fonction de visualisation UMAP pour annotation
    output$spatial_annotation_umap <- renderPlot({
      if ("umap" %in% names(obj@reductions)) {
        DimPlot(obj, 
                reduction = "umap", 
                group.by = cluster_col,
                label = TRUE, 
                label.size = 5, 
                repel = TRUE,
                pt.size = 1,
                raster = FALSE) +  # Ajoutez cette ligne
          ggtitle("UMAP - Annotated Clusters") +
          theme(plot.title = element_text(hjust = 0.5, size = 16),
                legend.text = element_text(size = 12))
      }
    })
    
    # Spatial plot with current names
    output$spatial_annotation_tissue <- renderPlot({
      if (length(obj@images) > 0) {
        tryCatch({
          SpatialDimPlot(obj, 
                         group.by = cluster_col,
                         label = TRUE,
                         label.size = 4,
                         pt.size.factor = 1.6,
                         alpha = c(0.3, 1)) +
            ggtitle("Spatial View - Annotated Clusters") +
            theme(plot.title = element_text(hjust = 0.5, size = 16),
                  legend.text = element_text(size = 12))
        }, error = function(e) {
          ggplot() + 
            geom_text(aes(x = 0.5, y = 0.5, label = "Spatial plot unavailable"), 
                      size = 6) +
            theme_void()
        })
      } else {
        ggplot() + 
          geom_text(aes(x = 0.5, y = 0.5, label = "No spatial coordinates found"), 
                    size = 6) +
          theme_void()
      }
    })
    
  })
  
  
  # Display selected plot in annotation tab
  output$annotation_selected_plot_display <- renderUI({
    selected_plot <- input$annotation_selected_plot
    
    if (is.null(selected_plot) || selected_plot == "none") {
      return(
        div(style = "text-align: center; padding: 50px; color: #666;",
            h4("No plot selected"),
            p("Choose a plot type above to display gene expression data")
        )
      )
    }
    
    # Trigger refresh when button is clicked
    input$refresh_annotation_plot
    
    plot_height <- "400px"
    
    switch(selected_plot,
           "violin" = {
             if (exists("spatial_violin_plot_reactive") && !is.null(spatial_violin_plot_reactive())) {
               plotOutput("annotation_violin_display", height = plot_height)
             } else {
               div(style = "text-align: center; padding: 30px; color: #orange;",
                   h5("No violin plot available"),
                   p("Generate a violin plot in the 'Gene Expression Visualization' tab first"))
             }
           },
           "dot" = {
             if (exists("spatial_dot_plot_reactive") && !is.null(spatial_dot_plot_reactive())) {
               plotOutput("annotation_dot_display", height = plot_height)
             } else {
               div(style = "text-align: center; padding: 30px; color: #orange;",
                   h5("No dot plot available"), 
                   p("Generate a dot plot in the 'Gene Expression Visualization' tab first"))
             }
           },
           "feature" = {
             if (exists("spatial_feature_plot_reactive") && !is.null(spatial_feature_plot_reactive())) {
               plotOutput("annotation_feature_display", height = plot_height)
             } else {
               div(style = "text-align: center; padding: 30px; color: #orange;",
                   h5("No feature plot available"),
                   p("Generate a feature plot in the 'Gene Expression Visualization' tab first"))
             }
           },
           div(style = "text-align: center; padding: 30px; color: #red;",
               h5("Unknown plot type selected"))
    )
  })
  
  # Render plots in annotation tab
  output$annotation_violin_display <- renderPlot({
    if (!is.null(spatial_violin_plot_reactive())) {
      spatial_violin_plot_reactive()
    }
  })
  
  output$annotation_dot_display <- renderPlot({
    if (!is.null(spatial_dot_plot_reactive())) {
      spatial_dot_plot_reactive()
    }
  })
  
  output$annotation_feature_display <- renderPlot({
    if (!is.null(spatial_feature_plot_reactive())) {
      spatial_feature_plot_reactive()
    }
  })
  
  # Save annotated object
  output$save_annotated_spatial_object <- downloadHandler(
    filename = function() {
      paste0("spatial_annotated_", Sys.Date(), ".rds")
    },
    content = function(file) {
      req(spatial_obj())
      
      showModal(modalDialog(
        title = "Saving Annotated Object",
        "Saving your annotated spatial object...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      saveRDS(spatial_obj(), file)
      
      removeModal()
      
    }
  )
  ############################## Enhanced Spatial Tissue Viewer Server ##############################
  ############################## Spatial Tissue Viewer - Separate Plots ##############################
  
  # Reactive values for H&E viewer
  hne_viewer_state <- reactiveValues(
    image = NULL,
    zoom_level = 1,
    center_x = 0.5,
    center_y = 0.5,
    view_width = 1,
    view_height = 1
  )
  
  # Reactive values for transcriptomics viewer
  trans_viewer_state <- reactiveValues(
    zoom_level = 1,
    center_x = NULL,
    center_y = NULL,
    view_width = NULL,
    view_height = NULL
  )
  
  
  # Realistic image reading function for R
  read_image_file <- function(file_path, file_name) {
    img_ext <- tolower(tools::file_ext(file_name))
    
    tryCatch({
      if (img_ext == "mrxs") {
        # MRXS is problematic in R - provide clear guidance
        stop(paste("MRXS files are not well supported in R.",
                   "\n\nRecommended solutions:",
                   "\n1. Convert MRXS to TIFF using ImageJ/Fiji:",
                   "\n   - Open MRXS in ImageJ",
                   "\n   - File > Export > Bio-Formats Exporter",
                   "\n   - Choose TIFF format",
                   "\n2. Use QuPath to export as PNG/TIFF",
                   "\n3. Use VIPS command line: vips copy file.mrxs file.tiff",
                   "\n4. Convert online with bioformats2raw tools",
                   "\n\nFor now, please use PNG, JPEG, or TIFF formats."))
        
      } else if (img_ext == "png") {
        img <- png::readPNG(file_path)
        return(list(image = img, type = "standard"))
        
      } else if (img_ext %in% c("jpg", "jpeg")) {
        img <- jpeg::readJPEG(file_path)
        return(list(image = img, type = "standard"))
        
      } else if (img_ext %in% c("tiff", "tif")) {
        # Handle multi-page TIFF files
        if (requireNamespace("tiff", quietly = TRUE)) {
          img <- tiff::readTIFF(file_path, all = FALSE)  # Read first page only
          return(list(image = img, type = "standard"))
        } else {
          stop("TIFF support requires the 'tiff' package. Install with: install.packages('tiff')")
        }
        
      } else if (img_ext %in% c("bmp")) {
        if (requireNamespace("bmp", quietly = TRUE)) {
          img <- bmp::read.bmp(file_path)
          return(list(image = img, type = "standard"))
        } else {
          stop("BMP support requires the 'bmp' package. Install with: install.packages('bmp')")
        }
        
      } else {
        stop(paste("Unsupported image format:", img_ext, 
                   ". Supported formats: PNG, JPEG, TIFF, BMP",
                   "\nFor MRXS files, please convert to TIFF using ImageJ or QuPath first."))
      }
      
    }, error = function(e) {
      stop(paste("Error reading image file:", e$message))
    })
  }
  
  # Enhanced load H&E image event with realistic expectations
  observeEvent(input$load_hne_image, {
    req(input$load_hne_image)
    
    tryCatch({
      showModal(modalDialog(
        title = "Loading Image",
        div(
          "Processing image...",
          br(),
          tags$small("Supported: PNG, JPEG, TIFF, BMP"),
          br(),
          tags$small("For MRXS files: Convert to TIFF using ImageJ first")
        ),
        easyClose = FALSE,
        footer = NULL
      ))
      
      img_path <- input$load_hne_image$datapath
      img_name <- input$load_hne_image$name
      img_ext <- tolower(tools::file_ext(img_name))
      file_size <- file.info(img_path)$size / (1024^2)  # Size in MB
      
      message(paste("Loading image:", img_name, "Type:", img_ext, "Size:", round(file_size, 1), "MB"))
      

      
      # Use enhanced reading function
      img_result <- read_image_file(img_path, img_name)
      
      # Store image data
      hne_viewer_state$image <- img_result$image
      hne_viewer_state$image_type <- img_result$type
      hne_viewer_state$file_name <- img_name
      hne_viewer_state$file_size <- file_size
      
      # Reset zoom
      hne_viewer_state$zoom_level <- 1
      hne_viewer_state$center_x <- 0.5
      hne_viewer_state$center_y <- 0.5
      
      # Display image information
      img_dims <- dim(hne_viewer_state$image)
      success_message <- paste("Image loaded successfully:",
                               "\nDimensions:", img_dims[2], "x", img_dims[1], "pixels",
                               "\nSize:", round(file_size, 1), "MB",
                               "\nFormat:", toupper(img_ext))
      
      showNotification(success_message, type = "success", duration = 8)
      removeModal()
      
    }, error = function(e) {
      removeModal()
      
      # Provide specific error messages and practical solutions
      error_msg <- e$message
      
      if (grepl("MRXS", error_msg, ignore.case = TRUE)) {
        showModal(modalDialog(
          title = "MRXS File Not Supported",
          div(
            h4("MRXS files require conversion", style = "color: #d9534f;"),
            p("R cannot directly read MRXS files. Please convert your file first:"),
            
            h5("Option 1: ImageJ/Fiji (Recommended)"),
            tags$ol(
              tags$li("Download ImageJ/Fiji (free)"),
              tags$li("Install Bio-Formats plugin"),
              tags$li("Open your MRXS file"),
              tags$li("File → Export → Bio-Formats Exporter"),
              tags$li("Choose TIFF format"),
              tags$li("Upload the converted TIFF")
            ),
            
            h5("Option 2: QuPath (Alternative)"),
            tags$ol(
              tags$li("Open MRXS in QuPath"),
              tags$li("File → Export snapshot"),
              tags$li("Choose PNG or TIFF"),
              tags$li("Upload the exported file")
            ),
            
            h5("Option 3: Online Conversion"),
            p("Search for 'MRXS to TIFF converter online' for web-based tools.")
          ),
          footer = modalButton("Close"),
          size = "l"
        ))
      } else {
        showNotification(paste("Error loading image:", error_msg), type = "error", duration = 10)
      }
    })
  })
  
  # Handle large file proceed action
  observeEvent(input$proceed_large_file, {
    removeModal()
    
    # Proceed with loading
    img_path <- input$load_hne_image$datapath
    img_name <- input$load_hne_image$name
    
    tryCatch({
      showModal(modalDialog(
        title = "Loading Large Image",
        "This may take several minutes and could affect browser performance...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      img_result <- read_image_file(img_path, img_name)
      
      hne_viewer_state$image <- img_result$image
      hne_viewer_state$image_type <- img_result$type
      hne_viewer_state$file_name <- img_name
      
      # Reset zoom
      hne_viewer_state$zoom_level <- 1
      hne_viewer_state$center_x <- 0.5
      hne_viewer_state$center_y <- 0.5
      
      removeModal()
      showNotification("Large image loaded successfully", type = "success")
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Failed to load large image:", e$message), type = "error")
    })
  })
  
  # Enhanced H&E plot rendering with better zoom handling
  output$hne_plot <- renderPlot({
    if (is.null(hne_viewer_state$image)) {
      # Instructions plot
      ggplot() + 
        geom_text(aes(x = 0.5, y = 0.5, 
                      label = "Load H&E image to begin\n\nSupported formats:\nPNG, JPEG, TIFF, MRXS*\n\n*MRXS requires ROpenSlide"), 
                  size = 6, color = "darkblue") +
        theme_void() +
        theme(panel.border = element_rect(color = "darkblue", fill = NA, size = 2))
    } else {
      
      tryCatch({
        img <- hne_viewer_state$image
        
        # Validate image dimensions
        if (is.null(dim(img)) || length(dim(img)) < 2) {
          stop("Invalid image dimensions")
        }
        
        # Calculate view window based on zoom
        zoom <- hne_viewer_state$zoom_level
        view_width <- 1 / zoom
        view_height <- 1 / zoom
        
        # Calculate bounds with safety checks
        center_x <- max(0, min(1, hne_viewer_state$center_x))
        center_y <- max(0, min(1, hne_viewer_state$center_y))
        
        x_min <- max(0, center_x - view_width/2)
        x_max <- min(1, center_x + view_width/2)
        y_min <- max(0, center_y - view_height/2)
        y_max <- min(1, center_y + view_height/2)
        
        # Convert to pixel coordinates
        img_height <- dim(img)[1]
        img_width <- dim(img)[2]
        
        x_min_px <- max(1, round(x_min * img_width))
        x_max_px <- min(img_width, round(x_max * img_width))
        y_min_px <- max(1, round(y_min * img_height))
        y_max_px <- min(img_height, round(y_max * img_height))
        
        # Ensure valid pixel ranges
        if (x_max_px <= x_min_px) x_max_px <- x_min_px + 1
        if (y_max_px <= y_min_px) y_max_px <- y_min_px + 1
        
        # Extract zoomed region
        if (length(dim(img)) == 3) {
          # Color image
          img_zoomed <- img[y_min_px:y_max_px, x_min_px:x_max_px, ]
        } else {
          # Grayscale
          img_zoomed <- img[y_min_px:y_max_px, x_min_px:x_max_px]
        }
        
        # Render using grid graphics for better performance
        grid::grid.newpage()
        grid::grid.raster(img_zoomed, interpolate = TRUE)
        
        # Add informative overlay
        info_text <- paste("Zoom:", round(zoom, 1), "x")
        if (hne_viewer_state$image_type == "mrxs") {
          info_text <- paste(info_text, "| Level:", hne_viewer_state$current_level)
        }
        
        # White outline for text visibility
        grid::grid.text(info_text, 
                        x = 0.02, y = 0.98, just = c("left", "top"),
                        gp = grid::gpar(col = "white", fontsize = 14, fontface = "bold"))
        grid::grid.text(info_text, 
                        x = 0.022, y = 0.978, just = c("left", "top"),
                        gp = grid::gpar(col = "black", fontsize = 14, fontface = "bold"))
        
      }, error = function(e) {
        message(paste("H&E plot rendering error:", e$message))
        # Fallback error plot
        grid::grid.newpage()
        grid::grid.text(paste("Error displaying image:", substr(e$message, 1, 50)), 
                        x = 0.5, y = 0.5, just = "center",
                        gp = grid::gpar(col = "red", fontsize = 16))
      })
    }
  })
  
  # Enhanced status info with better image information
  output$viewer_status_info <- renderText({
    status_parts <- c()
    
    # H&E status with more details
    if (!is.null(hne_viewer_state$image)) {
      img_dim <- dim(hne_viewer_state$image)
      img_type <- hne_viewer_state$image_type
      
      if (img_type == "mrxs") {
        status_parts <- c(status_parts, 
                          paste("H&E: MRXS file"),
                          paste("Level:", hne_viewer_state$current_level, "/", hne_viewer_state$max_level),
                          paste("Size:", img_dim[2], "x", img_dim[1], "pixels"),
                          paste("Zoom:", round(hne_viewer_state$zoom_level, 1), "x"))
      } else {
        status_parts <- c(status_parts, 
                          paste("H&E:", toupper(img_type)),
                          paste("Size:", img_dim[2], "x", img_dim[1], "pixels"),
                          paste("Zoom:", round(hne_viewer_state$zoom_level, 1), "x"))
      }
    } else {
      status_parts <- c(status_parts, "H&E: Not loaded")
    }
    
    # Transcriptomics status (unchanged)
    if (!is.null(spatial_obj()) && processing_states$clustered) {
      n_spots <- ncol(spatial_obj())
      n_clusters <- length(unique(spatial_obj()$seurat_clusters))
      status_parts <- c(status_parts,
                        "",  # Empty line for separation
                        paste("Transcriptomics:", format(n_spots, big.mark = ","), "spots"),
                        paste("Clusters:", n_clusters),
                        paste("Trans Zoom:", round(trans_viewer_state$zoom_level, 1), "x"))
    } else {
      status_parts <- c(status_parts, "", "Transcriptomics: Not ready")
    }
    
    paste(status_parts, collapse = "\n")
  })
  
  # Transcriptomics Plot rendering
  output$transcriptomics_plot <- renderPlot({
    if (is.null(spatial_obj()) || !processing_states$clustered) {
      # Empty plot
      ggplot() + 
        geom_text(aes(x = 0.5, y = 0.5, label = "Run clustering analysis first"), 
                  size = 8, color = "darkgreen") +
        theme_void() +
        theme(panel.border = element_rect(color = "darkgreen", fill = NA, size = 2))
    } else {
      # Get spatial data
      obj <- spatial_obj()
      coords <- GetTissueCoordinates(obj@images[[1]])
      
      if (is.null(coords)) {
        ggplot() + 
          geom_text(aes(x = 0.5, y = 0.5, label = "No spatial coordinates found"), 
                    size = 6, color = "red") +
          theme_void()
      } else {
        # Prepare plot data
        feature <- input$transcriptomics_feature
        feature_data <- FetchData(obj, vars = c(feature, "seurat_clusters"))
        
        plot_data <- data.frame(
          x = coords[,1],
          y = coords[,2],
          feature = feature_data[[feature]],
          cluster = feature_data$seurat_clusters
        )
        
        # Initialize view if needed
        if (is.null(trans_viewer_state$view_width)) {
          x_range <- range(plot_data$x, na.rm = TRUE)
          y_range <- range(plot_data$y, na.rm = TRUE)
          trans_viewer_state$center_x <- mean(x_range)
          trans_viewer_state$center_y <- mean(y_range)
          trans_viewer_state$view_width <- diff(x_range)
          trans_viewer_state$view_height <- diff(y_range)
        }
        
        # Calculate current view based on zoom
        zoom <- trans_viewer_state$zoom_level
        current_width <- trans_viewer_state$view_width / zoom
        current_height <- trans_viewer_state$view_height / zoom
        
        x_limits <- c(
          trans_viewer_state$center_x - current_width/2,
          trans_viewer_state$center_x + current_width/2
        )
        y_limits <- c(
          trans_viewer_state$center_y - current_height/2,
          trans_viewer_state$center_y + current_height/2
        )
        
        # Create plot
        p <- ggplot(plot_data, aes(x = x, y = y))
        
        if (feature == "seurat_clusters") {
          p <- p + geom_point(aes(color = as.factor(feature)), size = 1.5, alpha = 0.8) +
            scale_color_discrete(name = "Cluster")
        } else {
          p <- p + geom_point(aes(color = feature), size = 1.5, alpha = 0.8) +
            scale_color_viridis_c(name = feature)
        }
        
        p <- p +
          coord_fixed(xlim = x_limits, ylim = y_limits, expand = FALSE) +
          theme_minimal() +
          theme(
            panel.border = element_rect(color = "darkgreen", fill = NA, size = 2),
            legend.position = "right"
          ) +
          labs(title = paste("Zoom:", round(zoom, 1), "x"))
        
        return(p)
      }
    }
  })
  
  # H&E zoom controls
  observeEvent(input$hne_plot_click, {
    req(hne_viewer_state$image)
    
    click <- input$hne_plot_click
    zoom_factor <- input$hne_zoom_factor
    
    # Update center to click location
    hne_viewer_state$center_x <- click$x
    hne_viewer_state$center_y <- click$y
    hne_viewer_state$zoom_level <- hne_viewer_state$zoom_level * zoom_factor
    
    showNotification(paste("H&E zoomed to", round(hne_viewer_state$zoom_level, 1), "x"), 
                     type = "message", duration = 2)
  })
  
  observeEvent(input$hne_zoom_in, {
    req(hne_viewer_state$image)
    hne_viewer_state$zoom_level <- hne_viewer_state$zoom_level * input$hne_zoom_factor
    showNotification(paste("H&E zoom:", round(hne_viewer_state$zoom_level, 1), "x"), 
                     type = "message", duration = 2)
  })
  
  observeEvent(input$hne_zoom_out, {
    req(hne_viewer_state$image)
    hne_viewer_state$zoom_level <- max(1, hne_viewer_state$zoom_level / input$hne_zoom_factor)
    showNotification(paste("H&E zoom:", round(hne_viewer_state$zoom_level, 1), "x"), 
                     type = "message", duration = 2)
  })
  
  observeEvent(input$hne_reset_zoom, {
    hne_viewer_state$zoom_level <- 1
    hne_viewer_state$center_x <- 0.5
    hne_viewer_state$center_y <- 0.5
    showNotification("H&E view reset", type = "message", duration = 2)
  })
  
  # Transcriptomics zoom controls
  observeEvent(input$trans_plot_click, {
    req(spatial_obj())
    
    click <- input$trans_plot_click
    zoom_factor <- input$trans_zoom_factor
    
    trans_viewer_state$center_x <- click$x
    trans_viewer_state$center_y <- click$y
    trans_viewer_state$zoom_level <- trans_viewer_state$zoom_level * zoom_factor
    
    showNotification(paste("Transcriptomics zoomed to", round(trans_viewer_state$zoom_level, 1), "x"), 
                     type = "message", duration = 2)
  })
  
  observeEvent(input$trans_zoom_in, {
    req(spatial_obj())
    trans_viewer_state$zoom_level <- trans_viewer_state$zoom_level * input$trans_zoom_factor
    showNotification(paste("Transcriptomics zoom:", round(trans_viewer_state$zoom_level, 1), "x"), 
                     type = "message", duration = 2)
  })
  
  observeEvent(input$trans_zoom_out, {
    req(spatial_obj())
    trans_viewer_state$zoom_level <- max(1, trans_viewer_state$zoom_level / input$trans_zoom_factor)
    showNotification(paste("Transcriptomics zoom:", round(trans_viewer_state$zoom_level, 1), "x"), 
                     type = "message", duration = 2)
  })
  
  observeEvent(input$trans_reset_zoom, {
    trans_viewer_state$zoom_level <- 1
    if (!is.null(spatial_obj())) {
      # Reset to full view
      coords <- GetTissueCoordinates(spatial_obj()@images[[1]])
      if (!is.null(coords)) {
        x_range <- range(coords[,1], na.rm = TRUE)
        y_range <- range(coords[,2], na.rm = TRUE)
        trans_viewer_state$center_x <- mean(x_range)
        trans_viewer_state$center_y <- mean(y_range)
      }
    }
    showNotification("Transcriptomics view reset", type = "message", duration = 2)
  })
  
  # Status info display
  output$viewer_status_info <- renderText({
    status_parts <- c()
    
    # H&E status
    if (!is.null(hne_viewer_state$image)) {
      img_dim <- dim(hne_viewer_state$image)
      status_parts <- c(status_parts, 
                        paste("H&E:", img_dim[2], "x", img_dim[1], "pixels"),
                        paste("H&E Zoom:", round(hne_viewer_state$zoom_level, 1), "x"))
    } else {
      status_parts <- c(status_parts, "H&E: Not loaded")
    }
    
    # Transcriptomics status
    if (!is.null(spatial_obj()) && processing_states$clustered) {
      n_spots <- ncol(spatial_obj())
      n_clusters <- length(unique(spatial_obj()$seurat_clusters))
      status_parts <- c(status_parts,
                        paste("Transcriptomics:", n_spots, "spots,", n_clusters, "clusters"),
                        paste("Trans Zoom:", round(trans_viewer_state$zoom_level, 1), "x"))
    } else {
      status_parts <- c(status_parts, "Transcriptomics: Not ready")
    }
    
    paste(status_parts, collapse = "\n")
  })
  ############################## Marker Analysis ##############################
  
  # Update cluster choices when clusters are available
  observe({
    req(spatial_obj())
    
    if (processing_states$clustered && "seurat_clusters" %in% colnames(spatial_obj()@meta.data)) {
      obj <- spatial_obj()
      
      # Get cluster choices - use annotated names if available
      if ("annotated_clusters" %in% colnames(obj@meta.data)) {
        cluster_choices <- levels(obj$annotated_clusters)
      } else {
        cluster_choices <- levels(obj$seurat_clusters)
      }
      
      updateSelectInput(session, "spatial_marker_cluster", 
                        choices = cluster_choices,
                        selected = cluster_choices[1])
      
      updateSelectInput(session, "spatial_comparison_clusters", 
                        choices = cluster_choices,
                        selected = NULL)
    }
  })
  
  # Update available clusters for comparison based on selected target
  observeEvent(input$spatial_marker_cluster, {
    req(spatial_obj(), input$spatial_marker_cluster)
    
    obj <- spatial_obj()
    
    # Get all clusters
    if ("annotated_clusters" %in% colnames(obj@meta.data)) {
      all_clusters <- levels(obj$annotated_clusters)
    } else {
      all_clusters <- levels(obj$seurat_clusters)
    }
    
    # Remove the selected cluster from comparison options
    available_clusters <- setdiff(all_clusters, input$spatial_marker_cluster)
    
    updateSelectInput(session, "spatial_comparison_clusters", 
                      choices = available_clusters,
                      selected = available_clusters[1])
  })
  
  # Store marker results in reactive value
  spatial_marker_results_reactive <- reactiveVal(NULL)
  
  # Dynamic title for results
  output$marker_results_title <- renderText({
    if (input$spatial_comparison_type == "one_vs_all") {
      paste("Markers for", input$spatial_marker_cluster, "vs All Other Clusters")
    } else {
      comparison_clusters <- input$spatial_comparison_clusters
      if (length(comparison_clusters) > 0) {
        paste("Markers for", input$spatial_marker_cluster, "vs", 
              paste(comparison_clusters, collapse = ", "))
      } else {
        "Select clusters to compare"
      }
    }
  })
  
  # Find markers for selected cluster
  observeEvent(input$run_spatial_markers, {
    message("=== FIND SPATIAL MARKERS STARTED ===")
    
    req(spatial_obj())
    req(input$spatial_marker_cluster)
    
    # Check if comparison clusters are selected when needed
    if (input$spatial_comparison_type == "one_vs_selected" && 
        (is.null(input$spatial_comparison_clusters) || length(input$spatial_comparison_clusters) == 0)) {
      showNotification("Please select at least one cluster to compare against", type = "warning")
      return()
    }
    
    tryCatch({
      obj <- spatial_obj()
      
      # Check if clusters exist
      cluster_col <- if ("annotated_clusters" %in% colnames(obj@meta.data)) {
        "annotated_clusters"
      } else if ("seurat_clusters" %in% colnames(obj@meta.data)) {
        "seurat_clusters"
      } else {
        showNotification("No clusters found! Please run clustering first.", type = "error")
        return()
      }
      
      selected_cluster <- input$spatial_marker_cluster
      comparison_type <- input$spatial_comparison_type
      
      message(paste("Finding markers for cluster:", selected_cluster))
      message(paste("Comparison type:", comparison_type))
      
      # Prepare modal message
      modal_message <- if (comparison_type == "one_vs_all") {
        paste("Analyzing cluster", selected_cluster, "vs all other clusters...")
      } else {
        paste("Analyzing cluster", selected_cluster, "vs clusters:", 
              paste(input$spatial_comparison_clusters, collapse = ", "))
      }
      
      showModal(modalDialog(
        title = "Finding Marker Genes",
        div(
          modal_message,
          br(),
          tags$small("This may take a few minutes")
        ),
        easyClose = FALSE,
        footer = NULL
      ))
      
      # Set identity to the cluster column
      Idents(obj) <- cluster_col
      
      # Find markers based on comparison type
      if (comparison_type == "one_vs_all") {
        # One vs all comparison
        markers <- FindMarkers(
          obj,
          ident.1 = selected_cluster,
          logfc.threshold = input$spatial_marker_logfc,
          min.pct = input$spatial_marker_min_pct,
          only.pos = input$spatial_marker_only_pos,
          verbose = FALSE
        )
      } else {
        # One vs selected comparison
        markers <- FindMarkers(
          obj,
          ident.1 = selected_cluster,
          ident.2 = input$spatial_comparison_clusters,
          logfc.threshold = input$spatial_marker_logfc,
          min.pct = input$spatial_marker_min_pct,
          only.pos = input$spatial_marker_only_pos,
          verbose = FALSE
        )
      }
      
      # Add gene names as a column
      markers$gene <- rownames(markers)
      markers <- markers[, c("gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
      
      # Sort by fold change
      markers <- markers[order(-markers$avg_log2FC), ]
      
      # Round numeric columns
      markers$avg_log2FC <- round(markers$avg_log2FC, 3)
      markers$pct.1 <- round(markers$pct.1, 3)
      markers$pct.2 <- round(markers$pct.2, 3)
      markers$p_val <- formatC(markers$p_val, format = "e", digits = 2)
      markers$p_val_adj <- formatC(markers$p_val_adj, format = "e", digits = 2)
      
      # Store results
      spatial_marker_results_reactive(markers)
      
      message(paste("Found", nrow(markers), "marker genes"))
      
      removeModal()
      
      
      message("=== FIND SPATIAL MARKERS COMPLETED ===")
      
    }, error = function(e) {
      message(paste("=== FIND SPATIAL MARKERS ERROR ===", e$message))
      removeModal()
      showNotification(paste("Error finding markers:", e$message), type = "error")
    })
  })
  
  
  # Display marker results
  output$spatial_marker_results <- renderDT({
    markers <- spatial_marker_results_reactive()
    
    if (is.null(markers)) {
      return(datatable(data.frame(Message = "Click 'Find Markers' to analyze a cluster"),
                       options = list(dom = 't'),
                       rownames = FALSE))
    }
    
    datatable(
      markers,
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel'),
        columnDefs = list(
          list(className = 'dt-center', targets = 1:5)
        )
      ),
      rownames = FALSE,
      caption = paste("Marker genes for cluster", input$spatial_marker_cluster)
    ) %>%
      formatStyle('avg_log2FC',
                  backgroundColor = styleInterval(c(0.5, 1), c('white', 'lightblue', 'darkblue')),
                  color = styleInterval(c(1), c('black', 'white')))
  })
  
  # Download marker results - CORRECTED VERSION
  output$download_spatial_markers <- downloadHandler(
    filename = function() {
      cluster_name <- gsub(" ", "_", input$spatial_marker_cluster)
      paste0("markers_cluster_", cluster_name, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      markers <- spatial_marker_results_reactive()
      
      # Simply write the file, no notifications here
      if (!is.null(markers)) {
        write.csv(markers, file, row.names = FALSE)
      } else {
        # Write empty file if no data
        write.csv(data.frame(Message = "No results available"), file, row.names = FALSE)
      }
    }
  )
}