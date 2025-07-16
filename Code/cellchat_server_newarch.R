############################## Script Cell Chat Analysis ##############################
# cellchat_server.R
cellchat_server <- function(input, output, session) {
  ############################## Loading Data ##############################

  # Initialize reactive values
  seurat_object_cellchat <- reactiveVal(NULL)
  db_cellchat <- reactiveVal(NULL)
  all_logs <- reactiveVal(character())  # Pour accumuler les logs
  db_loaded <- reactiveVal(FALSE)       # Pour tracker si DB déjà chargée

  # Loading modal function
  loadingModal <- function(text) {
    showModal(modalDialog(
      title = "Processing...",
      p(text),
      footer = NULL,
      easyClose = FALSE
    ))
  }

  # Add logs function
  add_logs <- function(new_logs) {
    current_logs <- all_logs()
    all_logs(c(current_logs, new_logs))
  }

  ############################## Load Seurat Object ##############################
  observeEvent(input$seurat_file_cellchat, {
    logs <- character()

    tryCatch({
      loadingModal("Loading Seurat object...")
      logs <- c(logs, "Starting Seurat object loading process...")

      seurat_data <- readRDS(input$seurat_file_cellchat$datapath)
      logs <- c(logs, sprintf("Seurat object loaded with %d cells", ncol(seurat_data)))
      seurat_object_cellchat(seurat_data)

      # Update outputs with safe printing
      output$seurat_info_cellchat <- renderPrint({
        cat(print_seurat_info(seurat_data))
      })

      add_logs(logs)

      removeModal()
    }, error = function(e) {
      logs <- c(logs, sprintf("ERROR loading Seurat: %s", e$message))
      add_logs(logs)
      removeModal()
      showModal(modalDialog(
        title = "Error",
        tags$div(
          p(sprintf("Error loading Seurat: %s", e$message))
        ),
        easyClose = TRUE
      ))
    })
  })

  ############################## Load Database ##############################
  observeEvent(input$species_cellchat, {
    db_loaded(FALSE)
  })


  observeEvent(input$load_db_cellchat, {
    logs <- character()
    tryCatch({
      loadingModal("Loading GaspouDB...")
      logs <- c(logs, paste("Loading GaspouDB for", input$species_cellchat))
      db_path <- if(input$species_cellchat == "mouse") {
        "databases/GaspouDB_mouse.rds"
      } else {
        "databases/GaspouDB_human.rds"
      }
      GaspouDB <- readRDS(db_path)
      db_cellchat(GaspouDB)
      db_loaded(TRUE)
      logs <- c(logs, "GaspouDB loaded successfully")
      logs <- c(logs, paste("Database contains", nrow(GaspouDB$interaction), "interactions"))
      output$db_info_cellchat <- renderPrint({
        print_db_info(GaspouDB)
      })
      add_logs(logs)
      removeModal()
    }, error = function(e) {
      logs <- c(logs, paste("ERROR loading database:", e$message))
      add_logs(logs)
      removeModal()
      showModal(modalDialog(
        title = "Error",
        tags$div(
          p(paste("Error loading GaspouDB:", e$message))
        ),
        easyClose = TRUE
      ))
    })
  })
  # Update the accumulated logs display
  output$loading_logs <- renderPrint({
    cat(paste(all_logs(), collapse = "\n"))
  })

  # Helper functions
  print_seurat_info <- function(seurat_obj) {
    # Use proper methods for Seurat object
    tryCatch({
      info <- list(
        ncells = ncol(seurat_obj),
        nfeatures = nrow(seurat_obj),
        assays = names(seurat_obj@assays)
      )

      return(sprintf(
        "Dataset loaded:\nNumber of cells: %d\nNumber of features: %d\nAvailable assays: %s",
        info$ncells,
        info$nfeatures,
        paste(info$assays, collapse = ", ")
      ))
    }, error = function(e) {
      return("Error reading Seurat object properties")
    })
  }

  print_db_info <- function(db) {
    paste(
      "Database loaded:",
      "\nNumber of interactions:", nrow(db$interaction),
      "\nSignaling pathways:", length(unique(db$interaction$pathway_name))
    )
  }


  ###################### TAB 2 ###################
  
  # Update choices when Seurat object is loaded
  observe({
    req(seurat_object_cellchat())
    
    tryCatch({
      seurat_obj <- seurat_object_cellchat()
      
      # Create cluster_names column if it doesn't exist
      if(!"cluster_names" %in% colnames(seurat_obj@meta.data)) {
        seurat_obj <- create_cluster_names_column(seurat_obj)
        seurat_object_cellchat(seurat_obj)  # Update reactive object
      }
      
      # Define patterns to exclude
      excluded_patterns <- c(
        "percent.mt", "nCount_ATAC", "nFeature_ATAC", "nFeature_RNA", "nCount_RNA",
        "^RNA_snn_res", "^RNA_snn", "^pANN", "^DF", "^integrated", "^integrated_snn"
      )
      
      pattern <- paste(excluded_patterns, collapse = "|")
      meta_cols <- colnames(seurat_obj@meta.data)
      filtered_cols <- meta_cols[!grepl(pattern, meta_cols)]
      
      # Update choices
      updateSelectInput(session, "group_by_cellchat", choices = filtered_cols, selected = "cluster_names")
      updateSelectInput(session, "subset_column_cellchat", choices = filtered_cols)
      
    }, error = function(e) {
      print(paste("Error updating choices:", e$message))
    })
  })
  
  # Create and analyze CellChat objects
  cellchat_objects <- reactiveVal(list())
  analysis_logs <- reactiveVal(character())
  
  add_analysis_log <- function(new_logs) {
    current_logs <- analysis_logs()
    analysis_logs(c(current_logs, new_logs))
  }
  
  # Function to create cluster names column
  create_cluster_names_column <- function(seurat_obj) {
    tryCatch({
      # If cluster_names already exists with real names, keep them
      if("cluster_names" %in% colnames(seurat_obj@meta.data)) {
        existing_names <- unique(seurat_obj@meta.data$cluster_names)
        
        # Check if names are meaningful (not just generic)
        if(!all(grepl("^Cluster_[0-9]+$", existing_names))) {
          message("Preserving existing cluster names")
          Idents(seurat_obj) <- "cluster_names"
          return(seurat_obj)
        }
      }
      
      # Create cluster names based on available data
      if("cell_type" %in% colnames(seurat_obj@meta.data)) {
        # Use cell_type if available
        seurat_obj$cluster_names <- as.character(seurat_obj$cell_type)
      } else if("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
        # Use seurat clusters
        seurat_obj$cluster_names <- paste0("Cluster_", seurat_obj$seurat_clusters)
      } else {
        # Use current identities
        seurat_obj$cluster_names <- as.character(Idents(seurat_obj))
      }
      
      # Set identities to cluster_names
      Idents(seurat_obj) <- "cluster_names"
      return(seurat_obj)
      
    }, error = function(e) {
      warning(paste("Could not create cluster names:", e$message))
      return(seurat_obj)
    })
  }
  
 
  # Function to clean CellChat object factors
  clean_cellchat_factors <- function(cellchat_obj) {
    tryCatch({
      message("=== Cleaning CellChat factors ===")
      
      # Clean identities
      if(is.factor(cellchat_obj@idents)) {
        cellchat_obj@idents <- droplevels(cellchat_obj@idents)
      }
      
      # Clean meta$labels
      if(!is.null(cellchat_obj@meta$labels)) {
        if(is.factor(cellchat_obj@meta$labels)) {
          cellchat_obj@meta$labels <- droplevels(cellchat_obj@meta$labels)
        }
      } else {
        # Create labels from idents if missing
        cellchat_obj@meta$labels <- cellchat_obj@idents
      }
      
      # Ensure consistency between idents and labels
      if(!identical(levels(cellchat_obj@idents), levels(cellchat_obj@meta$labels))) {
        message("Synchronizing idents and labels...")
        cellchat_obj@meta$labels <- factor(as.character(cellchat_obj@meta$labels), 
                                           levels = levels(cellchat_obj@idents))
      }
      
      # Clean all factor columns in meta
      for(col in names(cellchat_obj@meta)) {
        if(is.factor(cellchat_obj@meta[[col]])) {
          cellchat_obj@meta[[col]] <- droplevels(cellchat_obj@meta[[col]])
        }
      }
      
      message("Unique identities: ", length(levels(cellchat_obj@idents)))
      
      return(cellchat_obj)
      
    }, error = function(e) {
      stop(paste("Error cleaning CellChat factors:", e$message))
    })
  }
  
  # Function to prepare Seurat object for CellChat WITHOUT removing cells
  prepare_seurat_for_cellchat <- function(seurat_obj, group_by_column) {
    tryCatch({
      message("=== Preparing Seurat object for CellChat ===")
      message("Original cells: ", ncol(seurat_obj))
      message("Grouping by: ", group_by_column)
      
      # Verify grouping column exists
      if(!group_by_column %in% colnames(seurat_obj@meta.data)) {
        stop(paste("Column", group_by_column, "not found in metadata"))
      }
      
      # Get grouping values
      group_values <- seurat_obj@meta.data[[group_by_column]]
      
      # Replace NA/empty values with "Unknown" instead of removing cells
      group_values <- as.character(group_values)
      group_values[is.na(group_values)] <- "Unknown"
      group_values[group_values == ""] <- "Unknown"
      group_values[group_values == "NA"] <- "Unknown"
      
      # Get unique values including "Unknown"
      unique_values <- unique(group_values)
      
      # Update the metadata with cleaned values
      seurat_obj@meta.data[[group_by_column]] <- factor(group_values, levels = unique_values)
      
      # Add sample column if it doesn't exist (to avoid CellChat warning)
      if(!"samples" %in% colnames(seurat_obj@meta.data)) {
        seurat_obj@meta.data$samples <- "sample1"
      }
      
      # Clean all factor columns to have only present levels
      for(col in colnames(seurat_obj@meta.data)) {
        if(is.factor(seurat_obj@meta.data[[col]])) {
          seurat_obj@meta.data[[col]] <- droplevels(seurat_obj@meta.data[[col]])
        }
      }
      
      message("Final cells: ", ncol(seurat_obj), " (NO CELLS REMOVED)")
      message("Groups found: ", paste(unique_values, collapse = ", "))
      
      return(seurat_obj)
      
    }, error = function(e) {
      stop(paste("Error preparing Seurat object:", e$message))
    })
  }
  
  # Updated observer to create CellChat objects
  observeEvent(input$create_cellchat, {
    req(seurat_object_cellchat(), db_cellchat())
    
    tryCatch({
      loadingModal("Creating CellChat objects...")
      seurat_obj <- seurat_object_cellchat()
      
      # Ensure cluster_names column exists
      if(!"cluster_names" %in% colnames(seurat_obj@meta.data)) {
        seurat_obj <- create_cluster_names_column(seurat_obj)
        seurat_object_cellchat(seurat_obj)  # Update reactive
      }
      
      # Join layers if using Seurat v5
      if("RNA" %in% names(seurat_obj@assays)) {
        seurat_obj[["RNA"]] <- JoinLayers(seurat_obj[["RNA"]])
      }
      
      if(input$use_subset_cellchat && !is.null(input$subset_column_cellchat)) {
        # Create subset CellChat objects
        subset_column <- input$subset_column_cellchat
        unique_values <- unique(seurat_obj@meta.data[[subset_column]])
        unique_values <- unique_values[!is.na(unique_values)]
        
        add_analysis_log(paste("Creating subsets by", subset_column))
        add_analysis_log(paste("Found values:", paste(unique_values, collapse = ", ")))
        
        for(value in unique_values) {
          add_analysis_log(sprintf("\nProcessing %s = %s", subset_column, value))
          
          # Create subset - keeping ALL cells from this subset
          subset_mask <- seurat_obj@meta.data[[subset_column]] == value
          subset_data <- seurat_obj[, subset_mask]
          
          add_analysis_log(paste("Subset has", ncol(subset_data), "cells"))
          
          # Prepare the subset WITHOUT removing cells
          subset_data <- prepare_seurat_for_cellchat(subset_data, input$group_by_cellchat)
          
          # Check if enough groups remain
          n_groups <- length(unique(subset_data@meta.data[[input$group_by_cellchat]]))
          if(n_groups < 2) {
            add_analysis_log(paste("Skipping - need 2+ groups, found", n_groups))
            next
          }
          
          # Extract data for CellChat
          data_input <- GetAssayData(subset_data, assay = "RNA", slot = "data")
          meta_data <- subset_data@meta.data
          
          # Ensure dimensions match
          if(ncol(data_input) != nrow(meta_data)) {
            stop(paste("Dimension mismatch:", ncol(data_input), "vs", nrow(meta_data)))
          }
          
          # Create CellChat object
          cellchat_obj <- createCellChat(
            object = data_input,
            meta = meta_data,
            group.by = input$group_by_cellchat
          )
          
          # Clean CellChat factors
          cellchat_obj <- clean_cellchat_factors(cellchat_obj)
          
          # Store the object
          objects <- cellchat_objects()
          objects[[as.character(value)]] <- cellchat_obj
          cellchat_objects(objects)
          
          add_analysis_log(paste("Created CellChat object for", value, "with", n_groups, "groups"))
        }
        
      } else {
        # Create single CellChat object for entire dataset
        add_analysis_log("Creating CellChat object for complete dataset")
        
        # Prepare the Seurat object WITHOUT removing cells
        seurat_obj <- prepare_seurat_for_cellchat(seurat_obj, input$group_by_cellchat)
        
        # Check number of groups
        n_groups <- length(unique(seurat_obj@meta.data[[input$group_by_cellchat]]))
        if(n_groups < 2) {
          stop(paste("Need at least 2 groups for CellChat. Found:", n_groups))
        }
        
        add_analysis_log(paste("Found", n_groups, "groups"))
        
        # Extract data for CellChat
        data_input <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
        meta_data <- seurat_obj@meta.data
        
        # Ensure dimensions match
        if(ncol(data_input) != nrow(meta_data)) {
          stop(paste("Dimension mismatch:", ncol(data_input), "vs", nrow(meta_data)))
        }
        
        # Create CellChat object
        cellchat_obj <- createCellChat(
          object = data_input,
          meta = meta_data,
          group.by = input$group_by_cellchat
        )
        
        # Clean CellChat factors
        cellchat_obj <- clean_cellchat_factors(cellchat_obj)
        
        # Store the object
        objects <- cellchat_objects()
        objects[["cellchat"]] <- cellchat_obj
        cellchat_objects(objects)
        
        add_analysis_log(paste("Created CellChat object with", n_groups, "groups"))
      }
      
      # Update UI
      updateSelectInput(session, "select_objects_cellchat", 
                        choices = names(cellchat_objects()))
      
      removeModal()
      add_analysis_log("\nCellChat objects created successfully!")
      
    }, error = function(e) {
      add_analysis_log(sprintf("ERROR: %s", e$message))
      removeModal()
      showModal(modalDialog(
        title = "Error",
        p(e$message),
        easyClose = TRUE
      ))
    })
  })
  # Observer to analyze CellChat objects
  observeEvent(input$analyze_cellchat, {
    req(cellchat_objects(), input$select_objects_cellchat)
    
    tryCatch({
      showModal(modalDialog(
        title = "Please Wait",
        "Analyzing Ligand-Receptor communications...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      objects <- cellchat_objects()
      selected_objects <- input$select_objects_cellchat
      
      for(obj_name in selected_objects) {
        add_analysis_log(sprintf("\nAnalyzing %s", obj_name))
        cellchat_obj <- objects[[obj_name]]
        
        # Final factor cleaning before analysis
        cellchat_obj <- clean_cellchat_factors(cellchat_obj)
        
        # Verify object validity
        if(length(levels(cellchat_obj@idents)) < 2) {
          add_analysis_log(paste("ERROR: Object", obj_name, "has less than 2 groups"))
          next
        }
        
        # Set database
        add_analysis_log("Setting database...")
        cellchat_obj@DB <- db_cellchat()
        
        # Run CellChat pipeline
        add_analysis_log("Subsetting data...")
        cellchat_obj <- subsetData(cellchat_obj)
        
        add_analysis_log("Identifying over-expressed genes...")
        cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
        
        add_analysis_log("Identifying over-expressed interactions...")
        cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)
        
        add_analysis_log("Computing communication probability...")
        cellchat_obj <- computeCommunProb(cellchat_obj)
        
        add_analysis_log("Aggregating network...")
        cellchat_obj <- aggregateNet(cellchat_obj)
        
        # Update stored object
        objects[[obj_name]] <- cellchat_obj
        add_analysis_log(sprintf("Analysis completed for %s", obj_name))
      }
      
      # Update reactive
      cellchat_objects(objects)
      
      removeModal()
      add_analysis_log("\nAll analyses completed successfully!")
      
    }, error = function(e) {
      removeModal()
      add_analysis_log(sprintf("ERROR in analysis: %s", e$message))
      showModal(modalDialog(
        title = "Error",
        p(sprintf("Error: %s", e$message)),
        easyClose = TRUE
      ))
    })
  })
  # Display analysis logs
  output$analysis_logs_cellchat <- renderPrint({
    cat(paste(analysis_logs(), collapse = "\n"))
  })





  observe({
    req(cellchat_objects())
    req(length(cellchat_objects()) > 0)
    print("Starting initialization...")

    object_choices <- names(cellchat_objects())

    # Mise à jour unique de tous les sélecteurs d'objets
    updateSelectInput(session, "subset_choose_cellchat", choices = object_choices)
    updateSelectInput(session, "cellchat_obj_chord", choices = object_choices)
    updateSelectInput(session, "cellchat_obj_global", choices = object_choices)
    updateSelectInput(session, "cellchat_obj_specific", choices = object_choices)
  })

  # Observer pour mettre à jour les listes déroulantes
  observe({
    req(input$subset_choose_cellchat)
    print("Updating bubble plot inputs...")

    tryCatch({
      cellchat_obj <- cellchat_objects()[[input$subset_choose_cellchat]]
      req(!is.null(cellchat_obj@idents))
      cell_types <- levels(cellchat_obj@idents)

      # Préserver les sélections existantes
      selectedSources <- isolate(input$sources_use_cellchat)
      selectedTargets <- isolate(input$targets_use_cellchat)

      # Mise à jour des cell types
      updateSelectizeInput(session, "sources_use_cellchat",
                           choices = cell_types,
                           selected = intersect(selectedSources, cell_types),
                           options = list(plugins = list('remove_button'))
      )
      updateSelectizeInput(session, "targets_use_cellchat",
                           choices = cell_types,
                           selected = intersect(selectedTargets, cell_types),
                           options = list(plugins = list('remove_button'))
      )

      # Mise à jour des LR et pathways
      if(!is.null(cellchat_obj@net) && !is.null(cellchat_obj@net$prob)) {
        # Get L-R interactions
        lr_interactions <- NULL
        tryCatch({
          communications <- subsetCommunication(cellchat_obj)
          lr_interactions <- unique(communications$interaction_name)
        }, error = function(e) {
          message("Error getting L-R interactions: ", e$message)
          # Fallback method
          if (!is.null(cellchat_obj@LR$LRsig)) {
            lr_interactions <- unique(cellchat_obj@LR$LRsig$interaction_name)
          }
        })

        # Get pathways
        pathways <- NULL
        tryCatch({
          pathway_data <- subsetCommunication(cellchat_obj, slot.name = "netP")
          pathways <- unique(pathway_data$pathway_name)
        }, error = function(e) {
          message("Error getting pathways: ", e$message)
          # Fallback method
          if (!is.null(cellchat_obj@netP$pathways)) {
            pathways <- unique(unlist(lapply(cellchat_obj@netP$pathways, function(x) names(x))))
          }
        })

        # Update UI dropdowns
        if (!is.null(lr_interactions)) {
          updateSelectizeInput(session, "selected_lr_cellchat",
                               choices = lr_interactions,
                               selected = NULL
          )
        }

        if (!is.null(pathways)) {
          updateSelectizeInput(session, "signaling_pathways_bubble",
                               choices = pathways,
                               selected = NULL
          )
        }
      }
    }, error = function(e) {
      print(paste("Error updating bubble plot inputs:", e$message))
      showNotification(paste("Error updating inputs:", e$message), type = "error")
    })
  })

  # Reactive value for plot
  bubble_plot <- reactiveVal()

  


  # Function to extract types only for displayed interactions
  get_interaction_type_info <- function(cellchat_obj, communications) {
    # Check if the database contains the interaction_type field
    if (!is.null(cellchat_obj@DB$interaction) &&
        "interaction_type" %in% colnames(cellchat_obj@DB$interaction) &&
        "interaction_name" %in% colnames(cellchat_obj@DB$interaction)) {

      # Get only the displayed interaction names from the filtered communications
      displayed_interactions <- unique(communications$interaction_name)

      # Get all interaction types from the database
      all_interaction_info <- data.frame(
        interaction_name = cellchat_obj@DB$interaction$interaction_name,
        interaction_type = cellchat_obj@DB$interaction$interaction_type,
        stringsAsFactors = FALSE
      )

      # Filter to only the displayed interactions
      interaction_info <- all_interaction_info[all_interaction_info$interaction_name %in% displayed_interactions, ]

      # Handle any NA values
      interaction_info$interaction_type[is.na(interaction_info$interaction_type)] <- "Unknown"

      # Sort by type and name
      interaction_info <- interaction_info[order(interaction_info$interaction_type,
                                                 interaction_info$interaction_name), ]

      return(interaction_info)
    }

    return(NULL)
  }


  # Modification of your observer to fix the issue
  observeEvent(input$generate_plot_cellchat, {
    req(input$subset_choose_cellchat)

    # Show loading indicator
    showNotification("Generating bubble plot...", type = "message", duration = 3)

    tryCatch({
      cellchat_obj <- cellchat_objects()[[input$subset_choose_cellchat]]

      # Check for required selections
      if (length(input$sources_use_cellchat) == 0 || length(input$targets_use_cellchat) == 0) {
        showNotification("Please select at least one source and one target cell type", type = "error")
        return()
      }

      # Generate cellchat communications
      communications <- subsetCommunication(cellchat_obj)

      # Filter by sources and targets
      communications <- communications[communications$source %in% input$sources_use_cellchat &
                                         communications$target %in% input$targets_use_cellchat, ]

      # Filter by threshold
      communications <- communications[communications$prob >= input$threshold_cellchat, ]

      # Exclude intra-cellular communications if requested
      if (input$exclude_intra) {
        communications <- communications[communications$source != communications$target, ]
      }

      # If specific L-R pairs are selected, use them
      pairLR <- NULL
      if (!is.null(input$selected_lr_cellchat) && length(input$selected_lr_cellchat) > 0) {
        pairLR <- data.frame(
          interaction_name = input$selected_lr_cellchat,
          stringsAsFactors = FALSE
        )

        # Also filter communications to match pairLR
        communications <- communications[communications$interaction_name %in% pairLR$interaction_name, ]
      }

      # If specific signaling pathways are selected, use them
      if (!is.null(input$signaling_pathways_bubble) && length(input$signaling_pathways_bubble) > 0) {
        pathway_comms <- subsetCommunication(cellchat_obj, signaling = input$signaling_pathways_bubble)
        pathway_lr <- unique(pathway_comms$interaction_name)

        if (is.null(pairLR)) {
          pairLR <- data.frame(
            interaction_name = pathway_lr,
            stringsAsFactors = FALSE
          )

          # Filter communications to match pathway LR
          communications <- communications[communications$interaction_name %in% pathway_lr, ]
        } else {
          # Intersection if both are specified
          pairLR <- pairLR[pairLR$interaction_name %in% pathway_lr, ]

          # Update communications filtering
          communications <- communications[communications$interaction_name %in% pairLR$interaction_name, ]
        }
      }

      # Create the bubble plot with correct options
      plot_result <- list()

      # Keep track of actual displayed LR pairs
      displayed_lr_names <- unique(communications$interaction_name)
      lr_count <- length(displayed_lr_names)

      # Use netVisual_bubble with pairLR.use if specified
      base_plot <- netVisual_bubble(
        object = cellchat_obj,
        sources.use = input$sources_use_cellchat,
        targets.use = input$targets_use_cellchat,
        pairLR.use = pairLR,  # This option controls the display of L-R pairs
        thresh = input$threshold_cellchat,
        title = input$subset_choose_cellchat,
        remove.isolate = input$exclude_intra,
        angle.x = 45,
        font.size = 10
      )

      # Apply coord_flip if necessary
      if (input$flip_axes) {
        plot_result$plot <- base_plot +
          coord_flip() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            axis.text.y = element_text(size = 10)
          )
      } else {
        plot_result$plot <- base_plot +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            axis.text.y = element_text(size = 10)
          )
      }

      # Get interaction type information - only for the displayed interactions
      interaction_info <- get_interaction_type_info(cellchat_obj, communications)
      if (!is.null(interaction_info)) {
        plot_result$interaction_info <- interaction_info

        # Update the count based on the filtered results
        lr_count <- nrow(interaction_info)
      }

      # Store the plot and info
      bubble_plot(plot_result)

      # Show notification with accurate count
      if (lr_count > 0) {
        # Count by type if available
        if (!is.null(interaction_info) && "interaction_type" %in% colnames(interaction_info)) {
          type_counts <- table(interaction_info$interaction_type)
          type_summary <- paste(names(type_counts), type_counts, sep=": ", collapse=", ")
          showNotification(paste("Plot generated with", lr_count, "ligand-receptor pairs.",
                                 "By type:", type_summary),
                           type = "message", duration = 5)
        } else {
          showNotification(paste("Plot generated with", lr_count, "ligand-receptor pairs."),
                           type = "message")
        }
      } else {
        showNotification("Warning: No ligand-receptor pairs found with current filters.",
                         type = "warning")
      }
    }, error = function(e) {
      message("Processing error: ", e$message)
      showNotification(paste("Error processing data:", e$message), type = "error")
    })
  })
      # Get interaction type information
      output$interaction_info_table <- renderTable({
        req(bubble_plot())

        info_df <- bubble_plot()$interaction_info

        if (is.null(info_df) || nrow(info_df) == 0) {
          return(data.frame(
            Message = "No interaction type information available."
          ))
        }

        # Vérifier que les deux colonnes existent bien
        if (!("interaction_name" %in% colnames(info_df)) ||
            !("interaction_type" %in% colnames(info_df))) {
          return(data.frame(
            Message = "Expected columns not found in the data."
          ))
        }

        # Format the table for display - IMPORTANT: inclure explicitement les deux colonnes
        display_table <- data.frame(
          Interaction = info_df$interaction_name,
          Type = info_df$interaction_type
        )

        return(display_table)
      }, striped = TRUE, bordered = TRUE, hover = TRUE, spacing = 'xs')

  # Render the plot without infinite refresh
  output$bubble_plot_cellchat <- renderPlot({
    req(bubble_plot())

    # Get the plot
    plot_obj <- bubble_plot()$plot

    # Return the plot
    plot_obj
  }, width = function() {
    # Use a default value if nothing is specified
    if(is.null(input$plot_width) || is.na(as.numeric(input$plot_width))) {
      return(1000)  # Default value
    }
    return(as.numeric(input$plot_width))
  }, height = function() {
    # Use a default value if nothing is specified
    if(is.null(input$plot_height) || is.na(as.numeric(input$plot_height))) {
      return(800)  # Default value
    }
    return(as.numeric(input$plot_height))
  })

  # Render the interaction info table
  output$interaction_info_table <- renderTable({
    req(bubble_plot())

    info_df <- bubble_plot()$interaction_info

    if (is.null(info_df) || nrow(info_df) == 0) {
      return(data.frame(
        Message = "No interaction type information available."
      ))
    }

    # Format the table for display
    display_table <- info_df %>%
      dplyr::rename(
        `Interaction` = interaction_name,
        `Type` = interaction_type
      )

    return(display_table)
  }, striped = TRUE, bordered = TRUE, hover = TRUE, spacing = 'xs')

  # Download handler for the plot
  output$download_bubble_plot <- downloadHandler(
    filename = function() {
      obj_name <- ifelse(!is.null(input$subset_choose_cellchat), input$subset_choose_cellchat, "cellchat")
      group_info <- if (!is.null(input$sources_use_cellchat) && !is.null(input$targets_use_cellchat)) {
        paste0("_", paste(head(input$sources_use_cellchat, 2), collapse="_"),
               "_to_", paste(head(input$targets_use_cellchat, 2), collapse="_"))
      } else ""
      paste0("BubblePlot_", obj_name, group_info, "_", format(Sys.time(), "%Y%m%d"), ".", input$bubble_plot_format)
    },
    content = function(file) {
      req(bubble_plot())

      tryCatch({
        # Get the plot
        p <- bubble_plot()$plot

        # Save with appropriate dimensions
        ggsave(
          file,
          plot = p,
          width = input$plot_width/72,
          height = input$plot_height/72,
          dpi = 300,
          device = input$bubble_plot_format
        )

        showNotification("Plot downloaded successfully!", type = "message")
      }, error = function(e) {
        message("Download error: ", e$message)
        showNotification(paste("Error saving plot:", e$message), type = "error")
      })
    }
  )

  # Download handler for interaction info
  output$download_info_table <- downloadHandler(
    filename = function() {
      obj_name <- ifelse(!is.null(input$subset_choose_cellchat), input$subset_choose_cellchat, "cellchat")
      paste0("InteractionTypes_", obj_name, "_", format(Sys.time(), "%Y%m%d"), ".csv")
    },
    content = function(file) {
      req(bubble_plot())

      info_df <- bubble_plot()$interaction_info

      if (!is.null(info_df) && nrow(info_df) > 0) {
        # Save as CSV
        write.csv(info_df, file, row.names = FALSE)
        showNotification("Interaction type information downloaded successfully!", type = "message")
      } else {
        # Create empty file with message
        write.csv(data.frame(
          Message = "No interaction type information available."
        ), file, row.names = FALSE)
        showNotification("Empty information table downloaded.", type = "warning")
      }
    }
  )

#################Tab 3 : Circle plot##############
  observe({
    req(input$cellchat_obj_global)
    print("Updating circle plot inputs...")
    tryCatch({
      cellchat_obj <- cellchat_objects()[[input$cellchat_obj_global]]
      req(!is.null(cellchat_obj@idents))
      cell_types <- levels(cellchat_obj@idents)
      updateSelectizeInput(session, "cell_types_to_show", choices = cell_types)
      output$cellchat_obj_title_global <- renderText({
        paste("Selected object:", input$cellchat_obj_global)
      })
    }, error = function(e) {
      print(paste("Error updating circle plot inputs:", e$message))
    })
  })

  # Créer des reactiveVal pour stocker les plots
  global_circle_plot <- reactiveVal()
  cell_type_plots <- reactiveVal()
  FONT_CONFIG <- list(
    base = list(
      display = 72
    )
  )

  # Fonction pour calculer la taille de police
  get_font_scaling <- function(dpi, for_export = FALSE) {
    if(for_export) {
      return(FONT_CONFIG$base$display / dpi)
    }
    return(1)
  }

  # Fonction pour créer le plot global
  create_global_plot <- function(cellchat_obj, plot_type, font_size = 0.5, margin_size = 0.2,
                                 title = NULL, for_export = FALSE, dpi = NULL) {
    groupSize <- as.numeric(table(cellchat_obj@idents))
    font_scale <- get_font_scaling(dpi, for_export)
    adjusted_font_size <- font_size * font_scale
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    par(mar = c(2, 2, 4, 2))
    par(cex = adjusted_font_size, cex.main = adjusted_font_size * 1.2)
    netVisual_circle(
      if(plot_type == "count") cellchat_obj@net$count else cellchat_obj@net$weight,
      vertex.weight = groupSize,
      weight.scale = TRUE,
      label.edge = FALSE,
      title.name = ifelse(!is.null(title),
                          paste(title, ifelse(plot_type == "count", "- Count", "- Weight")),
                          ifelse(plot_type == "count",
                                 "Number of interactions",
                                 "Interaction weights/strength")),
      margin = margin_size
    )
  }

  create_cell_type_plots <- function(cellchat_obj, cell_types, font_size = 0.5, margin_size = 0.15, title = NULL) {
    groupSize <- as.numeric(table(cellchat_obj@idents))
    mat <- cellchat_obj@net$weight

    n_plots <- length(cell_types)
    if(n_plots == 1) {
      par(mfrow=c(1,1), mar=c(2, 2, 4, 2), cex = font_size * 0.8, cex.main = font_size)
      mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
      mat2[cell_types, ] <- mat[cell_types, ]
      netVisual_circle(mat2,
                       vertex.weight = groupSize,
                       weight.scale = TRUE,
                       edge.weight.max = max(mat),
                       title.name = ifelse(!is.null(title), paste(title, "-", cell_types), cell_types),
                       margin = margin_size)
    } else {
      grid_dim <- ceiling(sqrt(n_plots))
      layout(matrix(1:(grid_dim*grid_dim), grid_dim, grid_dim))
      for(cell_type in cell_types) {
        par(mar=c(1, 1, 3, 1), oma=c(2, 2, 2, 2), cex = font_size * 0.6, cex.main = font_size * 0.8)
        mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
        mat2[cell_type, ] <- mat[cell_type, ]
        netVisual_circle(mat2,
                         vertex.weight = groupSize,
                         weight.scale = TRUE,
                         edge.weight.max = max(mat),
                         title.name = ifelse(!is.null(title), paste(title, "-", cell_type), cell_type),
                         margin = margin_size)
      }
    }
  }
  # Global circle plot with improved layout
  observeEvent(input$generate_global_plot, {
    req(input$cellchat_obj_global, input$plot_type_global)
    
    tryCatch({
      cellchat_obj <- cellchat_objects()[[input$cellchat_obj_global]]
      
      # Store plot function
      global_plot_function <- function() {
        groupSize <- as.numeric(table(cellchat_obj@idents))
        
        netVisual_circle(
          if(input$plot_type_global == "count") cellchat_obj@net$count else cellchat_obj@net$weight,
          vertex.weight = groupSize,
          weight.scale = TRUE,
          vertex.size = input$global_vertex_size,
          edge.width.max = input$global_edge_width,
          vertex.label.cex = input$global_label_size,
          label.edge = FALSE,
          title.name = paste(input$cellchat_obj_global, "-",
                             ifelse(input$plot_type_global == "count", 
                                    "Interaction Count", 
                                    "Interaction Strength")),
          margin = input$global_margin
        )
      }
      
      global_circle_plot(global_plot_function)
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # Dynamic UI for global plot
  output$global_plot_ui <- renderUI({
    if(!is.null(global_circle_plot())) {
      plotOutput("global_plot_display", 
                 height = paste0(input$global_plot_size, "px"),
                 width = paste0(input$global_plot_size, "px"))
    } else {
      div(class = "text-center text-muted", 
          style = "padding: 50px;",
          h4("No plot generated yet"))
    }
  })
  
  # Render global plot
  output$global_plot_display <- renderPlot({
    req(global_circle_plot())
    par(mar = c(0, 0, 2, 0))
    global_circle_plot()()
  })
  
  # Cell type specific plots - FIXED
  observeEvent(input$generate_specific_plots, {
    req(input$cellchat_obj_specific, input$cell_types_to_show)
    
    tryCatch({
      cellchat_obj <- cellchat_objects()[[input$cellchat_obj_specific]]
      groupSize <- as.numeric(table(cellchat_obj@idents))
      
      # Store plot function
      specific_plots_function <- function() {
        n_plots <- length(input$cell_types_to_show)
        n_cols <- min(input$n_cols_specific, n_plots)
        n_rows <- ceiling(n_plots / n_cols)
        
        # Set up layout
        par(mfrow = c(n_rows, n_cols), 
            mar = c(1, 1, 3, 1), 
            oma = c(0, 0, 2, 0))
        
        # Get the appropriate matrix
        mat <- if(input$signal_direction == "outgoing") {
          cellchat_obj@net$weight
        } else {
          t(cellchat_obj@net$weight)  # Transpose for incoming
        }
        
        # Plot for each cell type
        for(cell_type in input$cell_types_to_show) {
          # Create matrix with only selected cell type's communications
          mat_subset <- matrix(0, 
                               nrow = nrow(mat), 
                               ncol = ncol(mat), 
                               dimnames = dimnames(mat))
          
          if(cell_type %in% rownames(mat)) {
            mat_subset[cell_type, ] <- mat[cell_type, ]
          }
          
          # Make the plot
          netVisual_circle(
            mat_subset,
            vertex.weight = groupSize,
            weight.scale = TRUE,
            edge.width.max = max(mat),
            vertex.label.cex = input$specific_label_size,
            title.name = paste(cell_type, 
                               ifelse(input$signal_direction == "outgoing",
                                      "(Outgoing)", "(Incoming)")),
            margin = input$specific_margin
          )
        }
      }
      
      cell_type_plots(specific_plots_function)
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # Dynamic UI for specific plots
  output$specific_plots_ui <- renderUI({
    if(!is.null(cell_type_plots())) {
      n_plots <- length(input$cell_types_to_show)
      n_cols <- min(input$n_cols_specific, n_plots)
      n_rows <- ceiling(n_plots / n_cols)
      
      total_height <- input$specific_plot_height * n_rows
      
      plotOutput("specific_plots_display", 
                 height = paste0(total_height, "px"))
    } else {
      div(class = "text-center text-muted", 
          style = "padding: 50px;",
          h4("No plots generated yet"))
    }
  })
  
  # Render specific plots
  output$specific_plots_display <- renderPlot({
    req(cell_type_plots())
    cell_type_plots()()
  })
  
  # Download handlers
  output$download_global_plot <- downloadHandler(
    filename = function() {
      paste0("GlobalCirclePlot_", input$cellchat_obj_global, "_", 
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      png(file, width = 2400, height = 2400, res = 300)
      par(mar = c(0, 0, 2, 0))
      global_circle_plot()()
      dev.off()
    }
  )
  
  output$download_specific_plots <- downloadHandler(
    filename = function() {
      paste0("CellTypeSpecificPlots_", input$cellchat_obj_specific, "_", 
             format(Sys.time(), "%Y%m%d_%H%M%S"), ".png")
    },
    content = function(file) {
      n_plots <- length(input$cell_types_to_show)
      n_cols <- min(input$n_cols_specific, n_plots)
      n_rows <- ceiling(n_plots / n_cols)
      
      png(file, 
          width = 800 * n_cols, 
          height = 800 * n_rows, 
          res = 300)
      cell_type_plots()()
      dev.off()
    }
  )

  # Reactive value pour stocker le plot
  chord_plot <- reactiveVal()

  # Observer pour mettre à jour les inputs du chord plot
  observe({
    req(input$cellchat_obj_chord)
    print("Updating chord plot inputs...")
    tryCatch({
      cellchat_obj <- cellchat_objects()[[input$cellchat_obj_chord]]
      req(!is.null(cellchat_obj@idents))
      cell_types <- levels(cellchat_obj@idents)
      
      # Update selectizeInput instead of selectInput
      updateSelectizeInput(session, "sender_groups", 
                           choices = cell_types,
                           selected = NULL)
      updateSelectizeInput(session, "receiver_groups", 
                           choices = cell_types,
                           selected = NULL)
      
      if(!is.null(cellchat_obj@net) && !is.null(cellchat_obj@net$prob)) {
        all_comm <- subsetCommunication(cellchat_obj)
        if(!is.null(all_comm) && nrow(all_comm) > 0) {
          updateSelectizeInput(session, "selected_lr_chord",
                               choices = unique(all_comm$interaction_name))
        }
      }
      
      output$cellchat_obj_title_chord <- renderText({
        paste("Selected object:", input$cellchat_obj_chord)
      })
    }, error = function(e) {
      print(paste("Error updating chord plot inputs:", e$message))
    })
  })
  # Observer pour mettre à jour cell_types_to_show dans l'onglet specific
  observe({
    req(input$cellchat_obj_specific)
    print("Updating specific plot inputs...")
    tryCatch({
      cellchat_obj <- cellchat_objects()[[input$cellchat_obj_specific]]
      req(!is.null(cellchat_obj@idents))
      cell_types <- levels(cellchat_obj@idents)
      
      updateSelectizeInput(session, "cell_types_to_show", 
                           choices = cell_types,
                           selected = NULL)
      
      output$cellchat_obj_title_specific <- renderText({
        paste("Selected object:", input$cellchat_obj_specific)
      })
    }, error = function(e) {
      print(paste("Error updating specific plot inputs:", e$message))
    })
  })
  
  # UI dynamique pour les couleurs - version corrigée
  output$color_inputs_chord <- renderUI({
    tryCatch({
      req(input$cellchat_obj_chord)
      req(input$use_custom_colors)
      req(input$sender_groups)
      req(input$receiver_groups)

      # Groupes sélectionnés par l'utilisateur
      all_groups <- unique(c(input$sender_groups, input$receiver_groups))
      req(length(all_groups) > 0)

      # Générer les inputs de couleur
      color_inputs <- lapply(seq_along(all_groups), function(i) {
        group_name <- all_groups[i]
        colourInput(
          inputId = paste0("color_group_", i),
          label = paste("Color for", group_name),
          value = rainbow(length(all_groups))[i]
        )
      })

      # Retourner les inputs dans un div
      div(class = "color-inputs", color_inputs)

    }, error = function(e) {
      print(paste("Error in color inputs generation:", e$message))
      return(NULL)
    })
  })


  # Corrected chord plot generation function
  observeEvent(input$generate_chord_plot, {
    tryCatch({
      req(input$cellchat_obj_chord)
      req(length(input$sender_groups) > 0)
      req(length(input$receiver_groups) > 0)
      
      showNotification("Generating chord plot...", type = "message", duration = 2)
      
      cellchat_obj <- cellchat_objects()[[input$cellchat_obj_chord]]
      
      # Get custom colors if enabled
      custom_colors <- NULL
      if(isTRUE(input$use_custom_colors)) {
        all_groups <- unique(c(input$sender_groups, input$receiver_groups))
        custom_colors <- sapply(seq_along(all_groups), function(i) {
          color_input <- input[[paste0("color_group_", i)]]
          if(is.null(color_input)) return(scales::hue_pal()(length(all_groups))[i])
          return(color_input)
        })
        names(custom_colors) <- all_groups
      }
      
      # Store plot function for reuse
      chord_plot_function <- function() {
        # Use netVisual_chord_gene with correct parameters
        netVisual_chord_gene(
          cellchat_obj,
          sources.use = input$sender_groups,
          targets.use = input$receiver_groups,
          thresh = input$prob_threshold_chord,
          color.use = custom_colors,
          show.legend = TRUE,
          small.gap = 3,
          big.gap = 10,
          lab.cex = input$chord_label_size,
          # Remove problematic parameters
          directional = 1,
          link.target.prop = FALSE
          # Don't use annotationTrack and annotationTrackHeight with netVisual_chord_gene
        )
      }
      
      # Store the plot function
      chord_plot(chord_plot_function)
      
      # Get interaction count for notification
      communications <- tryCatch({
        subsetCommunication(cellchat_obj,
                            sources.use = input$sender_groups,
                            targets.use = input$receiver_groups,
                            thresh = input$prob_threshold_chord)
      }, error = function(e) {
        NULL
      })
      
      n_interactions <- if(!is.null(communications)) nrow(communications) else "unknown"
      
      # Show notification
      showNotification(
        paste("Chord plot generated with", n_interactions, "interactions"),
        type = "message"
      )
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
      print(paste("Chord plot error details:", e$message))
    })
  })
  
  # Dynamic UI for chord plot with adjustable size
  output$chord_plot_ui <- renderUI({
    if(!is.null(chord_plot())) {
      plotOutput("chord_plot_display", 
                 height = paste0(input$chord_plot_dims, "px"),
                 width = paste0(input$chord_plot_dims, "px"))
    } else {
      div(class = "text-center text-muted", 
          style = "padding: 50px;",
          h4("No plot generated yet"),
          p("Select parameters and click 'Generate Plot'"))
    }
  })
  
  # Render the chord plot
  output$chord_plot_display <- renderPlot({
    req(chord_plot())
    par(mar = c(0, 0, 0, 0))
    chord_plot()()  # Execute the stored function
  })
  
  # Improved download handler with better spacing
  output$download_chord_plot <- downloadHandler(
    filename = function() {
      obj_name <- input$cellchat_obj_chord
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      paste0("ChordPlot_", obj_name, "_", timestamp, ".", input$chord_plot_format)
    },
    content = function(file) {
      tryCatch({
        width_px <- input$chord_export_width
        height_px <- width_px  # Square plot
        
        if (input$chord_plot_format == "pdf") {
          pdf(file, width = width_px/72, height = height_px/72)
        } else if (input$chord_plot_format == "svg") {
          svg(file, width = width_px/72, height = height_px/72)
        } else {
          # PNG, JPEG, TIFF
          switch(input$chord_plot_format,
                 "png" = png(file, width = width_px, height = height_px, res = 300),
                 "jpeg" = jpeg(file, width = width_px, height = height_px, res = 300, quality = 100),
                 "tiff" = tiff(file, width = width_px, height = height_px, res = 300))
        }
        
        # Use larger margins for export to prevent label cutoff
        par(mar = c(2, 2, 2, 2))
        chord_plot()()  # Execute the stored function
        
        dev.off()
        showNotification("Plot downloaded successfully!", type = "message")
        
      }, error = function(e) {
        showNotification(paste("Download error:", e$message), type = "error")
      })
    }
  )
  
  
}
