
############################## Script Multiple Datasets Analysis ##############################

# multiple_datasets_server.R

multiple_datasets_server <-  function(input, output, session) {
  
  
  # Dans multiple_datasets_server_newarch.R - Section Loading Data
  
  ############################## Loading Data ##############################
  
  # Variables for the merge dataset part
  multiple_datasets_object <- reactiveVal()
  seurat_objects <- reactiveValues()
  data_loaded <- reactiveValues()
  merged_gene_tables <- reactiveValues()
  rv_metadata <- reactiveValues(num_fields = 1)
  
  shinyjs::disable("add_field")
  shinyjs::disable("add_metadata")
  
  output$datasets_loaded <- reactive({
    !is.null(multiple_datasets_object())
  })
  
  outputOptions(output, 'datasets_loaded', suspendWhenHidden = FALSE)
  
  # Load pre-integrated Seurat object
  observeEvent(input$load_seurat_file_merge, {
    tryCatch({
      loaded_seurat <- loadSeuratObject(
        rds_path = input$load_seurat_file_merge$datapath,
        add_dataset_column = TRUE,
        dataset_name = tools::file_path_sans_ext(basename(input$load_seurat_file_merge$name)),
        clean_before = TRUE,
        module_type = "multiple"
      )
      
      plot_result <- generateInitialPlot(
        seurat_obj = loaded_seurat,
        remove_axes = FALSE,
        remove_legend = FALSE
      )
      
      if (!is.null(plot_result)) {
        multiple_datasets_object(plot_result$seurat_obj)
      } else {
        multiple_datasets_object(loaded_seurat)
      }
      
      shinyjs::enable(c("add_field", "add_metadata", "runScalePCA"))
      showNotification("Seurat object loaded successfully!")
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # Processing individual datasets using modular functions
  observe_file_input <- function(index) {
    observeEvent(input[[paste0("merge", index)]], {
      # Check if dataset already processed
      if (!is.null(data_loaded[[paste0("loaded", index)]]) && 
          data_loaded[[paste0("loaded", index)]] == TRUE) {
        message(paste("Dataset", index, "already processed, skipping..."))
        return()
      }
      
      req(input[[paste0("merge", index)]])
      
      withProgress(message = paste("Processing dataset", index), value = 0, {
        tryCatch({
          # Extract parameters from UI
          qc_params <- extractQCParams(input)
          dataset_name <- input[[paste0("dataset_name", index)]]
          if (is.null(dataset_name) || dataset_name == "") {
            dataset_name <- paste0("Dataset_", index)
          }
          dataset_type <- input[[paste0("dataset_type_merge", index)]]
          
          # Use modular preprocessing function
          seurat_object <- preprocessRawDataset(
            file_path = input[[paste0("merge", index)]]$datapath,
            dataset_type = dataset_type,
            species = input$species_choice_merge,
            dataset_name = dataset_name,
            qc_params = qc_params
          )
          
          # Store processed object
          seurat_objects[[paste0("seurat_object", index)]] <- seurat_object
          data_loaded[[paste0("loaded", index)]] <- TRUE
          
          message(paste("Dataset", index, "processed and stored successfully"))
          showNotification(paste("Dataset", index, "processed successfully!"), type = "message")
          
          # Enable metadata controls
          shinyjs::enable("add_field")
          shinyjs::enable("add_metadata")
          
        }, error = function(e) {
          message(paste("Error in dataset", index, ":", e$message))
          showNotification(paste("Error processing dataset", index, ":", e$message), type = "error")
        })
      })
    }, ignoreInit = TRUE, once = TRUE)
  }
  
  # Simple merge using modular functions
  observeEvent(input$simple_merge, {
    tryCatch({
      showModal(modalDialog(
        title = "Processing Datasets",
        "Merging datasets without integration...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      # Collect all processed Seurat objects
      seurat_list <- list()
      for (i in 1:input$num_datasets) {
        seurat_obj <- seurat_objects[[paste0("seurat_object", i)]]
        if (is.null(seurat_obj)) {
          stop(paste("Dataset", i, "is not processed yet"))
        }
        seurat_list[[i]] <- seurat_obj
      }
      
      # Use modular integration function
      merged_object <- performDataIntegration(
        seurat_list = seurat_list,
        integration_method = "simple"
      )
      
      # Store integrated object
      multiple_datasets_object(merged_object)
      
      removeModal()
      showNotification("Datasets merged successfully without integration!", type = "success")
      updateUIElements()
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error during merge:", e$message), type = "error")
    })
  })
  
  # Modal for dataset loading
  observeEvent(input$open_file_input_modal, {
    showModal(modalDialog(
      title = "Load Datasets",
      easyClose = TRUE,
      footer = NULL,
      tagList(
        tags$button(
          type = "button",
          class = "close",
          "data-dismiss" = "modal",
          "Close window"
        ),
        
        # Dataset configuration
        fluidRow(
          column(6, numericInput('num_datasets', 'Number of datasets to upload:', value = 2, min = 1)),
          column(6, selectInput("species_choice_merge", "Select species:",
                                choices = c("Mouse" = "mouse", "Human" = "human", "Rat" = "rat"),
                                selected = "mouse"))
        ),
        
        tags$hr(style = "border-top: 1px solid #ddd;"),
        
        # Dynamic file inputs
        uiOutput("fileInputs"),
        
        tags$hr(style = "border-top: 1px solid #ddd;"),
        
        # QC Parameters
        h4("Quality Control Parameters"),
        fluidRow(
          column(4, numericInput("min_features_merge", "Minimum features per cell:",
                                 value = 200, min = 0)),
          column(4, numericInput("max_features_merge", "Maximum features per cell:",
                                 value = 4000, min = 0)),
          column(4, numericInput("max_mt_percent_merge", "RNA Mitochondrial %:",
                                 value = 5, min = 0, max = 100))
        ),
        helpText("These parameters will be applied during data processing to filter cells."),
        
        tags$hr(style = "border-top: 1px solid #ddd;"),
        
        # Processing strategy selection
        h4("Processing Strategy", style = "color: #2c3e50;"),
        radioButtons("merge_strategy", "Choose processing method:",
                     choices = list(
                       "Standard Integration" = "integration",
                       "Simple Merge (Preserve original clusters)" = "simple_merge"
                     ),
                     selected = "integration",
                     inline = FALSE),
        
        # Information panel for simple merge
        conditionalPanel(
          condition = "input.merge_strategy == 'simple_merge'",
          div(style = "background-color: #fff3e0; padding: 10px; border-radius: 5px; margin: 10px 0;",
              icon("info-circle"), 
              strong(" Simple Merge Mode:"),
              p("Datasets will be merged without integration. Original cluster identities 
              and annotations will be preserved.",
                style = "margin: 5px 0 0 0; font-size: 0.9em; color: #666;"),
              checkboxInput("preserve_annotations", "Preserve custom cluster names", 
                            value = TRUE, width = "100%")
          )
        ),
        
        tags$hr(style = "border-top: 1px solid #ddd;"),
        
        # Action buttons
        div(style = "text-align: center; margin-top: 15px;",
            conditionalPanel(
              condition = "input.merge_strategy == 'integration'",
              actionButton('integrate', 'Integrate Datasets', 
                           class = 'btn-primary btn-lg', 
                           disabled = TRUE)
            ),
            conditionalPanel(
              condition = "input.merge_strategy == 'simple_merge'",
              actionButton('simple_merge', 'Merge Datasets', 
                           class = 'btn-warning btn-lg', 
                           disabled = TRUE)
            )
        )
      )
    ))
    
    # Generate dynamic file inputs
    observeEvent(input$num_datasets, {
      req(input$num_datasets > 0)
      output$fileInputs <- renderUI({
        lapply(1:input$num_datasets, function(i) {
          fluidRow(
            column(4, textInput(paste0('dataset_name', i), paste0('Dataset name ', i), 
                                value = paste0("Dataset ", i))),
            column(4, selectInput(paste0("dataset_type_merge", i), "Data type:",
                                  choices = list("snRNA-seq" = "snRNA_merge", 
                                                 "Multiome" = "multiome_merge", 
                                                 "Seurat Object" = "seurat_object_merge"))),
            column(4, fileInput(paste0('merge', i), paste0('Choose file ', i), 
                                accept = c('.rds', '.zip')))
          )
        })
      })
    })
    
    # Setup file input observers
    observeEvent(input$num_datasets, {
      lapply(1:input$num_datasets, function(i) {
        observe_file_input(i)
      })
    })
    
    # Enable/disable action buttons based on file status
    observe({
      req(input$num_datasets)
      
      # Check file processing status
      files_status <- sapply(1:input$num_datasets, function(i) {
        dataset_type <- input[[paste0("dataset_type_merge", i)]]
        file_input <- input[[paste0("merge", i)]]
        is_processed <- !is.null(data_loaded[[paste0("loaded", i)]]) && 
          data_loaded[[paste0("loaded", i)]]
        
        if (!is.null(dataset_type)) {
          if (dataset_type == "seurat_object_merge") {
            return(!is.null(file_input))
          } else {
            return(is_processed)
          }
        }
        return(FALSE)
      })
      
      all_loaded <- all(files_status)
      enough_files <- length(files_status) >= 2
      should_enable <- all_loaded && enough_files
      
      shinyjs::toggleState("simple_merge", condition = should_enable)
      shinyjs::toggleState("integrate", condition = should_enable)
    })
  })
  
  # Standard integration using modular functions
  observeEvent(input$integrate, {
    tryCatch({
      showModal(modalDialog(
        title = "Please Wait",
        "Integrating datasets...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      # Collect all processed Seurat objects
      seurat_list <- list()
      for (i in 1:input$num_datasets) {
        seurat_obj <- seurat_objects[[paste0("seurat_object", i)]]
        if (is.null(seurat_obj)) {
          stop(paste("Dataset", i, "is not processed yet"))
        }
        seurat_list[[i]] <- seurat_obj
      }
      
      # Use modular integration function
      integrated_object <- performDataIntegration(
        seurat_list = seurat_list,
        integration_method = "standard"
      )
      
      # Store integrated object
      multiple_datasets_object(integrated_object)
      
      # Enable UI controls
      shinyjs::enable("add_field")
      shinyjs::enable("add_metadata")
      
      removeModal()
      showNotification("Datasets integrated successfully!", type = "success")
      updateUIElements()
      
      # Clean up memory
      cleanupIntegrationMemory()
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error during integration:", e$message), type = "error")
    })
  })
  
  # Metadata management using modular functions
  observeEvent(input$add_field, {
    rv_metadata$num_fields <- rv_metadata$num_fields + 1
  })
  
  reactive_metadata_fields <- reactive({
    req(multiple_datasets_object())
    all_metadata_fields <- colnames(multiple_datasets_object()@meta.data)
    
    # Fields to exclude
    exclude_fields <- c("orig.ident", "nCount_RNA", "nCount_ATAC", "nFeature_RNA",
                        "nFeature_ATAC", "percent.mt")
    
    # Pattern-based exclusions
    exclude_pattern <- "^snn_res|^pANN|^PC_|^RNA_snn|^ATAC_snn|^integrated_snn"
    exclude_fields <- c(exclude_fields, all_metadata_fields[grepl(exclude_pattern, all_metadata_fields)])
    
    include_fields <- setdiff(all_metadata_fields, exclude_fields)
    return(include_fields)
  })
  
  output$metadata_inputs <- renderUI({
    req(multiple_datasets_object())
    datasets <- unique(multiple_datasets_object()@meta.data$dataset)
    lapply(1:rv_metadata$num_fields, function(j) {
      fluidRow(
        column(6, textInput(paste0("metadata_name_", j), paste0("Metadata Field Name ", j), value = "")),
        column(6, lapply(datasets, function(dataset_name) {
          dataset_index <- which(datasets == dataset_name)
          display_name <- if(!is.null(input[[paste0("dataset_name", dataset_index)]]) && 
                             input[[paste0("dataset_name", dataset_index)]] != "") {
            input[[paste0("dataset_name", dataset_index)]]
          } else {
            dataset_name
          }
          textInput(paste0("metadata_value_", dataset_name, "_", j), 
                    paste0("Value for ", display_name, " ", j), value = "")
        }))
      )
    })
  })
  
  # Process metadata using modular function
  observeEvent(input$add_metadata, {
    tryCatch({
      req(multiple_datasets_object())
      
      # Use modular metadata processing function
      updated_seurat <- processMetadataFromUI(
        seurat_object = multiple_datasets_object(),
        input = input,
        num_fields = rv_metadata$num_fields
      )
      
      # Update the reactive object
      multiple_datasets_object(updated_seurat)
      
      showNotification("Metadata added successfully!", type = "success")
      
    }, error = function(e) {
      showNotification(paste("Error adding metadata:", e$message), type = "error")
    })
  })
  
  # Function to update UI elements after data loading
  updateUIElements <- function() {
    req(multiple_datasets_object())
    tryCatch({
      updateSelectInput(session, "group_by_select", choices = unique(multiple_datasets_object()@meta.data$dataset))
      
      # Update gene choices
      unique_genes <- rownames(LayerData(multiple_datasets_object(), assay = "RNA", layer = 'counts'))
      updatePickerInput(session, "geneInput_merge", choices = c("", unique_genes))
      
      message("UI elements updated successfully")
    }, error = function(e) {
      message(paste("Error updating UI elements:", e$message))
    })
  }
  
  # Reactive gene list for multiple datasets
  reactive_gene_list_merge <- reactive({
    req(multiple_datasets_object())
    unique_genes <- rownames(LayerData(multiple_datasets_object(), assay = "RNA", layer = 'counts'))
    c("", unique_genes)
  })
  
  # Function to update UI elements after data loading
  updateUIElements <- function() {
    updateSelectInput(session, "group_by_select", choices = unique(multiple_datasets_object()@meta.data$dataset))
    updatePickerInput(session, "geneInput_merge", choices = reactive_gene_list_merge())
  }
  
  
  ############################## Scaling and PCA reduction ##############################
  
  
  # Scaling, PCA and Elbowplot
  observeEvent(input$runScalePCA, {
    req(multiple_datasets_object())
    
    tryCatch({
      # Show modal dialog
      showModal(modalDialog(
        title = "Please Wait",
        "Scaling data and running PCA...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      showNotification("Scaling and PCA have begun...", type = "message")
      all_genes <- rownames(multiple_datasets_object())
      seurat_object_temp <- multiple_datasets_object()
      seurat_object_temp <- FindVariableFeatures(seurat_object_temp, selection.method = "vst", nfeatures = 3000)
      seurat_object_temp <- ScaleData(seurat_object_temp, features = all_genes)
      seurat_object_temp <- RunPCA(seurat_object_temp, npcs = 50)
      multiple_datasets_object(seurat_object_temp)
      req(multiple_datasets_object()[["pca"]])
      output$elbow_plot2 <- renderPlot({
        ElbowPlot(multiple_datasets_object())
      })
      showNotification("Scaling and PCA are complete.", type = "message")
      
      # Close the modal dialog
      removeModal()
    }, error = function(e) {
      removeModal() # Close the modal dialog in case of an error
      showNotification(paste0("Scaling and PCA errors:", e$message), type = "error")
    })
  })
  # Update Harmony variables choices
  observe({
    req(multiple_datasets_object())
    meta_cols <- colnames(multiple_datasets_object()@meta.data)
    updateSelectInput(session, "harmony_vars",
                      choices = meta_cols,
                      selected = "dataset")
  })
  
  # Run Harmony integration
  observeEvent(input$runHarmony, {
    req(multiple_datasets_object())
    
    tryCatch({
      showModal(modalDialog(
        title = "Please Wait",
        "Running Harmony integration...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      seurat_object_temp <- multiple_datasets_object()
      
      # Run Harmony
      seurat_object_temp <- RunHarmony(
        object = seurat_object_temp,
        group.by.vars = input$harmony_vars,
        dims.use = 1:input$harmony_dims,
        plot_convergence = TRUE
      )
      
      # Set Harmony as default reduction
      DefaultAssay(seurat_object_temp) <- "integrated"
      
      # Update the Seurat object
      multiple_datasets_object(seurat_object_temp)
      
      showNotification("Harmony integration completed successfully!", type = "message")
      removeModal()
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error in Harmony integration:", e$message), type = "error")
    })
  })
  
  # Observer forfind Clusters and display UMAP
  clustering_plot_merge <- reactiveVal()
  
  
  # Neighbors calculation and UMAP plotting with logging and error handling
  observeEvent(input$runFindNeighbors, {
    tryCatch({
      # Show modal dialog to indicate that the process is running
      showModal(modalDialog(
        title = "Please Wait",
        "Finding neighbors and running UMAP...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      # Start logging the process
      print("Finding neighbors and UMAP started...")
      
      # Ensure the Seurat object is loaded and available
      req(multiple_datasets_object())
      
      # Retrieve the Seurat object from the reactive value
      seurat_object_temp <- multiple_datasets_object()
      
      # Log the dimensions of the Seurat object
      print(paste("Seurat object retrieved. Dimensions:", dim(seurat_object_temp)[1], "cells,", dim(seurat_object_temp)[2], "features"))
      
      # Get the number of dimensions to use for the UMAP and neighbors calculation
      umap_dims <- input$dimension_2
      print(paste("Dimensions for FindNeighbors and UMAP:", umap_dims))
      
      # Step 1: Finding Neighbors
      print("Running FindNeighbors...")
      seurat_object_temp <- FindNeighbors(seurat_object_temp, dims = 1:umap_dims)
      print("FindNeighbors completed.")
      
      # Step 2: Running UMAP
      print("Running UMAP...")
      seurat_object_temp <- RunUMAP(seurat_object_temp, dims = 1:umap_dims)
      print("UMAP completed.")
      
      # Step 3: Store the updated Seurat object back to the reactive value
      multiple_datasets_object(seurat_object_temp)
      
      # Step 4: Check if UMAP embeddings are available
      umap_embeddings <- Embeddings(multiple_datasets_object(), reduction = "umap")
      print(paste("UMAP embedding dimensions:", dim(umap_embeddings)))
      
      # Step 5: Check the dimensions of metadata
      meta_data <- multiple_datasets_object()@meta.data
      print(paste("Metadata dimensions:", dim(meta_data)))
      
      # Step 6: Ensure UMAP embeddings and metadata are of the same size
      if (nrow(umap_embeddings) != nrow(meta_data)) {
        stop("Mismatch between UMAP embeddings and metadata. They must have the same number of rows.")
      }
      
      # Step 7: Check if 'orig.ident' column is present in metadata and print unique values
      if (!"orig.ident" %in% colnames(meta_data)) {
        stop("The 'orig.ident' column is missing in the metadata.")
      }
      print(paste("Unique 'orig.ident' values:", unique(meta_data$orig.ident)))
      
      # Step 8: Generating clustering plot
      print("Generating clustering plot...")
      
      # Log the number of unique values in the 'orig.ident' column
      unique_groups <- unique(meta_data$orig.ident)
      print(paste("Unique 'orig.ident' groups:", length(unique_groups), "values."))
      
      # Attempt to generate the clustering plot using `DimPlot`
      clustering_plot <- DimPlot(multiple_datasets_object(), group.by = "orig.ident") + ggtitle(NULL)
      
      # Log the plot data for debugging (make sure to confirm the object types and structure)
      plot_data <- as.data.frame(seurat_object_temp@meta.data)
      print(paste("Number of rows in meta.data for plotting:", nrow(plot_data)))
      print(paste("Number of rows in UMAP embeddings for plotting:", nrow(umap_embeddings)))
      
      # If the plot data and UMAP embeddings have different lengths, it will cause issues
      if (nrow(plot_data) != nrow(umap_embeddings)) {
        stop("Mismatch between plot data and UMAP embeddings. Cannot generate the plot.")
      }
      
      # Set the clustering plot reactive value
      clustering_plot_merge(clustering_plot)
      print("Clustering plot successfully generated.")
      
      # Step 9: Notify the user that the process is complete
      showNotification("Finding neighbors and UMAP completed.", type = "message")
      
      # Close the modal dialog after the process is done
      removeModal()
      
    }, error = function(e) {
      # Close the modal dialog and show an error notification
      removeModal()
      showNotification(paste("Error during find_neighbors: ", e$message), type = "error")
      print(paste("Error during find_neighbors: ", e$message))  # Log the full error
    })
  })
  
  
  
  
  # Observer for finding clusters and rendering UMAP
  observeEvent(input$runFindClusters, {
    tryCatch({
      # Show modal dialog
      showModal(modalDialog(
        title = "Please Wait",
        "Finding clusters and rendering UMAP...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      req(multiple_datasets_object())
      showNotification("Clustering begins...", type = "message")
      
      seurat_integrated_temp <- FindClusters(multiple_datasets_object(), resolution = input$resolution_step2, algorithm = as.integer(input$algorithm_select))
      seurat_integrated_temp <- JoinLayers(seurat_integrated_temp, assay = "RNA")
      multiple_datasets_object(seurat_integrated_temp)
      
      plot <- DimPlot(multiple_datasets_object(), reduction = "umap", repel = TRUE, label = FALSE) + ggtitle(NULL)
      
      if(input$remove_axes_umap_merge) { plot <- plot + NoAxes() }
      if(input$remove_legend_umap_merge) { plot <- plot + NoLegend() }
      
      clustering_plot_merge(plot)
      
      showNotification("Clustering and UMAP rendering completed.", type = "message")
      
      # Close the modal dialog
      removeModal()
    }, error = function(e) {
      removeModal() # Close the modal dialog in case of an error
      showNotification(paste0("Error during clustering: ", e$message), type = "error")
    })
  })
  
  
  
  output$UMAPPlot_cluster_merge <- renderPlot({
    req(clustering_plot_merge())
    clustering_plot_merge()
  })
  
  
  # Handler to download UMAP using modular functions
  output$downloadUMAP_merge <- createDownloadHandler(
    reactive_data = clustering_plot_merge,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(multiple_datasets_object(), default_name = "MergedUMAP") 
    }),
    data_name = "UMAP_plot",
    download_type = "plot",
    plot_params = list(
      file_type = input$umap_merge_format,
      width = if(input$umap_merge_format == "pdf") 11 else 10,
      height = if(input$umap_merge_format == "pdf") 8 else 6,
      dpi = input$dpi_umap_merge
    )
  )
  

  
  ############################## Visualize genes expressions ##############################
  
  
  # Observer for multiple datasets gene selection synchronization
  observeEvent(input$geneInput_merge, {
    selected_genes <- input$geneInput_merge
    
    if (!is.null(selected_genes) && length(selected_genes) > 0) {
      genes_text <- paste(selected_genes, collapse = ", ")
      updateTextInput(session, "gene_list_vln_merge", value = genes_text)
      updateTextInput(session, "gene_list_feature_merge", value = genes_text)
      updateTextInput(session, "gene_list_dot_merge", value = genes_text)
      updateTextInput(session, "gene_list_ridge_merge", value = genes_text)
      updateTextInput(session, "gene_list_genes_expression_merge", value = genes_text)
    }
  })
  
  
  # Function to safely get cluster choices from Seurat object
  getClusterChoices <- function(seurat_object) {
    if (is.null(seurat_object)) {
      return(NULL)
    }
    
    tryCatch({
      # Get cluster identities and convert to character vector
      cluster_idents <- Idents(seurat_object)
      cluster_choices <- as.character(unique(cluster_idents))
      
      # Remove any NA values and ensure it's a proper character vector
      cluster_choices <- cluster_choices[!is.na(cluster_choices)]
      
      return(cluster_choices)
    }, error = function(e) {
      message(paste("Error getting cluster choices:", e$message))
      return(NULL)
    })
  }
  
  # Corrected updateClusterChoices function
  updateClusterChoices <- function(session, seurat_object, cluster_input_ids = NULL) {
    if (is.null(seurat_object)) {
      return(NULL)
    }
    
    # Default cluster input IDs if not specified
    if (is.null(cluster_input_ids)) {
      cluster_input_ids <- list(
        select = c("select_cluster", "ident_1", "select_ident_subset"),
        checkbox = c("ident_2")
      )
    }
    
    tryCatch({
      # Get cluster choices safely
      cluster_choices <- getClusterChoices(seurat_object)
      
      if (is.null(cluster_choices) || length(cluster_choices) == 0) {
        message("No valid cluster choices found")
        return(NULL)
      }
      
      message(paste("Available clusters:", paste(cluster_choices, collapse = ", ")))
      
      # Update selectInput types
      if ("select" %in% names(cluster_input_ids)) {
        for (input_id in cluster_input_ids$select) {
          updateSelectInput(session, input_id, 
                            choices = cluster_choices,
                            selected = cluster_choices[1])
        }
      }
      
      # Update checkboxGroupInput types
      if ("checkbox" %in% names(cluster_input_ids)) {
        for (input_id in cluster_input_ids$checkbox) {
          updateCheckboxGroupInput(session, input_id, 
                                   choices = cluster_choices, 
                                   selected = character(0))
        }
      }
      
    }, error = function(e) {
      message(paste("Error updating cluster choices:", e$message))
    })
  }
  
  # Function to safely handle metadata columns
  getMetadataChoices <- function(seurat_object, column_name) {
    if (is.null(seurat_object) || is.null(column_name)) {
      return(NULL)
    }
    
    tryCatch({
      # Check if column exists
      if (!column_name %in% colnames(seurat_object@meta.data)) {
        message(paste("Column", column_name, "not found in metadata"))
        return(NULL)
      }
      
      # Get unique values and convert to character
      values <- seurat_object@meta.data[[column_name]]
      unique_values <- as.character(unique(values))
      
      # Remove NA values
      unique_values <- unique_values[!is.na(unique_values)]
      
      return(unique_values)
    }, error = function(e) {
      message(paste("Error getting metadata choices for", column_name, ":", e$message))
      return(NULL)
    })
  }
  
  
  
  # Observer pour la gestion des métadonnées et clusters
  observe({
    req(multiple_datasets_object())
    seurat_object <- multiple_datasets_object()
    
    # Patterns à exclure
    excluded_patterns <- c("orig.ident","percent.mt", "nCount_ATAC", "nFeature_ATAC",
                           "nFeature_RNA", "nCount_RNA", "^RNA_snn_",
                           "^pANN", "^DF", "^integrated", "^integrated_snn")
    pattern <- paste(excluded_patterns, collapse = "|")
    
    # Obtenir les colonnes de métadonnées filtrées
    metadata_fields <- colnames(seurat_object@meta.data)
    display_fields <- metadata_fields[!grepl(pattern, metadata_fields)]
    
    # Vérifier que display_fields n'est pas vide
    if (length(display_fields) == 0) {
      message("No valid metadata fields found")
      return()
    }
    
    tryCatch({
      # Mise à jour du sélecteur de groupement principal
      current_selection <- if (!is.null(input$group_by_select) && input$group_by_select %in% display_fields) {
        input$group_by_select
      } else {
        display_fields[1]
      }
      
      updateSelectInput(session, "group_by_select",
                        choices = display_fields,
                        selected = current_selection
      )
      
      # Attendre que group_by_select soit défini
      if (is.null(input$group_by_select)) {
        return()
      }
      
      # Gestion selon le groupement
      if (input$group_by_select == "dataset") {
        # Mettre à jour les choix de sous-groupement (avec option vide)
        sub_fields <- display_fields[display_fields != "dataset"]
        updateSelectInput(session, "metadata_to_compare",
                          choices = c("None" = "", sub_fields),
                          selected = input$metadata_to_compare
        )
        
        # Attendre que metadata_to_compare soit défini ou vide
        metadata_compare <- input$metadata_to_compare
        if (is.null(metadata_compare)) {
          return()
        }
        
        # Mise à jour des clusters selon le sous-groupement
        if (metadata_compare == "" || is.na(metadata_compare)) {
          # Si pas de sous-groupement, utiliser uniquement les datasets
          if ("dataset" %in% colnames(seurat_object@meta.data)) {
            datasets <- as.character(unique(seurat_object@meta.data[["dataset"]]))
            datasets <- datasets[!is.na(datasets) & datasets != ""]
            
            if (length(datasets) > 0) {
              # Mettre à jour tous les sélecteurs de clusters
              updateSelectInput(session, "cluster_order_dotplot_merge",
                                choices = datasets,
                                selected = datasets
              )
              updateSelectInput(session, "cluster_order_vln_merge",
                                choices = datasets,
                                selected = datasets
              )
              
              # Debug
              print("Updated with datasets:")
              print(datasets)
            }
          }
        } else {
          # Si sous-groupement sélectionné
          if (metadata_compare %in% colnames(seurat_object@meta.data)) {
            subclusters <- as.character(unique(seurat_object@meta.data[[metadata_compare]]))
            subclusters <- subclusters[!is.na(subclusters) & subclusters != ""]
            
            if (length(subclusters) > 0) {
              # Mettre à jour tous les sélecteurs de clusters
              updateSelectInput(session, "cluster_order_dotplot_merge",
                                choices = subclusters,
                                selected = subclusters
              )
              updateSelectInput(session, "cluster_order_vln_merge",
                                choices = subclusters,
                                selected = subclusters
              )
              
              # Debug
              print("Updated with subclusters:")
              print(subclusters)
            }
          }
        }
      } else {
        # Pour les autres groupements, comportement standard
        clusters <- character(0)
        
        if (input$group_by_select == "seurat_clusters") {
          # Cas spécial pour seurat_clusters - utiliser Idents()
          cluster_idents <- Idents(seurat_object)
          clusters <- as.character(unique(cluster_idents))
        } else if (input$group_by_select %in% colnames(seurat_object@meta.data)) {
          # Autres colonnes de métadonnées
          clusters <- as.character(unique(seurat_object@meta.data[[input$group_by_select]]))
        }
        
        # Nettoyer les clusters
        clusters <- clusters[!is.na(clusters) & clusters != ""]
        
        if (length(clusters) > 0) {
          # Debug
          print("Before updating with clusters:")
          print(clusters)
          
          # Mettre à jour tous les sélecteurs de clusters
          updateSelectInput(session, "cluster_order_dotplot_merge",
                            choices = clusters,
                            selected = clusters
          )
          updateSelectInput(session, "cluster_order_vln_merge",
                            choices = clusters,
                            selected = clusters
          )
          
          # Debug
          print("After updating inputs")
        }
      }
    }, error = function(e) {
      message(paste("Error in metadata observer:", e$message))
      # Ne pas afficher de notification d'erreur car cela peut perturber l'utilisateur
      # showNotification(paste("Error updating metadata choices:", e$message), type = "error")
    })
  }, priority = 1)
  
  # Updated selectInput to choose the genes to visualize
  observe({
    updatePickerInput(session, "geneInput_merge", choices = reactive_gene_list_merge())
  })
  
  observeEvent(input$geneInput_merge, {
    selected_genes <- input$geneInput_merge
    
    # Update text fields with selected genes
    updateTextInput(session, "gene_list_vln_merge", value = paste(selected_genes, collapse = ", "))
    updateTextInput(session, "gene_list_feature_merge", value = paste(selected_genes, collapse = ", "))
    updateTextInput(session, "gene_list_dot_merge", value = paste(selected_genes, collapse = ", "))
    updateTextInput(session, "gene_list_ridge_merge", value = paste(selected_genes, collapse = ", "))
    updateTextInput(session, "gene_list_genes_expression_merge", value = paste(selected_genes, collapse = ", "))
    
  })
  
  # FeaturePlot of selected gene
  feature_plot_merge <- reactiveVal()
  
  observeEvent(input$runFeaturePlot, {
    tryCatch({
      req(input$gene_list_feature_merge, multiple_datasets_object())
      genes <- unique(trimws(strsplit(input$gene_list_feature_merge, ",")[[1]]))
      seurat_object <- multiple_datasets_object()
      req(seurat_object)
      DefaultAssay(seurat_object) <- input$viz_assay_merge
      
      min_cut <- if (!is.na(input$min_cutoff_feature_merge)) input$min_cutoff_feature_merge else NA
      max_cut <- if (!is.na(input$max_cutoff_feature_merge)) input$max_cutoff_feature_merge else NA
      
      if(input$group_by_select == "dataset") {
        datasets <- unique(seurat_object$dataset)
        all_plots <- list()
        
        for(gene in genes) {
          gene_plots <- lapply(datasets, function(ds) {
            subset_obj <- subset(seurat_object, subset = dataset == ds)
            p <- FeaturePlot(object = subset_obj, features = gene, reduction = "umap", min.cutoff = min_cut, max.cutoff = max_cut) + ggtitle(paste(ds, "-", gene)) + theme(plot.title = element_text(size = 14, face = "bold"))
            
            if(input$add_nolegend_feature_merge) p <- p + NoLegend()
            if(input$add_noaxes_feature_merge) p <- p + NoAxes()
            return(p)
          })
          all_plots <- c(all_plots, gene_plots)
        }
        
        n_plots <- length(all_plots)
        n_rows <- ceiling(n_plots / 2)
        combined_plot <- wrap_plots(plots = all_plots, ncol = 2, nrow = n_rows)
        attr(combined_plot, "n_rows") <- n_rows
        
      } else {
        combined_plot <- FeaturePlot(object = seurat_object, features = genes, ncol = 2, reduction = "umap", min.cutoff = min_cut, max.cutoff = max_cut)
        if(input$add_nolegend_feature_merge) combined_plot <- combined_plot + NoLegend()
        if(input$add_noaxes_feature_merge) combined_plot <- combined_plot + NoAxes()
        attr(combined_plot, "n_rows") <- ceiling(length(genes) / 2)
      }
      
      feature_plot_merge(combined_plot)
      
    }, error = function(e) {
      showNotification(paste("Error in FeaturePlot:", e$message), type = "error")
    })
  })
  
  output$FeaturePlot2 <- renderPlot({
    req(feature_plot_merge())
    plot_obj <- feature_plot_merge()
    print(plot_obj)
  }, height = function() {
    plot_obj <- feature_plot_merge()
    n_rows <- attr(plot_obj, "n_rows") %||% 1
    base_height <- 400
    return(base_height * n_rows)
  })
  
  # VlnPlot of the selected gene
  vln_plot_merge <- reactiveVal()
  
  # VlnPlot
  observeEvent(input$runVlnPlot, {
    tryCatch({
      req(input$gene_list_vln_merge, multiple_datasets_object())
      genes <- unique(trimws(strsplit(input$gene_list_vln_merge, ",")[[1]]))
      seurat_object <- multiple_datasets_object()
      req(seurat_object)
      DefaultAssay(seurat_object) <- input$viz_assay_merge
      
      if(input$group_by_select == "dataset") {
        if(!is.null(input$metadata_to_compare) && input$metadata_to_compare != "") {
          if(input$comparison_mode == "split") {
            # Mode ventilation
            seurat_object$dataset_category <- paste(seurat_object$dataset,
                                                    seurat_object@meta.data[[input$metadata_to_compare]], sep = "_")
            
            # Filtrer selon les clusters sélectionnés
            if(!is.null(input$cluster_order_vln_merge) && length(input$cluster_order_vln_merge) > 0) {
              valid_combinations <- unlist(lapply(unique(seurat_object$dataset), function(ds) {
                paste(ds, input$cluster_order_vln_merge, sep = "_")
              }))
              cells_to_keep <- seurat_object$dataset_category %in% valid_combinations
              seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[cells_to_keep])
              seurat_object$dataset_category <- factor(seurat_object$dataset_category,
                                                       levels = valid_combinations)
            }
            
            plot <- VlnPlot(seurat_object,
                            features = genes,
                            group.by = "dataset_category",
                            pt.size = if(input$hide_vln_points_merge) 0 else 0.1,
                            ncol = 2,
                            raster = FALSE) +
              ggtitle(paste("Expression of", paste(genes, collapse=", "))) +
              theme(
                plot.title = element_text(hjust = 0.5, size = 14),
                axis.text.x = element_text(angle = 45, hjust = 1),
                plot.margin = margin(t = 20, r = 20, b = 60, l = 20)
              )
          } else {
            # Mode subset
            cells_subset <- seurat_object@meta.data[[input$metadata_to_compare]] %in% input$cluster_order_vln_merge
            seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[cells_subset])
            plot <- VlnPlot(seurat_object,
                            features = genes,
                            group.by = "dataset",
                            pt.size = if(input$hide_vln_points_merge) 0 else 0.1,
                            ncol = 2,
                            raster = FALSE) +
              ggtitle(paste("Expression of", paste(genes, collapse=", "))) +
              theme(
                plot.title = element_text(hjust = 0.5, size = 14),
                axis.text.x = element_text(angle = 45, hjust = 1),
                plot.margin = margin(t = 20, r = 20, b = 60, l = 20)
              )
          }
        } else {
          # Mode simple par dataset
          # Filtrer par dataset si des datasets sont sélectionnés
          if(!is.null(input$cluster_order_vln_merge) && length(input$cluster_order_vln_merge) > 0) {
            cells_to_keep <- seurat_object$dataset %in% input$cluster_order_vln_merge
            seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[cells_to_keep])
            seurat_object$dataset <- factor(seurat_object$dataset, levels = input$cluster_order_vln_merge)
          }
          
          plot <- VlnPlot(seurat_object,
                          features = genes,
                          group.by = "dataset",
                          pt.size = if(input$hide_vln_points_merge) 0 else 0.1,
                          ncol = 2,
                          raster = FALSE) +
            ggtitle(paste("Expression of", paste(genes, collapse=", "))) +
            theme(
              plot.title = element_text(hjust = 0.5, size = 14),
              axis.text.x = element_text(angle = 45, hjust = 1),
              plot.margin = margin(t = 20, r = 20, b = 60, l = 20)
            )
        }
      } else {
        # Mode standard pour les autres groupements
        # Filtrer selon les clusters sélectionnés
        if(!is.null(input$cluster_order_vln_merge) && length(input$cluster_order_vln_merge) > 0) {
          cells_to_keep <- seurat_object@meta.data[[input$group_by_select]] %in% input$cluster_order_vln_merge
          seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[cells_to_keep])
          seurat_object@meta.data[[input$group_by_select]] <- factor(
            seurat_object@meta.data[[input$group_by_select]],
            levels = input$cluster_order_vln_merge
          )
        }
        
        plot <- VlnPlot(seurat_object,
                        features = genes,
                        group.by = input$group_by_select,
                        pt.size = if(input$hide_vln_points_merge) 0 else 0.1,
                        ncol = 2,
                        raster = FALSE) +
          ggtitle(paste("Expression of", paste(genes, collapse=", "))) +
          theme(
            plot.title = element_text(hjust = 0.5, size = 14),
            axis.text.x = element_text(angle = 45, hjust = 1),
            plot.margin = margin(t = 20, r = 20, b = 60, l = 20)
          )
      }
      
      # Options de visualisation
      if(input$add_nolegend_vln_merge) plot <- plot + NoLegend()
      if(input$add_noaxes_vln_merge) plot <- plot + NoAxes()
      
      # Si plusieurs gènes, créer des plots séparés
      if(length(genes) > 1) {
        # Créer une liste de plots, un pour chaque gène
        plot_list <- lapply(genes, function(gene) {
          if(input$group_by_select == "dataset" &&
             !is.null(input$metadata_to_compare) &&
             input$metadata_to_compare != "" &&
             input$comparison_mode == "split") {
            p <- VlnPlot(seurat_object,
                         features = gene,
                         group.by = "dataset_category",
                         pt.size = if(input$hide_vln_points_merge) 0 else 0.1,
                         raster = FALSE) +
              ggtitle(gene) +
              theme(
                plot.title = element_text(hjust = 0.5, size = 14),
                axis.text.x = element_text(angle = 45, hjust = 1),
                plot.margin = margin(t = 20, r = 20, b = 60, l = 20)
              )
          } else if(input$group_by_select == "dataset") {
            p <- VlnPlot(seurat_object,
                         features = gene,
                         group.by = "dataset",
                         pt.size = if(input$hide_vln_points_merge) 0 else 0.1,
                         raster = FALSE) +
              ggtitle(gene) +
              theme(
                plot.title = element_text(hjust = 0.5, size = 14),
                axis.text.x = element_text(angle = 45, hjust = 1),
                plot.margin = margin(t = 20, r = 20, b = 60, l = 20)
              )
          } else {
            p <- VlnPlot(seurat_object,
                         features = gene,
                         group.by = input$group_by_select,
                         pt.size = if(input$hide_vln_points_merge) 0 else 0.1,
                         raster = FALSE) +
              ggtitle(gene) +
              theme(
                plot.title = element_text(hjust = 0.5, size = 14),
                axis.text.x = element_text(angle = 45, hjust = 1),
                plot.margin = margin(t = 20, r = 20, b = 60, l = 20)
              )
          }
          
          if(input$add_nolegend_vln_merge) p <- p + NoLegend()
          if(input$add_noaxes_vln_merge) p <- p + NoAxes()
          return(p)
        })
        
        # Combiner les plots en une grille 2×N
        plot <- wrap_plots(plot_list, ncol = 2)
      }
      
      attr(plot, "n_rows") <- ceiling(length(genes) / 2)
      vln_plot_merge(plot)
      
    }, error = function(e) {
      showNotification(paste("Error in VlnPlot:", e$message), type = "error")
    })
  })
  # Ajuster le rendu du plot
  output$VlnPlot2 <- renderPlot({
    req(vln_plot_merge())
    plot_obj <- vln_plot_merge()
    n_rows <- attr(plot_obj, "n_rows") %||% 1
    print(plot_obj)
  }, height = function() {
    plot_obj <- vln_plot_merge()
    n_rows <- attr(plot_obj, "n_rows") %||% 1
    return(400 * n_rows)  # Hauteur par ligne
  })
  
  
  # DotPlot reactive value
  dot_plot_merge <- reactiveVal()
  
  
  
  # DotPlot
  observeEvent(input$runDotPlot, {
    tryCatch({
      req(input$gene_list_dot_merge, multiple_datasets_object())
      genes <- unique(trimws(strsplit(input$gene_list_dot_merge, ",")[[1]]))
      seurat_object <- multiple_datasets_object()
      req(seurat_object)
      DefaultAssay(seurat_object) <- input$viz_assay_merge
      
      if(input$group_by_select == "dataset") {
        if(!is.null(input$metadata_to_compare) && input$metadata_to_compare != "") {
          if (input$comparison_mode == "split") {
            # Mode ventilation
            seurat_object$dataset_category <- paste(seurat_object$dataset,
                                                    seurat_object@meta.data[[input$metadata_to_compare]],
                                                    sep = "_")
            
            if (!is.null(input$cluster_order_dotplot_merge) && length(input$cluster_order_dotplot_merge) > 0) {
              valid_combinations <- unlist(lapply(unique(seurat_object$dataset), function(ds) {
                paste(ds, input$cluster_order_dotplot_merge, sep = "_")
              }))
              
              cells_to_keep <- seurat_object$dataset_category %in% valid_combinations
              seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[cells_to_keep])
              seurat_object$dataset_category <- factor(seurat_object$dataset_category,
                                                       levels = valid_combinations)
            }
            
            plot <- DotPlot(seurat_object, features = genes, group.by = "dataset_category")
          } else {
            # Mode subset : comparison directe des datasets pour les clusters sélectionnés
            cells_subset <- seurat_object@meta.data[[input$metadata_to_compare]] %in% input$cluster_order_dotplot_merge
            seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[cells_subset])
            plot <- DotPlot(seurat_object, features = genes, group.by = "dataset")
          }
        } else {
          # Mode simple par dataset sans ventilation
          plot <- DotPlot(seurat_object, features = genes, group.by = "dataset")
        }
      } else {
        # Mode standard pour les autres groupements
        if (!is.null(input$cluster_order_dotplot_merge) && length(input$cluster_order_dotplot_merge) > 0) {
          cells_to_keep <- seurat_object@meta.data[[input$group_by_select]] %in% input$cluster_order_dotplot_merge
          seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[cells_to_keep])
          seurat_object@meta.data[[input$group_by_select]] <- factor(
            seurat_object@meta.data[[input$group_by_select]],
            levels = input$cluster_order_dotplot_merge
          )
        }
        plot <- DotPlot(seurat_object, features = genes, group.by = input$group_by_select)
      }
      
      # Options de visualisation
      if (input$invert_axes) {
        plot <- plot + coord_flip()
      } else {
        plot <- plot + RotatedAxis()
      }
      if (input$add_noaxes_dot_merge) plot <- plot + NoAxes()
      if (input$add_nolegend_dot_merge) plot <- plot + NoLegend()
      
      dot_plot_merge(plot)
    }, error = function(e) {
      showNotification(paste("Error in DotPlot: ", e$message), type = "error")
    })
  })
  
  # Render DotPlot
  output$DotPlot2 <- renderPlot({
    req(dot_plot_merge())
    dot_plot_merge()
  })
  
  
  # Ridge Plot
  ridge_plot_merge <- reactiveVal()
  
  # Ridge Plot
  observeEvent(input$runRidgePlot, {
    tryCatch({
      req(input$gene_list_ridge_merge, multiple_datasets_object())
      genes <- unique(trimws(strsplit(input$gene_list_ridge_merge, ",")[[1]]))
      seurat_object <- multiple_datasets_object()
      req(seurat_object)
      DefaultAssay(seurat_object) <- input$viz_assay_merge
      
      print(paste("Genes for RidgePlot:", paste(genes, collapse=", ")))
      print(paste("Group by:", input$group_by_select))
      
      plot <- RidgePlot(
        seurat_object,
        features = genes,
        group.by = input$group_by_select
      )
      
      if (input$add_noaxes_ridge_merge) plot <- plot + NoAxes()
      if (input$add_nolegend_ridge_merge) plot <- plot + NoLegend()
      
      ridge_plot_merge(plot)
    }, error = function(e) {
      showNotification(paste("Error in RidgePlot: ", e$message), type = "error")
      print(paste("RidgePlot error details:", e$message))
    })
  })
  
  # Mettre à jour le rendu du RidgePlot
  output$Ridge_plot_merge <- renderPlot({
    req(ridge_plot_merge())
    ridge_plot_merge()
  })
  
  # Download handler for VlnPlot Merge using modular functions
  output$downloadVlnPlotMerge <- createDownloadHandler(
    reactive_data = vln_plot_merge,
    object_name_reactive = reactive({ getObjectNameForDownload(multiple_datasets_object(), default_name = "MergedDatasets") }),
    data_name = reactive({ paste0("VlnPlot_", paste(trimws(strsplit(input$gene_list_vln_merge, ",")[[1]]), collapse="_")) }),
    download_type = "plot",
    plot_params = list(file_type = input$plot_format_merge, width = 15, height = reactive({ (attr(vln_plot_merge(), "n_rows") %||% 1) * 8 }), dpi = input$dpi_input_merge)
  )
  
  # Download handler for FeaturePlot Merge using modular functions
  output$downloadFeaturePlotMerge <- createDownloadHandler(
    reactive_data = feature_plot_merge,
    object_name_reactive = reactive({ getObjectNameForDownload(multiple_datasets_object(), default_name = "MergedDatasets") }),
    data_name = reactive({ paste0("FeaturePlot_", paste(trimws(strsplit(input$gene_list_feature_merge, ",")[[1]]), collapse="_")) }),
    download_type = "plot",
    plot_params = list(file_type = input$plot_format_merge, width = reactive({ min(12 * length(trimws(strsplit(input$gene_list_feature_merge, ",")[[1]])), 24) }), height = reactive({ (attr(feature_plot_merge(), "n_rows") %||% 1) * 8 }), dpi = input$dpi_input_merge)
  )
  
  # Download handler for DotPlot Merge using modular functions
  output$downloadDotPlotMerge <- createDownloadHandler(
    reactive_data = dot_plot_merge,
    object_name_reactive = reactive({ getObjectNameForDownload(multiple_datasets_object(), default_name = "MergedDatasets") }),
    data_name = reactive({ paste0("DotPlot_", paste(trimws(strsplit(input$gene_list_dot_merge, ",")[[1]]), collapse="_")) }),
    download_type = "plot",
    plot_params = list(file_type = input$plot_format_merge, width = reactive({ min(12 + length(trimws(strsplit(input$gene_list_dot_merge, ",")[[1]])) * 0.5, 24) }), height = 8, dpi = input$dpi_input_merge)
  )
  
  # Download handler for RidgePlot Merge using modular functions
  output$downloadRidgePlotMerge <- createDownloadHandler(
    reactive_data = ridge_plot_merge,
    object_name_reactive = reactive({ getObjectNameForDownload(multiple_datasets_object(), default_name = "MergedDatasets") }),
    data_name = reactive({ paste0("RidgePlot_", paste(trimws(strsplit(input$gene_list_ridge_merge, ",")[[1]]), collapse="_")) }),
    download_type = "plot",
    plot_params = list(file_type = input$plot_format_merge, width = reactive({ min(10 + length(trimws(strsplit(input$gene_list_ridge_merge, ",")[[1]])) * 0.5, 20) }), height = 8, dpi = input$dpi_input_merge)
  )
  
  
  ################# Number of nuclei expressing a gene ###################
  number_of_nuclei_merge <- reactiveVal(NULL)
  
  
  # Observer pour l'analyse d'expression dans les datasets intégrés
  observeEvent(input$analyze_btn_genes_expression_merge, {
    tryCatch({
      req(input$gene_list_genes_expression_merge, multiple_datasets_object())
      
      # Utiliser la fonction modulaire pour les datasets intégrés
      result <- analyze_gene_expression(
        seurat_obj = multiple_datasets_object(),
        selected_genes = input$gene_list_genes_expression_merge,
        assay_name = input$viz_assay_merge,
        expression_threshold = input$logfc_threshold_genes_expression_merge %||% 0.1,
        is_integrated = TRUE  # Important: spécifier que c'est un dataset intégré
      )
      
      # Stocker les résultats
      number_of_nuclei_merge(result$data)
      
      # Afficher le tableau avec la fonction de rendu modulaire
      output$expression_summary_merge <- renderDT({
        render_expression_table(result, "expression_summary_merge")
      })
      
      # Notification de succès
      showNotification(
        paste("Expression analysis completed for", length(result$data$Gene), "genes"),
        type = "success"
      )
      
    }, error = function(e) {
      showNotification(paste("Error processing expression data:", e$message), type = "error")
    })
  })
  
  # Pour le téléchargement des données
  output$download_genes_number_expression_merge <- createDownloadHandler(
    reactive_data = number_of_nuclei_merge,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(multiple_datasets_object(), default_name = "IntegratedExpression") 
    }),
    data_name = "integrated_expression_analysis",
    download_type = "csv"
  )
  
  
  output$save_seurat_merge_2 <- createDownloadHandler(
    reactive_data = multiple_datasets_object,
    object_name_reactive = reactive({ getObjectNameForDownload(multiple_datasets_object(), default_name = "IntegratedSeurat") }),
    data_name = "integrated_object",
    download_type = "seurat",
    show_modal = TRUE
  )
  
  ############################## Heatmap and dual expression multi dataset ##############################
  
  # Tab 5: Heatmap and dual expression multi dataset
  
  heatmap_plot_multidataset <- reactiveVal()
  
  
  
  # Updated selectInput to choose how to group data
  observe({
    req(multiple_datasets_object())
    seurat_object <- multiple_datasets_object()
    
    # Patterns à exclure
    excluded_patterns <- c("orig.ident","percent.mt", "nCount_ATAC", "nFeature_ATAC",
                           "nFeature_RNA", "nCount_RNA", "^RNA_snn_",
                           "^pANN", "^DF", "^integrated", "^integrated_snn")
    pattern <- paste(excluded_patterns, collapse = "|")
    
    # Récupérer tous les champs des métadonnées et filtrer
    meta_cols <- colnames(seurat_object@meta.data)
    filtered_cols <- meta_cols[!grepl(pattern, meta_cols)]
    
    # Mettre à jour le selectInput
    updateSelectInput(session, "dataset_select_heatmap", choices = filtered_cols)
  })
  
  # Updated selectInput to choose the genes to visualize
  observe({
    updatePickerInput(session, "gene_select_heatmap_multi", choices = reactive_gene_list_merge())
  })
  
  observe({
    req(multiple_datasets_object())
    clusters <- levels(Idents(multiple_datasets_object()))
    updateSelectInput(session, "cluster_selector_merge", choices = clusters, selected = clusters)
    updateTextInput(session, "text_clusters_merge", value = paste(clusters, collapse = ","))
  })
  
  
  observeEvent(input$select_all_clusters_merge, {
    req(multiple_datasets_object())
    tryCatch({
      if (input$select_all_clusters_merge) {
        shinyjs::disable("text_clusters_merge")
      } else {
        shinyjs::enable("text_clusters_merge")
      }
    }, error = function(e) {
      showNotification(paste("Erreur lors de la sélection des clusters : ", e$message), type = "error")
    })
  })
  
  
  
  # Générer la heatmap en prenant en compte les clusters spécifiés ou tous les clusters
  observeEvent(input$generateHeatmapMulti, {
    showModal(modalDialog(
      title = "Processing",
      div(
        h4("Generating Multi-dataset Heatmap...", style = "text-align: center;"),
        p("Please wait while we process data across your datasets.", style = "text-align: center;"),
        p("This operation might take longer for multiple datasets and large gene sets.",
          style = "text-align: center; color: #666;")
      ),
      footer = NULL,
      easyClose = FALSE
    ))
    
    req(multiple_datasets_object())
    seurat_object <- multiple_datasets_object()
    selected_group_by <- input$dataset_select_heatmap
    selected_assay <- input$assay_select_heatmap
    
    tryCatch({
      # Vérifier si l'assay sélectionné existe
      if (!(selected_assay %in% names(seurat_object@assays))) {
        showNotification(paste("Assay", selected_assay, "not found in object. Using RNA instead."), type = "warning")
        selected_assay <- "RNA"
      }
      
      # Définir l'assay actif
      DefaultAssay(seurat_object) <- selected_assay
      
      # Filtrer les clusters si nécessaire
      if (!input$select_all_clusters_merge) {
        if (nchar(input$text_clusters_merge) > 0) {
          specified_clusters <- unlist(strsplit(trimws(input$text_clusters_merge), ",\\s*"))
          valid_clusters <- specified_clusters %in% levels(Idents(seurat_object))
          if (sum(valid_clusters) == 0) {
            showNotification("None of the specified clusters is valid.", type = "error")
            removeModal()
            return()
          }
          Idents(seurat_object) <- factor(Idents(seurat_object), levels = specified_clusters)
          seurat_object <- subset(seurat_object, idents = specified_clusters)
        } else {
          selected_clusters <- input$cluster_selector_merge
          Idents(seurat_object) <- factor(Idents(seurat_object), levels = selected_clusters)
          seurat_object <- subset(seurat_object, idents = selected_clusters)
        }
      }
      
      # Sélectionner les gènes
      if(input$use_top10_genes_merge) {
        # Utiliser l'assay sélectionné pour trouver les marqueurs
        markers <- FindAllMarkers(seurat_object, min.pct = 0.25, assay = selected_assay)
        selected_genes <- markers %>%
          group_by(cluster) %>%
          dplyr::filter(avg_log2FC > 1) %>%
          slice_head(n = 10) %>%
          ungroup() %>%
          pull(gene)
        
        # Afficher les gènes sélectionnés pour plus de transparence
        gene_message <- paste("Selected top genes:", paste(head(selected_genes, 10), collapse=", "),
                              ifelse(length(selected_genes) > 10, "...", ""))
        showNotification(gene_message, type = "message", duration = 10)
      } else {
        selected_genes <- input$gene_select_heatmap_multi
      }
      
      # Vérifier la validité des gènes dans l'assay sélectionné
      final_valid_genes <- selected_genes %in% rownames(seurat_object[[selected_assay]])
      
      if(sum(final_valid_genes) > 0) {
        valid_gene_list <- selected_genes[final_valid_genes]
        
        # Créer la heatmap
        plot <- DoHeatmap(seurat_object,
                          features = valid_gene_list,
                          group.colors = c("lightgrey", "blue", "red"),
                          group.by = selected_group_by,
                          assay = selected_assay) +
          theme(axis.text.y = element_text(size = 8, face = "bold")) +  # Améliorer la visibilité des noms de gènes
          labs(title = paste("Heatmap using", selected_assay, "assay"))
        
        # Stocker le plot et les gènes valides pour référence
        heatmap_plot_multidataset(plot)
        
        # Afficher le nombre de gènes inclus
        showNotification(paste(sum(final_valid_genes), "genes included in heatmap"), type = "message")
      } else {
        showNotification("No valid genes found in the selected assay.", type = "error")
      }
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
    
    removeModal()
  })
  # Heatmap Output
  output$heatmap_plot_multi <- renderPlot({
    req(heatmap_plot_multidataset())
    heatmap_plot_multidataset()
  })
  
  
  
  # Handler to download heatmap using modular functions
  output$download_heatmap_multi <- createDownloadHandler(
    reactive_data = heatmap_plot_multidataset,
    object_name_reactive = reactive({ getObjectNameForDownload(multiple_datasets_object(), default_name = "MultiDatasets") }),
    data_name = "heatmap",
    download_type = "plot",
    plot_params = list(file_type = "tiff", width = 10, height = 8, dpi = input$dpi_heatmap_multi)
  )
  
  scatter_plot_multidataset <- reactiveVal()
  
  
  # Updated selectInput to choose how to group data
  observe({
    updateSelectInput(session, "dataset_select_scatter", choices = reactive_metadata_fields())
  })
  
  
  # Mettez à jour les choix de selectInput pour les gènes/features et les clusters
  observe({
    updatePickerInput(session, "gene_select_scatter_multi1", choices = reactive_gene_list_merge())
    updatePickerInput(session, "gene_select_scatter_multi2", choices = reactive_gene_list_merge())
  })
  
  observe({
    req(multiple_datasets_object())
    clusters <- levels(Idents(multiple_datasets_object()))
    updateSelectInput(session, "cluster_selector_merge_scatter", choices = clusters, selected = clusters)
    updateTextInput(session, "text_clusters_merge_scatter", value = paste(clusters, collapse = ","))
  })
  
  
  observeEvent(input$select_all_clusters_merge_scatter, {
    req(multiple_datasets_object())
    tryCatch({
      if (input$select_all_clusters_merge_scatter) {
        shinyjs::disable("text_clusters_merge_scatter")
      } else {
        shinyjs::enable("text_clusters_merge_scatter")
      }
    }, error = function(e) {
      showNotification(paste("Erreur lors de la sélection des clusters : ", e$message), type = "error")
    })
  })
  
  
  
  observeEvent(input$generateScatterMulti, {
    req(multiple_datasets_object(), input$gene_select_scatter_multi1, input$gene_select_scatter_multi2)
    seurat_object <- multiple_datasets_object()
    specified_clusters <- NULL
    if (!input$select_all_clusters_merge_scatter && nchar(input$text_clusters_merge_scatter) > 0) {
      specified_clusters <- unlist(strsplit(trimws(input$text_clusters_merge_scatter), ",\\s*"))
    }
    tryCatch({
      if (length(specified_clusters) > 0) {
        valid_clusters <- specified_clusters %in% levels(Idents(seurat_object))
        if (sum(valid_clusters) == 0) {
          showNotification("None of the specified clusters is valid.", type = "error")
          return()
        }
        seurat_object <- subset(seurat_object, idents = specified_clusters)
      }
      plot <- FeatureScatter(
        object = seurat_object,
        feature1 = input$gene_select_scatter_multi1,
        feature2 = input$gene_select_scatter_multi2,
        group.by = input$dataset_select_scatter
      )
      scatter_plot_multidataset(plot)
    }, error = function(e) {
      showNotification(paste("Error in generating scatter plot: ", e$message), type = "error")
    })
  })
  
  
  # Rendre la heatmap dans l'interface utilisateur
  output$scatter_plot_multi <- renderPlot({
    req(scatter_plot_multidataset())
    scatter_plot_multidataset()
  })
  
  # Handler to download scatter plot using modular functions
  output$download_scatter_multi <- createDownloadHandler(
    reactive_data = scatter_plot_multidataset,
    object_name_reactive = reactive({ getObjectNameForDownload(multiple_datasets_object(), default_name = "MultiDatasets") }),
    data_name = "scatter_plot",
    download_type = "plot",
    plot_params = list(file_type = "tiff", width = 10, height = 8, dpi = input$dpi_scatter_multi)
  )
  
  # Variable réactive pour stocker le plot actuel
  current_plot_merge <- reactiveVal()
  
  # Observer les changements dans la sélection du type de plot et mettre à jour le plot
  observe({
    plot_type <- input$plot_type_select_merge
    if (plot_type == "FeaturePlot") {
      current_plot_merge(feature_plot_merge())
    } else if (plot_type == "VlnPlot") {
      current_plot_merge(vln_plot_merge())
    } else if (plot_type == "DotPlot") {
      current_plot_merge(dot_plot_merge())
    } else if (plot_type == "RidgePlot") {
      current_plot_merge(ridge_plot_merge())
    }
  })
  
  # Afficher le plot sélectionné
  output$selected_plot_display_merge <- renderPlot({
    req(current_plot_merge())
    current_plot_merge()
  })
  
  
  ############################## Assigning cell identity merge ##############################
  
  
  
  
  observe({
    if (!is.null(multiple_datasets_object())) {
      updateSelectInput(session, "select_cluster_merge", choices = unique(Idents(multiple_datasets_object())))
      updateSelectInput(session, "select_color_merge", choices = levels(Idents(multiple_datasets_object())))
      get_cluster_colors_merge(multiple_datasets_object())
    }
  })
  
  
  
  # Renaming clusters
  observeEvent(input$rename_single_cluster_merge_button, {
    req(input$select_cluster_merge, input$rename_single_cluster_merge, multiple_datasets_object())
    updated_seurat <- multiple_datasets_object()
    
    if (input$rename_single_cluster_merge %in% unique(Idents(updated_seurat))) {
      cells_to_merge <- which(Idents(updated_seurat) %in% c(input$select_cluster_merge, input$rename_single_cluster_merge))
      Idents(updated_seurat, cells = cells_to_merge) <- input$rename_single_cluster_merge
      showNotification(paste("Clusters merged under the name:", input$rename_single_cluster_merge), type = "message")
    } else {
      Idents(updated_seurat, cells = which(Idents(updated_seurat) == input$select_cluster_merge)) <- input$rename_single_cluster_merge
      showNotification(paste("Cluster renamed to:", input$rename_single_cluster_merge), type = "message")
    }
    
    # Utilisation de `setNames` pour renommer les identifiants des clusters
    new_ident <- setNames(levels(updated_seurat), levels(Idents(updated_seurat)))
    updated_seurat <- RenameIdents(updated_seurat, new_ident)
    
    multiple_datasets_object(updated_seurat)
    get_cluster_colors_merge(updated_seurat)
    updateSelectInput(session, "select_cluster_merge", choices = unique(Idents(updated_seurat)))
    
    # Mise à jour de 'ClusterIdents' pour refléter les nouveaux noms
    updated_seurat$ClusterIdents <- Idents(updated_seurat)
    multiple_datasets_object(updated_seurat)
    
    # Mise à jour des choix du sélecteur 'group_by_select' pour inclure 'ClusterIdents'
  })
  
  
  # Reactive variable for that tab
  cluster_colours_merge <- reactiveVal()
  
  # Fonction centralisée pour récupérer ou initialiser les couleurs des clusters pour multiple_datasets_object
  get_cluster_colors_merge <- function(seurat_object) {
    if (!is.null(seurat_object@misc$cluster_colors)) {
      cluster_colors <- seurat_object@misc$cluster_colors
    } else {
      cluster_colors <- scales::hue_pal()(length(unique(Idents(seurat_object))))
      names(cluster_colors) <- sort(unique(Idents(seurat_object)))
    }
    return(cluster_colors)
  }
  
  
  # Function to apply updated colors to Seurat object and update reactive value
  apply_cluster_colours <- function(seurat_object) {
    cluster_colors <- cluster_colours_merge()
    if (!is.null(cluster_colors)) {
      seurat_object@misc$cluster_colors <- cluster_colors  # Store the colors in the misc slot
    }
    return(seurat_object)
  }
  
  # Section de mise à jour des couleurs
  observeEvent(input$update_colour_merge_button, {
    message("Update the color of the selected cluster for multiple datasets")
    updated_seurat <- multiple_datasets_object()
    cluster_colors <- get_cluster_colors_merge(updated_seurat)
    
    if (!is.null(input$select_color_merge) && input$select_color_merge %in% names(cluster_colors)) {
      cluster_colors[input$select_color_merge] <- input$select_cluster_merge_color
      message(paste("Updating color for cluster", input$select_color_merge, "to", input$select_cluster_merge_color))
    } else {
      showNotification("Selected cluster is not valid.", type = "error")
      return()
    }
    
    updated_seurat@misc$cluster_colors <- cluster_colors
    multiple_datasets_object(updated_seurat)
    message("Current cluster colors:")
    print(updated_seurat@misc$cluster_colors)
    showNotification("Cluster colors saved in Seurat object.", type = "message")
  })
  
  #Création du plot UMAP directement avec plotly
  output$umap_finale_merge <- renderPlotly({
    req(multiple_datasets_object())
    message("Generating the final UMAP for multiple datasets")
    updated_seurat <- multiple_datasets_object()
    
    # Récupérer les coordonnées UMAP avec vérification des noms de colonnes
    umap_data <- as.data.frame(Embeddings(updated_seurat, reduction = "umap"))
    
    # Vérifier les noms des colonnes et les renommer si nécessaire
    message("UMAP column names: ", paste(colnames(umap_data), collapse=", "))
    dim_names <- colnames(umap_data)
    if (length(dim_names) >= 2) {
      # Renommer explicitement les colonnes pour être sûr
      colnames(umap_data)[1:2] <- c("UMAP_1", "UMAP_2")
    } else {
      showNotification("UMAP reduction not found or has incorrect dimensions", type = "error")
      return(NULL)
    }
    
    # Ajouter les informations de cluster
    umap_data$cluster <- Idents(updated_seurat)
    
    # Récupérer les couleurs personnalisées
    cluster_colors <- get_cluster_colors_merge(updated_seurat)
    
    # Créer un vecteur nommé qui mappe les IDs de cluster aux couleurs
    if (!all(levels(as.factor(umap_data$cluster)) %in% names(cluster_colors))) {
      missing_clusters <- setdiff(levels(as.factor(umap_data$cluster)), names(cluster_colors))
      default_colors <- scales::hue_pal()(length(missing_clusters))
      names(default_colors) <- missing_clusters
      cluster_colors <- c(cluster_colors, default_colors)
    }
    
    # Créer les données pour les centres des clusters (utilisés pour les étiquettes)
    centers <- data.frame(
      cluster = levels(as.factor(umap_data$cluster)),
      UMAP_1 = vapply(levels(as.factor(umap_data$cluster)), function(cl) {
        median(umap_data$UMAP_1[umap_data$cluster == cl], na.rm = TRUE)
      }, numeric(1)),
      UMAP_2 = vapply(levels(as.factor(umap_data$cluster)), function(cl) {
        median(umap_data$UMAP_2[umap_data$cluster == cl], na.rm = TRUE)
      }, numeric(1))
    )
    
    # Créer le graphique plotly directement avec les arguments en ligne
    p <- plot_ly() %>%
      add_trace(data = umap_data, x = ~UMAP_1, y = ~UMAP_2, color = ~cluster, colors = cluster_colors,
                type = "scatter", mode = "markers", marker = list(size = input$pt_size_merge, opacity = 0.7),
                hoverinfo = "text", text = ~paste("Cluster:", cluster)) %>%
      layout(title = list(text = input$plot_title_merge, font = list(size = 24)),
             xaxis = list(title = "UMAP 1", zeroline = FALSE),
             yaxis = list(title = "UMAP 2", zeroline = FALSE),
             hovermode = "closest", showlegend = FALSE)
    
    # Ajouter les annotations de manière sécurisée
    for (i in 1:nrow(centers)) {
      p <- p %>% add_annotations(x = centers$UMAP_1[i], y = centers$UMAP_2[i],
                                 text = as.character(centers$cluster[i]),
                                 showarrow = FALSE,
                                 font = list(size = input$label_font_size_merge * 4, color = "black"))
    }
    
    return(p)
  })
  
  output$save_seurat_merge_3 <- createDownloadHandler(
    reactive_data = multiple_datasets_object,
    object_name_reactive = reactive({ getObjectNameForDownload(multiple_datasets_object(), default_name = "IntegratedSeurat") }),
    data_name = "integrated_object",
    download_type = "seurat",
    show_modal = TRUE
  )
  
  ############################# Differential Gene Expression Analysis  #############################
  
  # Reactive variables for this tab
  diff_genes_compare <- reactiveVal(NULL)          # Stores DE genes for a cluster vs all others
  diff_genes_compare_cluster <- reactiveVal(NULL)  # Stores DE genes for group1 vs group2
  diff_genes_compare_datasets <- reactiveVal(NULL) # Stores DE genes for dataset1 vs dataset2
  filtered_umap_plot <- reactiveVal(NULL)          # Stores the filtered UMAP plot
  gene_table_storage <- reactiveVal(list())        # Stores all DE gene tables for Venn diagrams
  current_gene_lists <- reactiveVal(NULL)          # Stores current lists for Venn diagrams
  
  # Disable download buttons at startup
  shinyjs::disable("download_markers_single_cluster_merge")
  shinyjs::disable("download_markers_multiple_clusters_merge")
  shinyjs::disable("download_diff_dataset_cluster")
  shinyjs::disable("download_venn_diagram")
  
  ############################# UMAP Filtering and Display #############################
  
  # Dropdown menu to filter by Dataset
  output$dataset_filter_ui <- renderUI({
    req(multiple_datasets_object())
    selectInput("dataset_filter",
                label = "Filter by Dataset",
                choices = c(unique(multiple_datasets_object()@meta.data$dataset)),
                selected = unique(multiple_datasets_object()@meta.data$dataset),
                multiple = TRUE)
  })
  
  # Display the filtered UMAP plot
  output$filtered_umap_plot <- renderPlot({
    req(multiple_datasets_object(), input$dataset_filter)
    
    if (!"dataset" %in% colnames(multiple_datasets_object()@meta.data)) {
      stop("Dataset column not found in Seurat object metadata.")
    }
    
    cluster_colors <- get_cluster_colors_merge(multiple_datasets_object())
    show_labels <- input$show_labels_merge
    label_size <- ifelse(input$bold_labels_merge, 8, 5)
    
    plot <- if ("All Datasets" %in% input$dataset_filter || is.null(input$dataset_filter)) {
      DimPlot(multiple_datasets_object(),
              group.by = "ident",
              label = show_labels,
              label.size = label_size) +
        scale_color_manual(values = cluster_colors) +
        NoAxes() +
        NoLegend() +
        ggtitle(NULL)
    } else {
      valid_datasets <- input$dataset_filter %in% unique(multiple_datasets_object()@meta.data$dataset)
      subset_seurat <- subset(multiple_datasets_object(), subset = dataset %in% input$dataset_filter[valid_datasets])
      
      DimPlot(subset_seurat,
              group.by = "ident",
              label = show_labels,
              label.size = label_size) +
        scale_color_manual(values = cluster_colors) +
        NoAxes() +
        NoLegend() +
        ggtitle(NULL)
    }
    
    filtered_umap_plot(plot)
    plot
  })
  
  output$download_filtered_umap_plot <- createDownloadHandler(
    reactive_data = filtered_umap_plot,
    object_name_reactive = reactive({ getObjectNameForDownload(multiple_datasets_object(), default_name = "MergedUMAP") }),
    data_name = reactive({ paste0("filtered_UMAP", if(!is.null(input$dataset_filter)) paste0("_", paste(input$dataset_filter, collapse = "_")) else "") }),
    download_type = "plot",
    plot_params = list(file_type = input$filtered_umap_format, width = if(input$filtered_umap_format == "pdf") 11 else 10, height = if(input$filtered_umap_format == "pdf") 8 else 6, dpi = input$filtered_umap_plot_dpi)
  )
  
  output$save_seurat_merge_4 <- createDownloadHandler(
    reactive_data = multiple_datasets_object,
    object_name_reactive = reactive({ getObjectNameForDownload(multiple_datasets_object(), default_name = "IntegratedSeurat") }),
    data_name = "integrated_object",
    download_type = "seurat",
    show_modal = TRUE
  )
  
  
  ############################# Cluster vs All Comparison #############################
  
  # Reactive update of cluster dropdown
  observe({
    req(multiple_datasets_object())
    cluster_choices <- unique(Idents(multiple_datasets_object()))
    updateSelectInput(session, "selected_cluster", choices = cluster_choices, selected = cluster_choices[1])
  })
  
  # Reactive function for markers
  observeEvent(input$calculate_DE, {
    tryCatch({
      # Show modal
      showModal(modalDialog(
        title = "Please wait",
        "Calculating differentially expressed genes...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      showNotification("Calculating differentially expressed genes...", type = "message")
      
      subset_seurat <- multiple_datasets_object()
      
      # Calculate markers
      markers <- FindMarkers(subset_seurat,
                             ident.1 = input$selected_cluster,
                             min.pct = input$min_pct_merge,
                             logfc.threshold = input$logfc_threshold_merge, 
                             slot = 'data', 
                             assay = 'RNA')
      
      if (nrow(markers) > 0) {
        markers$p_val <- format(markers$p_val, scientific = TRUE, digits = 3)
        markers$p_val_adj <- format(markers$p_val_adj, scientific = TRUE, digits = 3)
        cleaned_gene_names <- clean_gene_names_for_html(rownames(markers))
        
        # Format for display
        markers$gene <- paste0('<a href="#" class="gene-name" data-gene="', cleaned_gene_names, '">', cleaned_gene_names, '</a>')
        markers_df <- as.data.frame(markers)
        
        # Store in standard reactive variable
        diff_genes_compare(markers_df)
        shinyjs::enable("download_markers_single_cluster_merge")
        
        # Store in global list for Venn diagrams
        table_name <- paste0("Cluster_", input$selected_cluster, "_vs_All_", format(Sys.time(), "%H%M%S"))
        description <- paste0("Cluster ", input$selected_cluster, " vs All Others")
        parameters <- list(
          min_pct = input$min_pct_merge,
          logfc_threshold = input$logfc_threshold_merge
        )
        
        store_de_table(
          table_data = markers_df,
          table_name = table_name,
          description = description,
          type = "single_cluster",
          parameters = parameters
        )
        
      } else {
        showNotification("No differentially expressed genes found for the selected cluster.", type = "message")
      }
      
      removeModal()
    }, error = function(e) {
      removeModal()
      showNotification(paste0("Error calculating differentially expressed genes: ", e$message), type = "error")
    })
  })
  
  # Display DE genes table
  output$DE_genes_table <- renderDataTable({
    tryCatch({
      datatable(diff_genes_compare(), escape = FALSE)
    }, error = function(e) {
      showNotification(paste0("Error displaying gene table: ", e$message), type = "error")
    })
  })
  
  output$download_venn_gene_lists <- createDownloadHandler(
    reactive_data = current_gene_lists,
    object_name_reactive = reactive({ getObjectNameForDownload(multiple_datasets_object(), default_name = "VennDiagram") }),
    data_name = "gene_lists",
    download_type = "ods"
  )

  ############################# Group vs Group Comparison #############################
  
  # UI to select clusters
  output$cluster1_compare_ui <- renderUI({
    req(multiple_datasets_object())
    selectInput("cluster1_compare",
                label = "Select clusters for group 1",
                choices = unique(Idents(multiple_datasets_object())),
                selected = NULL,
                multiple = TRUE)
  })
  
  output$cluster2_compare_ui <- renderUI({
    req(multiple_datasets_object())
    # Exclude clusters already selected in group 1
    remaining_clusters <- setdiff(unique(Idents(multiple_datasets_object())), input$cluster1_compare)
    selectInput("cluster2_compare",
                label = "Select clusters for group 2",
                choices = remaining_clusters,
                selected = NULL,
                multiple = TRUE)
  })
  
  # Observer to update second group choices
  observeEvent(input$cluster1_compare, {
    req(multiple_datasets_object())
    remaining_clusters <- setdiff(unique(Idents(multiple_datasets_object())), input$cluster1_compare)
    updateSelectInput(session, "cluster2_compare", choices = remaining_clusters)
  })
  
  # Cluster comparison
  observeEvent(input$compare_clusters_button, {
    tryCatch({
      showModal(modalDialog(
        title = "Please wait",
        "Finding differentially expressed genes...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      req(input$cluster1_compare, input$cluster2_compare, multiple_datasets_object())
      
      # Check that there's no overlap between groups
      if (any(input$cluster1_compare %in% input$cluster2_compare)) {
        showNotification("Groups must not have overlapping clusters!", type = "error")
        removeModal()
        return()
      }
      
      seurat_obj <- multiple_datasets_object()
      
      # Create new identity for comparison groups
      new_idents <- as.character(Idents(seurat_obj))
      new_idents[new_idents %in% input$cluster1_compare] <- "group1"
      new_idents[new_idents %in% input$cluster2_compare] <- "group2"
      Idents(seurat_obj) <- new_idents
      
      # Find markers
      temp_res <- FindMarkers(seurat_obj,
                              ident.1 = "group1",
                              ident.2 = "group2",
                              min.pct = input$min_pct_compare_merge,
                              logfc.threshold = input$logfc_threshold_compare_merge)
      
      # Format results
      temp_res$p_val <- format(temp_res$p_val, scientific = TRUE, digits = 3)
      temp_res$p_val_adj <- format(temp_res$p_val_adj, scientific = TRUE, digits = 3)
      temp_res$gene <- paste0('<a href="#" class="gene-name" data-gene="', rownames(temp_res), '">', rownames(temp_res), '</a>')
      
      # Add information about compared groups
      temp_res$comparison <- sprintf("Group1(%s) vs Group2(%s)",
                                     paste(input$cluster1_compare, collapse = ","),
                                     paste(input$cluster2_compare, collapse = ","))
      
      shinyjs::enable("download_markers_multiple_clusters_merge")
      diff_genes_compare_cluster(temp_res)
      
      # Store in global list for Venn diagrams
      group1_text <- paste(input$cluster1_compare, collapse = "_")
      group2_text <- paste(input$cluster2_compare, collapse = "_")
      table_name <- paste0("Clusters_", group1_text, "_vs_", group2_text, "_", format(Sys.time(), "%H%M%S"))
      description <- paste0("Clusters [", group1_text, "] vs [", group2_text, "]")
      parameters <- list(
        min_pct = input$min_pct_compare_merge,
        logfc.threshold = input$logfc_threshold_compare_merge,
        group1 = input$cluster1_compare,
        group2 = input$cluster2_compare
      )
      
      store_de_table(
        table_data = temp_res,
        table_name = table_name,
        description = description,
        type = "cluster_group",
        parameters = parameters
      )
      
      removeModal()
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Cluster comparison error:", e$message), type = "error")
    })
  })
  
  # Display DE genes table between groups
  output$diff_genes_table_compare <- renderDataTable({
    req(diff_genes_compare_cluster())
    tryCatch({
      datatable(diff_genes_compare_cluster(), escape = FALSE)
    }, error = function(e) {
      showNotification(paste0("Error displaying gene table: ", e$message), type = "error")
    })
  })
  
  ## Download gene comparison table between groups using modular functions
  output$download_markers_multiple_clusters_merge <- createDownloadHandler(
    reactive_data = reactive({ diff_genes_compare_cluster()[order(diff_genes_compare_cluster()$p_val), ] }),
    object_name_reactive = reactive({ getObjectNameForDownload(multiple_datasets_object(), default_name = "MergedDatasets") }),
    data_name = reactive({ paste0("diff_genes_comparison_", paste(input$cluster1_compare, collapse = "-"), "_VS_", paste(input$cluster2_compare, collapse = "-")) }),
    download_type = "csv"
  )
  
  ############################# Dataset vs Dataset Comparison #############################
  
  # UI to select datasets
  output$dataset1_compare_ui <- renderUI({
    req(multiple_datasets_object())
    selectInput("dataset1_compare",
                label = "Select datasets for group 1",
                choices = unique(multiple_datasets_object()@meta.data$dataset),
                selected = NULL,
                multiple = TRUE)
  })
  
  output$dataset2_compare_ui <- renderUI({
    req(multiple_datasets_object())
    remaining_datasets <- setdiff(unique(multiple_datasets_object()@meta.data$dataset),
                                  input$dataset1_compare)
    selectInput("dataset2_compare",
                label = "Select datasets for group 2",
                choices = remaining_datasets,
                selected = NULL,
                multiple = TRUE)
  })
  
  output$cluster_compare_ui <- renderUI({
    req(multiple_datasets_object())
    if (!input$all_clusters) {
      selectizeInput("cluster_compare",
                     label = "Select clusters to analyze",
                     choices = unique(Idents(multiple_datasets_object())),
                     selected = NULL,
                     multiple = TRUE,
                     options = list(
                       plugins = list('remove_button'),
                       placeholder = 'Select multiple clusters'
                     ))
    }
  })
  
  # Observer to update second group dataset choices
  observeEvent(input$dataset1_compare, {
    req(multiple_datasets_object())
    remaining_datasets <- setdiff(unique(multiple_datasets_object()@meta.data$dataset),
                                  input$dataset1_compare)
    updateSelectInput(session, "dataset2_compare", choices = remaining_datasets)
  })
  
  # Observer for comparison
  observeEvent(input$compare_datasets_button, {
    tryCatch({
      showModal(modalDialog(
        title = "Please wait",
        "Comparing datasets and clusters...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      req(input$dataset1_compare, input$dataset2_compare, multiple_datasets_object())
      
      # Check that there's no overlap between dataset groups
      if (any(input$dataset1_compare %in% input$dataset2_compare)) {
        showNotification("Groups must not have overlapping datasets!", type = "error")
        removeModal()
        return()
      }
      
      seurat_obj <- multiple_datasets_object()
      
      # Create a new column to identify dataset-cluster combinations
      seurat_obj$dataset_cluster <- paste(seurat_obj@meta.data$dataset,
                                          Idents(seurat_obj),
                                          sep = "_")
      
      # Create group identifiers
      group1_combinations <- expand.grid(
        dataset = input$dataset1_compare,
        cluster = if (input$all_clusters) unique(Idents(seurat_obj)) else input$cluster_compare,
        stringsAsFactors = FALSE
      )
      group2_combinations <- expand.grid(
        dataset = input$dataset2_compare,
        cluster = if (input$all_clusters) unique(Idents(seurat_obj)) else input$cluster_compare,
        stringsAsFactors = FALSE
      )
      
      # Create group IDs
      group1_ids <- apply(group1_combinations, 1, function(x) paste(x[1], x[2], sep = "_"))
      group2_ids <- apply(group2_combinations, 1, function(x) paste(x[1], x[2], sep = "_"))
      
      # Find markers
      temp_res <- FindMarkers(
        object = seurat_obj,
        ident.1 = group1_ids,
        ident.2 = group2_ids,
        group.by = "dataset_cluster",
        min.pct = input$min_pct_compare_dataset_merge,
        logfc.threshold = input$logfc_threshold_datasets
      )
      
      # Format results
      temp_res$p_val <- format(temp_res$p_val, scientific = TRUE, digits = 3)
      temp_res$p_val_adj <- format(temp_res$p_val_adj, scientific = TRUE, digits = 3)
      temp_res$gene <- paste0('<a href="#" class="gene-name" data-gene="', rownames(temp_res), '">', rownames(temp_res), '</a>')
      
      # Add comparison info
      temp_res$comparison <- sprintf(
        "Group1(Datasets: %s, Clusters: %s) vs Group2(Datasets: %s, Clusters: %s)",
        paste(input$dataset1_compare, collapse = ","),
        paste(if (input$all_clusters) "all" else input$cluster_compare, collapse = ","),
        paste(input$dataset2_compare, collapse = ","),
        paste(if (input$all_clusters) "all" else input$cluster_compare, collapse = ",")
      )
      
      diff_genes_compare_datasets(temp_res)
      
      if (nrow(temp_res) == 0) {
        showNotification("No differentially expressed genes found!", type = "error")
        removeModal()
        return()
      }
      
      # Store in global list for Venn diagrams
      dataset1_text <- paste(input$dataset1_compare, collapse = "_")
      dataset2_text <- paste(input$dataset2_compare, collapse = "_")
      cluster_text <- if (input$all_clusters) "AllClusters" else paste(input$cluster_compare, collapse = "_")
      
      table_name <- paste0("Datasets_", dataset1_text, "_vs_", dataset2_text, "_", cluster_text, "_", format(Sys.time(), "%H%M%S"))
      description <- paste0("Datasets [", dataset1_text, "] vs [", dataset2_text, "] for clusters [",
                            if (input$all_clusters) "All" else paste(input$cluster_compare, collapse = ","), "]")
      parameters <- list(
        min_pct = input$min_pct_compare_dataset_merge,
        logfc.threshold = input$logfc_threshold_datasets,
        dataset1 = input$dataset1_compare,
        dataset2 = input$dataset2_compare,
        clusters = if (input$all_clusters) "all" else input$cluster_compare
      )
      
      store_de_table(
        table_data = temp_res,
        table_name = table_name,
        description = description,
        type = "dataset_compare",
        parameters = parameters
      )
      
      shinyjs::enable("download_diff_dataset_cluster")
      removeModal()
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Comparison error:", e$message), type = "error")
    })
  })
  
  # Display DE genes table between datasets
  output$diff_dataset_cluster <- renderDataTable({
    req(diff_genes_compare_datasets())
    datatable(diff_genes_compare_datasets(), escape = FALSE)
  })
  
  # Download gene comparison table between datasets
  output$download_diff_dataset_cluster <- downloadHandler(
    filename = function() {
      dataset1_text <- paste(input$dataset1_compare, collapse = "-")
      dataset2_text <- paste(input$dataset2_compare, collapse = "-")
      cluster_text <- if (input$all_clusters) "AllClusters" else paste(input$cluster_compare, collapse = "-")
      paste0("diff-datasets-", dataset1_text, "-vs-", dataset2_text, "-", cluster_text, "-", Sys.Date(), ".csv")
    },
    content = function(file) {
      tryCatch({
        diff_genes_compare_datasets_DL <- diff_genes_compare_datasets()
        req(!is.null(diff_genes_compare_datasets_DL))
        diff_genes_compare_datasets_DL_sorted <- diff_genes_compare_datasets_DL[order(diff_genes_compare_datasets_DL$p_val), ]
        write.csv(diff_genes_compare_datasets_DL_sorted, file)
      }, error = function(e) {
        showNotification(paste0("Error downloading gene comparison table: ", e$message), type = "error")
      })
    },
    contentType = "text/csv"
  )
  
  ############################## Cluster Composition Multiple ##############################
  
  # Reactive value to store cluster composition data
  cluster_composition_multiple <- reactiveVal(NULL)
  
  # Observer for generating cluster composition table for multiple datasets
  observeEvent(input$generate_cluster_table_multiple, {
    tryCatch({
      req(multiple_datasets_object())
      
      showModal(modalDialog(
        title = "Generating Cluster Composition Table",
        "Analyzing cluster composition across datasets...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      seurat_obj <- multiple_datasets_object()
      
      # Use existing modular function from cells_genes_expressions_newarch.R
      cluster_composition <- create_cluster_composition_table(seurat_obj, is_integrated = TRUE)
      
      # Store in reactive value
      cluster_composition_multiple(cluster_composition)
      
      removeModal()
      showNotification("Cluster composition table generated successfully!", type = "message")
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error generating cluster composition table:", e$message), type = "error")
    })
  })
  
  # Render cluster composition table for multiple datasets using existing function
  output$cluster_table_multiple <- renderDT({
    req(cluster_composition_multiple())
    render_cluster_composition_table(cluster_composition_multiple(), is_integrated = TRUE)
  })
  
  # Download handler for multiple datasets cluster composition using modular functions
  output$download_cluster_composition_multiple <- createDownloadHandler(
    reactive_data = reactive({ data <- cluster_composition_multiple(); data$Size_Bar <- NULL; return(data) }),
    object_name_reactive = reactive({ getObjectNameForDownload(multiple_datasets_object(), default_name = "MultipleDatasets") }),
    data_name = "cluster_composition",
    download_type = "csv"
  )
  ###############Co-expression calculation#############

  
  ###############Co-expression calculation using modular functions#############
  
  # Reactive values for storing results
  gene_coexpression_data_multiple <- reactiveVal(NULL)
  coexpression_plot_multiple <- reactiveVal(NULL)
  
  # Update assay choices when object is loaded
  observe({
    req(multiple_datasets_object())
    updateSelectInput(session, "coexpr_assay_multiple", choices = names(multiple_datasets_object()@assays), selected = DefaultAssay(multiple_datasets_object()))
  })
  
  # Observer for co-expression analysis using modular function
  observeEvent(input$analyze_coexpression_multiple, {
    req(multiple_datasets_object(), input$gene_text_coexpression_multiple)
    
    tryCatch({
      # Use modular function from cells_genes_expressions_newarch.R
      coexpr_results <- analyze_gene_coexpression(seurat_obj = multiple_datasets_object(), genes = input$gene_text_coexpression_multiple, assay_name = input$coexpr_assay_multiple, expression_threshold = input$coexpr_threshold_multiple %||% 0, is_integrated = TRUE)
      
      # Store results and create plot using modular functions
      gene_coexpression_data_multiple(coexpr_results)
      coexpression_plot_multiple(create_coexpression_plot(coexpr_results$data, coexpr_results$genes_analyzed))
      
      showNotification(paste("Co-expression analysis completed for", length(coexpr_results$genes_analyzed), "genes"), type = "success")
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # Render results using modular functions
  output$gene_coexpression_table_multiple <- renderDT({
    req(gene_coexpression_data_multiple())
    render_coexpression_table(gene_coexpression_data_multiple(), "gene_coexpression_table_multiple")
  })
  
  output$gene_coexpression_plot_multiple <- renderPlot({
    req(coexpression_plot_multiple())
    coexpression_plot_multiple()
  }, height = 700)
  
  # Render summary statistics using modular approach
  output$coexpression_summary_stats_multiple <- renderTable({
    req(gene_coexpression_data_multiple())
    create_coexpression_summary_stats(gene_coexpression_data_multiple()$data, gene_coexpression_data_multiple()$genes_analyzed)
  }, striped = TRUE, hover = TRUE, spacing = 'l')
  
  # Download handlers using modular functions
  output$download_coexpression_table_multiple <- createDownloadHandler(
    reactive_data = reactive({ data <- gene_coexpression_data_multiple()$data; data$Coexpression_Visual <- NULL; return(data) }),
    object_name_reactive = reactive({ getObjectNameForDownload(multiple_datasets_object(), default_name = "MultipleDatasets") }),
    data_name = "coexpression_analysis",
    download_type = "csv"
  )
  
  output$download_coexpression_plot_multiple <- createDownloadHandler(
    reactive_data = coexpression_plot_multiple,
    object_name_reactive = reactive({ getObjectNameForDownload(multiple_datasets_object(), default_name = "MultipleDatasets") }),
    data_name = "coexpression_plot",
    download_type = "plot",
    plot_params = list(file_type = "pdf", width = 12, height = 8, dpi = 300)
  )

  ############################# Venn Diagram Comparison #############################
  
  # Initialize reactive storage for Venn diagrams
  venn_plot_rendered <- reactiveVal(NULL)
  
  # Update select inputs with available gene tables
  observe({
    updateVennSelectInputs(session, gene_table_storage(), 
                           c("venn_table_1", "venn_table_2", "venn_table_3"))
  })
  
  # Enable/disable generate button based on available tables
  observe({
    tables <- gene_table_storage()
    if (length(tables) < 2) {
      shinyjs::disable("generate_venn_btn")
    } else {
      shinyjs::enable("generate_venn_btn")
    }
  })
  
  # Generate Venn diagram using modular functions
  observeEvent(input$generate_venn_btn, {
    showModal(modalDialog(title = "Generating Venn Diagram", "Processing...", 
                          easyClose = FALSE, footer = NULL))
    
    # Use modular function for complete Venn generation
    result <- processVennGeneration(
      table_storage = gene_table_storage(),
      selected_tables = c(input$venn_table_1, input$venn_table_2, input$venn_table_3),
      filter_params = list(
        significant_only = c(input$significant_only_venn_1, input$significant_only_venn_2, 
                             input$significant_only_venn_3),
        log_fc_threshold = c(input$log_fc_threshold_venn_1, input$log_fc_threshold_venn_2, 
                             input$log_fc_threshold_venn_3),
        p_val_threshold = input$p_val_threshold_venn,
        use_adjusted_p = input$use_adjusted_p_venn,
        direction = input$venn_direction
      ),
      colors = c(input$venn_color_1, input$venn_color_2, input$venn_color_3)
    )
    
    removeModal()
    
    if (result$success) {
      # Store results using modular structure
      venn_plot_rendered(result$venn_plot)
      current_gene_lists(result$overlaps)
      updateSelectInput(session, "selected_gene_set", choices = names(result$overlaps))
      shinyjs::enable("download_venn_diagram")
    } else {
      showNotification(result$message, type = "error")
    }
  })
  
  # Render Venn diagram
  output$venn_plot <- renderPlot({
    req(venn_plot_rendered())
    grid.draw(venn_plot_rendered())
  })
  
  # Display gene table for selected overlap
  output$venn_gene_table <- renderDT({
    req(current_gene_lists(), input$selected_gene_set)
    overlaps <- current_gene_lists()
    selected_genes <- overlaps[[input$selected_gene_set]]
    
    if (length(selected_genes) == 0) {
      return(data.frame(Gene = character(0)))
    }
    
    gene_df <- data.frame(Gene = selected_genes)
    datatable(gene_df, 
              options = list(pageLength = 15, scrollX = TRUE, dom = 'Bfrtip', 
                             buttons = c('copy', 'csv', 'excel')), 
              rownames = FALSE)
  })
  
  # Download handlers using modular functions
  output$download_venn_diagram <- createDownloadHandler(
    reactive_data = venn_plot_rendered,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(multiple_datasets_object(), default_name = "VennDiagram") 
    }),
    data_name = "venn_comparison",
    download_type = "plot",
    plot_params = list(
      file_type = input$venn_diagram_format,
      width = 8,
      height = 6,
      dpi = input$venn_diagram_dpi
    )
  )
  
  output$download_venn_gene_lists <- createDownloadHandler(
    reactive_data = reactive({
      gene_lists <- current_gene_lists()
      if (is.null(gene_lists)) return(NULL)
      
      # Convert gene lists to data frame
      all_genes <- data.frame()
      for (list_name in names(gene_lists)) {
        genes <- gene_lists[[list_name]]
        if (length(genes) > 0) {
          df <- data.frame(Gene = genes, Gene_Set = list_name, stringsAsFactors = FALSE)
          all_genes <- rbind(all_genes, df)
        } else {
          df <- data.frame(Gene = "No genes", Gene_Set = list_name, stringsAsFactors = FALSE)
          all_genes <- rbind(all_genes, df)
        }
      }
      return(all_genes)
    }),
    object_name_reactive = reactive({ getObjectNameForDownload(multiple_datasets_object(), default_name = "VennDiagram") }),
    data_name = "gene_lists",
    download_type = "csv"
  )
  
  ############################# Helper Functions for DE Table Storage #############################
  
  # Function to store DE tables using modular approach
  store_de_table <- function(table_data, table_name, description, type, parameters) {
    # Ensure table is non-NULL and has rows
    if (is.null(table_data) || nrow(table_data) == 0) {
      showNotification("Cannot store an empty table.", type = "warning")
      return()
    }
    
    # Use modular storage function
    storeDETable(gene_table_storage(), table_data, table_name, description, type, parameters)
    
    # Update reactive storage
    current_storage <- gene_table_storage()
    current_storage[[table_name]] <- list(
      data = table_data,
      description = description,
      type = type,
      parameters = parameters,
      timestamp = Sys.time()
    )
    gene_table_storage(current_storage)
    
    showNotification(paste("Table", table_name, "stored successfully for Venn analysis"), 
                     type = "message")
  }
  
  
  ############################## Subseting seurat object ##############################
  
  
  # Variable réactive pour cet onglet
  subset_seurat_merge <- reactiveVal(NULL)
  shinyjs::disable("download_subset_merge")
  
  # Mise à jour des choix d'identité de cellules pour le sous-ensemble
  observe({
    req(multiple_datasets_object())
    seurat_object <- multiple_datasets_object()
    
    tryCatch({
      # Convertir les identifiants en caractères
      cluster_idents <- Idents(seurat_object)
      cluster_choices <- as.character(unique(cluster_idents))
      cluster_choices <- cluster_choices[!is.na(cluster_choices)]
      
      if (length(cluster_choices) > 0) {
        updateSelectInput(session, "select_ident_subset_merge", 
                          choices = cluster_choices,
                          selected = cluster_choices[1])
      }
    }, error = function(e) {
      message(paste("Error updating subset choices:", e$message))
    })
  })
  
  # Réinitialisation de l'objet subset_seurat_merge lors du changement de multiple_datasets_object
  observe({
    if (!is.null(multiple_datasets_object())) {
      subset_seurat_merge(multiple_datasets_object())
    }
  })
  
  # Affichage de l'UMAP global
  output$global_umap_merge <- renderPlot({
    req(multiple_datasets_object())
    plot_data <- DimPlot(multiple_datasets_object(), group.by = "ident", label = TRUE) +
      theme(axis.line = element_line(size = 0.5)) +
      NoLegend()+ggtitle(NULL)
    print(plot_data)
  })
  
  # Affichage de l'UMAP du sous-ensemble
  output$subset_umap_merge <- renderPlot({
    req(subset_seurat_merge())
    plot_data <- DimPlot(subset_seurat_merge(), group.by = "ident", label = TRUE) +
      theme(axis.line = element_line(size = 0.5)) +
      NoLegend()+ggtitle(NULL)
    print(plot_data)
  })
  
  # Téléchargement de l'objet Seurat sous-ensemble
  output$download_subset_merge <- createDownloadHandler(
    reactive_data = subset_seurat_merge,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(multiple_datasets_object(), default_name = "SubsetSeurat") 
    }),
    data_name = "subset",
    download_type = "seurat",
    show_modal = TRUE
  )
  
  # Mise à jour de la sélection de sous-ensemble
  observeEvent(input$apply_subset_merge, {
    tryCatch({
      req(multiple_datasets_object())
      
      # Si aucun identifiant n'est sélectionné, retourne l'objet original
      if (length(input$select_ident_subset_merge) == 0) {
        subset_seurat_merge(multiple_datasets_object())
        shinyjs::enable("download_subset_merge")
        return()
      }
      
      subsetted_seurat <- subset(multiple_datasets_object(), idents = input$select_ident_subset_merge)
      
      # Vérifiez s'il y a des cellules dans le sous-ensemble
      if (nrow(subsetted_seurat@meta.data) == 0) {
        showNotification("No cells found with the selected identities.", type = "error")
        return()
      }
      
      subset_seurat_merge(subsetted_seurat)
      shinyjs::enable("download_subset_merge")
    }, error = function(e) {
      showNotification(paste0("Error updating subset selection: ", e$message), type = "error")
    })
  })
  
  observeEvent(input$apply_gene_subset_merge, {
    tryCatch({
      req(multiple_datasets_object())
      gene_list <- unlist(strsplit(input$gene_list_merge, ","))
      gene_list <- trimws(gene_list)
      
      if (length(gene_list) == 0) {
        showNotification("Please enter valid gene names.", type = "error")
        return()
      }
      
      expression_matrix <- sapply(gene_list, function(gene) {
        FetchData(multiple_datasets_object(), vars = gene) >= input$expression_threshold_merge
      })
      
      cells_to_keep <- which(rowSums(expression_matrix) >= input$num_genes_to_express_merge)
      
      if (length(cells_to_keep) == 0) {
        showNotification("No cells found with the specified gene expression criteria.", type = "error")
        return()
      }
      
      subsetted_seurat <- subset(multiple_datasets_object(), cells = cells_to_keep)
      subset_seurat_merge(subsetted_seurat)
      shinyjs::enable("download_subset_merge")
    }, error = function(e) {
      showNotification(paste0("Error during gene-based subset: ", e$message), type = "error")
    })
  })
  
  # Mise à jour des choix de colonnes de métadonnées avec exclusions
  observe({
    req(multiple_datasets_object())
    
    # Liste complète des colonnes de métadonnées
    all_metadata_fields <- colnames(multiple_datasets_object()@meta.data)
    
    # Champs à exclure
    exclude_fields <- c("orig.ident", "nCount_RNA", "nCount_ATAC", "nFeature_RNA",
                        "nFeature_ATAC", "percent.mt")
    
    # Exclusion des champs qui contiennent certains motifs (optionnel)
    exclude_pattern <- "^snn_res|^pANN|^PC_|^RNA_snn|^ATAC_snn|^integrated_snn"
    exclude_fields <- c(exclude_fields, all_metadata_fields[grepl(exclude_pattern, all_metadata_fields)])
    
    # Filtrer les colonnes à afficher
    available_columns <- setdiff(all_metadata_fields, exclude_fields)
    
    # Mettre à jour le sélecteur
    updateSelectInput(session, "metadata_column_subset", choices = available_columns)
  })
  
  # Interface dynamique pour les valeurs de métadonnées
  output$metadata_values_ui <- renderUI({
    req(multiple_datasets_object(), input$metadata_column_subset)
    
    # Obtenir les valeurs uniques pour la colonne sélectionnée
    col_values <- unique(multiple_datasets_object()@meta.data[[input$metadata_column_subset]])
    
    # Détecter le type de données
    if (is.numeric(col_values)) {
      # Pour les données numériques, offrir une plage
      min_val <- min(col_values, na.rm = TRUE)
      max_val <- max(col_values, na.rm = TRUE)
      
      tagList(
        sliderInput("metadata_num_range", "Value range:",
                    min = min_val, max = max_val,
                    value = c(min_val, max_val)),
        checkboxInput("invert_metadata_selection", "Invert selection", value = FALSE)
      )
    } else {
      # Pour les données catégorielles, offrir une liste de sélection
      selectInput("metadata_cat_values", "Select values:",
                  choices = sort(as.character(col_values)),
                  multiple = TRUE,
                  selected = sort(as.character(col_values))[1])
    }
  })
  
  # Gestion du bouton de subset par métadonnées
  observeEvent(input$apply_metadata_subset, {
    tryCatch({
      req(multiple_datasets_object(), input$metadata_column_subset)
      
      meta_col <- input$metadata_column_subset
      col_values <- multiple_datasets_object()@meta.data[[meta_col]]
      
      # Logique de filtrage différente selon le type de données
      if (is.numeric(col_values)) {
        # Pour les données numériques, utiliser la plage
        req(input$metadata_num_range)
        min_val <- input$metadata_num_range[1]
        max_val <- input$metadata_num_range[2]
        
        if (input$invert_metadata_selection) {
          cells_to_keep <- rownames(multiple_datasets_object()@meta.data)[col_values < min_val | col_values > max_val]
        } else {
          cells_to_keep <- rownames(multiple_datasets_object()@meta.data)[col_values >= min_val & col_values <= max_val]
        }
      } else {
        # Pour les données catégorielles, utiliser les valeurs sélectionnées
        req(input$metadata_cat_values)
        selected_values <- input$metadata_cat_values
        
        cells_to_keep <- rownames(multiple_datasets_object()@meta.data)[col_values %in% selected_values]
      }
      
      # Vérifier si des cellules correspondent au critère
      if (length(cells_to_keep) == 0) {
        showNotification("No cells found matching the metadata criteria.", type = "error")
        return()
      }
      
      # Créer le sous-ensemble
      subsetted_seurat <- subset(multiple_datasets_object(), cells = cells_to_keep)
      subset_seurat_merge(subsetted_seurat)
      
      # Activer le bouton de téléchargement
      shinyjs::enable("download_subset_merge")
      
      showNotification(paste0("Subset created with ", length(cells_to_keep), " cells."), type = "message")
    }, error = function(e) {
      showNotification(paste0("Error during metadata-based subset: ", e$message), type = "error")
    })
  })
  
  
}
