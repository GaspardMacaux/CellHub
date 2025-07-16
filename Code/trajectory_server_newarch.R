############################## Differential expressed genes along the trajectory ##############################

# trajectory_server.R


trajectory_server <- function(input, output, session) {

  shinyjs::useShinyjs() # button deactivation

  ############################## Monocle Conversion and trajectory ##############################

  # Reactive values
  seurat_monocle <- reactiveVal()
  monocle_object <- reactiveVal()
  selected_cell_id <- reactiveVal(NULL)
  trajectory_plot <- reactiveVal()
  pseudotime_plot <- reactiveVal()

  # Disable buttons initially
  shinyjs::disable("convertToMonocle")
  shinyjs::disable("constructGraph")
  shinyjs::disable("startRootSelection")


  output$monocle_clusters_info <- renderInfoBox({
    if (!is.null(seurat_monocle())) {
      n_clusters <- 0
      all_clusters <- ""
      display_text <- ""
      
      # Check for renamed clusters first
      if (!is.null(Idents(seurat_monocle()))) {
        clusters <- unique(Idents(seurat_monocle()))
        n_clusters <- length(clusters)
        all_clusters <- paste(as.character(clusters), collapse = ", ")
        
        # Pour l'affichage, tronquer si trop long
        if (nchar(all_clusters) > 30) {
          display_text <- paste0(substr(all_clusters, 1, 27), "...")
        } else {
          display_text <- all_clusters
        }
        
      } else if ("seurat_clusters" %in% colnames(seurat_monocle()@meta.data)) {
        n_clusters <- length(unique(seurat_monocle()$seurat_clusters))
        display_text <- "Original clusters"
      }
      
      infoBox(
        "Clusters",
        value = n_clusters,
        subtitle = tags$span(
          display_text,
          title = all_clusters,  # Tooltip avec la liste complète
          style = "cursor: help;"
        ),
        icon = icon("layer-group"),
        color = "green",
        fill = TRUE
      )
    } else {
      infoBox(
        "Clusters",
        "No data",
        icon = icon("layer-group"),
        color = "light-blue"
      )
    }
  })
  
  # Info boxes for Monocle data
  output$monocle_cells_info <- renderInfoBox({
    if (!is.null(seurat_monocle())) {
      infoBox(
        "Cells",
        format(ncol(seurat_monocle()), big.mark = ","),
        icon = icon("circle"),
        color = "blue",
        fill = TRUE
      )
    } else {
      infoBox(
        "Cells",
        "No data",
        icon = icon("circle"),
        color = "light-blue"
      )
    }
  })
  
  # For Seurat loading
  observeEvent(input$load_seurat_file_monocle, {
    message("Attempting to read file at: ", input$load_seurat_file_monocle$datapath)

    tryCatch({
      temp_seurat <- readRDS(input$load_seurat_file_monocle$datapath)
      message("File successfully read.")
      seurat_monocle(temp_seurat)
      showNotification("Seurat object has been successfully loaded!")
      shinyjs::enable("convertToMonocle")
    }, error = function(e) {
      message("Error loading Seurat object: ", e$message)
      showNotification(paste("An error occurred: ", e$message), type = "error")
    })
  })

  # Initialize gene loadings
  gene_loadings <- NULL


  initialize_gene_loadings <- function(monocle_obj, gene_loadings_matrix) {
    gene_names <- rowData(monocle_obj)$gene_short_name
    if (!is.null(gene_names)) {
      common_genes <- intersect(gene_names, rownames(gene_loadings_matrix))
      if (length(common_genes) == nrow(gene_loadings_matrix)) {
        rownames(gene_loadings_matrix) <- common_genes
        gene_loadings <<- gene_loadings_matrix
      } else {
        showNotification("Mismatch between gene names and gene loadings matrix.", type = "error")
      }
    } else {
      showNotification("Gene names not found in Monocle object.", type = "error")
    }
  }

  observeEvent(input$convertToMonocle, {
    req(seurat_monocle())
    
    showModal(modalDialog(
      title = "Converting to Monocle",
      "This may take a moment...",
      footer = NULL,
      easyClose = FALSE
    ))
    
    tryCatch({
      seurat_obj <- seurat_monocle()
      message("Starting Seurat to Monocle conversion...")
      message(paste("Seurat object dimensions:", nrow(seurat_obj), "genes x", ncol(seurat_obj), "cells"))
      
      # Get current cell identities
      current_idents <- Idents(seurat_obj)
      message(paste("Current Idents levels:", paste(levels(current_idents), collapse = ", ")))
      
      # Create Monocle object using counts
      message("Creating Monocle object from Seurat...")
      
      # Extract the count matrix
      count_matrix <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
      
      # Create cell metadata with all Seurat metadata
      cell_metadata <- seurat_obj@meta.data
      cell_metadata$ClusterIdents <- as.character(current_idents)
      
      # Create gene metadata
      gene_metadata <- data.frame(
        gene_short_name = rownames(count_matrix),
        row.names = rownames(count_matrix)
      )
      
      # Create the cell_data_set
      monocle_temp <- new_cell_data_set(
        count_matrix,
        cell_metadata = cell_metadata,
        gene_metadata = gene_metadata
      )
      
      message("Monocle object created. Transferring reductions...")
      
      # Transfer PCA - CRITICAL: Ensure dimensions match
      if ("pca" %in% names(seurat_obj@reductions)) {
        pca_embeddings <- Embeddings(seurat_obj, "pca")
        # Verify dimensions match
        if (nrow(pca_embeddings) == ncol(monocle_temp)) {
          reducedDims(monocle_temp)[["PCA"]] <- pca_embeddings
          message(paste("PCA transferred:", ncol(pca_embeddings), "components"))
        } else {
          message(paste("PCA dimension mismatch: embeddings have", nrow(pca_embeddings), 
                        "cells but Monocle has", ncol(monocle_temp)))
        }
      }
      
      # Transfer UMAP - CRITICAL: Ensure dimensions match
      if ("umap" %in% names(seurat_obj@reductions)) {
        umap_embeddings <- Embeddings(seurat_obj, "umap")
        # Verify dimensions match
        if (nrow(umap_embeddings) == ncol(monocle_temp)) {
          reducedDims(monocle_temp)[["UMAP"]] <- umap_embeddings
          message(paste("UMAP transferred:", nrow(umap_embeddings), "cells"))
        } else {
          message(paste("UMAP dimension mismatch: embeddings have", nrow(umap_embeddings), 
                        "cells but Monocle has", ncol(monocle_temp)))
          # This is a critical error - we need UMAP for trajectory
          showNotification("UMAP dimensions don't match cell count. Cannot proceed.", type = "error")
          removeModal()
          return()
        }
      }
      
      # Calculate size factors
      message("Calculating size factors...")
      monocle_temp <- estimate_size_factors(monocle_temp)
      
      # Set up cluster partitions for Monocle
      message("Setting up cluster partitions...")
      
      # CRITICAL: Run cluster_cells to set up the internal structure properly
      monocle_temp <- cluster_cells(
        monocle_temp,
        reduction_method = "UMAP",
        resolution = NULL,  # Don't re-cluster, just set up structure
        k = 20
      )
      
      # Now override the clusters with our Seurat clusters
      monocle_temp@clusters[["UMAP"]]$clusters <- factor(cell_metadata$ClusterIdents)
      monocle_temp@clusters[["UMAP"]]$partitions <- factor(cell_metadata$ClusterIdents)
      
      # Store cluster info in colData as well
      colData(monocle_temp)$monocle_clusters <- factor(cell_metadata$ClusterIdents)
      
      # Validation
      message("Validating converted object...")
      message(paste("Monocle cells:", ncol(monocle_temp)))
      message(paste("UMAP dimensions:", nrow(reducedDims(monocle_temp)[["UMAP"]]), "x", 
                    ncol(reducedDims(monocle_temp)[["UMAP"]])))
      message(paste("Clusters:", paste(unique(colData(monocle_temp)$ClusterIdents), collapse = ", ")))
      
      # Store the object
      monocle_object(monocle_temp)
      
      # Enable next button
      shinyjs::enable("constructGraph")
      removeModal()
      
      showNotification("Conversion successful! Ready to construct trajectory.", type = "message")
      
    }, error = function(e) {
      removeModal()
      message("Error during conversion:", e$message)
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  

  # Construct trajectory graph
  observeEvent(input$constructGraph, {
    req(monocle_object())
    
    showModal(modalDialog(
      title = "Constructing Graph",
      "Creating trajectory using existing Seurat clusters...",
      footer = NULL,
      easyClose = FALSE
    ))
    
    tryCatch({
      # Get the monocle object
      monocle_temp <- monocle_object()
      
      message("Starting graph construction...")
      message(paste("Cells:", ncol(monocle_temp), "Genes:", nrow(monocle_temp)))
      
      # Check which cluster column to use
      cluster_column <- NULL
      if ("ClusterIdents" %in% colnames(colData(monocle_temp))) {
        cluster_column <- "ClusterIdents"
      } else if ("seurat_clusters" %in% colnames(colData(monocle_temp))) {
        cluster_column <- "seurat_clusters"
      } else {
        showNotification("No cluster information found", type = "error")
        removeModal()
        return()
      }
      
      message(paste("Using cluster column:", cluster_column))
      cluster_data <- colData(monocle_temp)[[cluster_column]]
      message(paste("Unique clusters:", paste(unique(cluster_data), collapse = ", ")))
      
      # IMPORTANT: Set up clusters structure properly before graph construction
      cell_names <- colnames(monocle_temp)
      
      # Create single partition
      single_partition <- factor(rep("1", ncol(monocle_temp)))
      names(single_partition) <- cell_names
      
      # Set up clusters with names
      cluster_factor <- factor(cluster_data)
      names(cluster_factor) <- cell_names
      
      # Initialize clusters slot
      monocle_temp@clusters <- SimpleList()
      monocle_temp@clusters[["UMAP"]] <- SimpleList(
        partitions = single_partition,
        clusters = cluster_factor
      )
      
      message("Cluster structure initialized")
      
      # Construct the trajectory graph
      monocle_temp <- learn_graph(
        monocle_temp,
        use_partition = FALSE,
        close_loop = FALSE
      )
      
      message("Graph construction completed")
      
      # Update the monocle object
      monocle_object(monocle_temp)
      
      # Create trajectory plot
      base_plot <- plot_cells(
        monocle_temp,
        color_cells_by = cluster_column,
        label_groups_by_cluster = TRUE,
        label_leaves = FALSE,
        label_branch_points = TRUE,
        cell_size = 0.5
      ) + 
        ggtitle("Trajectory on Clusters") +
        theme_minimal()
      
      # Store and display plot
      trajectory_plot(base_plot)
      output$trajectoryPlot <- renderPlot({
        base_plot
      })
      
      # Update cluster choices for root selection
      clusters <- sort(unique(colData(monocle_temp)[[cluster_column]]))
      updateSelectInput(session, "root_cluster_select", 
                        choices = clusters,
                        selected = clusters[1])
      
      # Enable root selection button
      shinyjs::enable("set_root_cell")
      
      removeModal()
      showNotification("Graph construction completed! Now select a root cluster.", type = "message")
      
    }, error = function(e) {
      removeModal()
      message("Error during graph construction:", e$message)
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # Set root cell based on selected cluster (partie qui calcule le pseudotime)
  observeEvent(input$set_root_cell, {
    req(monocle_object(), input$root_cluster_select)
    
    monocle_data <- monocle_object()
    
    # Determine cluster column
    cluster_column <- if("ClusterIdents" %in% colnames(colData(monocle_data))) {
      "ClusterIdents"
    } else {
      "seurat_clusters"
    }
    
    # Get cells in selected cluster
    cells_in_cluster <- which(colData(monocle_data)[[cluster_column]] == input$root_cluster_select)
    cell_barcodes <- colnames(monocle_data)[cells_in_cluster]
    
    message(paste("Selected cluster:", input$root_cluster_select, "with", length(cells_in_cluster), "cells"))
    
    # Find center cell of the cluster
    umap_coords <- reducedDims(monocle_data)[["UMAP"]][cells_in_cluster, ]
    center <- colMeans(umap_coords)
    distances <- sqrt(rowSums((umap_coords - center)^2))
    center_idx <- which.min(distances)
    root_cell <- cell_barcodes[center_idx]
    
    # Update selected cell
    selected_cell_id(root_cell)
    
    # Update info
    output$root_info <- renderText({
      paste("Root cluster:", input$root_cluster_select,
            "\nNumber of cells:", length(cells_in_cluster),
            "\nRoot cell selected (center of cluster)")
    })
    
    # Order cells with pseudotime
    showModal(modalDialog(
      title = "Ordering Cells",
      "Computing pseudotime trajectories...",
      footer = NULL,
      easyClose = FALSE
    ))
    
    tryCatch({
      message(paste("Ordering cells with root:", root_cell))
      
      # Order cells
      monocle_temp <- order_cells(
        monocle_data,
        reduction_method = "UMAP",
        root_cells = root_cell
      )
      
      message("Cell ordering completed")
      
      # IMPORTANT: Extract pseudotime using the pseudotime() function
      pt_values <- pseudotime(monocle_temp)
      
      if (is.null(pt_values) || all(is.na(pt_values))) {
        message("ERROR: Pseudotime calculation failed!")
        removeModal()
        showNotification("Error: Could not calculate pseudotime", type = "error")
        return()
      }
      
      # Add pseudotime to colData if not already there
      colData(monocle_temp)$pseudotime <- pt_values
      
      # Check pseudotime values
      message(paste("Pseudotime range:", round(min(pt_values, na.rm = TRUE), 2), 
                    "to", round(max(pt_values, na.rm = TRUE), 2)))
      message(paste("Cells with valid pseudotime:", sum(!is.na(pt_values))))
      
      # Verify it's in colData now
      if ("pseudotime" %in% colnames(colData(monocle_temp))) {
        message("Pseudotime successfully added to colData")
      }
      
      # CRITICAL: Update the monocle object
      monocle_object(monocle_temp)
      
      # Create updated trajectory plot with clusters
      traj_plot <- plot_cells(
        monocle_temp,
        color_cells_by = cluster_column,
        cell_size = 0.5,
        label_groups_by_cluster = TRUE,
        label_leaves = FALSE,
        label_branch_points = TRUE,
        graph_label_size = 1.5
      ) +
        ggtitle(paste("Trajectory - Root:", input$root_cluster_select)) +
        theme_minimal()
      
      # Create pseudotime plot - Use the pseudotime values directly
      pseudo_plot <- plot_cells(
        monocle_temp,
        color_cells_by = "pseudotime",
        cell_size = 0.5,
        label_cell_groups = FALSE,
        label_leaves = FALSE,
        label_branch_points = TRUE,
        graph_label_size = 1.5
      ) +
        ggtitle("Pseudotime Trajectory") +
        theme_minimal() +
        viridis::scale_color_viridis(
          name = "Pseudotime",
          option = "plasma"
        )
      
      # Store plots
      trajectory_plot(traj_plot)
      pseudotime_plot(pseudo_plot)
      
      # Update plot outputs
      output$trajectoryPlot <- renderPlot({
        trajectory_plot()
      })
      
      output$pseudotimePlot <- renderPlot({
        pseudotime_plot()
      })
      
      removeModal()
      showNotification("Pseudotime calculated successfully!", type = "message")
      
    }, error = function(e) {
      removeModal()
      message("Error in cell ordering:", e$message)
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  
  

  # Download handlers
  output$download_trajectory_umap <- downloadHandler(
    filename = function() {
      paste("trajectory_plot_", Sys.Date(), ".", input$trajectory_download_format, sep = "")
    },
    content = function(file) {
      tryCatch({
        req(trajectory_plot())
        ggsave(
          file,
          plot = trajectory_plot(),
          device = input$trajectory_download_format,
          dpi = input$trajectory_download_dpi,
          width = 10,
          height = 8
        )
      }, error = function(e) {
        showNotification(paste("Error saving plot:", e$message), type = "error")
      })
    }
  )

  output$download_pseudotime_umap <- downloadHandler(
    filename = function() {
      paste("pseudotime_plot_", Sys.Date(), ".", input$trajectory_download_format, sep = "")
    },
    content = function(file) {
      tryCatch({
        req(pseudotime_plot())
        ggsave(
          file,
          plot = pseudotime_plot(),
          device = input$trajectory_download_format,
          dpi = input$trajectory_download_dpi,
          width = 10,
          height = 8
        )
      }, error = function(e) {
        showNotification(paste("Error saving plot:", e$message), type = "error")
      })
    }
  )


  # Reactive values to store results
  pseudotime_diff_gene_results <- reactiveVal(NULL)
  gene_trajectory_plot <- reactiveVal(NULL)
  pr_deg_ids <- reactiveVal(NULL)
  
  # Function to identify genes that change along pseudotime
  CalculateDiffGenesPseudotime <- function() {
    req(monocle_object())
    
    monocle_obj <- monocle_object()
    
    # Check if trajectory is calculated
    if (is.null(monocle_obj@principal_graph) || length(monocle_obj@principal_graph) == 0) {
      showNotification("Trajectory not calculated. Complete trajectory analysis first.", type = "error")
      return(NULL)
    }
    
    # IMPORTANT: Check if pseudotime exists
    if (!"pseudotime" %in% colnames(colData(monocle_obj))) {
      showNotification("Pseudotime not calculated. Please select a root cell first.", type = "error")
      return(NULL)
    }
    
    # Additional check: verify pseudotime has valid values
    pt_values <- pseudotime(monocle_obj)
    if (all(is.na(pt_values))) {
      showNotification("Pseudotime values are all NA. Please re-run root cell selection.", type = "error")
      return(NULL)
    }
    
    tryCatch({
      # Get q-value threshold
      q_value_threshold <- input$sig_genes_cutoff
      
      # Log analysis parameters
      message("Starting differential gene analysis along trajectory...")
      message(paste("Q-value threshold:", q_value_threshold))
      message(paste("Number of cells with pseudotime:", sum(!is.na(pt_values))))
      
      # Show progress
      showNotification("Calculating differential genes... This may take a few minutes.", 
                       type = "message", duration = NULL)
      
      # Run differential test
      n_cores <- min(4, parallel::detectCores() - 1)
      
      markers <- graph_test(
        monocle_obj, 
        neighbor_graph = "principal_graph", 
        cores = n_cores
      )
      # Check if results are valid
      if (is.null(markers) || nrow(markers) == 0) {
        showNotification("No results from differential test.", type = "error")
        return(NULL)
      }
      
      # Log initial results
      message(paste("Total genes tested:", nrow(markers)))
      message(paste("Genes with q < 0.05:", sum(markers$q_value < 0.05)))
      message(paste("Genes with q < 0.01:", sum(markers$q_value < 0.01)))
      
      # Filter by q-value
      sig_genes <- markers[markers$q_value < q_value_threshold, ]
      
      # Check if we found genes
      if (nrow(sig_genes) == 0) {
        # Show distribution of q-values to help user
        min_q <- min(markers$q_value)
        showNotification(
          paste0("No genes with q-value < ", q_value_threshold,
                 ". Minimum q-value found: ", format(min_q, scientific = TRUE),
                 ". Try increasing threshold."), 
          type = "warning"
        )
        return(NULL)
      }
      
      # Sort by q-value
      sig_genes <- sig_genes[order(sig_genes$q_value), ]
      
      # Add gene name column and rank
      sig_genes$gene_name <- rownames(sig_genes)
      sig_genes$rank <- 1:nrow(sig_genes)
      
      # Calculate effect size (Moran's I gives us a measure of spatial autocorrelation)
      sig_genes$effect_size <- abs(sig_genes$morans_I)
      
      # Log results
      message(paste("Found", nrow(sig_genes), "differentially expressed genes"))
      message("Top 5 genes:")
      for (i in 1:min(5, nrow(sig_genes))) {
        message(paste("  ", sig_genes$gene_name[i], 
                      "- q-value:", format(sig_genes$q_value[i], scientific = TRUE),
                      "- Moran's I:", round(sig_genes$morans_I[i], 3)))
      }
      
      # Store results
      pseudotime_diff_gene_results(sig_genes)
      pr_deg_ids(rownames(sig_genes))
      
      # Update gene picker with better organization
      # Show top genes first
      gene_choices <- rownames(sig_genes)
      names(gene_choices) <- paste0(
        sig_genes$gene_name, 
        " (q=", format(sig_genes$q_value, digits = 2, scientific = TRUE), ")"
      )
      
      updatePickerInput(
        session, "gene_picker",
        choices = gene_choices,
        selected = gene_choices[1:min(5, length(gene_choices))],  # Pre-select top 5
        options = list(
          `live-search` = TRUE,
          `actions-box` = TRUE,
          `selected-text-format` = "count > 3"
        )
      )
      
      
      
      return(sig_genes)
      
    }, error = function(e) {
      message(paste("Error in differential analysis:", e$message))
      showNotification(paste("Error:", e$message), type = "error")
      return(NULL)
    })
  }
  
  # Run differential analysis when button clicked
  observeEvent(input$run_diff_gene_pseudotime, {
    req(monocle_object())
    
    showModal(modalDialog(
      title = "Analyzing Gene Expression",
      "Calculating differential genes along pseudotime...",
      "This analysis may take several minutes depending on dataset size.",
      footer = NULL,
      easyClose = FALSE
    ))
    
    result <- CalculateDiffGenesPseudotime()
    removeModal()
  })
  
  # Display differential expression results table - IMPROVED VERSION
  output$diffGeneTable <- renderDT({
    req(pseudotime_diff_gene_results())
    
    # Get data
    display_data <- pseudotime_diff_gene_results()
    
    # Select relevant columns and rename for clarity
    display_data <- display_data[, c("rank", "gene_name", "status", "p_value", "q_value", 
                                     "morans_I", "morans_test_statistic")]
    
    # Rename columns for better understanding
    colnames(display_data) <- c("Rank", "Gene", "Status", "P-value", "Q-value", 
                                "Moran's I", "Test Statistic")
    
    # Format numeric columns
    display_data$`P-value` <- formatC(display_data$`P-value`, format = "e", digits = 2)
    display_data$`Q-value` <- formatC(display_data$`Q-value`, format = "e", digits = 2)
    display_data$`Moran's I` <- round(display_data$`Moran's I`, 4)
    display_data$`Test Statistic` <- round(display_data$`Test Statistic`, 2)
    
    # Create interactive table
    datatable(
      display_data,
      options = list(
        pageLength = 25,  # Show more genes by default
        scrollX = TRUE,
        order = list(list(0, 'asc')),  # Order by rank
        columnDefs = list(
          list(className = 'dt-center', targets = c(0, 2)),  # Center rank and status
          list(className = 'dt-right', targets = c(3, 4, 5, 6))  # Right-align numbers
        ),
        dom = 'Bfrtip',  # Add buttons
        buttons = c('copy', 'csv', 'excel')  # Export options
      ),
      rownames = FALSE,
      filter = 'top',
      caption = "Genes differentially expressed along pseudotime trajectory"
    ) %>%
      formatStyle(
        'Q-value',
        backgroundColor = styleInterval(c(0.01, 0.05), c('lightgreen', 'lightyellow', 'white'))
      ) %>%
      formatStyle(
        'Moran\'s I',
        color = styleInterval(c(0.5), c('black', 'red')),
        fontWeight = styleInterval(c(0.5), c('normal', 'bold'))
      )
  })
  
  # Download handler for results - IMPROVED
  output$download_pseudotime_diff_genes <- downloadHandler(
    filename = function() {
      paste0("pseudotime_differential_genes_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(pseudotime_diff_gene_results())
      
      # Add metadata to the file
      data_to_save <- pseudotime_diff_gene_results()
      
      # Create a summary header
      summary_info <- data.frame(
        Analysis = "Monocle3 Pseudotime Differential Expression",
        Date = Sys.Date(),
        Total_Genes_Tested = nrow(monocle_object()),
        Significant_Genes = nrow(data_to_save),
        Q_value_Threshold = input$sig_genes_cutoff
      )
      
      # Write both summary and results
      write.csv(data_to_save, file, row.names = FALSE)
    }
  )

  # Plot selected genes along trajectory - IMPROVED VERSION
  observeEvent(input$visualize_gene_trajectory, {
    req(monocle_object(), input$gene_picker)
    
    if (length(input$gene_picker) == 0) {
      showNotification("Please select at least one gene", type = "warning")
      return()
    }
    
    # Check if too many genes selected
    if (length(input$gene_picker) > 20) {
      showNotification("Warning: Plotting many genes may be slow. Consider selecting fewer genes.", 
                       type = "warning", duration = 5)
    }
    
    showModal(modalDialog(
      title = "Generating Gene Expression Plot",
      "Creating visualization of gene expression along pseudotime...",
      footer = NULL,
      easyClose = FALSE
    ))
    
    tryCatch({
      # Get the monocle object
      cds <- monocle_object()
      
      # Verify genes exist in the dataset
      available_genes <- rownames(cds)
      valid_genes <- input$gene_picker[input$gene_picker %in% available_genes]
      
      if (length(valid_genes) == 0) {
        removeModal()
        showNotification("Selected genes not found in dataset", type = "error")
        return()
      }
      
      if (length(valid_genes) < length(input$gene_picker)) {
        missing_genes <- setdiff(input$gene_picker, valid_genes)
        showNotification(paste("Warning: Some genes not found:", paste(missing_genes, collapse = ", ")), 
                         type = "warning")
      }
      
      # Log the analysis
      message(paste("Plotting", length(valid_genes), "genes along pseudotime"))
      
      # Create the plot with better parameters
      plot <- plot_genes_in_pseudotime(
        cds,
        valid_genes,
        min_expr = 0.1,  # Minimum expression threshold
        cell_size = 0.75,  # Size of points
        nrow = NULL,  # Let it auto-arrange
        ncol = min(3, length(valid_genes)),  # Max 3 columns
        panel_order = valid_genes,  # Keep gene order
        color_cells_by = "pseudotime",  # Color by pseudotime
        trend_formula = "~ splines::ns(pseudotime, df=3)"  # Smooth trend line
      )
      
      # Enhance plot appearance
      plot <- plot +
        theme_bw(base_size = 12) +
        theme(
          legend.position = "bottom",
          legend.title = element_text(size = 10, face = "bold"),
          strip.text = element_text(size = 11, face = "bold"),
          strip.background = element_rect(fill = "#E8F4FD", color = "black"),
          panel.spacing = unit(0.5, "lines"),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          axis.text = element_text(size = 9),
          axis.title = element_text(size = 10, face = "bold")
        ) +
        labs(
          title = "Gene Expression Dynamics Along Pseudotime",
          x = "Pseudotime",
          y = "Expression",
          color = "Pseudotime"
        ) +
        scale_color_viridis_c(option = "plasma")  # Match the pseudotime plot colors
      
      # If only one gene, make plot larger
      if (length(valid_genes) == 1) {
        plot <- plot + theme(aspect.ratio = 0.6)
      }
      
      # Store and display
      gene_trajectory_plot(plot)
      
      # Render with appropriate height based on number of genes
      plot_height <- 300 + (ceiling(length(valid_genes) / 3) * 200)
      
      output$geneTrajectoryPlot <- renderPlot({
        gene_trajectory_plot()
      }, height = plot_height)
      
      removeModal()
      
      # Show completion message
      showNotification(
        paste("Successfully plotted", length(valid_genes), "genes along pseudotime"), 
        type = "success"
      )
      
    }, error = function(e) {
      removeModal()
      message(paste("Error in gene trajectory plot:", e$message))
      
      # Try alternative plotting method if plot_genes_in_pseudotime fails
      if (grepl("plot_genes_in_pseudotime", e$message)) {
        showNotification("Trying alternative visualization method...", type = "message")
        
        tryCatch({
          # Alternative: Create manual plot
          cds <- monocle_object()
          
          # Extract expression data
          expr_data <- as.matrix(exprs(cds)[valid_genes, , drop = FALSE])
          pseudotime_data <- pseudotime(cds)
          
          # Create data frame for plotting
          plot_data <- data.frame(
            pseudotime = rep(pseudotime_data, each = length(valid_genes)),
            expression = as.vector(expr_data),
            gene = rep(valid_genes, ncol(expr_data))
          )
          
          # Remove NA pseudotime values
          plot_data <- plot_data[!is.na(plot_data$pseudotime), ]
          
          # Create plot manually
          plot <- ggplot(plot_data, aes(x = pseudotime, y = expression)) +
            geom_point(aes(color = pseudotime), size = 0.5, alpha = 0.6) +
            geom_smooth(method = "loess", se = TRUE, color = "red", size = 1) +
            facet_wrap(~ gene, scales = "free_y", ncol = min(3, length(valid_genes))) +
            theme_bw() +
            scale_color_viridis_c(option = "plasma") +
            labs(title = "Gene Expression Along Pseudotime",
                 x = "Pseudotime", 
                 y = "Expression") +
            theme(legend.position = "bottom",
                  strip.text = element_text(face = "bold"))
          
          gene_trajectory_plot(plot)
          output$geneTrajectoryPlot <- renderPlot(plot)
          
          showNotification("Created alternative visualization", type = "message")
          
        }, error = function(e2) {
          showNotification(paste("Error:", e2$message), type = "error")
        })
      } else {
        showNotification(paste("Error:", e$message), type = "error")
      }
    })
  })
  
  # Download handler for trajectory plot - IMPROVED
  output$download_trajectory_plot <- downloadHandler(
    filename = function() {
      # Include gene info in filename
      n_genes <- length(input$gene_picker)
      gene_part <- if(n_genes <= 3) {
        paste(input$gene_picker, collapse = "_")
      } else {
        paste0(n_genes, "_genes")
      }
      paste0("gene_trajectory_", gene_part, "_", Sys.Date(), ".", input$trajectory_download_format)
    },
    content = function(file) {
      req(gene_trajectory_plot())
      
      # Calculate dimensions based on number of genes
      n_genes <- length(input$gene_picker)
      plot_height <- 6 + (ceiling(n_genes / 3) * 2)
      plot_width <- min(12, 4 * min(3, n_genes))
      
      ggsave(
        file, 
        plot = gene_trajectory_plot(),
        width = plot_width, 
        height = plot_height, 
        dpi = input$dpi_selection,
        device = input$trajectory_download_format
      )
    }
  )
  # Reactive value to store the plot
  gene_path_plot <- reactiveVal()
  
  # Update gene list when monocle object changes
  observeEvent(monocle_object(), {
    req(monocle_object())
    
    # Get gene names from monocle object
    gene_names <- rowData(monocle_object())$gene_short_name
    if (is.null(gene_names) || length(gene_names) == 0) {
      gene_names <- rownames(rowData(monocle_object()))
    } else {
      # Replace NA values with rownames
      na_indices <- which(is.na(gene_names) | gene_names == "")
      if (length(na_indices) > 0) {
        gene_names[na_indices] <- rownames(rowData(monocle_object()))[na_indices]
      }
    }
    
    # Log gene count
    message(paste("Updated gene list with", length(gene_names), "genes"))
    
    # Update gene picker with better organization
    updatePickerInput(
      session, "gene_selection", 
      choices = gene_names, 
      selected = NULL,
      options = list(
        `live-search` = TRUE,
        `actions-box` = TRUE,
        title = "Select genes to visualize"
      )
    )
  })
  
  # Create reactive value to store selected genes
  selected_genes_monocle <- reactiveVal("")
  
  # Get Seurat object name for file naming
  get_seurat_object_name_monocle <- function() {
    obj_name <- tryCatch({
      if (!is.null(input$load_seurat_file_monocle)) {
        base_name <- tools::file_path_sans_ext(basename(input$load_seurat_file_monocle$name))
        if (is.null(base_name) || base_name == "") {
          "SeuratObject"
        } else {
          # Clean filename - remove special characters
          gsub("[^A-Za-z0-9_-]", "_", base_name)
        }
      } else {
        "SeuratObject"
      }
    }, error = function(e) {
      "SeuratObject"
    })
    return(obj_name)
  }
  
  # Generate gene visualization plots - VERSION ADAPTÉE
  observeEvent(input$generate_gene_path, {
    req(monocle_object(), input$gene_selection)
    
    if (length(input$gene_selection) == 0) {
      showNotification("Please select at least one gene to visualize", type = "warning")
      return()
    }
    
    if (length(input$gene_selection) > 12) {
      showNotification("Warning: Plotting more than 12 genes may be hard to read", 
                       type = "warning", duration = 5)
    }
    
    showModal(modalDialog(
      title = "Generating Plot",
      paste("Creating visualization for", length(input$gene_selection), "gene(s)..."),
      footer = NULL,
      easyClose = FALSE
    ))
    
    tryCatch({
      cds <- monocle_object()
      
      # Verify genes exist in dataset
      valid_genes <- input$gene_selection[input$gene_selection %in% rownames(rowData(cds))]
      
      if (length(valid_genes) == 0) {
        removeModal()
        showNotification("None of the selected genes were found in the dataset", type = "error")
        return()
      }
      
      if (length(valid_genes) < length(input$gene_selection)) {
        missing <- setdiff(input$gene_selection, valid_genes)
        showNotification(paste("Warning:", length(missing), "gene(s) not found:", 
                               paste(head(missing, 3), collapse = ", ")), 
                         type = "warning")
      }
      
      message(paste("Plotting", length(valid_genes), "genes on trajectory"))
      
      # Store selected genes for download handler
      selected_genes_str <- paste(valid_genes, collapse = "_")
      selected_genes_monocle(selected_genes_str)
      
      legend_setting <- theme(legend.position = "right")  
      
      # Check if trajectory graph should be shown
      show_trajectory <- !is.null(cds@principal_graph) && length(cds@principal_graph) > 0
      
      # CELL SIZE FIXE (plus de input$cell_size)
      fixed_cell_size <- 0.5
      
      # Create plot based on number of genes
      if (length(valid_genes) == 1) {
        # Single gene plot
        gene <- valid_genes[1]
        
        plot <- plot_cells(
          cds,
          genes = gene,
          show_trajectory_graph = show_trajectory,
          label_cell_groups = FALSE,
          label_leaves = FALSE,
          label_branch_points = show_trajectory,
          cell_size = fixed_cell_size * 1.5,  # Un peu plus grand pour un seul gène
          cell_stroke = 0.1,
          trajectory_graph_color = "black",
          trajectory_graph_segment_size = 0.75
        ) +
          ggtitle(paste("Expression of", gene)) +
          theme_minimal(base_size = 12) +
          theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
            axis.text = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank(),
            axis.ticks = element_blank()
          ) +
          legend_setting
        
      } else {
        # Multi-gene visualization
        plots <- list()
        
        for (i in seq_along(valid_genes)) {
          gene <- valid_genes[i]
          
          p <- plot_cells(
            cds,
            genes = gene,
            show_trajectory_graph = show_trajectory,
            label_cell_groups = FALSE,
            label_leaves = FALSE,
            label_branch_points = FALSE,
            cell_size = fixed_cell_size,  # Taille fixe
            cell_stroke = 0,
            trajectory_graph_color = "black",
            trajectory_graph_segment_size = 0.5
          ) +
            ggtitle(gene) +
            theme_minimal(base_size = 10) +
            theme(
              plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
              axis.text = element_blank(),
              axis.title = element_blank(),
              panel.grid = element_blank(),
              axis.ticks = element_blank(),
              plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm")
            ) +
            legend_setting
          
          plots[[i]] <- p
        }
        
        # Calculate optimal layout
        n_plots <- length(plots)
        n_cols <- ifelse(n_plots <= 3, n_plots, ceiling(sqrt(n_plots)))
        
        # Combine plots using patchwork
        plot <- wrap_plots(plots, ncol = n_cols) +
          plot_annotation(
            title = "Gene Expression on Trajectory",
            subtitle = paste(length(valid_genes), "genes visualized"),
            theme = theme(
              plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
              plot.subtitle = element_text(size = 11, hjust = 0.5)
            )
          )
      }
      
      # Store and display the plot
      gene_path_plot(plot)
      
      removeModal()
      showNotification(paste("Successfully visualized", length(valid_genes), "gene(s)"), 
                       type = "message")
      
    }, error = function(e) {
      removeModal()
      message(paste("Error in gene visualization:", e$message))
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # Render the plot with dynamic height
  output$gene_path_plot_output <- renderPlot({
    req(gene_path_plot())
    gene_path_plot()
  }, height = function() {
    # Dynamic height based on number of genes
    if (!is.null(input$gene_selection)) {
      n_genes <- length(input$gene_selection)
      if (n_genes == 1) {
        return(500)
      } else {
        n_rows <- ceiling(n_genes / ifelse(n_genes <= 3, n_genes, ifelse(n_genes <= 6, 3, 4)))
        return(200 + n_rows * 250)
      }
    }
    return(500)
  })
  
  # Download handler with gene names and format choice - IMPROVED
  output$download_gene_path_plot <- downloadHandler(
    filename = function() {
      # Get gene info for filename
      n_genes <- length(strsplit(selected_genes_monocle(), "_")[[1]])
      genes_part <- if(n_genes <= 3) {
        gsub("_", "-", selected_genes_monocle())
      } else {
        paste0(n_genes, "genes")
      }
      
      # Truncate if too long
      if(nchar(genes_part) > 50) {
        genes_part <- paste0(substr(genes_part, 1, 47), "...")
      }
      
      # Get the Seurat object name
      seurat_part <- get_seurat_object_name_monocle()
      
      # Create filename
      paste0(seurat_part, "_GeneTrajectory_", genes_part, "_",
             format(Sys.time(), "%Y%m%d"), ".", input$plot_format_monocle)
    },
    content = function(file) {
      req(gene_path_plot())
      
      # DIMENSIONS FIXES (plus de input$plot_width_monocle et plot_height_monocle)
      n_genes <- length(strsplit(selected_genes_monocle(), "_")[[1]])
      
      # Dimensions automatiques selon le nombre de gènes
      if (n_genes == 1) {
        width <- 8
        height <- 6
      } else {
        n_cols <- ifelse(n_genes <= 3, n_genes, ceiling(sqrt(n_genes)))
        n_rows <- ceiling(n_genes / n_cols)
        width <- min(16, 4 * n_cols)
        height <- min(20, 3 + 3 * n_rows)
      }
      
      tryCatch({
        ggsave(
          file, 
          plot = gene_path_plot(),
          dpi = input$plot_dpi_monocle,
          width = width, 
          height = height,
          device = input$plot_format_monocle,
          limitsize = FALSE
        )
        showNotification("Plot successfully downloaded!", type = "message")
      }, error = function(e) {
        showNotification(paste("Error saving plot:", e$message), type = "error")
      })
    }
  )
  ############################## Gene modules along trajectory ##############################


  # Reactive values
  gene_module_df <- reactiveVal(NULL)
  module_heatmap_plot <- reactiveVal(NULL)
  module_visualization_plot <- reactiveVal(NULL)
  
  # Find gene modules along trajectory - VERSION SIMPLIFIÉE
  observeEvent(input$find_gene_modules, {
    req(monocle_object(), pr_deg_ids())
    
    if (length(pr_deg_ids()) == 0) {
      showNotification("No differentially expressed genes found. Run differential gene test first.", type = "error")
      return()
    }
    
    showModal(modalDialog(
      title = "Finding Gene Modules",
      "Grouping genes into modules based on expression patterns...",
      footer = NULL,
      easyClose = FALSE
    ))
    
    tryCatch({
      # Get monocle object and DE genes
      cds <- monocle_object()
      deg_genes <- pr_deg_ids()
      
      message("Starting module analysis...")
      message(paste("Total DE genes:", length(deg_genes)))
      
      # Ensure valid genes
      all_genes <- rownames(cds)
      deg_genes <- deg_genes[deg_genes %in% all_genes]
      
      # Limit genes
      max_genes <- min(2000, length(deg_genes))
      if (length(deg_genes) > max_genes) {
        message(paste("Limiting to", max_genes, "genes"))
        deg_genes <- deg_genes[1:max_genes]
      }
      
      # Try Monocle's find_gene_modules first
      modules <- NULL
      
      tryCatch({
        # Subset to DE genes
        cds_subset <- cds[deg_genes, ]
        
        # Use find_gene_modules with simple parameters
        modules <- find_gene_modules(
          cds_subset,
          resolution = input$resolution_value,
          random_seed = 123
        )
        
      }, error = function(e) {
        message(paste("Monocle module finding failed:", e$message))
        message("Using alternative clustering approach...")
        
        # Alternative: Simple clustering based on expression patterns
        # Extract normalized expression matrix
        norm_counts <- normalized_counts(cds)
        gene_expr <- as.matrix(norm_counts[deg_genes, ])
        
        # Remove genes with no variation
        gene_vars <- apply(gene_expr, 1, var)
        keep_genes <- gene_vars > 0
        gene_expr <- gene_expr[keep_genes, ]
        kept_genes <- deg_genes[keep_genes]
        
        if (nrow(gene_expr) < 10) {
          showNotification("Too few variable genes for module analysis", type = "error")
          removeModal()
          return()
        }
        
        # Log transform
        gene_expr <- log1p(gene_expr)
        
        # Calculate correlation matrix (subset if too many genes)
        if (nrow(gene_expr) > 500) {
          set.seed(123)
          sample_idx <- sample(1:nrow(gene_expr), 500)
          cor_mat <- cor(t(gene_expr[sample_idx, ]), use = "complete.obs")
        } else {
          cor_mat <- cor(t(gene_expr), use = "complete.obs")
        }
        
        # Replace NA with 0
        cor_mat[is.na(cor_mat)] <- 0
        
        # Hierarchical clustering
        dist_mat <- as.dist(1 - cor_mat)
        hc <- hclust(dist_mat, method = "ward.D2")
        
        # Cut tree to get modules
        n_modules <- min(10, max(3, floor(length(kept_genes) / 20)))
        clusters <- cutree(hc, k = n_modules)
        
        # Create module dataframe
        modules <- data.frame(
          id = names(clusters),
          module = as.character(clusters),
          stringsAsFactors = FALSE
        )
      })
      
      if (is.null(modules) || nrow(modules) == 0) {
        removeModal()
        showNotification("Failed to find gene modules", type = "error")
        return()
      }
      
      # Process results
      modules$module <- as.character(modules$module)
      n_modules <- length(unique(modules$module))
      
      message(paste("Found", n_modules, "gene modules"))
      
      # Store results
      gene_module_df(modules)
      
      # Update UI
      module_choices <- sort(unique(modules$module))
      updateSelectInput(session, "selected_modules",
                        choices = module_choices,
                        selected = module_choices[1:min(3, length(module_choices))])
      
      # Create summary
      module_summary <- data.frame(
        Module = paste("Module", sort(unique(modules$module))),
        Genes = as.numeric(table(modules$module))
      )
      
      output$module_summary_table <- renderDT({
        datatable(
          module_summary,
          options = list(pageLength = 10, dom = 'tip'),
          rownames = FALSE
        )
      })
      
      removeModal()
      showNotification(paste("Found", n_modules, "gene modules!"), type = "message")
      
    }, error = function(e) {
      removeModal()
      message(paste("Error in module finding:", e$message))
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  # Create module heatmap - IMPROVED
  observeEvent(input$generate_module_heatmap, {
    req(monocle_object(), gene_module_df())
    
    showModal(modalDialog(
      title = "Generating Module Heatmap",
      "Creating heatmap of module expression patterns...",
      footer = NULL,
      easyClose = FALSE
    ))
    
    tryCatch({
      cds <- monocle_object()
      modules <- gene_module_df()
      
      # Determine grouping column (prefer ClusterIdents)
      group_col <- NULL
      if ("ClusterIdents" %in% colnames(colData(cds))) {
        group_col <- "ClusterIdents"
      } else if ("seurat_clusters" %in% colnames(colData(cds))) {
        group_col <- "seurat_clusters"
      } else {
        group_col <- colnames(colData(cds))[1]
      }
      
      message(paste("Using grouping column:", group_col))
      
      # Create cell grouping dataframe
      cell_group_df <- data.frame(
        cell = colnames(cds),
        cell_group = as.character(colData(cds)[[group_col]])
      )
      
      # Aggregate expression by module
      agg_mat <- aggregate_gene_expression(
        cds,
        modules,
        cell_group_df
      )
      
      # Ensure we have a matrix
      if (!is.matrix(agg_mat)) {
        agg_mat <- as.matrix(agg_mat)
      }
      
      # Add module labels
      rownames(agg_mat) <- paste("Module", rownames(agg_mat))
      
      # Create color annotation for clusters
      n_groups <- ncol(agg_mat)
      cluster_colors <- colorRampPalette(c("#3366CC", "#DC3912", "#FF9900", "#109618", 
                                           "#990099", "#0099C6", "#DD4477"))(n_groups)
      names(cluster_colors) <- colnames(agg_mat)
      
      # Create heatmap with better colors
      heatmap <- pheatmap::pheatmap(
        agg_mat,
        scale = "row",  # Scale by row (module)
        clustering_method = "ward.D2",
        color = colorRampPalette(c("blue", "white", "red"))(100),
        main = paste("Module Expression by", group_col),
        fontsize = 10,
        fontsize_row = 9,
        fontsize_col = 9,
        annotation_colors = list(cluster = cluster_colors),
        border_color = NA,
        show_rownames = TRUE,
        show_colnames = TRUE,
        angle_col = 45
      )
      
      # Store plot
      module_heatmap_plot(heatmap)
      
      # Display plot
      output$module_heatmap <- renderPlot({
        grid::grid.newpage()
        grid::grid.draw(heatmap$gtable)
      })
      
      removeModal()
      showNotification("Module heatmap generated successfully!", type = "success")
      
    }, error = function(e) {
      removeModal()
      message(paste("Error generating heatmap:", e$message))
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # Visualize selected modules - IMPROVED
  observeEvent(input$visualize_modules, {
    req(monocle_object(), gene_module_df(), input$selected_modules)
    
    if (length(input$selected_modules) == 0) {
      showNotification("Please select at least one module to visualize", type = "warning")
      return()
    }
    
    showModal(modalDialog(
      title = "Visualizing Modules",
      paste("Creating visualization for", length(input$selected_modules), "module(s)..."),
      footer = NULL,
      easyClose = FALSE
    ))
    
    tryCatch({
      cds <- monocle_object()
      modules <- gene_module_df()
      
      # Get top genes from each selected module
      genes_per_module <- input$genes_per_module
      top_genes_list <- list()
      
      for (mod in input$selected_modules) {
        module_genes <- modules[modules$module == mod, "id"]
        
        # If we have expression data, sort by variance
        if (nrow(module_genes) > genes_per_module) {
          gene_vars <- apply(exprs(cds)[module_genes$id, ], 1, var)
          top_idx <- order(gene_vars, decreasing = TRUE)[1:genes_per_module]
          top_genes <- module_genes$id[top_idx]
        } else {
          top_genes <- module_genes$id
        }
        
        top_genes_list[[paste0("Module_", mod)]] <- top_genes
      }
      
      # Flatten the list
      all_genes <- unlist(top_genes_list)
      gene_labels <- rep(names(top_genes_list), sapply(top_genes_list, length))
      
      message(paste("Plotting", length(all_genes), "genes from", length(input$selected_modules), "modules"))
      
      # Create plot
      plot <- plot_genes_in_pseudotime(
        cds,
        all_genes,
        min_expr = 0.1,
        ncol = min(3, length(all_genes)),
        panel_order = all_genes,
        color_cells_by = "pseudotime",
        trend_formula = "~ splines::ns(pseudotime, df=3)"
      ) +
        ggtitle(paste("Top Genes from Selected Modules")) +
        theme_bw(base_size = 10) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
          strip.text = element_text(face = "bold", size = 9),
          strip.background = element_rect(fill = "#E8F4FD"),
          legend.position = "bottom"
        ) +
        scale_color_viridis_c(option = "plasma")
      
      # Add module labels as facet titles if possible
      
      # Store and display
      module_visualization_plot(plot)
      output$module_visualization <- renderPlot({
        plot
      })
      
      removeModal()
      showNotification("Module visualization created!", type = "success")
      
    }, error = function(e) {
      removeModal()
      message(paste("Error visualizing modules:", e$message))
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # Download handlers remain the same but with better error handling
  output$download_module_heatmap <- downloadHandler(
    filename = function() {
      paste0("module_heatmap_", Sys.Date(), ".", input$module_download_format)
    },
    content = function(file) {
      req(module_heatmap_plot())
      
      if (input$module_download_format == "pdf") {
        pdf(file, width = input$module_plot_width, height = input$module_plot_height)
      } else {
        png(file, width = input$module_plot_width, height = input$module_plot_height, 
            units = "in", res = input$module_plot_dpi)
      }
      
      grid::grid.newpage()
      grid::grid.draw(module_heatmap_plot()$gtable)
      dev.off()
    }
  )
  
  output$download_module_genes <- downloadHandler(
    filename = function() {
      paste0("gene_modules_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(gene_module_df())
      modules <- gene_module_df()
      modules$module <- paste("Module", modules$module)
      write.csv(modules, file, row.names = FALSE)
    }
  )
}
