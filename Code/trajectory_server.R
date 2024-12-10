############################## Differential expressed genes along the trajectory ############################## 

# trajectory_server.R


trajectory_server <- function(input, output, session) {
  
  shinyjs::useShinyjs() # button deactivation
  
  ############################## Monocle Conversion and trajectory ############################## 
  
  seurat_monocle <- reactiveVal() # Seurat object loaded by the user
  monocle_object <- reactiveVal() # Monocle object converted by the user from a seurat object
  selected_cell_id <- reactiveVal(NULL) # Cell id to connect the pointer and the id of cell when the user is clicking on the umap for the trajectory
  
  # Désactiver initialement les boutons
  shinyjs::disable("convertToMonocle")
  shinyjs::disable("constructGraph")
  shinyjs::disable("startRootSelection")
  
  # For Seurat loading
  observeEvent(input$load_seurat_file_monocle, {
    req(input$load_seurat_file_monocle) # Verify file input exists
    
    # Show loading modal
    showModal(modalDialog(
      title = "Loading Seurat Object",
      div(
        style = "text-align: center;",
        h4("Please wait while the Seurat object is being loaded..."),
        p("This process includes:", style = "text-align: left;"),
        tags$ul(
          style = "text-align: left;",
          tags$li("Validating file format"),
          tags$li("Loading object into memory"),
          tags$li("Checking required dimensional reductions (PCA, UMAP)"),
          tags$li("Preparing object for trajectory analysis")
        ),
        p("This may take a few moments depending on the size of your dataset.")
      ),
      footer = NULL,
      easyClose = FALSE
    ))
    
    # Validate file extension
    if(!grepl("\\.rds$", input$load_seurat_file_monocle$name, ignore.case = TRUE)) {
      removeModal()
      showNotification("Please select a .rds file", type = "error")
      return()
    }
    
    # Check file path exists
    if(!file.exists(input$load_seurat_file_monocle$datapath)) {
      removeModal()
      showNotification("File not found", type = "error")
      return()
    }
    
    message("Attempting to read file at: ", input$load_seurat_file_monocle$datapath)
    
    tryCatch({
      temp_seurat <- readRDS(input$load_seurat_file_monocle$datapath)
      
      # Validate it's a Seurat object
      if(!inherits(temp_seurat, "Seurat")) {
        removeModal()
        showNotification("File does not contain a valid Seurat object", type = "error")
        return()
      }
      
      # Check required reductions exist
      if(!all(c("pca", "umap") %in% names(temp_seurat@reductions))) {
        removeModal()
        showNotification("Seurat object must contain PCA and UMAP reductions", type = "error")
        return()
      }
      
      seurat_monocle(temp_seurat)
      removeModal()
      showNotification("Seurat object has been successfully loaded!")
      shinyjs::enable("convertToMonocle")
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Error loading file:", e$message), type = "error")
    })
  })
  
  # Stocker les chargements de gènes dans une variable globale
  gene_loadings <- NULL
  
  
  
  # Fonction pour initialiser les chargements de gènes
  initialize_gene_loadings <- function(monocle_obj, gene_loadings_matrix) {
    gene_names <- rowData(monocle_obj)$gene_short_name
    if (!is.null(gene_names)) {
      # Trouver les noms de gènes qui sont dans gene_loadings_matrix
      common_genes <- intersect(gene_names, rownames(gene_loadings_matrix))
      if (length(common_genes) == nrow(gene_loadings_matrix)) {
        rownames(gene_loadings_matrix) <- common_genes
        print("Initializing gene loadings matrix with gene names:")
        print(head(gene_loadings_matrix))
        gene_loadings <<- gene_loadings_matrix
      } else {
        showNotification("Mismatch between common gene names length and gene loadings matrix rows.", type = "error")
      }
    } else {
      showNotification("Gene names not found in Monocle object.", type = "error")
    }
  }
  
  # Fonction de conversion de Seurat à Monocle
  observeEvent(input$convertToMonocle, {
    req(seurat_monocle())
    
    tryCatch({
      message("Début de la conversion de l'objet Seurat en objet Monocle.")
      
      # Initialisation de l'objet Monocle à partir de Seurat
      sce <- as.cell_data_set(seurat_monocle(), default.reduction = "pca")
      message("Conversion en cell_data_set réussie.")
      
      # Conversion en cell_data_set pour manipulation
      monocle_temp <- as(sce, "cell_data_set")
      message("Conversion en cell_data_set réussie (vérification).")
      
      # Transfert des données PCA et UMAP
      reducedDims(monocle_temp)[["PCA"]] <- Embeddings(seurat_monocle(), "pca")
      reducedDims(monocle_temp)[["UMAP"]] <- Embeddings(seurat_monocle(), "umap")
      monocle_temp@reduce_dim_aux[['PCA']][['model']][['svd_v']] <- seurat_monocle()@reductions[["pca"]]@feature.loadings.projected
      monocle_temp@reduce_dim_aux@listData[["UMAP"]] <- Embeddings(seurat_monocle(), "umap")
      
      # Estimation des facteurs de taille avec Monocle
      monocle_temp <- estimate_size_factors(monocle_temp)
      message("Estimation des facteurs de taille réussie.")
      
      # Ajout des noms de gènes à l'objet Monocle
      rowData(monocle_temp)$gene_short_name <- rownames(rowData(monocle_temp))
      message("Noms de gènes ajoutés à l'objet Monocle.")
      
      # Vérification et initialisation des gene_loadings
      gene_loadings_matrix <- seurat_monocle()@reductions[["pca"]]@feature.loadings
      initialize_gene_loadings(monocle_temp, gene_loadings_matrix)
      
      # Mise à jour de l'objet Monocle global
      monocle_object(monocle_temp)
      message("Objet Monocle converti avec succès et stocké.")
      
      # Enable the buttons after conversion
      shinyjs::enable("constructGraph")
      showNotification("Converted to Monocle object successfully!")
    }, error = function(e) {
      message("Log: Error during conversion: ", e$message)
      showNotification(paste("An error occurred during conversion: ", e$message), type = "error")
    })
  })
  
  cell_data <- reactive({
    data <- as.data.frame(reducedDims(monocle_object())$UMAP)
    data$cell_id <- as.character(rownames(data))
    return(data)
  })
  
  get_cell_id_from_coordinates <- function(x, y, cell_data_func) {
    data <- cell_data_func()
    distances <- rowSums((data[,1:2] - c(x, y))^2)
    closest_cell <- which.min(distances)
    return(data$cell_id[closest_cell])
  }
  
  observeEvent(input$constructGraph, {
    req(monocle_object())
    
    showModal(modalDialog(
      title = "Please wait",
      "Constructing the graph. This might take a while...",
      easyClose = FALSE,
      footer = NULL
    ))
    
    tryCatch({
      monocle_temp <- monocle_object()
      
      # Clustering cells
      k_value <- min(as.integer(input$cluster_k), ncol(monocle_temp) - 3)
      monocle_temp <- cluster_cells(monocle_temp, reduction_method = "UMAP", k = k_value)
      
      # Check UMAP data
      if (is.null(reducedDims(monocle_temp)$UMAP)) {
        stop("UMAP reduction data is NULL")
      }
      
      # Learn graph
      monocle_temp <- learn_graph(monocle_temp, use_partition = FALSE)
      monocle_object(monocle_temp)
      
      # Create plotly plot directly
      umap_coords <- reducedDims(monocle_temp)$UMAP
      clusters <- monocle_temp@clusters$UMAP$clusters
      
      # Create plotly object directly
      plot_data <- plot_ly(
        x = umap_coords[,1],
        y = umap_coords[,2],
        type = "scatter",
        mode = "markers",
        marker = list(
          size = 3,  # Taille réduite des points
          opacity = 0.7  # Ajout de transparence
        ),        color = as.factor(clusters),
        source = "trajectoryPlot"
      ) %>%
        layout(
          title = "UMAP Trajectory",
          xaxis = list(title = "UMAP 1"),
          yaxis = list(title = "UMAP 2"),
          showlegend = TRUE,
          height = 500
        ) %>%
        event_register("plotly_click")
      
      output$trajectoryPlot <- renderPlotly({
        plot_data
      })
      
      showNotification("Graph construction completed!")
      shinyjs::enable("startRootSelection")
      
    }, error = function(e) {
      message("Log: Error during graph construction: ", e$message)
      showNotification(paste("An error occurred during graph construction: ", e$message), type = "error")
    }, finally = {
      removeModal()
    })
  })
  # Gestion améliorée de la sélection de la racine
  root_selection_started <- reactiveVal(FALSE)
  
  # Réinitialiser la sélection quand on clique sur le bouton
  observeEvent(input$startRootSelection, {
    root_selection_started(FALSE)  # Reset d'abord
    Sys.sleep(0.1)  # Petit délai pour assurer la réinitialisation
    root_selection_started(TRUE)   # Puis activation
    showNotification("Please click on the graphic to select the root cell.")
  })
  
  # Observer séparé pour les clics sur le plot
  observeEvent(event_data("plotly_click", source = "trajectoryPlot"), {
    req(root_selection_started())
    clicked_point <- event_data("plotly_click", source = "trajectoryPlot")
    
    if (!is.null(clicked_point)) {
      sel_id <- get_cell_id_from_coordinates(clicked_point$x, clicked_point$y, cell_data)
      selected_cell_id(sel_id)
      
      # Afficher la boîte de dialogue de confirmation
      showModal(modalDialog(
        title = "Confirm root selection",
        paste("Have you selected the right point with the cell ID:", sel_id, "?"),
        footer = tagList(
          actionButton("confirmRoot", "Confirm"),
          modalButton("Cancel")
        )
      ))
    }
  }, ignoreInit = TRUE)
  
  # Gestion de la confirmation
  observeEvent(input$confirmRoot, {
    removeModal()
    req(monocle_object(), selected_cell_id())
    
    # Créer le plot pseudotime
    monocle_temp <- order_cells(monocle_object(), 
                                reduction_method = "UMAP", 
                                root_cells = selected_cell_id())
    
    output$pseudotimePlot <- renderPlot({
      plot_cells(monocle_temp,
                 color_cells_by = "pseudotime",
                 cell_size = 0.5,
                 label_groups_by_cluster = FALSE,
                 label_leaves = FALSE,
                 label_branch_points = FALSE,
                 graph_label_size = 4)
    })
    
    showNotification("Cell ordering completed!")
    root_selection_started(FALSE)  # Reset après confirmation
  })
  
  # Reset de la sélection si on annule
  observeEvent(input$modalButton, {
    root_selection_started(FALSE)
    selected_cell_id(NULL)
  })
  
  ############################## Differential expressed genes along the trajectory ############################## 
  
  pseudotime_diff_gene_results <- reactiveVal()
  gene_modules_result <- reactiveVal()  
  pr_deg_ids <- reactiveVal() 
  
  # Fonction modifiée pour calculer tous les marqueurs différentiels une seule fois
  pr_deg_ids <- reactiveVal(NULL)
  
  # Ajout dans la fonction CalculateDiffGenesPseudotime pour mettre à jour gene_list
  CalculateDiffGenesPseudotime <- function() {
    if (is.null(monocle_object())) {
      showNotification("Monocle object is NULL, please check the previous steps.", type = "error")
      return(NULL)
    }
    tryCatch({
      showNotification("Calculation of differential markers in progress...", type = "message", duration = 10)
      markers <- graph_test(monocle_object(), neighbor_graph = "principal_graph", cores = 4)
      ids <- row.names(subset(markers, q_value < 0.05))
      filtered_results <- subset(markers, row.names(markers) %in% ids)
      
      print("IDs de gènes différentiels trouvés:")
      print(ids)
      
      # Stocker les résultats et les IDs de gènes différentiels dans des variables réactives
      pseudotime_diff_gene_results(filtered_results)
      pr_deg_ids(ids)
      updatePickerInput(session, "gene_picker", choices = row.names(filtered_results), selected = row.names(filtered_results)[1])
      updatePickerInput(session, "gene_list", choices = row.names(filtered_results), selected = row.names(filtered_results)[1]) # Ajout pour le nouveau pickerInput
      
    }, error = function(e) {
      showNotification(paste0("Error when calculating differential markers: ", e$message), type = "error")
    })
  }
  
  
  observeEvent(input$run_diff_gene_pseudotime, {
    tryCatch({
      showModal(modalDialog(
        title = "Please Wait",
        "Calculating differential genes based on pseudotime...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      CalculateDiffGenesPseudotime()
      
      removeModal()
    }, error = function(e) {
      removeModal()
      showNotification(paste0("Error during calculation of differential genes: ", e$message), type = "error")
    })
  })
  
  observe({
    output$diffGeneTable <- renderDataTable({
      datatable(pseudotime_diff_gene_results(), escape = FALSE, options = list(pageLength = 10))
    })
  })
  
  output$download_pseudotime_diff_genes <- downloadHandler(
    filename = function() {
      paste("pseudotime-diff-genes-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      diff_gene_results <- pseudotime_diff_gene_results()
      req(!is.null(diff_gene_results))
      
      write.csv(diff_gene_results, file, row.names = FALSE)
    },
    contentType = "text/csv"
  )
  
  gene_trajectory_plot <- reactiveVal()
  
  observeEvent(input$visualize_gene_trajectory, {
    tryCatch({
      req(monocle_object(), input$gene_picker)
      plot <- plot_cells(monocle_object(), genes = input$gene_picker, show_trajectory_graph = FALSE, label_cell_groups = FALSE, label_leaves = FALSE, cell_size = 1)
      gene_trajectory_plot(plot)
      print("Trajectory performed.")
    }, error = function(e) {
      showNotification(paste0("Trajectory error: ", e$message), type = "error")
    })
  })
  
  output$geneTrajectoryPlot <- renderPlot({
    req(gene_trajectory_plot())
    gene_trajectory_plot()
  })
  
  output$download_trajectory_plot <- downloadHandler(
    filename = function() {
      paste("gene_trajectory_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      req(gene_trajectory_plot())
      ggsave(file, plot = gene_trajectory_plot(), dpi = input$dpi_selection, width = 10, height = 6)
    }
  )
  ############################## Genes along path ##############################
  
  gene_path_plot <- reactiveVal()
  
  # Mettre à jour les pickerInput pour les gènes et les types de cellules
  observeEvent(monocle_object(), {
    req(monocle_object())
    update_picker_inputs(session, monocle_object())
  })
  
  # Function to update pickerInput for cell types using the correct columns
  update_picker_inputs <- function(session, monocle_obj) {
    if (!is.null(monocle_obj) && !is.null(rowData(monocle_obj)) && !is.null(colData(monocle_obj))) {
      # Extract unique values for "seurat_clusters" and "ident"
      seurat_clusters_choices <- if ("seurat_clusters" %in% colnames(colData(monocle_obj))) {
        unique(colData(monocle_obj)$seurat_clusters)
      } else {
        NULL
      }
      
      ident_choices <- if ("ident" %in% colnames(colData(monocle_obj))) {
        unique(colData(monocle_obj)$ident)
      } else {
        NULL
      }
      
      # Combine choices for the picker input
      cell_type_choices <- c(seurat_clusters_choices, ident_choices)
      
      # Update the picker input with the available choices
      updatePickerInput(session, "cell_type", choices = cell_type_choices)
    } else {
      showNotification("Monocle object is not properly loaded or missing required data.", type = "error")
    }
  }
  
  
  
  # Update cell types based on the user's choice between "seurat_clusters" and "ident"
  observeEvent(input$cell_type_choice, {
    req(monocle_object())
    monocle_obj <- monocle_object()
    
    if (input$cell_type_choice == "seurat_clusters" && "seurat_clusters" %in% colnames(colData(monocle_obj))) {
      cell_type_choices <- unique(colData(monocle_obj)$seurat_clusters)
    } else if (input$cell_type_choice == "ident" && "ident" %in% colnames(colData(monocle_obj))) {
      cell_type_choices <- unique(colData(monocle_obj)$ident)
    } else {
      cell_type_choices <- NULL
      showNotification("No valid cell type column found in Monocle object.", type = "error")
    }
    
    updatePickerInput(session, "cell_type", choices = cell_type_choices)
  })
  
  
  

  
  
  




 
}