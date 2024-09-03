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
    message("Attempting to read file at: ", input$load_seurat_file_monocle$datapath)
    
    tryCatch({
      temp_seurat <- readRDS(input$load_seurat_file_monocle$datapath)
      message("File successfully read.")
      message("Log: Seurat object loaded.")
      
      seurat_monocle(temp_seurat)
      showNotification("Seurat object has been successfully loaded!")
      
      shinyjs::enable("convertToMonocle")
    }, error = function(e) {
      message("Log: Error loading Seurat object: ", e$message)
      showNotification(paste("An error occurred: ", e$message), type = "error")
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
      
      # Clustering des cellules pour Monocle avec méthode de réduction UMAP
      k_value <- as.integer(input$cluster_k)
      monocle_temp <- cluster_cells(monocle_temp, reduction_method = "UMAP", k = k_value)
      message("Clustering des cellules avec méthode UMAP réussi.")
      
      # Vérifier les données nécessaires pour `learn_graph`
      if (is.null(reducedDims(monocle_temp)$UMAP)) {
        stop("UMAP reduction data is NULL")
      }
      
      # Apprendre le graphe
      monocle_temp <- learn_graph(monocle_temp)
      
      # Vérifier si le graphe a été construit
      if (is.null(monocle_temp@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex)) {
        stop("Graph construction failed, adjacency matrix is NULL")
      }
      
      monocle_object(monocle_temp)
      
      p <- plot_cells(monocle_temp, color_cells_by = "cluster", cell_size = 0.5, group_label_size = 3, label_groups_by_cluster = F, label_leaves = FALSE, label_branch_points = FALSE, graph_label_size = 5)
      p <- p + theme_void() # Pour retirer les éléments inutiles
      
      p_plotly <- ggplotly(p, source = "trajectoryPlot") %>% event_register("plotly_click")
      
      output$trajectoryPlot <- renderPlotly({
        p_plotly
      })
      
      showNotification("Graph construction completed!")
      cat("Log: Graph construction completed.\n")
      shinyjs::enable("startRootSelection")
    }, error = function(e) {
      message("Log: Error during graph construction: ", e$message)
      showNotification(paste("An error occurred during graph construction: ", e$message), type = "error")
    }, finally = {
      removeModal()
    })
  })
  
  root_selection_started <- reactiveVal(FALSE)
  
  observeEvent(input$startRootSelection, {
    showNotification("Please click on the graphic to select the root cell.")
    root_selection_started(TRUE)
  })
  
  observeEvent(event_data("plotly_click", source = "trajectoryPlot"), {
    if(root_selection_started()) {
      clicked_point <- event_data("plotly_click", source = "trajectoryPlot")
      cat("Clicked coordinates: x = ", clicked_point$x, ", y = ", clicked_point$y, "\n")
      print(head(cell_data()))
      
      sel_id <- get_cell_id_from_coordinates(clicked_point$x, clicked_point$y, cell_data)
      cat("Selected ID from coordinates: ", sel_id, "\n")
      
      selected_cell_id(sel_id)
      cat("Log: Selected root cell ID: ", sel_id, "\n")
      
      if (is.null(sel_id)) {
        cat("Log: No root cell selected.\n")
      } else {
        showModal(modalDialog(
          title = "Confirm root selection",
          paste("Have you selected the right point with the cell ID:", sel_id, "?"),
          footer = tagList(
            actionButton("confirmRoot", "Confirm"),
            modalButton("Annuler")
          )
        ))
      }
    }
  })
  
  observeEvent(input$confirmRoot, {
    removeModal()
    root_selection_started(FALSE)
    
    req(monocle_object())
    if (is.null(selected_cell_id())) {
      showNotification("No root cell selected!", type = "error")
      cat("Log: No root cell selected.\n")
      return()
    }
    
    monocle_temp <- order_cells(monocle_object(), reduction_method = "UMAP", root_cells = selected_cell_id())
    output$pseudotimePlot <- renderPlot({
      plot_cells(monocle_temp, color_cells_by = "pseudotime", cell_size = 0.5, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE, graph_label_size = 4)
    })
    
    showNotification("Cell ordering completed!")
    cat("Log: Cell ordering completed.\n")
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
  
  # Fonction pour mettre à jour les pickerInput
  update_picker_inputs <- function(session, monocle_obj) {
    if (!is.null(monocle_obj) && !is.null(rowData(monocle_obj)) && !is.null(colData(monocle_obj))) {
      dataset_choices <- unique(colData(monocle_obj)$dataset)
      cluster_choices <- unique(colData(monocle_obj)$cluster)
      
      updatePickerInput(session, "cell_type", choices = c(dataset_choices, cluster_choices))
    } else {
      showNotification("Monocle object is not properly loaded.", type = "error")
    }
  }
  
  # Mettre à jour les types de cellules en fonction du choix (dataset ou cluster)
  observeEvent(input$cell_type_choice, {
    req(monocle_object())
    if (input$cell_type_choice == "Dataset") {
      cell_type_choices <- unique(colData(monocle_object())$dataset)
    } else {
      cell_type_choices <- unique(colData(monocle_object())$cluster)
    }
    updatePickerInput(session, "cell_type", choices = cell_type_choices)
  })
  
  
  # Observer le bouton pour tracer les gènes le long du chemin
  observeEvent(input$plot_gene_path, {
    req(input$gene_list, input$cell_type)
    tryCatch({
      plot <- plot_genes_dynamics_along_path(monocle_object(), input$gene_list, input$cell_type)
      gene_path_plot(plot)
      showNotification("Gene dynamics plot completed!")
    }, error = function(e) {
      showNotification(paste("Error during gene dynamics plotting:", e$message), type = "error")
    })
  })
  
  # Rendu du plot
  output$genePathPlot <- renderPlot({
    req(gene_path_plot())
    gene_path_plot()
  })
  
  # Télécharger le plot
  output$download_gene_path_plot <- downloadHandler(
    filename = function() {
      paste("gene_path_", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      req(gene_path_plot())
      ggsave(file, plot = gene_path_plot(), dpi = input$dpi_selection_path, width = 10, height = 6)
    }
  )
  
  # Fonction pour tracer les gènes le long du chemin
  plot_genes_dynamics_along_path <- function(cds, genes, cell_type) {
    require(monocle3)
    
    # Filtrer le cds pour les gènes d'intérêt et le type de cellule (dataset ou cluster)
    lineage_cds <- cds[rowData(cds)$gene_short_name %in% genes & 
                         (colData(cds)$dataset %in% cell_type | colData(cds)$cluster %in% cell_type), ]
    
    # Ordonner les cellules le long du pseudotime
    lineage_cds <- order_cells(lineage_cds)
    
    # Tracer les gènes le long du pseudotime
    plot <- plot_genes_in_pseudotime(lineage_cds,
                                     color_cells_by = "pseudotime",
                                     min_expr = 0.5)
    
    return(plot)
  }
  
  



  calculateAndFindModules <- function() {
    print("IDs de gènes différentiels (pr_deg_ids):")
    print(pr_deg_ids())
    
    if (is.null(pr_deg_ids()) || is.null(monocle_object())) {
      showNotification("Identifiants de gènes différentiels ou objet Monocle manquant.", type = "error")
      return(NULL)
    }
    
    cds <- monocle_object()
    deg_ids <- pr_deg_ids()
    
    # Vérifier l'existence de gene_loadings
    if (is.null(gene_loadings)) {
      showNotification("Les chargements de gènes ne sont pas disponibles.", type = "error")
      return(NULL)
    }
    
    print("Noms des lignes dans gene_loadings:")
    print(rownames(gene_loadings))
    
    # Filtrer les identifiants de gènes différentiels qui sont dans gene_loadings
    deg_ids <- deg_ids[deg_ids %in% rownames(gene_loadings)]
    
    print("IDs de gènes après filtrage avec gene_loadings:")
    print(deg_ids)
    
    if (length(deg_ids) == 0) {
      showNotification("Aucun identifiant de gène différentiel valide pour l'analyse des modules.", type = "error")
      return(NULL)
    }
    
    # Assurez-vous que les noms de gènes correspondent
    gene_names <- rownames(rowData(cds))
    loading_gene_names <- rownames(gene_loadings)
    common_genes <- intersect(gene_names, loading_gene_names)
    
    if (length(common_genes) != length(loading_gene_names)) {
      showNotification("Les noms de gènes dans 'cds' et 'gene_loadings' ne correspondent pas.", type = "error")
      return(NULL)
    }
    
    filtered_cds <- cds[deg_ids, , drop = FALSE]
    print("CDS filtré pour l'analyse des modules:")
    print(filtered_cds)
    
    # Ajuster n_neighbors en fonction de la taille de l'ensemble de données
    dataset_size <- ncol(filtered_cds)
    n_neighbors_value <- min(15, dataset_size - 1)
    
    print(paste("Nombre de voisins (n_neighbors):", n_neighbors_value))
    print(paste("Taille de l'ensemble de données:", dataset_size))
    
    # Vérifier si n_neighbors_value est valide
    if (n_neighbors_value <= 0) {
      showNotification("La taille de l'ensemble de données est trop petite pour le nombre de voisins spécifié.", type = "error")
      return(NULL)
    }
    
    # Conversion en matrix pour utilisation avec le package umap
    exprs_matrix <- as.matrix(reducedDims(filtered_cds)$PCA)
    if (ncol(exprs_matrix) < 2) {
      showNotification("Nombre insuffisant de dimensions dans les données pour exécuter UMAP.", type = "error")
      return(NULL)
    }
    
    # Appel à la fonction umap du package umap
    umap_results <- umap::umap(exprs_matrix, n_neighbors = n_neighbors_value, n_components = 2)
    umap_embedding <- umap_results$layout
    
    reducedDims(filtered_cds)$UMAP <- umap_embedding
    
    # Appel correct à find_gene_modules avec les données UMAP
    gene_modules_df <- find_gene_modules(
      cds = filtered_cds, 
      reduction_method = "UMAP",
      max_components = 2,
      umap.metric = "cosine",
      umap.min_dist = 0.1,
      umap.n_neighbors = n_neighbors_value, # Utilisation de n_neighbors_value dynamique
      umap.fast_sgd = FALSE,
      umap.nn_method = "annoy",
      k = 20,
      leiden_iter = 1,
      partition_qval = 0.05,
      weight = FALSE,
      resolution = 1e-2,
      random_seed = 0L,
      cores = 1,
      verbose = FALSE,
      preprocess_method = "PCA",
      nn_control = list(method = "annoy", metric = "euclidean", n_trees = 50, search_k = 2000, cores = 1, grain_size = 1)
    )
    
    return(gene_modules_df)
  }
  
  
  observeEvent(input$run_analysis_button, {
    tryCatch({
      showModal(modalDialog(
        title = "Please Wait",
        "Calculating gene modules...",
        easyClose = FALSE,
        footer = NULL
      ))
      
      modules_df <- calculateAndFindModules()
      if (!is.null(modules_df)) {
        gene_modules_result(modules_df)
        updatePickerInput(session, "module_picker", choices = unique(modules_df$module))
      }
      
      removeModal()
    }, error = function(e) {
      removeModal()
      showNotification(paste0("Error during module calculation: ", e$message), type = "error")
    })
  })
  
  observeEvent(input$plot_module_button, {
    req(input$module_picker)
    selected_module <- input$module_picker
    modules_df <- gene_modules_result()
    selected_genes <- modules_df[modules_df$module == selected_module, "gene"]
    
    output$module_plot <- renderPlot({
      plotGeneModule(monocle_object(), selected_genes)
    })
  })
  
  observeEvent(input$module_picker, {
    selected_module <- input$module_picker
    modules_df <- gene_modules_result()
    selected_genes <- modules_df[modules_df$module == selected_module, "gene"]
    print(paste("Gènes sélectionnés pour le module:", selected_genes))
  })
}