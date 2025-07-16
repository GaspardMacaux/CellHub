############################## Script Single Dataset Analysis ##############################

# single_dataset_server.R

single_dataset_server <- function(input, output, session) {

  ############################## Loading Data ##############################
  # Tab 1 : Loading Data
  
  shinyjs::useShinyjs() # button deactivation
  
  # Reactives objects for the single dataset part
  single_dataset_object <- reactiveVal(NULL)
  gene_list <- reactiveValues(features = NULL)
  
  # Observer pour le chargement de données
  observeEvent(input$file, {
    showNotification("Uploading and processing data...", type = "message")
    
    tryCatch({
      # Get file extension
      file_extension <- tolower(tools::file_ext(input$file$name))
      message("File extension detected: ", file_extension)
      
      # Validate file type
      if (!file_extension %in% c("zip", "h5", "hdf5")) {
        stop("Please upload a .zip or .h5/.hdf5 file.")
      }
      
      # Clean workspace before loading
      cleanWorkspace("single")
      
      # Load raw data using our general function
      seuratObj <- loadRawData(
        file_path = input$file$datapath,
        dataset_type = input$dataset_type,
        species = input$species_choice
      )
      
      # Get dataset name
      dataset_name <- getDatasetFileName(
        list(input$file),
        "SingleDataset"
      )
      
      # Update project name
      seuratObj@project.name <- dataset_name
      
      # Store Seurat object
      single_dataset_object(seuratObj)
      
      # Success notification
      showNotification(
        paste0("Data loaded successfully! Cells: ", ncol(seuratObj), 
               ", Genes: ", nrow(seuratObj)), 
        type = "message"
      )
      
      # Update UI elements
      updateUIElements()
      
    }, error = function(e) {
      showNotification(paste("Error processing files:", e$message), type = "error")
      message("Error during file processing: ", e$message)
    })
  })
  
  # Load and restore colors after loading Seurat object for single dataset
  observeEvent(input$load_seurat_file, {
    tryCatch({
      loaded_seurat <- loadSeuratObject(
        rds_path = input$load_seurat_file$datapath,
        clean_before = TRUE,
        module_type = "single"
      )
      
      single_dataset_object(loaded_seurat)
      
      # Generate plot if UMAP exists
      plot_result <- generateInitialPlot(
        seurat_obj = loaded_seurat,
        remove_axes = input$remove_axes %||% FALSE,
        remove_legend = input$remove_legend %||% FALSE
      )
      
      if (!is.null(plot_result)) {
        single_dataset_object(plot_result$seurat_obj)  # Update with colors
        clustering_plot(plot_result$plot)
      }
      
      showNotification("Seurat object loaded successfully!")
      updateUIElements()
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  # Fonction réactive pour récupérer les clusters depuis l'objet Seurat
  get_clusters <- reactive({
    req(single_dataset_object())
    return(getClusters(single_dataset_object()))
  })



 

  # Function to update UI elements after data loading or dataset changes
  updateUIElements <- function() {
    req(single_dataset_object())  # Ensure the object is not null

    # Update cluster-related inputs
    cluster_choices <- unique(Idents(single_dataset_object()))
    updateSelectInput(session, "ident_1", choices = cluster_choices)
    updateCheckboxGroupInput(session, "ident_2", choices = cluster_choices, selected = "")
    updateSelectInput(session, "cluster_order_vln", choices = cluster_choices)
    updateSelectInput(session, "cluster_order_dot", choices = cluster_choices)


    # Update gene-related inputs
    gene_list <- sort(rownames(LayerData(single_dataset_object(), assay = "RNA", layer = 'counts')))
    updatePickerInput(session, "gene_select", choices = gene_list)
    updatePickerInput(session, "gene_select_heatmap", choices = gene_list)
    updatePickerInput(session, "gene_select_genes_analysis", choices = gene_list)
  }

  # Observer for Seurat object changes
  observeEvent(single_dataset_object(), {
    updateUIElements()
  })




  ############################## QC metrics and normalization ##############################
  # Count nuclei before and after QC
  nuclei_count <- reactive({
    req(single_dataset_object())
    seuratRNA_subset <- subset(single_dataset_object(),
                               subset = nFeature_RNA > input$nFeature_range[1] &
                                 nFeature_RNA < input$nFeature_range[2] &
                                 percent.mt < input$percent.mt_max
    )
    if (exists("seuratRNA_subset")) {
      # Log QC parameters for verification
      message(sprintf("QC Parameters used:"))
      message(sprintf("- nFeature range: %d to %d", input$nFeature_range[1], input$nFeature_range[2]))
      message(sprintf("- Max MT%%: %f", input$percent.mt_max))
      message(sprintf("Original nuclei: %d", dim(single_dataset_object())[2]))
      message(sprintf("Filtered nuclei: %d", dim(seuratRNA_subset)[2]))

      # Verify filters were actually applied
      if(any(seuratRNA_subset$nFeature_RNA < input$nFeature_range[1]) ||
         any(seuratRNA_subset$nFeature_RNA > input$nFeature_range[2]) ||
         any(seuratRNA_subset$percent.mt > input$percent.mt_max)) {
        warning("Some nuclei outside QC parameters remain!")
      }

      return(dim(seuratRNA_subset)[2])
    } else {
      warning("Subsetting failed!")
      return(0)
    }
  })

  # Text output number of nuclei with verification
  output$nuclei_count <- renderInfoBox({
    tryCatch({
      count_value <- nuclei_count()
      total_nuclei <- dim(single_dataset_object())[2]

      # Additional verification
      if(count_value > total_nuclei) {
        warning("Filtered count larger than total count!")
      }

      message("Number of Nuclei in infoBox: ", count_value)

      infoBox(
        title = "Number of Nuclei",
        value = paste0(
          count_value,
          " (",
          round((count_value/total_nuclei)*100, 2),
          "% retained)"
        ),
        color = "blue",
        icon = icon("dna"),
        width = 1
      )
    }, error = function(e) {
      warning("Error in nuclei count: ", e$message)
      infoBox(
        title = "Number of Nuclei",
        value = "Error",
        icon = icon("dna"),
        color = "red"
      )
    })
  })

  # Observer pour la normalisation
  observeEvent(input$normalize_data, {
    req(single_dataset_object())
    
    showModal(modalDialog(
      title = "Traitement",
      "Normalisation des données et identification des caractéristiques variables...",
      footer = NULL,
      easyClose = FALSE
    ))
    
    tryCatch({
      # Récupérer l'objet Seurat
      seurat_obj <- single_dataset_object()
      
      # Définir l'assay actif sur RNA et normaliser
      message("Définition de l'assay par défaut sur RNA")
      DefaultAssay(seurat_obj) <- "RNA"
      
      seurat_obj <- normalizeData(seurat_obj)
      
      # Mettre à jour l'objet Seurat
      single_dataset_object(seurat_obj)
      
      # Afficher les caractéristiques variables
      output$variable_feature_plot <- renderPlot({
        plot1 <- VariableFeaturePlot(seurat_obj)
        LabelPoints(plot = plot1, points = head(VariableFeatures(seurat_obj), 10), repel = TRUE)
      })
      
      removeModal()
      showNotification("Données normalisées avec succès", type = "message")
      
    }, error = function(e) {
      removeModal()
      showNotification(paste("Erreur pendant la normalisation:", e$message), type = "error")
      message(paste("Erreur détaillée:", e$message))
    })
  })
  

  
  # Display QC metrics on a VlnPlot chart with verification
  output$vlnplot <- renderPlot({
    req(input$QCmetrics, single_dataset_object())
    tryCatch({
      if (input$QCmetrics) {
        # Verify metrics exist
        required_metrics <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
        missing_metrics <- setdiff(required_metrics, colnames(single_dataset_object()@meta.data))
        if(length(missing_metrics) > 0) {
          stop(paste("Missing metrics:", paste(missing_metrics, collapse = ", ")))
        }

        VlnPlot(single_dataset_object(),
                features = required_metrics,
                ncol = 3)
      }
    }, error = function(e) {
      showNotification(paste0("Error displaying violin plot: ", e$message), type = "error")
    })
  })

  # Display scatter plots with verification
  output$scatter_plot1 <- renderPlot({
    req(input$show_plots, single_dataset_object())
    tryCatch({
      if (input$show_plots) {
        # Create subset with verification
        seuratRNA_subset <- subset(single_dataset_object(),
                                   subset = nFeature_RNA > input$nFeature_range[1] &
                                     nFeature_RNA < input$nFeature_range[2] &
                                     percent.mt < input$percent.mt_max
        )

        # Verify subsetting worked
        if(ncol(seuratRNA_subset) == 0) {
          stop("No cells pass the current filters!")
        }

        FeatureScatter(seuratRNA_subset,
                       feature1 = "nCount_RNA",
                       feature2 = "percent.mt")
      }
    }, error = function(e) {
      showNotification(paste0("Error in scatter plot 1: ", e$message), type = "error")
    })
  })

  # Display scatter plot 2 with same verifications
  output$scatter_plot2 <- renderPlot({
    req(input$show_plots, single_dataset_object())
    tryCatch({
      if (input$show_plots) {
        seuratRNA_subset <- subset(single_dataset_object(),
                                   subset = nFeature_RNA > input$nFeature_range[1] &
                                     nFeature_RNA < input$nFeature_range[2] &
                                     percent.mt < input$percent.mt_max
        )

        if(ncol(seuratRNA_subset) == 0) {
          stop("No cells pass the current filters!")
        }

        FeatureScatter(seuratRNA_subset,
                       feature1 = "nCount_RNA",
                       feature2 = "nFeature_RNA")
      }
    }, error = function(e) {
      showNotification(paste0("Error in scatter plot 2: ", e$message), type = "error")
    })
  })

  # Apply QC filters with thorough verification
  observeEvent(input$apply_qc, {
    req(single_dataset_object())

    tryCatch({
      showModal(modalDialog(
        title = "Please Wait",
        "Applying QC filters...",
        easyClose = FALSE,
        footer = NULL
      ))

      seurat_obj <- single_dataset_object()

      # Record original metrics
      orig_total <- ncol(seurat_obj)
      orig_mt_mean <- mean(seurat_obj$percent.mt)
      orig_feature_mean <- mean(seurat_obj$nFeature_RNA)

      # Apply filters with verification
      seurat_obj <- subset(seurat_obj,
                           subset = nFeature_RNA > input$nFeature_range[1] &
                             nFeature_RNA < input$nFeature_range[2] &
                             percent.mt < input$percent.mt_max
      )

      # Verify filtering
      if(any(seurat_obj$nFeature_RNA <= input$nFeature_range[1]) ||
         any(seurat_obj$nFeature_RNA >= input$nFeature_range[2]) ||
         any(seurat_obj$percent.mt >= input$percent.mt_max)) {
        stop("QC filtering failed - some nuclei outside parameters remain")
      }

      # Update object if verification passed
      single_dataset_object(seurat_obj)

      # Generate detailed report
      num_nuclei_after_qc <- ncol(seurat_obj)
      new_mt_mean <- mean(seurat_obj$percent.mt)
      new_feature_mean <- mean(seurat_obj$nFeature_RNA)

      msg <- paste0(
        "QC filters applied successfully:\n",
        sprintf("- Nuclei: %d → %d (%d%% retained)\n",
                orig_total, num_nuclei_after_qc,
                round(num_nuclei_after_qc/orig_total*100)),
        sprintf("- Mean MT%%: %.2f%% → %.2f%%\n",
                orig_mt_mean, new_mt_mean),
        sprintf("- Mean Features: %.1f → %.1f",
                orig_feature_mean, new_feature_mean)
      )

      showNotification(msg, type = "message", duration = 10)
      removeModal()

      # Update UI
      updateUIElements()

    }, error = function(e) {
      removeModal()
      showNotification(paste0("Error applying QC: ", e$message), type = "error")
    })
  })


  observeEvent(input$normalize_data, {
    req(input$num_var_features)  # S'assurer que l'utilisateur a bien fourni une valeur
    tryCatch({
      # Afficher la boîte de dialogue modale
      showModal(modalDialog(
        title = "Please Wait",
        "Normalizing data...",
        easyClose = FALSE,
        footer = NULL
      ))

      # Normalisation des données
      normalized_seurat <- NormalizeData(single_dataset_object(), normalization.method = "LogNormalize", scale.factor = input$scale_factor)
      print("Data normalization completed")

      # Identification des variable features basée sur l'entrée de l'utilisateur
      normalized_seurat <- FindVariableFeatures(normalized_seurat, selection.method = "vst", nfeatures = input$num_var_features)
      print("Variable features identification completed")
      single_dataset_object(normalized_seurat)

      # Fermer la boîte de dialogue modale
      removeModal()
    }, error = function(e) {
      removeModal() # Fermer la boîte de dialogue modale en cas d'erreur
      showNotification(paste0("Error during data normalization or variable features identification: ", e$message), type = "error")
    })
  })


  # Plot rendering of variable features
  output$variable_feature_plot <- renderPlot({
    req(input$show_plots, single_dataset_object())
    if (!is.null(single_dataset_object()) && length(VariableFeatures(single_dataset_object())) > 0) {
      plot1 <- VariableFeaturePlot(single_dataset_object())
      top10 <- head(VariableFeatures(single_dataset_object()), 10)
      plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
      plot2
    } else {
      return(NULL)
    }
  })

  ############################## Scaling, PCA and elbow plot ##############################
  #  Tab 3: Scaling and PCA and elbow plot


  # Scaling data
  observeEvent(input$scale_button, {
    showNotification("Scaling data...", type = "message")
    tryCatch({
      showModal(modalDialog(
        title = "Please Wait",
        "Running PCA on the dataset...",
        easyClose = FALSE,
        footer = NULL
      ))
      print("Scaling data for all genes")
      scaled_seurat <-ScaleData(single_dataset_object(),assay= "RNA", features = rownames(single_dataset_object()))
      single_dataset_object(scaled_seurat)
      print("Data scaling completed")
      shinyjs::enable("pca_button")
      showNotification("Data scaling completed successfully.", type = "message")
    }, error = function(e) {
      print(paste("Error during data scaling:", e$message))
      showNotification(paste("Error during data scaling:", e$message), type = "error")
    })

    tryCatch({
      pca_result <- RunPCA(single_dataset_object(), features = VariableFeatures(single_dataset_object()))
      single_dataset_object(pca_result)
      removeModal()
    }, error = function(e) {
      showNotification(paste0("Error running PCA:", e$message), type = "error")
      removeModal()
    })

  })


  #Printing PCA results
  output$pca_results <- renderPrint({
    req(!is.null(single_dataset_object()))
    if ("pca" %in% names(single_dataset_object())) {
      print(single_dataset_object()$pca, dims = 1:5, nfeatures = 5)
    }
  })

  #Vizdimloading plot
  output$loading_plot <- renderPlot({
    req(!is.null(single_dataset_object()))
    tryCatch(VizDimLoadings(single_dataset_object(), dims = 1:2, reduction = "pca"), error = function(e){})
  })

  # Dimplot
  output$dim_plot <- renderPlot({
    req(!is.null(single_dataset_object()))
    tryCatch(DimPlot(single_dataset_object(), reduction = "pca"), error = function(e){})
  })

  # Show the ElbowPlot
  observeEvent(input$run_elbow, {
    tryCatch({
      output$elbow_plot <- renderPlot({
        req(single_dataset_object())
        ElbowPlot(single_dataset_object(), ndims = 50)
      })}, error = function(e) {
        showNotification(paste("Error running elbow plot:", e$message), type = "error")
      })
  })

  ############################## Neighbors calculation and clustering ##############################

  clustering_plot <- reactiveVal()

  # Finding neighbors
  observeEvent(input$run_neighbors, {
    tryCatch({
      req(!is.null(single_dataset_object()))

      # Afficher la boîte de dialogue modale
      showModal(modalDialog(
        title = "Please Wait",
        "Looking for neighbors and calculating UMAP...",
        easyClose = FALSE,
        footer = NULL
      ))

      seurat_tmp <- single_dataset_object()
      seurat_tmp <- FindNeighbors(seurat_tmp, dims = 1:input$dimension_1)
      seurat_tmp <- RunUMAP(seurat_tmp, dims = 1:input$dimension_1)
      single_dataset_object(seurat_tmp)
      print("Neighbors found and UMAP calculated.")

      plot <- DimPlot(single_dataset_object(), group.by = "ident") + ggtitle(NULL)
      if(input$remove_axes) { plot <- plot + NoAxes() }
      if(input$remove_legend) { plot <- plot + NoLegend() }
      clustering_plot(plot)

      # Fermer la boîte de dialogue modale
      removeModal()
    }, error = function(e) {
      removeModal() # Fermer la boîte de dialogue modale en cas d'erreur
      showNotification(paste0("Neighbor search error: ", e$message), type = "error")
    })
  })

# Clustering
  observeEvent(input$run_clustering, {
    tryCatch({
      req(!is.null(single_dataset_object()))

      # Afficher la boîte de dialogue modale
      showModal(modalDialog(
        title = "Please Wait",
        "Clustering process started...",
        easyClose = FALSE,
        footer = NULL
      ))

      seurat_tmp <- single_dataset_object()
      seurat_tmp <- FindClusters(seurat_tmp, resolution = input$resolution_step1, algorithm = as.integer(input$algorithm_select))
      single_dataset_object(seurat_tmp)

      plot <- DimPlot(single_dataset_object(), group.by = "ident", label = FALSE) + ggtitle("")
      if(input$remove_axes) { plot <- plot + NoAxes() }
      if(input$remove_legend) { plot <- plot + NoLegend() }
      clustering_plot(plot)

      print("Clustering performed.")

      # Fermer la boîte de dialogue modale
      removeModal()
    }, error = function(e) {
      removeModal() # Fermer la boîte de dialogue modale en cas d'erreur
      showNotification(paste0("Clustering error: ", e$message), type = "error")
    })
  })



  output$clustering_plot <- renderPlot({
    req(clustering_plot())
    clustering_plot()
  })

  #SAVING UMAP PLOT
  output$downloadUMAP <- createDownloadHandler(
    reactive_data = clustering_plot,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(single_dataset_object(), default_name = "UMAP_plot")
    }),
    data_name = "UMAP",
    download_type = "plot",
    plot_params = list(
      file_type = reactive({ input$umap_export_format }),
      width = reactive({ ifelse(input$umap_export_format == "pdf", 11, 10) }),
      height = reactive({ ifelse(input$umap_export_format == "pdf", 8, 6) }),
      dpi = reactive({ input$dpi_umap })
    )
  )
  
######################Doublet Finder#########################
  
  observeEvent(input$run_doubletfinder, {
    req(single_dataset_object())

    showModal(modalDialog(
      title = "Running DoubletFinder",
      p("Processing your data..."),
      footer = NULL,
      easyClose = FALSE
    ))

    tryCatch({
      seu <- single_dataset_object()

      # Parameter sweep
      sweep.res.list <- paramSweep(seu, PCs = 1:input$pc_use, sct = FALSE)
      sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
      bcmvn <- find.pK(sweep.stats)

      # Calculate expected doublets
      nExp_poi <- round(input$doublet_rate/100 * ncol(seu))

      # Find optimal pK
      pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))

      # Run DoubletFinder
      seu <- doubletFinder(seu,
                           PCs = 1:input$pc_use,
                           pN = input$pN_value,
                           pK = pK,
                           nExp = nExp_poi,
                           reuse.pANN = FALSE,
                           sct = FALSE)

      # Get results column
      df_col <- grep("DF.classifications", colnames(seu@meta.data), value = TRUE)[1]

      # Create plots
      output$doublet_umap <- renderPlot({
        DimPlot(seu,
                reduction = "umap",
                group.by = df_col,
                cols = c("Singlet" = "grey", "Doublet" = "red")) +
          ggtitle("Doublets Identification")
      })

      output$doublet_stats <- renderDT({
        stats <- table(seu@meta.data[[df_col]])
        data.frame(
          Category = names(stats),
          Count = as.numeric(stats),
          Percentage = round(as.numeric(stats)/sum(stats)*100, 2)
        )
      })

      # Store results
      single_dataset_object(seu)

      removeModal()
      showNotification("DoubletFinder analysis completed!", type = "message")

    }, error = function(e) {
      removeModal()
      showNotification(paste("Error in DoubletFinder:", e$message), type = "error")
    })
  })


  # Add in server.R
  observeEvent(input$remove_doublets, {
    req(single_dataset_object())

    tryCatch({
      seu <- single_dataset_object()

      # Get the doublet classification column
      df_col <- grep("DF.classifications", colnames(seu@meta.data), value = TRUE)[1]

      if(is.null(df_col)) {
        showNotification("Please run doublet detection first", type = "warning")
        return()
      }

      # Keep only singlets
      seu_filtered <- subset(seu, cells = rownames(seu@meta.data)[seu@meta.data[[df_col]] == "Singlet"])

      # Update object
      single_dataset_object(seu_filtered)

      # Update plots
      output$doublet_umap <- renderPlot({
        DimPlot(seu_filtered, reduction = "umap", group.by = "seurat_clusters") +
          ggtitle("UMAP (Doublets Removed)")
      })

      # Show success message
      showNotification(
        sprintf("Removed %d doublets from dataset",
                ncol(seu) - ncol(seu_filtered)),
        type = "message"
      )

    }, error = function(e) {
      showNotification(paste("Error removing doublets:", e$message), type = "error")
    })
  })

  ############################## Visualization of expressed genes ##############################


  # Observer pour mettre à jour les choix d'assays
  observe({
    req(single_dataset_object())
    updateAssayChoices(session, single_dataset_object())
  })
  
  # Observer pour mettre à jour les choix de gènes
  observeEvent(c(single_dataset_object(), input$viz_assay), {
    req(single_dataset_object(), input$viz_assay)
    updateGeneChoices(session, single_dataset_object(), input$viz_assay)
  })
  
  # Observer pour mettre à jour les sélecteurs de cluster
  observe({
    req(single_dataset_object())
    updateClusterChoices(session, single_dataset_object())
  })
  
  
  # Observer to update text inputs when genes are selected in pickerInput
  observeEvent(input$gene_select, {
    selected_genes <- input$gene_select
    
    if (!is.null(selected_genes) && length(selected_genes) > 0) {
      # Convert selected genes to comma-separated string
      genes_text <- paste(selected_genes, collapse = ", ")
      
      # Update all gene text inputs with selected genes
      updateTextInput(session, "gene_list_feature", value = genes_text)
      updateTextInput(session, "gene_list_dotplot", value = genes_text)
      updateTextInput(session, "gene_list_ridge_plot", value = genes_text)
      updateTextInput(session, "gene_list_vln", value = genes_text)  # If you have this one too
    }
  })
  
  # FeaturePlot avec options
  feature_plot <- reactiveVal()

  observeEvent(input$show_feature, {
    tryCatch({
      req(input$gene_list_feature, single_dataset_object())
      genes <- unique(trimws(unlist(strsplit(input$gene_list_feature, ","))))
      seurat_object <- single_dataset_object()
      req(seurat_object)
      DefaultAssay(seurat_object) <- input$viz_assay
      present_genes <- genes[genes %in% rownames(LayerData(seurat_object, assay = "RNA", layer = "counts"))]

      if (length(present_genes) > 0) {
        print("Creating plot")
        min_cut <- if (!is.na(input$min_cutoff)) input$min_cutoff else NA
        max_cut <- if (!is.na(input$max_cutoff)) input$max_cutoff else NA

        # Configuration de base du thème
        base_theme <- theme(
          axis.text = element_text(size = input$axis_text_size),
          axis.title = element_text(size = input$title_text_size),
          plot.title = element_text(size = input$title_text_size),
          legend.text = element_text(size = input$axis_text_size),
          axis.line = element_line(linewidth = input$axis_line_width),
          axis.ticks = element_line(linewidth = input$axis_line_width)
        )

        # Créer le plot avec les paramètres appropriés
        plot <- FeaturePlot(
          seurat_object,
          features = present_genes,
          blend = input$show_coexpression && length(present_genes) > 1,
          blend.threshold = 1,
          order = TRUE,
          min.cutoff = min_cut,
          max.cutoff = max_cut
        ) + base_theme

        # Ajout des modifications conditionnelles
        if (input$add_noaxes_feature) {
          print("Adding NoAxes")
          plot <- plot + NoAxes()
        }
        if (input$add_nolegend_feature) {
          print("Adding NoLegend")
          plot <- plot + NoLegend()
        }

        feature_plot(plot)
      }
    }, error = function(e) {
      showNotification(paste("Error in FeaturePlot: ", e$message), type = "error")
      print(paste("Error details:", e$message))
    })
  })

  output$feature_plot <- renderPlot({
    req(feature_plot())
    feature_plot()
  })

  # VlnPlot
  vln_plot <- reactiveVal()
  observeEvent(input$show_vln, {
    tryCatch({
      req(input$gene_list_vln)
      print("Starting VlnPlot generation")

      genes <- unique(trimws(unlist(strsplit(input$gene_list_vln, ","))))
      print(paste("Processing genes:", paste(genes, collapse = ", ")))

      seurat_object <- single_dataset_object()
      req(seurat_object)
      DefaultAssay(seurat_object) <- input$viz_assay

      if (!is.null(input$cluster_order_vln) && length(input$cluster_order_vln) > 0) {
        print(paste("Selected clusters:", paste(input$cluster_order_vln, collapse=", ")))

        seurat_subset <- seurat_object
        Idents(seurat_subset) <- factor(
          Idents(seurat_subset),
          levels = input$cluster_order_vln
        )
        seurat_subset <- subset(seurat_subset, idents = input$cluster_order_vln)

        plot <- VlnPlot(
          object = seurat_subset,
          features = genes,
          pt.size = ifelse(input$hide_vln_points, 0, 1)
        ) +
          theme(
            axis.text = element_text(size = input$axis_text_size),
            axis.title = element_text(size = input$title_text_size),
            plot.title = element_text(size = input$title_text_size),
            legend.text = element_text(size = input$axis_text_size),
            axis.line = element_line(linewidth = input$axis_line_width),
            axis.ticks = element_line(linewidth = input$axis_line_width)
          )

      } else {
        plot <- VlnPlot(
          object = seurat_subset,
          features = genes,
          pt.size = ifelse(input$hide_vln_points, 0, 1)
        ) +
          theme(
            axis.text = element_text(size = input$axis_text_size),
            axis.title = element_text(size = input$title_text_size),
            plot.title = element_text(size = input$title_text_size),
            legend.text = element_text(size = input$axis_text_size),
            axis.line = element_line(linewidth = input$axis_line_width),
            axis.ticks = element_line(linewidth = input$axis_line_width)
          )
      }

      if (!is.null(input$add_noaxes_vln) && input$add_noaxes_vln) {
        plot <- plot + NoAxes()
      }
      if (!is.null(input$add_nolegend_vln) && input$add_nolegend_vln) {
        plot <- plot + NoLegend()
      }

      vln_plot(plot)
      print("Plot successfully stored")

    }, error = function(e) {
      print(paste("Error in VlnPlot:", e$message))
      print("Error details:")
      print(e)
    })
  })

  output$vln_plot <- renderPlot({
    req(vln_plot())
    vln_plot()
  })


  # DotPlot reactive value
  dot_plot <- reactiveVal()
  # DotPlot avec gestion correcte des clusters
  observeEvent(input$show_dot, {
    tryCatch({
      req(input$gene_list_dotplot, single_dataset_object())
      genes <- unique(trimws(unlist(strsplit(input$gene_list_dotplot, ","))))
      seurat_object <- single_dataset_object()
      req(seurat_object)
      DefaultAssay(seurat_object) <- input$viz_assay
   
      if (!is.null(input$cluster_order_dot) && length(input$cluster_order_dot) > 0) {
        print(paste("Selected clusters:", paste(input$cluster_order_dot, collapse=", ")))

        seurat_subset <- seurat_object
        Idents(seurat_subset) <- factor(
          Idents(seurat_subset),
          levels = input$cluster_order_dot
        )
        seurat_subset <- subset(seurat_subset, idents = input$cluster_order_dot)

        plot <- DotPlot(
          seurat_subset,
          features = genes
        ) +
          RotatedAxis() +
          theme(
            axis.text = element_text(size = input$axis_text_size),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title = element_text(size = input$title_text_size),
            plot.title = element_text(size = input$title_text_size),
            legend.text = element_text(size = input$axis_text_size),
            axis.line = element_line(linewidth = input$axis_line_width),
            axis.ticks = element_line(linewidth = input$axis_line_width)
          )
      } else {
        plot <- DotPlot(
          seurat_subset,
          features = genes
        ) +
          RotatedAxis() +
          theme(
            axis.text = element_text(size = input$axis_text_size),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title = element_text(size = input$title_text_size),
            plot.title = element_text(size = input$title_text_size),
            legend.text = element_text(size = input$axis_text_size),
            axis.line = element_line(linewidth = input$axis_line_width),
            axis.ticks = element_line(linewidth = input$axis_line_width)
          )
      }

      if (input$add_noaxes_dot) plot <- plot + NoAxes()
      if (input$add_nolegend_dot) plot <- plot + NoLegend()
      if (input$invert_axes) plot <- plot + coord_flip()

      dot_plot(plot)
    }, error = function(e) {
      showNotification(paste("Error in DotPlot: ", e$message), type = "error")
      print(paste("Error details:", e$message))
    })
  })

  # Observer pour mettre à jour les choix de clusters
  observe({
    req(single_dataset_object())

    # Récupérer les identifiants de clusters uniques
    cluster_choices <- unique(Idents(single_dataset_object()))

    # Mettre à jour les deux selectInput
    updateSelectInput(session,
                      "cluster_order_vln",  # Pour le VlnPlot
                      choices = cluster_choices,
                      selected = cluster_choices)

    updateSelectInput(session,
                      "cluster_order_dot",  # Pour le DotPlot
                      choices = cluster_choices,
                      selected = cluster_choices)
  })
  # Display the DotPlot
  output$dot_plot <- renderPlot({
    req(dot_plot())
    dot_plot()
  })

  # RidgePlot
  ridge_plot <- reactiveVal()

  observeEvent(input$show_ridge, {
    tryCatch({
      req(input$gene_list_ridge_plot, single_dataset_object())
      genes <- unique(trimws(unlist(strsplit(input$gene_list_ridge_plot, ","))))
      seurat_object <- single_dataset_object()
      req(seurat_object)
      DefaultAssay(seurat_object) <- input$viz_assay

      plot <- RidgePlot(
        seurat_object,
        features = genes
      ) +
        theme(
          axis.text = element_text(size = input$axis_text_size),
          axis.title = element_text(size = input$title_text_size),
          plot.title = element_text(size = input$title_text_size),
          legend.text = element_text(size = input$axis_text_size),
          axis.line = element_line(linewidth = input$axis_line_width),
          axis.ticks = element_line(linewidth = input$axis_line_width)
        )

      # Options de base uniquement
      if (input$add_noaxes_ridge) plot <- plot + NoAxes()
      if (input$add_nolegend_ridge) plot <- plot + NoLegend()

      ridge_plot(plot)
    }, error = function(e) {
      showNotification(paste("Error in RidgePlot: ", e$message), type = "error")
    })
  })

  output$ridge_plot <- renderPlot({
    req(ridge_plot())
    ridge_plot()
  })

  # Download handler for RidgePlot
  output$download_ridge_plot <- createDownloadHandler(
    reactive_data = ridge_plot,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(single_dataset_object(), default_name = "RidgePlot")
    }),
    data_name = "RidgePlot",
    download_type = "plot",
    plot_params = list(
      file_type = input$plot_format,
      width = 10,
      height = 8,
      dpi = input$dpi_plot
    )
  )
  
  # Download handler for DotPlot
  output$downloadDotPlot <- createDownloadHandler(
    reactive_data = dot_plot,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(single_dataset_object(), default_name = "DotPlot")
    }),
    data_name = "DotPlot",
    download_type = "plot",
    plot_params = list(
      file_type = input$plot_format,
      width = 10,
      height = 8,
      dpi = input$dpi_plot
    )
  )
  
  # Download handler for ViolinPlot
  output$downloadVlnPlot <- createDownloadHandler(
    reactive_data = vln_plot,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(single_dataset_object(), default_name = "VlnPlot")
    }),
    data_name = "VlnPlot",
    download_type = "plot",
    plot_params = list(
      file_type = input$plot_format,
      width = 10,
      height = 8,
      dpi = input$dpi_plot
    )
  )
  
  # Download handler for FeaturePlot
  output$downloadFeaturePlot <- createDownloadHandler(
    reactive_data = feature_plot,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(single_dataset_object(), default_name = "FeaturePlot")
    }),
    data_name = "FeaturePlot",
    download_type = "plot",
    plot_params = list(
      file_type = input$plot_format,
      width = 10,
      height = 8,
      dpi = input$dpi_plot
    )
  )
  
  # Analyse de l'expression des gènes
  number_of_nuclei <- reactiveVal(NULL)

  observeEvent(input$analyze_btn, {
    req(input$gene_select_genes_analysis, single_dataset_object())
    
    tryCatch({
      result <- analyze_gene_expression(
        seurat_obj = single_dataset_object(),
        selected_genes = input$gene_select_genes_analysis,
        assay_name = input$viz_assay,
        expression_threshold = input$expression_threshold %||% 0.1,
        is_integrated = FALSE
      )
      
      number_of_nuclei(result$data)
      
      output$expression_summary <- renderDT({
        render_expression_table(result, "expression_summary")
      })
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  
  # Rendre le plot actuel dans l'UI
  output$selected_plot_display <- renderPlot({
    req(current_plot())
    current_plot()
  })

  # Réactive pour stocker et mettre à jour le plot actuel basé sur la sélection de l'utilisateur
  current_plot <- reactiveVal()
  
  
  
  
  observe({
    # Mise à jour du plot actuel en fonction de la sélection
    plot_type <- input$plot_type_select
    if (plot_type == "FeaturePlot") {
      current_plot(feature_plot())
    } else if (plot_type == "VlnPlot") {
      current_plot(vln_plot())
    } else if (plot_type == "DotPlot") {
      current_plot(dot_plot())
    } else if (plot_type == "RidgePlot") {
      current_plot(ridge_plot())
    }
  })
  
  output$download_genes_number_expresion <- createDownloadHandler(
    reactive_data = number_of_nuclei,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(single_dataset_object(), default_name = "GeneExpression")
    }),
    data_name = "number_of_gene_expression",
    download_type = "csv"
  )
  
  # Download handler for Seurat object with modal
  output$save_seurat_object_2 <- createDownloadHandler(
    reactive_data = single_dataset_object,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(single_dataset_object(), default_name = "SingleDataset")
    }),
    data_name = "seurat_object",
    download_type = "seurat",
    show_modal = TRUE  
  )

  ############################## Heatmaps and scatter plots ##############################

  heatmap_plot <- reactiveVal()


  observe({
    req(single_dataset_object())
    clusters <- levels(Idents(single_dataset_object()))
    updateSelectInput(session, "cluster_selector", choices = clusters, selected = clusters)
  })


  observeEvent(input$select_all_clusters, {
    req(single_dataset_object())
    if(input$select_all_clusters) {
      shinyjs::disable("text_clusters")
    } else {
      shinyjs::enable("text_clusters")
    }
  })

  observeEvent(input$generateCustomHeatmap, {
    showModal(modalDialog(
      title = "Processing",
      div(
        h4("Generating Heatmap...", style = "text-align: center;"),
        p("Please wait while we process your data.", style = "text-align: center;"),
        p("This may take a few moments depending on the number of genes and clusters selected.",
          style = "text-align: center; color: #666;")
      ),
      footer = NULL,
      easyClose = FALSE
    ))
    req(single_dataset_object())
    seurat_object <- single_dataset_object()

    tryCatch({
      # Gestion des clusters (code existant inchangé)
      if (!input$select_all_clusters) {
        if (nchar(input$text_clusters) > 0) {
          specified_clusters <- unlist(strsplit(trimws(input$text_clusters), ",\\s*"))
          valid_clusters <- specified_clusters %in% levels(Idents(seurat_object))
          if (sum(valid_clusters) == 0) {
            showNotification("None of the specified clusters is valid.", type = "error")
            return()
          }
          Idents(seurat_object) <- factor(Idents(seurat_object), levels = specified_clusters)
          seurat_object <- subset(seurat_object, idents = specified_clusters)
        }
      }

      # Sélection des gènes - nouvelle logique
      if(input$use_top10_genes) {
        # Utiliser les top 10 gènes par cluster
        markers <- FindAllMarkers(seurat_object, min.pct = 0.25)
        selected_genes <- markers %>%
          group_by(cluster) %>%
          dplyr::filter(avg_log2FC > 1) %>%
          slice_head(n = 10) %>%
          ungroup() %>%
          pull(gene)
      } else {
        # Utiliser les gènes sélectionnés manuellement
        selected_genes <- input$gene_select_heatmap
      }

      valid_genes <- selected_genes %in% rownames(seurat_object[["RNA"]])
      if (sum(valid_genes) == 0) {
        showNotification("None of the selected genes is valid.", type = "error")
        return()
      }

      # Générer la heatmap
      heatmap_plot(
        DoHeatmap(seurat_object, features = selected_genes[valid_genes]) + NoLegend()
      )

    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
    removeModal()
  })

  # Rendre le plot stocké dans l'interface utilisateur
  output$heatmap_single <- renderPlot({
    req(heatmap_plot())
    heatmap_plot()
  })

  # Download handler for heatmap
  output$download_heatmap_single <- createDownloadHandler(
    reactive_data = heatmap_plot,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(single_dataset_object(), default_name = "Heatmap")
    }),
    data_name = "heatmap",
    download_type = "plot",
    plot_params = list(
      file_type = "tiff",
      width = 10,
      height = 8,
      dpi = input$dpi_heatmap_single
    )
  )


  # Mettez à jour les choix de selectInput pour les gènes/features et les clusters
  observe({
    req(single_dataset_object())
    updateGeneChoices(session, single_dataset_object(), "RNA", c("feature1_select", "feature2_select"))
    updateClusterTextInputs(session, single_dataset_object(), c("scatter_text_clusters", "text_clusters"))
  })

  scatter_plot <- reactiveVal()

  observeEvent(input$generateFeatureScatter, {
    req(input$feature1_select, input$feature2_select)

    tryCatch({
      seurat_object <- single_dataset_object()
      if (!input$scatter_select_all_clusters) {
        if (nchar(input$scatter_text_clusters) > 0) {
          specified_clusters <- unlist(strsplit(trimws(input$scatter_text_clusters), ",\\s*"))
          seurat_object <- subset(seurat_object, idents = specified_clusters)
        }
      }

      plot <- FeatureScatter(
        object = seurat_object,
        feature1 = input$feature1_select,
        feature2 = input$feature2_select,
        group.by = "ident"
      )
      scatter_plot(plot)
    }, error = function(e) {
      # Gérer l'erreur
      showNotification(paste("Erreur lors de la génération du scatter plot :", e$message), type = "error")
    })
  })

  output$scatterPlot <- renderPlot({
    req(scatter_plot())
    scatter_plot()
  })

  # Download handler for scatter plot
  output$download_scatter_plot <- createDownloadHandler(
    reactive_data = scatter_plot,
    object_name_reactive = reactive({ 
      getObjectNameForDownload(single_dataset_object(), default_name = "FeatureScatter")
    }),
    data_name = "FeatureScatter",
    download_type = "plot",
    plot_params = list(
      file_type = "tiff",
      width = 10,
      height = 8,
      dpi = input$dpi_scatter
    )
  )



  ############################## Final UMAP ##############################
  # Tab 8: Final UMAP

  #Reactive variable for that tab
  cluster_colours <- reactiveVal()  # Initialize a list to store cluster colors

  observe({
    req(single_dataset_object())
    updateClusterChoices(session, single_dataset_object(), list(select = c("select_cluster", "cluster_select")))
  })


  # Function for renaming each cluster
  observeEvent(input$rename_single_cluster_button, {
    req(input$select_cluster, input$rename_single_cluster, single_dataset_object())
    
    updated_seurat <- single_dataset_object()
    
    # Si le nouveau nom de cluster existe déjà, fusionnez les clusters
    if (input$rename_single_cluster %in% unique(Idents(updated_seurat))) {
      cells_to_merge <- which(Idents(updated_seurat) %in% c(input$select_cluster, input$rename_single_cluster))
      Idents(updated_seurat, cells = cells_to_merge) <- input$rename_single_cluster
      showNotification("Clusters merged under the name: ", input$rename_single_cluster, type = "message")
    } else {
      # Sinon, renommez simplement le cluster sélectionné
      Idents(updated_seurat, cells = which(Idents(updated_seurat) == input$select_cluster)) <- input$rename_single_cluster
      showNotification("Cluster renamed to: ", input$rename_single_cluster, type = "message")
    }
    
    # AJOUT: Mettre à jour la colonne cluster_names avec les identifiants actuels
    updated_seurat$cluster_names <- as.character(Idents(updated_seurat))
    
    # Vérifier que les noms sont correctement stockés
    if(!all(sort(unique(as.character(Idents(updated_seurat)))) == 
            sort(unique(as.character(updated_seurat$cluster_names))))) {
      showNotification("Warning: Cluster names mismatch between Idents and metadata", type = "warning")
    } else {
      showNotification("Cluster names synchronized in metadata", type = "message")
    }
    
    # Mettre à jour l'objet Seurat
    single_dataset_object(updated_seurat)
    
    # Mettre à jour les options de sélection des clusters avec les nouveaux noms de clusters
    updateSelectInput(session, "select_cluster", choices = unique(Idents(single_dataset_object())))
    updateSelectInput(session, "cluster_select", choices = unique(Idents(single_dataset_object())))
  })
  # Fonction centralisée pour récupérer ou initialiser les couleurs des clusters
  get_cluster_colors <- function(seurat_object) {
    if (!is.null(seurat_object@misc$cluster_colors)) {
      cluster_colors <- seurat_object@misc$cluster_colors
    } else {
      cluster_colors <- scales::hue_pal()(length(unique(Idents(seurat_object))))
      names(cluster_colors) <- sort(unique(Idents(seurat_object)))
    }
    return(cluster_colors)
  }

  # Sauvegarder les couleurs des clusters dans l'objet Seurat après mise à jour
  observeEvent(input$update_colour, {
    message("Update the color of the selected cluster")

    # Charger les couleurs actuelles à partir de l'objet Seurat
    updated_seurat <- single_dataset_object()

    # Récupérer les couleurs actuelles ou initialiser les couleurs
    cluster_colors <- get_cluster_colors(updated_seurat)

    # Mettre à jour la couleur du cluster sélectionné
    if (!is.null(input$cluster_select) && input$cluster_select %in% names(cluster_colors)) {
      cluster_colors[input$cluster_select] <- input$cluster_colour
      message(paste("Updating color for cluster", input$cluster_select, "to", input$cluster_colour))
    } else {
      showNotification("Selected cluster is not valid.", type = "error")
      return()
    }

    # Sauvegarder les nouvelles couleurs dans @misc
    updated_seurat@misc$cluster_colors <- cluster_colors
    single_dataset_object(updated_seurat)  # Mettre à jour l'objet Seurat avec les nouvelles couleurs

    # Vérifier la mise à jour des couleurs
    message("Current cluster colors:")
    print(updated_seurat@misc$cluster_colors)

    showNotification("Cluster colors saved in Seurat object.", type = "message")
  })



  # Fonction pour afficher le UMAP avec les couleurs mises à jour
  output$umap_finale <- renderPlotly({
    req(single_dataset_object())
    message("Generating the final UMAP")

    updated_seurat <- single_dataset_object()

    # Récupérer les couleurs depuis @misc
    cluster_colors <- get_cluster_colors(updated_seurat)

    # Générer le plot UMAP avec les couleurs mises à jour
    plot_data <- DimPlot(
      updated_seurat,
      group.by = "ident",
      pt.size = input$pt_size,
      label = TRUE,
      label.size = input$label_font_size
    ) +
      scale_color_manual(values = cluster_colors) +
      NoLegend() +
      theme(axis.line = element_line(size = 0.5))

    # Convertir en plot interactif avec plotly
    interactive_plot <- ggplotly(plot_data, tooltip = "text") %>%
      layout(
        title = list(text = input$plot_title, font = list(size = 24)),
        hovermode = "closest"
      )

    message("UMAP final généré")
    return(interactive_plot)
  })


  output$save_seurat_object_1 <- createDownloadHandler(
    reactive_data = single_dataset_object,
    object_name_reactive = reactive({ getObjectNameForDownload(single_dataset_object(), default_name = "SingleDataset") }),
    data_name = "seurat_object",
    download_type = "seurat",
    show_modal = TRUE
  )
  
  ############################## Find markers for a specific cluster ##############################

  # Variables réactives pour cet onglet
  gene_tables_10 <- reactiveValues()
  all_markers_10 <- reactiveVal()
  umap_plot <- reactiveVal()
  
  
  # Function to calculate all differential markers once only
  clean_gene_names_for_html <- function(gene_names) {
    tryCatch({
      # Obtenir les noms de gènes originaux depuis les rownames
      original_names <- rownames(single_dataset_object())
      
      # Ne supprimer les suffixes que si le nom de gène existe sans suffixe
      cleaned_names <- sapply(gene_names, function(gene) {
        base_name <- gsub("\\.\\d+$", "", gene)  # Supprime le .X à la fin
        if (base_name %in% original_names) {
          return(base_name)
        }
        return(gene)  # Garde le suffixe si nécessaire
      })
      
      return(cleaned_names)
    }, error = function(e) {
      showNotification(paste("Error cleaning gene names:", e$message), type = "error")
      return(gene_names)  # Retourne les noms non modifiés en cas d'erreur
    })
  }
  
  

  # UMAP pour un cluster spécifique avec couleurs mises à jour
  output$umap_plot <- renderPlot({
    req(single_dataset_object())
    message("UMAP generation for a specific cluster")

    # Utiliser la fonction centralisée pour obtenir les couleurs
    cluster_colors <- get_cluster_colors(single_dataset_object())

    # Gérer l'affichage des labels en fonction de la case cochée
    label_option <- input$show_labels  # TRUE si la case est cochée, FALSE sinon

    # Générer le plot UMAP avec ou sans les labels selon l'option choisie
    plot_data <- DimPlot(
      single_dataset_object(),
      group.by = "ident",
      pt.size = input$pt_size,
      label = label_option,  # Activer ou désactiver les labels
      label.size = input$label_font_size
    ) +
      theme(axis.line = element_line(size = 0.5)) +
      theme_void() +
      theme(legend.position = "none") +  # Forcer la suppression de la légende
      scale_color_manual(values = cluster_colors)  # Appliquer les couleurs des clusters

    plot_data <- plot_data + NoLegend() + ggtitle(NULL)
    umap_plot(plot_data)
    message("UMAP for a specific cluster generated")
    return(plot_data)
  })


  # Mettre à jour les choix pour le selectInput cluster
  observe({
    req(single_dataset_object())
    updateSelectInput(session, "cluster", choices = levels(single_dataset_object()))
  })

  # Download handler for UMAP cluster plot
  output$downloadUMAPCluster <- createDownloadHandler(
    reactive_data = umap_plot,
    object_name_reactive = reactive({ getObjectNameForDownload(single_dataset_object(), default_name = "UMAP_plot") }),
    data_name = "UMAP",
    download_type = "plot",
    plot_params = list(
      file_type = input$umap_cluster_format,
      width = ifelse(input$umap_cluster_format == "pdf", 11, 10),
      height = ifelse(input$umap_cluster_format == "pdf", 8, 6),
      dpi = input$dpi_umap_cluster
    )
  )
  
  # Download handler for Seurat object
  output$save_seurat_object_3 <- createDownloadHandler(
    reactive_data = single_dataset_object,
    object_name_reactive = reactive({ getObjectNameForDownload(single_dataset_object(), default_name = "SingleDataset") }),
    data_name = "seurat_object",
    download_type = "seurat",
    show_modal = TRUE
  )
  

  calculate_markers_for_cluster <- function() {
    # Vérifier que les objets et les entrées nécessaires sont disponibles
    req(single_dataset_object(), input$cluster, input$min_pct_single, input$logfc_threshold_single)
    
    tryCatch({
      # Calcul des marqueurs différentiels
      markers <- FindMarkers(single_dataset_object(),
                             ident.1 = input$cluster,
                             min.pct = input$min_pct_single,
                             logfc.threshold = input$logfc_threshold_single)
      
      # Formatage des valeurs p et ajustement des valeurs p pour l'affichage
      markers$p_val <- format(markers$p_val, scientific = TRUE, digits = 3)
      markers$p_val_adj <- format(markers$p_val_adj, scientific = TRUE, digits = 3)
      # Nettoyer les noms de gènes avant de générer le HTML
      cleaned_gene_names <- clean_gene_names_for_html(rownames(markers))
      markers$gene <- paste0('<a href="#" class="gene-name" data-gene="', cleaned_gene_names, '">', cleaned_gene_names, '</a>')
      
      # Stocker les résultats dans une valeur réactive
      all_markers_10(markers)
      
      # Stocker dans la liste globale pour les diagrammes de Venn
      table_name <- paste0("SingleDataset_Cluster_", input$cluster, "_vs_All_", format(Sys.time(), "%H%M%S"))
      description <- paste0("Cluster ", input$cluster, " vs All Others")
      parameters <- list(
        min_pct = input$min_pct_single,
        logfc_threshold = input$logfc_threshold_single
      )
      
      result <- storeDETable(single_gene_table_storage(), markers, table_name, description, "cluster_group", parameters)
      if (result$success) {
        single_gene_table_storage(result$storage)
      }
      
      # Retourner TRUE en cas de succès
      return(TRUE)
    }, error = function(e) {
      # Gestion des erreurs, par exemple, affichage d'une notification dans l'application
      showNotification(paste0("Error in calculate_markers_for_cluster: ", e$message), type = "error")
      # Retourner FALSE en cas d'erreur
      return(FALSE)
    })
  }
  
  # Observer for the find_markers button
  observeEvent(input$find_markers, {
    tryCatch({
      # Show modal dialog
      showModal(modalDialog(
        title = "Please Wait",
        "Finding markers for the selected cluster...",
        easyClose = FALSE,
        footer = NULL
      ))

      calculate_markers_for_cluster()
      update_gene_tables_display_10()

      # Close the modal dialog
      removeModal()
    }, error = function(e) {
      removeModal() # Close the modal dialog in case of an error
      showNotification(paste0("Error finding markers: ", e$message), type = "error")
    })
  })

  # Fonction pour mettre à jour le tableau de gènes pour l'onglet 10
  update_gene_tables_display_10 <- function() {
    markers <- all_markers_10()
    num_genes_to_display_10 <- input$number_genes_10

    gene_tables_10$table <- head(markers, n = num_genes_to_display_10)

    output$gene_tables_10 <- renderUI({
      tags$div(style = "width: 100%; font-size: 75%;",
               tagList(
                 DTOutput("table_10")
               )
      )
    })
  }

  # Observeur pour la mise à jour de l'affichage uniquement quand le nombre de gènes change
  observeEvent(input$number_genes_10, {
    if (!is.null(input$number_genes_10) && !is.null(all_markers_10())) {
      update_gene_tables_display_10()
    }
  })

  # Ce bloc s'exécute pour mettre à jour le tableau de sortie pour l'onglet 10
  observe({
    output$table_10 <- renderDataTable({
      datatable(gene_tables_10$table, escape = FALSE)
    })
  })

  observe({
    req(single_dataset_object())
    updateSelectInput(session, "cluster1", choices = levels(single_dataset_object()))
  })

  observeEvent(input$cluster1, {
    updateSelectInput(session, "cluster2",
                      choices = setdiff(levels(single_dataset_object()), input$cluster1))
  })

  # Variables réactives
  gene_tables_new <- reactiveValues()
  all_markers_new <- reactiveVal()

  calculate_markers_for_comparison <- function() {
    req(single_dataset_object(), input$cluster1, input$cluster2,
        input$min_pct_comparison, input$logfc_threshold_comparison)
    
    tryCatch({
      # Vérifier qu'il n'y a pas de chevauchement entre les groupes
      if(any(input$cluster1 %in% input$cluster2)) {
        showNotification("Les groupes de clusters doivent être distincts", type = "error")
        return(FALSE)
      }
      
      seurat_obj <- single_dataset_object()
      
      # Créer une nouvelle identité pour les groupes de comparaison
      new_idents <- as.character(Idents(seurat_obj))
      new_idents[new_idents %in% input$cluster1] <- "group1"
      new_idents[new_idents %in% input$cluster2] <- "group2"
      Idents(seurat_obj) <- new_idents
      
      # Calcul des marqueurs
      markers <- FindMarkers(seurat_obj,
                             ident.1 = "group1",
                             ident.2 = "group2",
                             min.pct = input$min_pct_comparison,
                             logfc.threshold = input$logfc_threshold_comparison)
      
      # Formater les résultats
      markers$gene <- rownames(markers)
      markers$p_val <- format(markers$p_val, scientific = TRUE, digits = 3)
      markers$p_val_adj <- format(markers$p_val_adj, scientific = TRUE, digits = 3)
      
      # Ajouter des informations sur les groupes comparés
      markers$comparison <- sprintf("Group1(%s) vs Group2(%s)",
                                    paste(input$cluster1, collapse=","),
                                    paste(input$cluster2, collapse=","))
      
      # Ajouter les liens HTML
      cleaned_gene_names <- clean_gene_names_for_html(markers$gene)
      markers$gene <- paste0('<a href="#" class="gene-name" data-gene="',
                             cleaned_gene_names, '">',
                             cleaned_gene_names, '</a>')
      
      # Stocker les résultats
      all_markers_new(markers)
      
      # Stocker dans la liste globale pour les diagrammes de Venn
      group1_text <- paste(input$cluster1, collapse = "_")
      group2_text <- paste(input$cluster2, collapse = "_")
      table_name <- paste0("SingleDataset_Clusters_", group1_text, "_vs_", group2_text, "_", format(Sys.time(), "%H%M%S"))
      description <- paste0("Clusters [", group1_text, "] vs [", group2_text, "]")
      parameters <- list(
        min_pct = input$min_pct_comparison,
        logfc.threshold = input$logfc_threshold_comparison,
        group1 = input$cluster1,
        group2 = input$cluster2
      )
      
      result <- storeDETable(single_gene_table_storage(), markers, table_name, description, "cluster_group", parameters)
      if (result$success) {
        single_gene_table_storage(result$storage)
      }
      
      return(TRUE)
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
      return(FALSE)
    })
  }


  # Fonction pour mettre à jour le tableau de gènes pour cet onglet
  update_gene_tables_display_new <- function() {
    markers <- all_markers_new()
    num_genes_to_display_new <- input$n_diff_markers

    gene_tables_new$table <- head(markers, n = num_genes_to_display_new)

    output$gene_tables_new <- renderUI({
      tags$div(style = "width: 100%; font-size: 75%;",
               tagList(
                 DTOutput("table_new")
               )
      )
    })
  }

  # Observer for the compare_markers button
  observeEvent(input$compare_markers, {
    tryCatch({
      # Show modal dialog
      showModal(modalDialog(
        title = "Please Wait",
        "Comparing markers...",
        easyClose = FALSE,
        footer = NULL
      ))

      calculate_markers_for_comparison()
      update_gene_tables_display_new()

      # Close the modal dialog
      removeModal()
    }, error = function(e) {
      removeModal() # Close the modal dialog in case of an error
      showNotification(paste0("Error comparing markers: ", e$message), type = "error")
    })
  })

  # Observer to update display only when gene count changes
  observeEvent(input$n_diff_markers, {
    if (!is.null(input$n_diff_markers) && !is.null(all_markers_new())) {
      update_gene_tables_display_new()
    }
  })

  # Observer to update the output table for this tab
  observe({
    output$table_new <- renderDataTable({
      datatable(gene_tables_new$table, escape = FALSE)
    })
  })



  # Download gene comparison table - single cluster
  output$download_markers_single_cluster <- createDownloadHandler(
    reactive_data = reactive({ all_markers_10()[order(all_markers_10()$p_val), ] }),
    object_name_reactive = reactive({ getObjectNameForDownload(single_dataset_object(), default_name = "DiffGenes") }),
    data_name = reactive({ paste("comparison", input$cluster, "VS", "all_others_clusters", sep = "_") }),
    download_type = "csv"
  )
  
  # Download gene comparison table - multiple clusters
  output$download_markers_multiple_clusters <- createDownloadHandler(
    reactive_data = reactive({ all_markers_new()[order(all_markers_new()$p_val), ] }),
    object_name_reactive = reactive({ getObjectNameForDownload(single_dataset_object(), default_name = "DiffGenes") }),
    data_name = reactive({ paste("comparison", input$cluster1, "VS", input$cluster2, sep = "_") }),
    download_type = "csv"
  )
  
  
  
######################################### Venn Diagramm ################################
  
  # Initialize reactive storage
  single_gene_table_storage <- reactiveVal(list())
  single_current_gene_lists <- reactiveVal(NULL)
  single_venn_plot_rendered <- reactiveVal(NULL)
  
  # Update select inputs
  observe({
    updateVennSelectInputs(session, single_gene_table_storage(), c("venn_table_1_single", "venn_table_2_single", "venn_table_3_single"))
  })
  
  # Enable/disable generate button
  observe({
    tables <- single_gene_table_storage()
    if (length(tables) < 2) {
      shinyjs::disable("generate_venn_btn_single")
    } else {
      shinyjs::enable("generate_venn_btn_single")
    }
  })
  
  # Generate Venn diagram
  observeEvent(input$generate_venn_btn_single, {
    showModal(modalDialog(title = "Generating Venn Diagram", "Processing...", easyClose = FALSE, footer = NULL))
    result <- processVennGeneration(
      table_storage = single_gene_table_storage(),
      selected_tables = c(input$venn_table_1_single, input$venn_table_2_single, input$venn_table_3_single),
      filter_params = list(
        significant_only = c(input$significant_only_venn_1_single, input$significant_only_venn_2_single, input$significant_only_venn_3_single),
        log_fc_threshold = c(input$log_fc_threshold_venn_1_single, input$log_fc_threshold_venn_2_single, input$log_fc_threshold_venn_3_single),
        p_val_threshold = input$p_val_threshold_venn_single,
        use_adjusted_p = input$use_adjusted_p_venn_single,
        direction = input$venn_direction_single
      ),
      colors = c(input$venn_color_1_single, input$venn_color_2_single, input$venn_color_3_single)
    )
    
    removeModal()
    
    if (result$success) {
      single_venn_plot_rendered(result$venn_plot)
      single_current_gene_lists(result$overlaps)
      updateSelectInput(session, "selected_gene_set_single", choices = names(result$overlaps))
      shinyjs::enable("download_venn_diagram_single")
    } else {
      showNotification(result$message, type = "error")
    }
  })
  
  # Render Venn diagram
  output$venn_plot_single <- renderPlot({
    req(single_venn_plot_rendered())
    grid.draw(single_venn_plot_rendered())
  })
  
  # Display gene table
  output$venn_gene_table_single <- renderDT({
    req(single_current_gene_lists(), input$selected_gene_set_single)
    overlaps <- single_current_gene_lists()
    selected_genes <- overlaps[[input$selected_gene_set_single]]
    if (length(selected_genes) == 0) {
      return(data.frame(Gene = character(0)))
    }
    gene_df <- data.frame(Gene = selected_genes)
    datatable(gene_df, options = list(pageLength = 15, scrollX = TRUE, dom = 'Bfrtip', buttons = c('copy', 'csv', 'excel')), rownames = FALSE)
  })
  
  

  
  
  output$download_venn_diagram_single <- createDownloadHandler(
    reactive_data = single_venn_plot_rendered,
    object_name_reactive = reactive({ getObjectNameForDownload(single_dataset_object(), default_name = "VennDiagram") }),
    data_name = "venn_comparison",
    download_type = "plot",
    plot_params = list(
      file_type = input$venn_diagram_format_single,
      width = 8,
      height = 6,
      dpi = input$venn_diagram_dpi_single
    )
  )
  

  output$download_venn_gene_lists_single <- createDownloadHandler(
    reactive_data = single_current_gene_lists,
    object_name_reactive = reactive({ getObjectNameForDownload(single_dataset_object(), default_name = "VennDiagram") }),
    data_name = "gene_lists",
    download_type = "ods"
  )
  
##########################CLuster Composition#########################
  
  # Add this reactive value with other reactive values in single dataset server
  cluster_composition_single <- reactiveVal(NULL)
  
  # Observer for generating cluster composition table for single dataset
  observeEvent(input$generate_cluster_composition_single, {
    req(single_dataset_object())
    
    tryCatch({
      cluster_composition <- create_cluster_composition_table(
        single_dataset_object(),
        is_integrated = FALSE
      )
      cluster_composition_single(cluster_composition)
      
    }, error = function(e) {
      showNotification(paste("Error:", e$message), type = "error")
    })
  })
  
  output$cluster_composition_single <- renderDT({
    req(cluster_composition_single())
    render_cluster_composition_table(cluster_composition_single(), is_integrated = FALSE)
  })
  
  # Download handler for single dataset cluster composition
  output$download_cluster_composition_single <- createDownloadHandler(
    reactive_data = reactive({ 
      data <- cluster_composition_single()
      data$Size_Bar <- NULL  # Remove HTML column
      return(data)
    }),
    object_name_reactive = reactive({ getObjectNameForDownload(single_dataset_object(), default_name = "ClusterComposition") }),
    data_name = "cluster_composition",
    download_type = "csv"
  )
  
  
  #####################Co-expression###############################
  # Add these reactive values at the beginning of your server function
  gene_coexpression_data_single <- reactiveVal(NULL)
  coexpression_plot_single <- reactiveVal(NULL)
  observeEvent(input$analyze_coexpression_single, {
    req(single_dataset_object(), input$gene_text_coexpression_single)
    
    tryCatch({
      # Validate input genes before analysis
      genes_input <- trimws(strsplit(input$gene_text_coexpression_single, ",")[[1]])
      genes_input <- genes_input[genes_input != ""]
      
      if (length(genes_input) != 2) {
        showNotification("Please enter exactly 2 gene names for co-expression analysis", type = "error")
        return()
      }
      
      # Create thresholds vector with the two specific thresholds
      thresholds <- c(input$gene1_threshold_single, input$gene2_threshold_single)
      names(thresholds) <- genes_input
      
      # Validate genes exist in the dataset
      available_genes <- rownames(single_dataset_object()[[DefaultAssay(single_dataset_object())]])
      missing_genes <- setdiff(genes_input, available_genes)
      if (length(missing_genes) > 0) {
        showNotification(paste("Genes not found in dataset:", paste(missing_genes, collapse = ", ")), type = "error")
        return()
      }
      
      # Analyze gene coexpression with proper error handling
      coexpr_results <- analyze_gene_coexpression(
        seurat_obj = single_dataset_object(),
        genes = genes_input,  # ✅ Ensure genes is properly passed
        assay_name = DefaultAssay(single_dataset_object()),
        expression_thresholds = thresholds,  # Use individual thresholds instead of single threshold
        is_integrated = FALSE
      )
      
      # Validate results before creating plot
      if (is.null(coexpr_results) || is.null(coexpr_results$data)) {
        showNotification("No coexpression data generated", type = "error")
        return()
      }
      
      # Store results
      gene_coexpression_data_single(coexpr_results)
      
      # Create plot with proper validation
      if (!is.null(coexpr_results$genes_analyzed) && length(coexpr_results$genes_analyzed) > 0) {
        plot <- create_coexpression_plot(
          coexpr_results$data,           # ✅ Remove 'data =' to fix the error
          coexpr_results$genes_analyzed  # ✅ Remove 'genes =' to fix the error
        )
        coexpression_plot_single(plot)
        

      } else {
        showNotification("No genes available for plotting", type = "error")
      }
      
    }, error = function(e) {
      showNotification(paste("Error in coexpression analysis:", e$message), type = "error")
      message("Coexpression error details: ", e$message)
    })
  })
  
  output$gene_coexpression_table_single <- renderDT({
    req(gene_coexpression_data_single())
    render_coexpression_table(gene_coexpression_data_single(), "gene_coexpression_table_single")
  })
  

  # Render plot
  output$gene_coexpression_plot_single <- renderPlot({
    req(coexpression_plot_single())
    coexpression_plot_single()
  }, height = 600)
  
  # Download handlers for coexpression analysis
  output$download_coexpression_table_single <- createDownloadHandler(
    reactive_data = gene_coexpression_data_single,
    object_name_reactive = reactive({ getObjectNameForDownload(single_dataset_object(), default_name = "Coexpression") }),
    data_name = "coexpression_analysis",
    download_type = "csv"
  )
  
  output$download_coexpression_plot_single <- createDownloadHandler(
    reactive_data = coexpression_plot_single,
    object_name_reactive = reactive({ getObjectNameForDownload(single_dataset_object(), default_name = "Coexpression") }),
    data_name = "coexpression_plot",
    download_type = "plot",
    plot_params = list(
      file_type = "pdf",
      width = 10,
      height = 8,
      dpi = 300
    )
  )
  ############################## Subseting ##############################

  # Variables wich stock the subset of single_dataset_object()
  subset_seurat <- reactiveVal()

  observe({
    { subset_seurat(single_dataset_object())
    }
  })

  observe({
    if (!is.null(single_dataset_object())) {
      updateSelectInput(session, "select_ident_subset", choices = unique(Idents(single_dataset_object())))
    }
  })

  #Run subset
  observeEvent(input$apply_subset, {
    req(single_dataset_object())
    tryCatch({
      subset_seurat_temp <- single_dataset_object()

      if (length(input$select_ident_subset) > 0) {
        subset_seurat_temp <- subset(x = subset_seurat_temp, idents = input$select_ident_subset)

        if (nrow(subset_seurat_temp@meta.data) == 0) {
          showNotification("No cells found with the selected identities.", type = "error")
          return()
        }
      }

      subset_seurat(subset_seurat_temp)



    }, error = function(e) {
      showNotification(paste("Error applying subset: ", e$message), type = "error")
    })
  })




  observeEvent(input$apply_gene_subset, {
    tryCatch({
      req(single_dataset_object())

      # Vérifiez que gene_list_subset est correctement récupéré
      print(paste("Raw gene list input:", input$gene_list_subset))

      # Récupérer et traiter la liste des gènes
      gene_list <- unlist(strsplit(input$gene_list_subset, ","))
      gene_list <- trimws(gene_list)
      gene_list <- as.character(gene_list)  # Convertir explicitement en caractères

      if (length(gene_list) == 0 || all(gene_list == "")) {
        showNotification("Please enter valid gene names.", type = "error")
        return()
      }

      # Vérifier si tous les gènes sont présents dans l'objet Seurat
      seurat_genes <- rownames(LayerData(single_dataset_object(), assay = "RNA", layer = 'counts'))
      print(paste("Seurat genes:", paste(head(seurat_genes), collapse = ", ")))  # Afficher les premiers gènes pour vérifier

      missing_genes <- setdiff(gene_list, seurat_genes)
      print(paste("Missing genes:", paste(missing_genes, collapse = ", ")))

      if (length(missing_genes) > 0) {
        showNotification(paste("The following genes are not found in the dataset:", paste(missing_genes, collapse = ", ")), type = "error")
        return()
      }

      # Créer une matrice logique pour l'expression génique
      expression_matrix <- sapply(gene_list, function(gene) {
        data <- FetchData(single_dataset_object(), vars = gene)
        if (is.list(data)) {
          data <- unlist(data)
        }
        print(paste("Data for gene", gene, ":", paste(data, collapse = ", ")))  # Afficher les données de chaque gène
        as.numeric(data) >= input$expression_threshold
      })

      # Convertir la matrice logique en dataframe pour éviter les erreurs de typage
      expression_matrix <- as.data.frame(expression_matrix)

      # Sélection des cellules exprimant le nombre minimum de gènes spécifiés
      cells_to_keep <- colnames(single_dataset_object())[which(rowSums(expression_matrix) >= input$num_genes_to_express)]

      # Vérifier si des cellules sont sélectionnées
      if (length(cells_to_keep) == 0) {
        showNotification("No cells found with the specified gene expression criteria.", type = "error")
        return()
      }

      subset_seurat_temp <- subset(single_dataset_object(), cells = cells_to_keep)
      subset_seurat(subset_seurat_temp)
      showNotification(paste("Subsetting applied. Number of cells retained:", length(cells_to_keep)), type = "message")
    }, error = function(e) {
      showNotification(paste0("Error during gene-based subset: ", e$message), type = "error")
    })
  })


  # Downloading Seurat subset object
  output$download_subset_seurat <- createDownloadHandler(
    reactive_data = subset_seurat,
    object_name_reactive = reactive({ getObjectNameForDownload(single_dataset_object(), default_name = "SeuratSubset") }),
    data_name = "seurat_subset",
    download_type = "seurat",
    show_modal = TRUE
  )
  
  # Umap plot with all clusters
  reactivePlotAll <- reactive({
    req(single_dataset_object())
    plot_data_all <- DimPlot(single_dataset_object(), group.by = "ident", label = TRUE) +
      theme(axis.line = element_line(size = 0.5)) +
      NoLegend()+ggtitle(NULL)
    return(plot_data_all)
  })

  # Umap plot with filtered data
  reactivePlotSubset <- reactive({
    req(subset_seurat())
    plot_data_subset <- DimPlot(subset_seurat(), group.by = "ident", label = TRUE) +
      theme(axis.line = element_line(size = 0.5)) +
      NoLegend()+ggtitle(NULL)
    return(plot_data_subset)
  })

  # Rendering Umap plot with all clusters
  output$global_umap <- renderPlot({
    plot_data_all <- reactivePlotAll()
    print(plot_data_all)
  })

  # Rendering Umap plot with filtered data
  output$subset_umap <- renderPlot({
    plot_data_subset <- reactivePlotSubset()
    print(plot_data_subset)
  })





}

