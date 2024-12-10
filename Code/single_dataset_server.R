############################## Script Single Dataset Analysis ##############################

# single_dataset_server.R

single_dataset_server <- function(input, output, session) {

  ############################## Loading Data ##############################
  # Tab 1 : Loading Data

  shinyjs::useShinyjs() # button deactivation

  # Reactives objects for the single dataset part
  single_dataset_object <- reactiveVal(NULL)
  gene_list <- reactiveValues(features = NULL)

  # Function to clean the workspace
  cleanWorkspaceSingleDataset <- function() {
    message("Cleaning workspace...")

    # Clean reactive objects
    single_dataset_object(NULL)
    gene_list$features <- NULL

    subset_seurat(NULL)
    clustering_plot(NULL)
    feature_plot(NULL)
    vln_plot(NULL)
    dot_plot(NULL)
    ridge_plot(NULL)
    heatmap_plot(NULL)
    scatter_plot(NULL)

    # Delete temp directories if they exist
    if (dir.exists("tempdir")) {
      unlink("tempdir", recursive = TRUE, force = TRUE)
      message("Deleted tempdir.")
    }
    if (dir.exists("dataDir")) {
      unlink("dataDir", recursive = TRUE, force = TRUE)
      message("Deleted dataDir.")
    }

    # Delete all directories starting with "unzipped"
    unzipped_dirs <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
    unzipped_dirs <- unzipped_dirs[grepl("unzipped", unzipped_dirs)]
    if (length(unzipped_dirs) > 0) {
      sapply(unzipped_dirs, unlink, recursive = TRUE, force = TRUE)
      message("Deleted unzipped directories.")
    }

    gc()  # Garbage collection
    message("Workspace cleaned.")
  }


  # Fonction réactive pour récupérer les clusters depuis l'objet Seurat
  get_clusters <- reactive({
    req(single_dataset_object())  # Assure que l'objet Seurat est disponible
    return(levels(Idents(single_dataset_object())))  # Retourne les clusters
  })



  # Function to unzip and analyze contents with handling nested ZIP files
  extract_and_process_zip <- function(zip_path, exdir) {
    message(paste("Extracting zip file:", zip_path, "to", exdir))
    unzip(zip_path, exdir = exdir)

    # Look for nested zip files
    nested_zips <- list.files(exdir, pattern = "\\.zip$", recursive = TRUE, full.names = TRUE)
    while (length(nested_zips) > 0) {
      for (nested_zip in nested_zips) {
        message(paste("Extracting nested zip file:", nested_zip))
        unzip(nested_zip, exdir = exdir)
        unlink(nested_zip)
      }
      nested_zips <- list.files(exdir, pattern = "\\.zip$", recursive = TRUE, full.names = TRUE)
    }

    # Decompress .gz files if necessary
    gz_files <- list.files(exdir, pattern = "\\.gz$", recursive = TRUE, full.names = TRUE)
    if (length(gz_files) > 0) {
      lapply(gz_files, function(file) {
        tryCatch({
          message(paste("Decompressing file:", file))
          gunzip(file, overwrite = TRUE)
        }, error = function(e) {
          message(paste("Error decompressing file:", file, "->", e$message))
        })
      })
    }

    return(exdir)
  }

  # Function to compress uncompressed files
  compress_files <- function(files) {
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

  # Function to process the dataset after loading
  processDataset <- function(dataDir, dataset_type, mt_pattern) {
    cleanWorkspaceSingleDataset()
    tryCatch({
      message("Processing dataset...")

      # List of files after extraction
      all_files <- list.files(dataDir, recursive = TRUE, full.names = TRUE)
      message("Files found in data directory:", paste(all_files, collapse = ", "))

      # Check for required files
      compressed_files <- grep("\\.gz$", all_files, value = TRUE)
      uncompressed_files <- grep("matrix.mtx$|features.tsv$|barcodes.tsv$", all_files, value = TRUE)

      if (length(compressed_files) == 3) {
        message("Using compressed files.")
        file_paths <- compressed_files
      } else if (length(uncompressed_files) == 3) {
        message("Found uncompressed files. Compressing...")
        file_paths <- compress_files(uncompressed_files)
      } else {
        stop("Necessary 10X files (matrix.mtx, features.tsv, barcodes.tsv) are not found in the specified directory.")
      }

      dataDir <- dirname(file_paths[1])
      message("Using data directory:", dataDir)

      if (dataset_type == "snRNA") {
        single_dataset_data <- Read10X(dataDir)
        if (is.null(single_dataset_data)) {
          stop("Unable to load 10X data. Please check the file format.")
        }
        seuratObj <- CreateSeuratObject(counts = single_dataset_data, project = "single_dataset", min.cells = 3, min.features = 200)
        seuratObj[["percent.mt"]] <- PercentageFeatureSet(seuratObj, pattern = mt_pattern)
        single_dataset_object(seuratObj)
        showNotification("snRNA-seq data processed successfully.", type = "message")
      } else if (dataset_type == "multiome") {
        rna.data <- Read10X(dataDir)$`Gene Expression`
        atac.data <- Read10X(dataDir)$`Peaks`
        if (is.null(rna.data) || is.null(atac.data)) {
          stop("Unable to load Multiome data. Please check the file format.")
        }
        rna.seurat <- CreateSeuratObject(counts = rna.data, project = "RNA", min.cells = 3, min.features = 200)
        atac.seurat <- CreateSeuratObject(counts = atac.data, project = "ATAC", min.cells = 3, min.features = 200)
        rna.seurat[["percent.mt"]] <- PercentageFeatureSet(rna.seurat, pattern = mt_pattern)
        single_dataset_object(rna.seurat)
        showNotification("Multiome data processed successfully.", type = "message")
      }
      updateUIElements()
    }, error = function(e) {
      showNotification(paste("Error processing data: ", e$message), type = "error")
      message("Error during dataset processing:", e$message)
    })
  }

  # Observer pour le chargement de données
  observeEvent(input$file, {
    cleanWorkspaceSingleDataset()

    showNotification("Uploading and processing data...", type = "message")
    mt_pattern <- ifelse(input$species_choice == "mouse", "^mt-", "^MT-")

    tryCatch({
      file_extension <- tools::file_ext(input$file$name)
      message("File extension detected:", file_extension)

      if (file_extension == "zip") {
        dataDir <- tempdir()
        extract_and_process_zip(input$file$datapath, exdir = dataDir)
      } else {
        stop("Please upload a .zip file.")
      }

      # Verify file integrity after extraction
      required_files <- c("matrix.mtx", "features.tsv", "barcodes.tsv")
      files_found <- list.files(dataDir, recursive = TRUE, full.names = TRUE)

      if (!all(required_files %in% basename(files_found))) {
        stop("Necessary 10X files (matrix.mtx, features.tsv, barcodes.tsv) are not found in the specified directory.")
      }

      processDataset(dataDir, input$dataset_type, mt_pattern)

      # Call to update the UI elements after the data is loaded
      updateUIElements()

    }, error = function(e) {
      showNotification(paste("Error processing files: ", e$message), type = "error")
      message("Error during file processing:", e$message)
    })
  })



  # Charger et restaurer les couleurs après chargement de l'objet Seurat
  observeEvent(input$load_seurat_file, {
    cleanWorkspaceSingleDataset()
    message("Attempting to read file at: ", input$load_seurat_file$datapath)
    tryCatch({
      loaded_seurat <- readRDS(input$load_seurat_file$datapath)
      message("File successfully read.")
      single_dataset_object(loaded_seurat)

      # Si des couleurs de clusters sont présentes, les restaurer
      if (!is.null(loaded_seurat@misc$cluster_colors)) {
        message("Restoring cluster colors from Seurat object.")
        print(loaded_seurat@misc$cluster_colors)
      } else {
        message("No cluster colors found in Seurat object. Applying default colors.")
      }

      showNotification("The Seurat object has been successfully loaded!")
      updateUIElements()

    }, error = function(e) {
      showNotification(paste("An error occurred: ", e), type = "error")
      message("An error occurred: ", e)
    })
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


  output$downloadUMAP <- downloadHandler(
    filename = function() {
      paste("UMAP_plot", Sys.Date(), ".tiff", sep = "")
    },
    content = function(file) {
      tryCatch({
        req(clustering_plot())
        ggsave(file, plot = clustering_plot(), dpi = input$dpi_umap, width = 10, height = 6)
      }, error = function(e) {
        showNotification(paste0("Download error: ", e$message), type = "error")
      })
    }
  )


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

  ############################## Calculate differential expressed genes for each cluster ##############################

  # Variables réactives pour cet onglet
  gene_tables <- reactiveValues()
  all_markers <- reactiveVal()

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

  # Fonction modifiée pour calculer tous les marqueurs différentiels une seule fois
  calculate_all_markers <- function() {
    if (is.null(single_dataset_object())) {
      showNotification("single_dataset_object is NULL, please check the previous steps.", type = "error")
      return(NULL)
    }
    tryCatch({
      showNotification("Calculation of differential markers in progress...", type = "message", duration = 10)
      markers <- FindAllMarkers(single_dataset_object(),
                                min.pct = input$min_pct_all_single,
                                logfc.threshold = input$logfc_threshold_all_single)

      # Examiner tous les gènes avec suffixes numériques
      all_genes <- rownames(GetAssayData(single_dataset_object()))
      genes_with_suffix <- grep("\\.\\d+$", all_genes, value=TRUE)

      # Afficher des détails
      print("Exemples de gènes avec suffixes:")
      print(head(sort(genes_with_suffix), 20))

      print("Nombre total de gènes avec suffixes:")
      print(length(genes_with_suffix))

      print("Distribution des suffixes:")
      suffixes <- as.numeric(gsub(".*\\.", "", genes_with_suffix))
      print(table(suffixes))

      print("Gènes avec plusieurs versions:")
      base_genes <- gsub("\\.\\d+$", "", genes_with_suffix)
      duplicated_genes <- unique(base_genes[duplicated(base_genes)])
      for(gene in head(duplicated_genes)) {
        print(paste("Gene:", gene))
        print(grep(paste0("^", gene, "\\."), genes_with_suffix, value=TRUE))
      }

      print("Distribution of suffixes:")
      suffixes <- gsub(".*\\.", "", grep("\\.\\d+$", rownames(GetAssayData(single_dataset_object())), value=TRUE))
      print(table(suffixes))
      markers <- as.data.frame(markers)
      cleaned_gene_names <- clean_gene_names_for_html(rownames(markers))
      markers$gene <- paste0('<span class="gene-name">', cleaned_gene_names, '</span>')
      markers$p_val <- format(markers$p_val, scientific = TRUE, digits = 20)
      markers$p_val_adj <- format(markers$p_val_adj, scientific = TRUE, digits =20)
      markers$gene <- paste0('<a href="#" class="gene-name" data-gene="', cleaned_gene_names, '">', cleaned_gene_names, '</a>')
      all_markers(markers)
    }, error = function(e) {
      showNotification(paste0("Error when calculating differential markers: ", e$message), type = "error")
    })
  }

  # Function to update gene tables according to the number selected
  update_gene_tables_display <- function() {
    tryCatch({
      markers <- all_markers()
      num_genes_to_display <- input$number_genes
      for (cluster in unique(markers$cluster)) {
        gene_tables[[paste0("table_", cluster)]] <- head(markers[markers$cluster == cluster, ], n = num_genes_to_display)
      }
      output$diff_genes_tables <- renderUI({
        tagList(
          lapply(names(gene_tables), function(name) {
            tags$div(style = "width: 100%; font-size: 75%;",
                     tagList(
                       h3(paste0("Cluster ", stringr::str_replace(name, "table_", ""))),
                       DTOutput(name),
                       hr()
                     )
            )
          })
        )
      })
    }, error = function(e) {
      showNotification(paste0("Error when updating gene tables: ", e$message), type = "error")
    })
  }

  # Observer for run_DE button to calculate all markers
  observeEvent(input$run_DE, {
    tryCatch({
      # Afficher la boîte de dialogue modale
      showModal(modalDialog(
        title = "Please Wait",
        "Calculating differential expression...",
        easyClose = FALSE,
        footer = NULL
      ))

      calculate_all_markers()
      update_gene_tables_display()

      # Fermer la boîte de dialogue modale
      removeModal()
    }, error = function(e) {
      removeModal() # Fermer la boîte de dialogue modale en cas d'erreur
      showNotification(paste0("Error calculating differential expression:", e$message), type = "error")
    })
  })

  # Observer to update display only when gene number changes
  observeEvent(input$number_genes, {
    tryCatch({
      if (!is.null(input$number_genes) && !is.null(all_markers())) {
        update_gene_tables_display()
      }
    }, error = function(e) {
      showNotification(paste0("Error:", e$message), type = "error")
    })  })

  #update the output genes table
  observe({
    lapply(names(gene_tables), function(name) {
      output[[name]] <- renderDataTable({
        datatable(gene_tables[[name]], escape = FALSE)
      })
    })
  })



  # Download seurat object
  output$save_seurat <- downloadHandler(
    tryCatch({

      filename = function() {
        showModal(modalDialog(
          title = "Please Wait",
          "Preparing the seurat object for download...",
          footer = NULL,
          easyClose = FALSE
        ))
        return("seurat_object.rds")
      }}, error = function(e) {
        showNotification(paste0("Error:", e$message), type = "error")
      }),
    content = function(file) {
      saveRDS(single_dataset_object(), file)
      removeModal()
    }
  )


  # Download markers as a CSV
  output$download_DE <- downloadHandler(
    tryCatch({
      filename = function() {
        paste("DE_genes_", Sys.Date(), ".csv", sep = "")
      }}, error = function(e) {
        showNotification(paste0("Error:", e$message), type = "error")
      }),
    content = function(file) {
      markers <- all_markers()
      write.csv(markers, file)
    }
  )



  ############################## Visualization of expressed genes ##############################
  # Observer pour mettre à jour les choix de gènes selon l'assay sélectionné 
  observeEvent(c(single_dataset_object(), input$viz_assay), {
    req(single_dataset_object(), input$viz_assay)
    gene_list <- sort(rownames(LayerData(single_dataset_object(), 
                                         assay = input$viz_assay, 
                                         layer = 'counts')))
    updatePickerInput(session, "gene_select", choices = gene_list)
    updatePickerInput(session, "gene_select_heatmap", choices = gene_list)
    updatePickerInput(session, "gene_select_genes_analysis", choices = gene_list)
  })

  # Observer pour mettre à jour les gènes sélectionnés
  observeEvent(input$gene_select, {
    selected_gene <- input$gene_select
    update_genes <- function(existing_genes, new_gene) {
      if (existing_genes == "") {
        return(new_gene)
      } else {
        genes <- strsplit(existing_genes, ",")[[1]]
        genes <- unique(c(trimws(genes), new_gene))
        return(paste(genes, collapse = ", "))
      }
    }
    updateTextInput(session, "gene_list_vln", value = update_genes(input$gene_list_vln, selected_gene))
    updateTextInput(session, "gene_list_feature", value = update_genes(input$gene_list_feature, selected_gene))
    updateTextInput(session, "gene_list_dotplot", value = update_genes(input$gene_list_dotplot, selected_gene))
    updateTextInput(session, "gene_list_ridge_plot", value = update_genes(input$gene_list_ridge_plot, selected_gene))
  })

  # Mise à jour des sélecteurs de cluster
  observe({
    if (!is.null(single_dataset_object())) {
      cluster_choices <- unique(Idents(single_dataset_object()))  # Utiliser Idents()
      updateSelectInput(session, "select_cluster", choices = cluster_choices)
      updateSelectInput(session, "ident_1", choices = cluster_choices)
      updateCheckboxGroupInput(session, "ident_2", choices = cluster_choices, selected = "")
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
      DefaultAssay(seurat_object) <-input$viz_assay
      present_genes <- genes[genes %in% rownames(LayerData(seurat_object, assay = "RNA", layer = "counts"))]
      
      if (length(present_genes) > 0) {
        print("Création du plot")
        
        # Préparer les paramètres de seuils
        min_cut <- if (!is.na(input$min_cutoff)) input$min_cutoff else NA
        max_cut <- if (!is.na(input$max_cutoff)) input$max_cutoff else NA
        
        if (input$show_coexpression && length(present_genes) > 1) {
          print("Mode coexpression activé")
          plot <- FeaturePlot(
            seurat_object,
            features = present_genes,
            blend = TRUE,
            blend.threshold = 1,
            order = TRUE,
            min.cutoff = min_cut,  # Ajout du seuil minimum
            max.cutoff = max_cut   # Ajout du seuil maximum
          )
        } else {
          print("Mode standard")
          plot <- FeaturePlot(
            seurat_object,
            features = present_genes,
            min.cutoff = min_cut,  # Ajout du seuil minimum
            max.cutoff = max_cut   # Ajout du seuil maximum
          )
        }
        
        if (input$add_noaxes_feature) {
          print("Ajout NoAxes")
          plot <- plot + NoAxes()
        }
        if (input$add_nolegend_feature) {
          print("Ajout NoLegend")
          plot <- plot + NoLegend()
        }
        
        print(paste("Seuils appliqués - Min:", min_cut, "Max:", max_cut))
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

  output$downloadFeaturePlot <- downloadHandler(
    filename = function() {
      paste("FeaturePlot", Sys.Date(), ".tiff", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = feature_plot(), dpi = input$dpi_plot, width = 10, height = 8)
    }
  )

  # VlnPlot
  vln_plot <- reactiveVal()
  
  observeEvent(input$show_vln, {
    tryCatch({
      req(input$gene_list_vln)
      print("Starting VlnPlot generation")
      
      seurat_object <- single_dataset_object()
      req(seurat_object)
      print("Got Seurat object")
      
      # Définir l'assay sélectionné
      DefaultAssay(seurat_object) <- input$viz_assay
      print(paste("Using assay:", input$viz_assay))
      
      plot <- VlnPlot(
        object = seurat_object,
        features = input$gene_list_vln,
        group.by = "seurat_clusters",  # utilisé seurat_clusters qui marchait
        pt.size = ifelse(input$hide_vln_points, 0, 1)
      )
      print("Plot created")
      
      # Ajouter les options
      if(input$add_noaxes_vln) plot <- plot + NoAxes()
      if(input$add_nolegend_vln) plot <- plot + NoLegend()
      
      vln_plot(plot)
      print("Plot stored")
      
    }, error = function(e) {
      print(paste("Basic error:", e$message))
    })
  })
  
  output$vln_plot <- renderPlot({
    req(vln_plot())
    vln_plot()
  })
  

  output$downloadVlnPlot <- downloadHandler(
    filename = function() {
      paste("VlnPlot", Sys.Date(), ".tiff", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = vln_plot(), dpi = input$dpi_plot, width = 10, height = 8)
    }
  )

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

      # Gestion des clusters sélectionnés
      if (!is.null(input$cluster_order_dot) && length(input$cluster_order_dot) > 0) {
        print(paste("Clusters sélectionnés:", paste(input$cluster_order_dot, collapse=", ")))

        # Créer une copie pour ne pas modifier l'objet original
        seurat_subset <- seurat_object

        # Redéfinir les identifiants avec le nouvel ordre
        Idents(seurat_subset) <- factor(
          Idents(seurat_subset),
          levels = input$cluster_order_dot
        )

        # Subset uniquement sur les clusters sélectionnés
        seurat_subset <- subset(seurat_subset, idents = input$cluster_order_dot)

        # Utiliser l'objet subset pour le plot
        plot <- DotPlot(
          seurat_subset,
          features = genes
        ) + RotatedAxis()
      } else {
        # Si aucun cluster n'est sélectionné, utiliser tous les clusters
        plot <- DotPlot(
          seurat_object,
          features = genes
        ) + RotatedAxis()
      }

      # Options de base
      if (input$add_noaxes_dot) plot <- plot + NoAxes()
      if (input$add_nolegend_dot) plot <- plot + NoLegend()

      # Inversion des axes si demandé
      if (input$invert_axes) {
        plot <- plot + coord_flip()
      }

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


  # Téléchargement du DotPlot
  output$downloadDotPlot <- downloadHandler(
    filename = function() {
      paste("DotPlot", Sys.Date(), ".tiff", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = dot_plot(), dpi = input$dpi_plot, width = 10, height = 8)
    }
  )


  # RidgePlot
  ridge_plot <- reactiveVal()

  observeEvent(input$show_ridge, {
    tryCatch({
      req(input$gene_list_ridge_plot, single_dataset_object())
      genes <- unique(trimws(unlist(strsplit(input$gene_list_ridge_plot, ","))))
      seurat_object <- single_dataset_object()
      req(seurat_object)
      DefaultAssay(seurat_object) <- input$viz_assay

      # Plot de base
      plot <- RidgePlot(
        seurat_object,
        features = genes
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

  output$download_ridge_plot <- downloadHandler(
    filename = function() {
      paste("RidgePlot", Sys.Date(), ".tiff", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = ridge_plot(), dpi = input$dpi_plot, width = 10, height = 8)
    }
  )

  # Analyse de l'expression des gènes
  number_of_nuclei <- reactiveVal(NULL)
  observeEvent(input$analyze_btn, {
    tryCatch({
      req(input$gene_select_genes_analysis, single_dataset_object())
      seurat_obj <- single_dataset_object()
      
      # Utiliser l'assay sélectionné
      DefaultAssay(seurat_obj) <- input$viz_assay
      # Récupérer les données d'expression des gènes dans l'assay sélectionné
      gene_data <- GetAssayData(seurat_obj, assay = input$viz_assay, slot = "counts")
      print(paste("Using assay:", input$viz_assay))
      
      available_genes <- rownames(gene_data)
      # Le reste du code reste identique...
      selected_genes <- intersect(input$gene_select_genes_analysis, available_genes)
      if (length(selected_genes) == 0) {
        showNotification("No selected genes are present in the dataset.", type = "error")
        return()
      }

      # Log pour les gènes sélectionnés
      message("Selected genes: ", paste(selected_genes, collapse = ", "))

      # Récupérer les identifiants des clusters pour chaque noyau
      cluster_info <- FetchData(seurat_obj, vars = "ident")

      # Log pour les informations des clusters
      message("Cluster information: ", paste(unique(cluster_info$ident), collapse = ", "))

      # Compter le nombre total de noyaux dans chaque cluster (tous les noyaux)
      cluster_counts_total <- table(cluster_info$ident)

      # Log pour le nombre total de noyaux dans chaque cluster
      message("Total nuclei per cluster: ", paste(names(cluster_counts_total), "=", cluster_counts_total, collapse = ", "))

      # Initialiser une liste pour stocker les résumés d'expression par gène
      expression_summary_list <- lapply(selected_genes, function(gene) {
        gene_data <- gene_data[gene, ]  # Extraire les données d'expression pour le gène

        # Déterminer les noyaux exprimant le gène (expression > 0)
        expressed_indices <- gene_data > 0

        # Log pour vérifier combien de noyaux expriment ce gène
        message("Nuclei expressing ", gene, ": ", sum(expressed_indices))

        # Obtenir les clusters pour les noyaux exprimant le gène
        clusters_expressing <- cluster_info[expressed_indices, "ident"]

        # Compter le nombre de noyaux exprimant le gène par cluster
        cluster_counts_expressed <- table(clusters_expressing)

        # Log pour le nombre de noyaux exprimant le gène par cluster
        message("Nuclei expressing ", gene, " per cluster: ", paste(names(cluster_counts_expressed), "=", cluster_counts_expressed, collapse = ", "))

        # Créer un dataframe avec les résultats
        df <- data.frame(
          Gene = gene,
          Cluster = names(cluster_counts_expressed),
          Cells_Expressed = as.numeric(cluster_counts_expressed),  # Noyaux exprimant le gène
          Total_Cells_in_Cluster = as.numeric(cluster_counts_total[names(cluster_counts_expressed)]),  # Nombre total de noyaux par cluster
          Percentage_of_Cells = (as.numeric(cluster_counts_expressed) / as.numeric(cluster_counts_total[names(cluster_counts_expressed)])) * 100  # Pourcentage
        )

        return(df)
      })

      # Filtrer les résultats valides et afficher
      valid_expression_summary <- Filter(Negate(is.null), expression_summary_list)
      if (length(valid_expression_summary) > 0) {
        expression_df <- do.call(rbind, valid_expression_summary)
        number_of_nuclei(expression_df)
      } else {
        expression_df <- data.frame()
        showNotification("No results meet the thresholds specified.", type = "warning")
      }

      # Afficher les résultats dans un tableau
      output$expression_summary <- renderDataTable({
        datatable(expression_df, options = list(pageLength = 10, scrollX = TRUE))
      })

    }, error = function(e) {
      showNotification(paste("Error processing expression data: ", e$message), type = "error")
    })
  })




  output$download_genes_number_expresion <- downloadHandler(
    filename = function() {
      paste("number_of_gene_expression", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(number_of_nuclei())
      write.csv(number_of_nuclei(), file, row.names = FALSE)
    },
    contentType = "text/csv"
  )

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

  # Rendre le plot actuel dans l'UI
  output$selected_plot_display <- renderPlot({
    req(current_plot())
    current_plot()
  })


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

  output$download_heatmap_single <- downloadHandler(
    filename = function() {
      paste("heatmap-", Sys.Date(), ".tiff", sep="")
    },
    content = function(file) {
      tryCatch({
        req(heatmap_plot())
        ggsave(file, plot = heatmap_plot(), dpi = input$dpi_heatmap_single, width = 10, height = 8)
      }, error = function(e) {
        error_message(paste("Error in Heatmap: ", e$message))
      })
    }
  )


  # Mettez à jour les choix de selectInput pour les gènes/features et les clusters
  observe({
    req(single_dataset_object())
    gene_list <- rownames(LayerData(single_dataset_object(), assay = "RNA", layer = 'counts'))
    updatePickerInput(session, "feature1_select", choices = gene_list)
    updatePickerInput(session, "feature2_select", choices = gene_list)

    cluster_list <- levels(Idents(single_dataset_object()))
    updateTextInput(session, "scatter_text_clusters", value = paste(cluster_list, collapse = ","))
    updateTextInput(session, "text_clusters", value = paste(cluster_list, collapse = ","))
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

  output$download_scatter_plot <- downloadHandler(
    filename = function() {
      paste0("FeatureScatter_", Sys.Date(), ".tiff")
    },
    content = function(file) {
      tryCatch({
        ggsave(file, plot = scatter_plot(), device = "tif", width = 10, height = 8, dpi = input$dpi_scatter)
      }, error = function(e) {
        showNotification(paste("Erreur lors du téléchargement :", e$message), type = "error")
      })
    }
  )



############################## Genes analysis ##############################

 number_of_nuclei <- reactiveVal(NULL)
  # Ajuster le seuil d'expression
  observeEvent(input$analyze_btn, {
    tryCatch({
      req(input$gene_select_genes_analysis, single_dataset_object(), input$expression_threshold)  # Ajouter le seuil d'expression

      seurat_obj <- single_dataset_object()

      # Assurez-vous que l'assay "RNA" est bien utilisé
      DefaultAssay(seurat_obj) <- "RNA"

      # Récupérer les données d'expression des gènes dans l'assay "RNA"
      gene_data <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
      available_genes <- rownames(gene_data)

      # Sélectionner les gènes présents dans le jeu de données
      selected_genes <- intersect(input$gene_select_genes_analysis, available_genes)
      if (length(selected_genes) == 0) {
        showNotification("No selected genes are present in the dataset.", type = "error")
        return()
      }

      # Log pour les gènes sélectionnés
      message("Selected genes: ", paste(selected_genes, collapse = ", "))

      # Récupérer les identifiants des clusters pour chaque noyau
      cluster_info <- FetchData(seurat_obj, vars = "ident")

      # Log pour les informations des clusters
      message("Cluster information: ", paste(unique(cluster_info$ident), collapse = ", "))

      # Compter le nombre total de noyaux dans chaque cluster (tous les noyaux)
      cluster_counts_total <- table(cluster_info$ident)

      # Log pour le nombre total de noyaux dans chaque cluster
      message("Total nuclei per cluster: ", paste(names(cluster_counts_total), "=", cluster_counts_total, collapse = ", "))

      # Initialiser une liste pour stocker les résumés d'expression par gène
      expression_summary_list <- lapply(selected_genes, function(gene) {
        gene_data <- gene_data[gene, ]  # Extraire les données d'expression pour le gène

        # Déterminer les noyaux exprimant le gène selon le seuil d'expression
        expressed_indices <- gene_data > input$expression_threshold  # Utiliser le seuil ici

        # Log pour vérifier combien de noyaux expriment ce gène
        message("Nuclei expressing ", gene, ": ", sum(expressed_indices))

        # Obtenir les clusters pour les noyaux exprimant le gène
        clusters_expressing <- cluster_info[expressed_indices, "ident"]

        # Compter le nombre de noyaux exprimant le gène par cluster
        cluster_counts_expressed <- table(clusters_expressing)

        # Log pour le nombre de noyaux exprimant le gène par cluster
        message("Nuclei expressing ", gene, " per cluster: ", paste(names(cluster_counts_expressed), "=", cluster_counts_expressed, collapse = ", "))

        # Créer un dataframe avec les résultats
        df <- data.frame(
          Gene = gene,
          Cluster = names(cluster_counts_expressed),
          Cells_Expressed = as.numeric(cluster_counts_expressed),  # Noyaux exprimant le gène
          Total_Cells_in_Cluster = as.numeric(cluster_counts_total[names(cluster_counts_expressed)]),  # Nombre total de noyaux par cluster
          Percentage_of_Cells = (as.numeric(cluster_counts_expressed) / as.numeric(cluster_counts_total[names(cluster_counts_expressed)])) * 100  # Pourcentage
        )

        return(df)
      })

      # Vérification des résultats valides et affichage
      valid_expression_summary <- Filter(Negate(is.null), expression_summary_list)

      # Si des résultats valides existent, les combiner en un seul dataframe
      if (length(valid_expression_summary) > 0) {
        expression_df <- do.call(rbind, valid_expression_summary)
        message("Data for expression summary successfully created:")
        print(expression_df)  # Log pour vérifier les données

        # Mettre à jour la variable réactive
        number_of_nuclei(expression_df)

        # Afficher les résultats dans un tableau si les données sont présentes
        output$expression_summary <- renderDataTable({
          datatable(expression_df, options = list(pageLength = 10, scrollX = TRUE))
        })
      } else {
        # Si aucun résultat n'est trouvé
        showNotification("No results meet the thresholds specified.", type = "warning")
        output$expression_summary <- renderDataTable({ NULL })  # Supprimer le tableau si aucun résultat n'est trouvé
      }

    }, error = function(e) {
      showNotification(paste("Error processing expression data: ", e$message), type = "error")
      print(paste("Error: ", e$message))  # Log complet de l'erreur pour debug
    })
  })



  # Download gene comparison table
  output$download_genes_number_expresion <- downloadHandler(
    filename = function() {
      paste("number_of_nuclei_expression", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      tryCatch({
        req(number_of_nuclei())
        write.csv(number_of_nuclei(), file, row.names = FALSE)
      }, error = function(e) {
        showNotification(paste0("Error downloading number of nuclei comparison table: ", e$message), type = "error")
      })
    },
    contentType = "text/csv"
  )


  ############################## Final UMAP ##############################
  # Tab 8: Final UMAP

  #Reactive variable for that tab
  cluster_colours <- reactiveVal()  # Initialize a list to store cluster colors

  observe({
    if (!is.null(single_dataset_object())) {
      updateSelectInput(session, "select_cluster", choices = unique(Idents(single_dataset_object())))
      updateSelectInput(session, "cluster_select", choices = unique(Idents(single_dataset_object())))

    }
  })

  #function for renaming each cluster
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

    single_dataset_object(updated_seurat)  # Mettez à jour l'objet Seurat renommé
    # Mettez à jour les options de sélection des clusters avec les nouveaux noms de clusters
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



  ############################## Find markers for a specific cluster ##############################

  # Variables réactives pour cet onglet
  gene_tables_10 <- reactiveValues()
  all_markers_10 <- reactiveVal()
  umap_plot <- reactiveVal()

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

  # Download handler for the UMAP plot
  output$downloadUMAPCluster <- downloadHandler(
    filename = function() {
      paste("UMAP_plot", Sys.Date(), ".tiff", sep = "")
    },
    content = function(file) {
      tryCatch({
        req(umap_plot())
        ggsave(file, plot = umap_plot(), dpi = input$dpi_umap_cluster, width = 10, height = 6)
      }, error = function(e) {
        # Log or notification of the error
        showNotification(paste("Error in UMAP Plot: ", e$message), type = "error")
      })
    }
  )


  # Fonction pour calculer tous les marqueurs différentiels pour le cluster choisi
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

  # Fonction améliorée pour calculer les marqueurs différentiels entre groupes de clusters
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



  # Download gene comparison table
  output$download_markers_single_cluster <- downloadHandler(
    filename = function() {
      paste("diff-genes-comparison-", input$cluster, "-VS-", "all_others_clusters", "-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      tryCatch({
        req(!is.null(all_markers_10))
        all_markers_10_temp <- all_markers_10()
        all_markers_10_temp_sorted <- all_markers_10_temp[order(all_markers_10_temp$p_val), ]
        write.csv(all_markers_10_temp_sorted, file)
      }, error = function(e) {
        showNotification(paste0("Error downloading gene comparison table: ", e$message), type = "error")
      })
    },
    contentType = "text/csv"
  )

  # Download gene comparison table
  output$download_markers_multiple_clusters <- downloadHandler(
    filename = function() {
      paste("diff-genes-comparison-", input$cluster1, "-VS-", input$cluster2, "-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      tryCatch({
        req(!is.null(all_markers_new))

        all_markers_new_temp <- all_markers_new()
        all_markers_new_temp_sorted <- all_markers_new_temp[order(all_markers_new_temp$p_val), ]
        write.csv(all_markers_new_temp_sorted, file)
      }, error = function(e) {
        showNotification(paste0("Error downloading gene comparison table: ", e$message), type = "error")
      })
    },
    contentType = "text/csv"
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
  output$download_subset_seurat <- downloadHandler(
    filename = function() {
      paste("seurat_subset_", Sys.Date(), ".rds", sep="")
    },
    content = function(file) {
      # Display the modal box
      showModal(modalDialog(
        title = "Please Wait",
        "Preparing the seurat object for download...",
        easyClose = FALSE,
        footer = NULL
      ))

      tryCatch({
        req(subset_seurat())
        saveRDS(subset_seurat(), file)

        # Remove the modal box after saving the Seurat subset
        removeModal()

      }, error = function(e) {
        removeModal()
        showNotification(paste0("Error while downloading Seurat subset: ", e$message), type = "error")
      })
    }
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

