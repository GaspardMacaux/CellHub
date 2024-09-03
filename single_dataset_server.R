############################## Script Single Dataset Analysis ##############################

# single_dataset_server.R

single_dataset_server <- function(input, output, session) {
  
  ############################## Loading Data ##############################
  # Tab 1 : Loading Data
  
  shinyjs::useShinyjs() # button deactivation
  
  # Reactives objects for the single dataset part
  single_dataset_object <- reactiveVal(NULL)
  gene_list <- reactiveValues(features = NULL)
  
  # Fonction pour nettoyer l'espace de travail
  cleanWorkspaceSingleDataset <- function() {
    message("Cleaning workspace...")
    
    # Nettoyage des objets réactifs
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
    
    # Suppression des dossiers si existants
    if (dir.exists("tempdir")) {
      unlink("tempdir", recursive = TRUE, force = TRUE)
      message("Deleted tempdir.")
    }
    if (dir.exists("dataDir")) {
      unlink("dataDir", recursive = TRUE, force = TRUE)
      message("Deleted dataDir.")
    }
    
    # Suppression de tous les dossiers commençant par "unzipped"
    unzipped_dirs <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
    unzipped_dirs <- unzipped_dirs[grepl("unzipped", unzipped_dirs)]
    if (length(unzipped_dirs) > 0) {
      sapply(unzipped_dirs, unlink, recursive = TRUE, force = TRUE)
      message("Deleted unzipped directories.")
    }
    
    gc()  # Garbage collection
    message("Workspace cleaned.")
  }
  
  # Fonction pour dézipper et analyser le contenu avec gestion des fichiers ZIP imbriqués
  extract_and_process_zip <- function(zip_path, exdir) {
    message(paste("Extracting zip file:", zip_path, "to", exdir))
    unzip(zip_path, exdir = exdir)
    
    # Rechercher des fichiers ZIP imbriqués
    nested_zips <- list.files(exdir, pattern = "\\.zip$", recursive = TRUE, full.names = TRUE)
    while (length(nested_zips) > 0) {
      for (nested_zip in nested_zips) {
        message(paste("Extracting nested zip file:", nested_zip))
        unzip(nested_zip, exdir = exdir)
        unlink(nested_zip)  # Supprimer le sous-fichier ZIP après extraction
      }
      # Vérifier s'il y a encore des fichiers ZIP après la première extraction
      nested_zips <- list.files(exdir, pattern = "\\.zip$", recursive = TRUE, full.names = TRUE)
    }
    
    # Déterminer s'il faut décompresser les fichiers .gz
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
  
  # Fonction pour compresser les fichiers non compressés
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
  
  # Fonction pour charger les données à partir du fichier
  processDataset <- function(dataDir, dataset_type, mt_pattern) {
    cleanWorkspaceSingleDataset()
    tryCatch({
      message("Processing dataset...")
      
      # Liste des fichiers après extraction
      all_files <- list.files(dataDir, recursive = TRUE, full.names = TRUE)
      message("Files found in data directory:", paste(all_files, collapse = ", "))
      
      # Vérifier la présence des fichiers nécessaires
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
    
    # Vérification de l'intégrité des fichiers après extraction
    required_files <- c("matrix.mtx", "features.tsv", "barcodes.tsv")
    files_found <- list.files(dataDir, recursive = TRUE, full.names = TRUE)
    
    if (!all(required_files %in% basename(files_found))) {
      stop("Necessary 10X files (matrix.mtx, features.tsv, barcodes.tsv) are not found in the specified directory.")
    }
    
    processDataset(dataDir, input$dataset_type, mt_pattern)
  }, error = function(e) {
    showNotification(paste("Error processing files: ", e$message), type = "error")
    message("Error during file processing:", e$message)
  })
})

  
  
  # Loading seurat object
  observeEvent(input$load_seurat_file, {
    cleanWorkspaceSingleDataset()
    message("Attempting to read file at: ", input$load_seurat_file$datapath)
    tryCatch({
      loaded_seurat <- readRDS(input$load_seurat_file$datapath)
      message("File successfully read.")
      single_dataset_object(loaded_seurat)
      showNotification("The Seurat object has been successfully loaded!")
      updateUIElements() # Call the function to update UI elements
      
    }, error = function(e) {
      showNotification(paste("An error occurred: ", e), type = "error")
      message("An error occurred: ", e)
    })
  })
  
  # Function to update UI elements after data loading
  updateUIElements <- function() {
    updateSelectInput(session, "ident_1", choices = unique(Idents(single_dataset_object())))
    updateCheckboxGroupInput(session, "ident_2", choices = unique(Idents(single_dataset_object())), selected = "")
    updatePickerInput(session, "gene_select", choices = sort(rownames(LayerData(single_dataset_object(), assay = "RNA", layer = 'counts'))))
    updatePickerInput(session, "gene_select_heatmap", choices = sort(rownames(LayerData(single_dataset_object(), assay = "RNA", layer = 'counts'))))
    updatePickerInput(session, "gene_select_genes_analysis", choices = sort(rownames(LayerData(single_dataset_object(), assay = "RNA", layer = 'counts'))))
  }
  
  # Observe event for gene selection
  observeEvent(single_dataset_object(), {
    gene_list <- sort(rownames(LayerData(single_dataset_object(), assay = "RNA", layer = 'counts')))
    updatePickerInput(session, "gene_select", choices = gene_list)
    updatePickerInput(session, "gene_select_heatmap", choices = gene_list)
    updatePickerInput(session, "gene_select_genes_analysis", choices = gene_list)
  })



  ############################## QC metrics and normalization ##############################
  # Tab 2 : QC metrics and normalization

  #number of nuclei:
  nuclei_count <- reactive({
    req(single_dataset_object())
    seuratRNA_subset <-subset(single_dataset_object(), subset = nFeature_RNA > input$nFeature_range[1] & nFeature_RNA < input$nFeature_range[2] & percent.mt < input$percent.mt_max)
    if (exists("seuratRNA_subset")) {
      return(dim(seuratRNA_subset)[2])
    } else {
      return(0)
    }
  })

  #text output number of nuclei
  output$nuclei_count <- renderInfoBox({
    # S'assurer que la fonction nuclei_count est exécutée et renvoie une valeur
    count <- tryCatch({
      count_value <- nuclei_count()  # Appel de la fonction pour obtenir la valeur
      infoBox(
        title = "Number of Nuclei",
        value = count_value,  # Utiliser la valeur retournée ici
        color = "blue",
        icon = icon("dna"),
        width = 1
      )
    }, error = function(e) {
      # Gérer l'erreur si nécessaire
      infoBox(
        title = "Number of Nuclei",
        value = "Error",
        icon = icon("dna"),
        color = "blue"
      )
    })
  })


  # Display QC metrics on a VlnPlot chart
  output$vlnplot <- renderPlot({
    tryCatch({
      if (input$QCmetrics) {
        VlnPlot(single_dataset_object(), features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
        # Vérification que les gènes personnalisés sont présents dans l'objet

      }
    }, error = function(e) {
      showNotification(paste0("Error displaying violin plot: ", e$message), type = "error")
    })
  })

  # Display the scatter plot 1
  output$scatter_plot1 <- renderPlot({
    tryCatch({
      if (input$show_plots) {
        seuratRNA_subset <- subset(single_dataset_object(), subset = nFeature_RNA > input$nFeature_range[1] & nFeature_RNA < input$nFeature_range[2] & percent.mt < input$percent.mt_max)
        FeatureScatter(seuratRNA_subset, feature1 = "nCount_RNA", feature2 = "percent.mt")

      }
    }, error = function(e) {
    })
  })

  # Display the scatter plot 2
  output$scatter_plot2 <- renderPlot({
    tryCatch({
      if (input$show_plots) {
        seuratRNA_subset <- subset(single_dataset_object(), subset = nFeature_RNA > input$nFeature_range[1] & nFeature_RNA < input$nFeature_range[2] & percent.mt < input$percent.mt_max)
        FeatureScatter(seuratRNA_subset, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
      }
    }, error = function(e) {
      showNotification(paste0("Error displaying scatter plot: ", e$message), type = "error")
    })
  })

  # Apllying filtering with the qc aprameter
  observeEvent(input$apply_qc, {
    tryCatch({
      # Show modal dialog
      showModal(modalDialog(
        title = "Please Wait",
        "Applying QC filters...",
        easyClose = FALSE,
        footer = NULL
      ))

      seurat_obj <- single_dataset_object()

      # Ajouter un tryCatch pour le sous-ensemble
      seurat_obj <- tryCatch({
        subset(seurat_obj, subset = nFeature_RNA > input$nFeature_range[1] &
                 nFeature_RNA < input$nFeature_range[2] &
                 percent.mt < input$percent.mt_max)
      }, error = function(e) {
        stop(paste("Error subsetting the Seurat object: ", e$message))
      })

      single_dataset_object(seurat_obj)
      num_nuclei_after_qc <- dim(seurat_obj)[2]
      msg <- paste("QC filters applied. Number of nuclei retained:", num_nuclei_after_qc)
      showNotification(msg, type = "message")

      # Close the modal dialog
      removeModal()

    }, error = function(e) {
      removeModal() # Close the modal dialog in case of an error
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

      # Vérification de la présence des gènes personnalisés
      custom_genes <- c("gene_id_EGFP_1", "gene_id_EGFP_2", "Myh1", "EGFP1")
      counts_matrix <- GetAssayData(object = normalized_seurat, slot = "counts")
      for (gene in custom_genes) {
        if (gene %in% rownames(counts_matrix)) {
          cat(gene, "is present in the data.\n")
        } else {
          cat(gene, "is NOT present in the data.\n")
        }
      }

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
      paste("UMAP_plot", Sys.Date(), ".png", sep = "")
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

  ############################## Calculate differential expressed genes for each cluster ##############################

  # Variables réactives pour cet onglet
  gene_tables <- reactiveValues()
  all_markers <- reactiveVal()

  # Function to calculate all differential markers once only
  clean_gene_names_for_html <- function(gene_names) {
    tryCatch({
      cleaned_names <- gsub("\\..*$", "", gene_names) # Enlève le suffixe après le point
      cleaned_names <- gsub("\\.", "", cleaned_names) # Enlève aussi les points restants
      return(cleaned_names)
    }, error = function(e) {
      showNotification(paste0("Error when preparing genes names", e$message), type = "error")
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
  # Observer pour mettre à jour les choix de gènes
  observeEvent(single_dataset_object(), {
    gene_list <- sort(rownames(LayerData(single_dataset_object(), assay = "RNA", layer = 'counts')))
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

  # Observer pour mettre à jour les choix d'identifiants uniques
  observe({
    req(!is.null(single_dataset_object()))
    unique_idents <- unique(Idents(single_dataset_object()))
    updateSelectInput(session, "ident_1", choices = unique_idents)
    updateCheckboxGroupInput(session, "ident_2", choices = unique_idents, selected = "")
  })

  # FeaturePlot
  feature_plot <- reactiveVal()

  observeEvent(input$show_feature, {
    tryCatch({
      req(input$gene_list_feature, single_dataset_object())
      genes <- unique(trimws(unlist(strsplit(input$gene_list_feature, ","))))

      seurat_object <- single_dataset_object()
      req(seurat_object)

      DefaultAssay(seurat_object) <- "RNA"

      present_genes <- genes[genes %in% rownames(LayerData(seurat_object, assay = "RNA", layer = "counts"))]
      missing_genes <- genes[!(genes %in% rownames(LayerData(seurat_object, assay = "RNA", layer = "counts")))]

      if (length(missing_genes) > 0) {
        showNotification(paste("The following genes are not present in the dataset:", paste(missing_genes, collapse = ", ")), type = "warning")
      }

      if (length(present_genes) > 0) {
        min_cutoff <- ifelse(input$min_cutoff == "None", NA, input$min_cutoff)
        max_cutoff <- ifelse(input$max_cutoff == "None", NA, input$max_cutoff)

        plot <- FeaturePlot(seurat_object, features = present_genes, min.cutoff = min_cutoff, max.cutoff = max_cutoff)

        if (input$show_coexpression && length(present_genes) > 1) {
          plot <- FeaturePlot(seurat_object, features = present_genes, blend = TRUE, blend.threshold = 1, order = TRUE, min.cutoff = min_cutoff, max.cutoff = max_cutoff)
        }

        if (input$add_noaxes_feature) { plot <- plot + NoAxes() }
        if (input$add_nolegend_feature) { plot <- plot + NoLegend() }

        plot <- plot + theme(
          text = element_text(size = input$text_size, face = if (input$bold_text) "bold" else "plain"),
          axis.text.x = element_text(size = input$text_size, face = if (input$bold_text) "bold" else "plain"),
          axis.text.y = element_text(size = input$text_size, face = if (input$bold_text) "bold" else "plain"),
          axis.title.x = element_text(size = input$text_size, face = if (input$bold_text) "bold" else "plain"),
          axis.title.y = element_text(size = input$text_size, face = if (input$bold_text) "bold" else "plain")
        )

        if (input$hide_x_labels) {
          plot <- plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
        }
        if (input$hide_y_labels) {
          plot <- plot + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
        }
        if (input$hide_x_title) {
          plot <- plot + xlab("")
        }
        if (input$hide_y_title) {
          plot <- plot + ylab("")
        }

        feature_plot(plot)
      } else {
        showNotification("None of the requested genes are present in the dataset.", type = "error")
      }
    }, error = function(e) {
      showNotification(paste("Error in FeaturePlot: ", e$message), type = "error")
    })
  })

  output$feature_plot <- renderPlot({
    req(feature_plot())
    feature_plot()
  })

  output$downloadFeaturePlot <- downloadHandler(
    filename = function() {
      paste("FeaturePlot", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = feature_plot(), dpi = input$dpi_plot, width = 10, height = 8)
    }
  )

  # VlnPlot
  vln_plot <- reactiveVal()

  observeEvent(input$show_vln, {
    tryCatch({
      req(input$gene_list_vln, single_dataset_object())
      genes <- unique(trimws(unlist(strsplit(input$gene_list_vln, ","))))

      seurat_object <- single_dataset_object()
      req(seurat_object)

      DefaultAssay(seurat_object) <- "RNA"

      plot <- VlnPlot(seurat_object, features = genes, log = TRUE, pt.size = ifelse(input$hide_vln_points, 0, 1))

      plot <- plot + theme(
        text = element_text(size = input$text_size, face = if (input$bold_text) "bold" else "plain"),
        axis.text.x = element_text(size = input$text_size, face = if (input$bold_text) "bold" else "plain"),
        axis.text.y = element_text(size = input$text_size, face = if (input$bold_text) "bold" else "plain"),
        axis.title.x = element_text(size = input$text_size, face = if (input$bold_text) "bold" else "plain"),
        axis.title.y = element_text(size = input$text_size, face = if (input$bold_text) "bold" else "plain")
      )

      if (input$hide_x_labels) {
        plot <- plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
      }
      if (input$hide_y_labels) {
        plot <- plot + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
      }
      if (input$hide_x_title) {
        plot <- plot + xlab("")
      }
      if (input$hide_y_title) {
        plot <- plot + ylab("")
      }

      if (input$add_nolegend_feature) {
        plot <- plot + NoLegend()
      }

      vln_plot(plot)
    }, error = function(e) {
      showNotification(paste("Error in VlnPlot: ", e$message), type = "error")
    })
  })

  output$vln_plot <- renderPlot({
    req(vln_plot())
    vln_plot()
  })

  output$downloadVlnPlot <- downloadHandler(
    filename = function() {
      paste("VlnPlot", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      ggsave(file, plot = vln_plot(), dpi = input$dpi_plot, width = 10, height = 8)
    }
  )

  # DotPlot
  dot_plot <- reactiveVal()

  observeEvent(input$show_dot, {
    tryCatch({
      req(input$gene_list_dotplot, single_dataset_object())
      genes <- unique(trimws(unlist(strsplit(input$gene_list_dotplot, ","))))

      seurat_object <- single_dataset_object()
      req(seurat_object)

      DefaultAssay(seurat_object) <- "RNA"

      plot <- DotPlot(seurat_object, features = genes) + RotatedAxis()

      plot <- plot + theme(
        text = element_text(size = input$text_size, face = if (input$bold_text) "bold" else "plain"),
        axis.text.x = element_text(size = input$text_size, face = if (input$bold_text) "bold" else "plain"),
        axis.text.y = element_text(size = input$text_size, face = if (input$bold_text) "bold" else "plain"),
        axis.title.x = element_text(size = input$text_size, face = if (input$bold_text) "bold" else "plain"),
        axis.title.y = element_text(size = input$text_size, face = if (input$bold_text) "bold" else "plain")
      )

      if (input$hide_x_labels) {
        plot <- plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
      }
      if (input$hide_y_labels) {
        plot <- plot + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
      }
      if (input$hide_x_title) {
        plot <- plot + xlab("")
      }
      if (input$hide_y_title) {
        plot <- plot + ylab("")
      }

      if (input$add_nolegend_feature) {
        plot <- plot + NoLegend()
      }

      dot_plot(plot)
    }, error = function(e) {
      showNotification(paste("Error in DotPlot: ", e$message), type = "error")
    })
  })

  output$dot_plot <- renderPlot({
    req(dot_plot())
    dot_plot()
  })

  output$downloadDotPlot <- downloadHandler(
    filename = function() {
      paste("DotPlot", Sys.Date(), ".png", sep = "")
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

      DefaultAssay(seurat_object) <- "RNA"

      plot <- RidgePlot(seurat_object, features = genes)

      plot <- plot + theme(
        text = element_text(size = input$text_size, face = if (input$bold_text) "bold" else "plain"),
        axis.text.x = element_text(size = input$text_size, face = if (input$bold_text) "bold" else "plain"),
        axis.text.y = element_text(size = input$text_size, face = if (input$bold_text) "bold" else "plain"),
        axis.title.x = element_text(size = input$text_size, face = if (input$bold_text) "bold" else "plain"),
        axis.title.y = element_text(size = input$text_size, face = if (input$bold_text) "bold" else "plain")
      )

      if (input$hide_x_labels) {
        plot <- plot + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
      }
      if (input$hide_y_labels) {
        plot <- plot + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
      }
      if (input$hide_x_title) {
        plot <- plot + xlab("")
      }
      if (input$hide_y_title) {
        plot <- plot + ylab("")
      }

      if (input$add_nolegend_feature) {
        plot <- plot + NoLegend()
      }

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
      paste("RidgePlot", Sys.Date(), ".png", sep = "")
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

      gene_data <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
      available_genes <- rownames(gene_data)
      selected_genes <- intersect(input$gene_select_genes_analysis, available_genes)
      if (length(selected_genes) == 0) {
        showNotification("No selected genes are present in the dataset.", type = "error")
        return()
      }

      data <- FetchData(seurat_obj, vars = c("ident", selected_genes))

      expression_summary_list <- lapply(selected_genes, function(gene) {
        gene_data <- data[[gene]]

        expressed_indices <- gene_data > input$logfc_threshold
        cells_expressed <- sum(expressed_indices, na.rm = TRUE)
        total_cells <- length(gene_data)

        if (cells_expressed > 0) {
          cluster_info <- data$ident[expressed_indices]
          cluster_counts <- table(cluster_info)
          pct_per_cluster <- prop.table(cluster_counts) * 100

          df <- data.frame(
            Gene = gene,
            Cluster = names(cluster_counts),
            Cells_Expressed = as.numeric(cluster_counts),
            Percentage_of_Cells = as.numeric(pct_per_cluster),
            Total_Cells = total_cells
          )
          return(df)
        } else {
          return(NULL)
        }
      })

      valid_expression_summary <- Filter(Negate(is.null), expression_summary_list)
      if (length(valid_expression_summary) > 0) {
        expression_df <- do.call(rbind, valid_expression_summary)
        number_of_nuclei(expression_df)
      } else {
        expression_df <- data.frame()
        showNotification("No results meet the thresholds specified.", type = "warning")
      }

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
    req(single_dataset_object(), input$gene_select_heatmap)

    selected_genes <- input$gene_select_heatmap
    selected_genes <- trimws(selected_genes)
    valid_genes <- selected_genes %in% rownames(single_dataset_object()[["RNA"]])

    if (sum(valid_genes) == 0) {
      showNotification("none of the genes selected is valid.", type = "error")
      return()
    }

    seurat_object <- single_dataset_object()

    tryCatch({
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
        } else {
          selected_clusters <- input$cluster_selector
          Idents(seurat_object) <- factor(Idents(seurat_object), levels = selected_clusters)
          seurat_object <- subset(seurat_object, idents = selected_clusters)
        }
      }

      # Générer la heatmap
      heatmap_plot(
        DoHeatmap(seurat_object, features = selected_genes[valid_genes]) + NoLegend()
      )
    }, error = function(e) {
      showNotification("Genes are not found", type = "error")
    })
  })


  # Rendre le plot stocké dans l'interface utilisateur
  output$heatmap_single <- renderPlot({
    req(heatmap_plot())
    heatmap_plot()
  })

  output$download_heatmap_single <- downloadHandler(
    filename = function() {
      paste("heatmap-", Sys.Date(), ".png", sep="")
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
      paste0("FeatureScatter_", Sys.Date(), ".png")
    },
    content = function(file) {
      tryCatch({
        ggsave(file, plot = scatter_plot(), device = "png", width = 10, height = 8, dpi = input$dpi_scatter)
      }, error = function(e) {
        showNotification(paste("Erreur lors du téléchargement :", e$message), type = "error")
      })
    }
  )



############################## Genes analysis ##############################

 number_of_nuclei <- reactiveVal(NULL)
  observeEvent(input$analyze_btn, {
    tryCatch({
      req(input$gene_select_genes_analysis, single_dataset_object())

      # Récupération de l'objet Seurat
      seurat_obj <- single_dataset_object()

      # Vérification de la présence des gènes dans l'objet Seurat
      gene_data <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
      available_genes <- rownames(gene_data)
      if (!all(input$gene_select_genes_analysis %in% available_genes)) {
        showNotification("One or more selected genes are not present in the dataset.", type = "error")
        return()
      }

      # Extraction des données y compris les identités des clusters
      data <- FetchData(seurat_obj, vars = c("ident", input$gene_select_genes_analysis))

      # Calcul des métriques pour chaque gène
      expression_summary_list <- lapply(input$gene_select_genes_analysis, function(gene) {
        gene_data <- data[[gene]]
        # Appliquer le seuil pour considérer un gène comme exprimé
        expressed_indices <- gene_data > input$logfc_threshold
        cells_expressed <- sum(expressed_indices)
        total_cells <- length(gene_data)

        if (cells_expressed > 0) {
          cluster_info <- data$ident[expressed_indices]
          cluster_counts <- table(cluster_info)
          pct_per_cluster <- prop.table(cluster_counts) * 100

          df <- data.frame(
            Gene = gene,
            Cluster = names(cluster_counts),
            Cells_Expressed = as.numeric(cluster_counts),
            Percentage_of_Cells = as.numeric(pct_per_cluster),
            Total_Cells = total_cells
          )
          return(df)
        } else {
          return(NULL)
        }
      })

      # Filtrer les résultats vides et combiner tous les dataframes
      valid_expression_summary <- Filter(Negate(is.null), expression_summary_list)
      if (length(valid_expression_summary) > 0) {
        expression_df <- do.call(rbind, valid_expression_summary)
      } else {
        expression_df <- data.frame()
        showNotification("No results meet the thresholds specified.", type = "warning")
      }

      # Affichage du tableau de résultats
      output$expression_summary <- renderDataTable({
        datatable(expression_df, options = list(pageLength = 10, scrollX = TRUE))
      })

    }, error = function(e) {
      showNotification("Error processing expression data: " %||% e$message, type = "error")
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


  observe({
    if (!is.null(single_dataset_object())) {
      message("Update cluster colors")
      current_colours <- scales::hue_pal()(length(unique(Idents(single_dataset_object()))))
      names(current_colours) <- sort(unique(Idents(single_dataset_object())))
      cluster_colours(current_colours)
    }
  })

  observeEvent(input$update_colour, {
    message("Update the color of the selected cluster")
    current_colours <- cluster_colours()
    current_colours[input$cluster_select] <- input$cluster_colour
    cluster_colours(current_colours)
  })

  output$umap_finale <- renderPlotly({
    req(single_dataset_object())
    message("Generating the final UMAP")


    plot_data <- DimPlot(single_dataset_object(), group.by = "ident", pt.size = input$pt_size, label = TRUE, label.size = input$label_font_size) +
      theme(axis.line = element_line(size = 0.5)) +
      scale_color_manual(values = cluster_colours())

    plot_data <- plot_data + NoLegend()
    interactive_plot <- ggplotly(plot_data, tooltip = "text")

    interactive_plot <- interactive_plot %>%
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

  # Ajout de l'UMAP en haut de l'onglet
  output$umap_cluster_10 <- renderPlot({
    req(single_dataset_object())
    message("UMAP generation for a specific cluster")
    plot_data <- DimPlot(single_dataset_object(), group.by = "ident", pt.size = input$pt_size, label = TRUE, label.size = input$label_font_size) +
      theme(axis.line = element_line(size = 0.5)) +
      scale_color_manual(values = cluster_colours())
    plot_data <- plot_data + NoLegend()+ ggtitle(NULL)
    message("UMAP for a specific cluster generated")
    return(plot_data)
  })

  # Mettre à jour les choix pour le selectInput cluster
  observe({
    req(single_dataset_object())
    updateSelectInput(session, "cluster", choices = levels(single_dataset_object()))
  })

  # Variable réactive pour stocker le graphique UMAP
  umap_plot <- reactive({
    req(single_dataset_object())
    DimPlot(single_dataset_object(), group.by = "ident", pt.size = input$pt_size, label = TRUE, label.size = input$label_font_size) +
      theme(axis.line = element_line(size = 0.5)) +
      scale_color_manual(values = cluster_colours()) +
      NoLegend() + ggtitle(NULL)
  })

  output$umap_plot <- renderPlot({
    req(umap_plot())
    umap_plot()
  })

  output$downloadUMAPCluster <- downloadHandler(
    filename = function() {
      paste("UMAP_plot", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      tryCatch({
        req(umap_plot())
        ggsave(file, plot = umap_plot(), dpi = input$dpi_umap_cluster, width = 10, height = 6)
      }, error = function(e) {
        # Log ou notification de l'erreur
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

  # Cluster selection
  observeEvent(input$cluster1, {
    updateSelectInput(session, "cluster2", choices = setdiff(levels(single_dataset_object()), input$cluster1))
  })

  # Variables réactives pour cet onglet
  gene_tables_new <- reactiveValues()
  all_markers_new <- reactiveVal()

  # Fonction pour calculer tous les marqueurs différentiels pour la comparaison de clusters
  calculate_markers_for_comparison <- function() {
    # Vérifier que les objets et les entrées nécessaires sont disponibles
    req(single_dataset_object(), input$cluster1, input$cluster2, input$min_pct_comparison, input$logfc_threshold_comparison)

    tryCatch({
      # Vérifier que les clusters choisis sont différents
      if (input$cluster1 == input$cluster2) {
        showNotification("Les clusters sélectionnés ne doivent pas être identiques.", type = "error")
        return(FALSE)
      }

      # Calcul des marqueurs différentiels
      markers <- FindMarkers(single_dataset_object(),
                             ident.1 = input$cluster1,
                             ident.2 = input$cluster2,
                             min.pct = input$min_pct_comparison,
                             logfc.threshold = input$logfc_threshold_comparison)

      # Formatage des valeurs p et ajustement des valeurs p pour l'affichage
      markers$p_val <- format(markers$p_val, scientific = TRUE, digits = 3)
      markers$p_val_adj <- format(markers$p_val_adj, scientific = TRUE, digits = 3)
      cleaned_gene_names <- clean_gene_names_for_html(rownames(markers))
      markers$gene <- paste0('<a href="#" class="gene-name" data-gene="', cleaned_gene_names, '">',cleaned_gene_names, '</a>')

      # Stocker les résultats dans une valeur réactive
      all_markers_new(markers)

      # Retourner TRUE en cas de succès
      return(TRUE)
    }, error = function(e) {
      # Gestion des erreurs, par exemple, affichage d'une notification dans l'application
      showNotification(paste0("Error in calculate_markers_for_comparison: ", e$message), type = "error")
      # Retourner FALSE en cas d'erreur
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
      # Notification de l'erreur
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

