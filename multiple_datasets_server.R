
############################## Script Multiple Datasets Analysis ##############################

# multiple_datasets_server.R

    multiple_datasets_server <-  function(input, output, session) {

      ############################## Loading Data ##############################
      
      # Variables for the merge dataset part
      multiple_datasets_object <- reactiveVal()
      seurat_objects <- reactiveValues()
      data_loaded <- reactiveValues()
      merged_gene_tables <- reactiveValues()
      rv_metadata <- reactiveValues(num_fields = 1)
      
      shinyjs::disable("add_field")
      shinyjs::disable("add_metadata")
      
      cleanWorkspaceMultipleDatasets <- function() {
        # Réinitialisation des objets réactifs
        multiple_datasets_object(NULL)
        seurat_objects <<- list()
        data_loaded <<- list()
        merged_gene_tables <<- list()
        rv_metadata$num_fields <- 1
        
        # Nettoyage des plots réactifs
        clustering_plot_merge(NULL)
        feature_plot_merge(NULL)
        vln_plot_merge(NULL)
        dot_plot_merge(NULL)
        ridge_plot_merge(NULL)
        heatmap_plot_multidataset(NULL)
        
        # Suppression des répertoires temporaires s'ils existent
        if (dir.exists("unzipped")) {
          unlink("unzipped", recursive = TRUE, force = TRUE)
        }
        if (dir.exists("tempData")) {
          unlink("tempData", recursive = TRUE, force = TRUE)
        }
        
        # Suppression de tous les dossiers commençant par "unzipped"
        unzipped_dirs <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
        unzipped_dirs <- unzipped_dirs[grepl("unzipped", unzipped_dirs)]
        if (length(unzipped_dirs) > 0) {
          sapply(unzipped_dirs, unlink, recursive = TRUE, force = TRUE)
        }
        rm(list=ls())
        print(list.files(recursive = TRUE))
      }
      
      # Load a seurat object
      observeEvent(input$load_seurat_file_merge, {
        cleanWorkspaceMultipleDatasets()
        message("Attempting to read file at: ", input$load_seurat_file_merge$datapath)
        tryCatch({
          loaded_seurat <- readRDS(input$load_seurat_file_merge$datapath)
          message("File successfully read.")
          loaded_seurat <- JoinLayers(loaded_seurat, assay = "RNA")
          
          # Initialiser la colonne "ClusterIdents" avec les identifiants actuels des clusters si elle n'existe pas
          if (!"ClusterIdents" %in% colnames(loaded_seurat@meta.data)) {
            loaded_seurat$ClusterIdents <- Idents(loaded_seurat)
          }
          
          multiple_datasets_object(loaded_seurat)
          showNotification("The Seurat object has been successfully loaded!")
          shinyjs::enable("add_field")
          shinyjs::enable("add_metadata")
          shinyjs::enable("runScalePCA")
          
          # Mise à jour des choix du sélecteur 'group_by_select' pour inclure 'ClusterIdents'
          update_group_by_choices()
        }, error = function(e) {
          showNotification(paste("An error occurred: ", e), type = "error")
        })
      })
      
      # Processing of loaded data
      observe_file_input <- function(index) {
        observeEvent(input[[paste0("merge", index)]], {
          if (is.null(input[[paste0("merge", index)]])) {
            return(NULL)
          }
          
          if (!is.null(data_loaded[[paste0("loaded", index)]]) && data_loaded[[paste0("loaded", index)]] == TRUE) {
            return(NULL)
          }
          
          tryCatch({
            # Show modal dialog
            showModal(modalDialog(
              title = "Please Wait",
              "Processing data...",
              easyClose = FALSE,
              footer = NULL
            ))
            
            seurat_object <- NULL
            mt_pattern_merge <- ifelse(input$species_choice_merge == "mouse", "^mt-", "^MT-")
            dataset_type <- input[[paste0("dataset_type_merge", index)]]
            if (is.null(dataset_type)) {
              return(NULL)
            }
            if(dataset_type == "seurat_object_merge") {
              seurat_object <- readRDS(input[[paste0("merge", index)]]$datapath)
              message("Seurat object loading complete")
            }
            else {
              if (dir.exists(paste0("unzipped", index))) {
                unlink(paste0("unzipped", index), recursive = TRUE)
              }
              unzip(input[[paste0("merge", index)]]$datapath, exdir = paste0("unzipped", index))
              message("Decompression complete")
              data <- Read10X(paste0("unzipped", index))
              message("Data read-out complete")
              if (dataset_type == "snRNA_merge") {
                seurat_object <- CreateSeuratObject(counts = data, project = input[[paste0("dataset_name", index)]], min.cells = 3, min.features = 200)
                seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = mt_pattern_merge)
                message("Seurat object creation for snRNA-seq completed")
              }
              else if (dataset_type == "multiome_merge") {
                rna.data <- data$`Gene Expression`
                atac.data <- data$`Peaks`
                rna.seurat <- CreateSeuratObject(counts = rna.data, project = input[[paste0("dataset_name", index)]], min.cells = 3, min.features = 200)
                atac.seurat <- CreateSeuratObject(counts = atac.data, project = "ATAC", min.cells = 3, min.features = 200)
                rna.seurat[["percent.mt"]] <- PercentageFeatureSet(rna.seurat, pattern = mt_pattern_merge)
                seurat_object <- rna.seurat
                message("Seurat object creation for Multiome completed")
              }
              seurat_object <- NormalizeData(seurat_object)
              seurat_object <- FindVariableFeatures(seurat_object, nfeatures = 4000)
              seurat_object <- ScaleData(seurat_object)
              seurat_object <- RunPCA(seurat_object)
              message("Seurat object processing complete")
            }
            seurat_object@meta.data$dataset <- input[[paste0("dataset_name", index)]]
            seurat_object$orig.ident <- seurat_object@meta.data$dataset
            seurat_objects[[paste0("seurat_object", index)]] <- seurat_object
            data_loaded[[paste0("loaded", index)]] <- TRUE
            showNotification("Data processed successfully!", type = "message")
            
            # Close the modal dialog
            removeModal()
            
            # Mise à jour des éléments de l'interface utilisateur
            updateUIElements()
          }, error = function(e) {
            removeModal() # Close the modal dialog in case of an error
          })
        })
      }
      
      output$fileInputs <- renderUI({
        lapply(1:input$num_datasets, function(i) {
          fluidRow(
            column(3, textInput(paste0('dataset_name', i), paste0('Dataset name', i), value = paste0("Dataset ", i))),
            column(3, selectInput(paste0("dataset_type_merge", i), "Data type:",
                                  choices = list("snRNA-seq" = "snRNA_merge",
                                                 "Multiome" = "multiome_merge",
                                                 "Seurat Object" = "seurat_object_merge"))),
            column(6, fileInput(paste0('merge', i), paste0('Choose a Seurat or .gz file', i), accept=c('rds', 'gz')))
          )
        })
      })
      
      # Observe the change in the number of datasets and create the necessary entries
      observeEvent(input$num_datasets, {
        lapply(1:input$num_datasets, function(i) {
          observe_file_input(i)
        })
      })
      
      # Count the number of dataset loaded
      observe({
        if(is.numeric(input$num_datasets) && !is.null(input$num_datasets)) {
          num_loaded_datasets <- sum(sapply(1:input$num_datasets, function(i) {
            !is.null(data_loaded[[paste0("loaded", i)]]) && data_loaded[[paste0("loaded", i)]] == TRUE
          }))
          if (num_loaded_datasets == input$num_datasets) {
            shinyjs::enable("integrate")
          } else {
            shinyjs::disable("integrate")
          }
        } else {
          shinyjs::disable("integrate")
          warning("input$num_datasets is not valid: ", input$num_datasets)
        }
      })
      
      # Integration function
      integrate_data <- function(seurat_list) {
        print("Starting integration")
        print("Selecting integration features")
        features <- SelectIntegrationFeatures(
          object.list = seurat_list,
          nfeatures = 2000
        )
        print("Finding Integration Anchors")
        anchors <- FindIntegrationAnchors(
          object.list = seurat_list,
          dims = 1:30
        )
        print("Integrating Data")
        seurat_integrated_temp <- IntegrateData(
          anchorset = anchors,
          dims = 1:30,
          features.to.integrate = features
        )
        DefaultAssay(seurat_integrated_temp) <- "integrated"
        new_metadata <- seurat_integrated_temp@meta.data
        colnames(new_metadata)[colnames(new_metadata) == "orig.ident"] <- "dataset"
        seurat_integrated_temp <- AddMetaData(seurat_integrated_temp, metadata = new_metadata)
        
        print("Finished integration")
        print(paste("Integrated object class: ", class(seurat_integrated_temp)))
        print(paste("Integrated object dimensions: ", dim(seurat_integrated_temp)))
        
        return(seurat_integrated_temp)
      }
      
      # Integration button
      observeEvent(input$integrate, {
        tryCatch({
          # Show modal dialog
          showModal(modalDialog(
            title = "Please Wait",
            "Integrating datasets...",
            easyClose = FALSE,
            footer = NULL
          ))
          
          print("Integrate button pressed")
          seurat_list <- list()
          for (i in 1:input$num_datasets) {
            seurat_list[[i]] <- seurat_objects[[paste0("seurat_object", i)]]
          }
          print("Copied Seurat objects from reactiveValues to a standard list")
          print(paste("seurat_list length:", length(seurat_list)))
          
          if (length(seurat_list) < 2) {
            stop("Please upload at least two datasets for integration")
          }
          
          integrated_object <- integrate_data(seurat_list)
          multiple_datasets_object(integrated_object)  # Update the integrated object
          print("Finished data integration")
          
          # Check if the object has been updated
          if (!is.null(multiple_datasets_object())) {
            print("The integrated object is updated in multiple_datasets_object.")
            shinyjs::enable("add_field")
            shinyjs::enable("add_metadata")
          } else {
            print("Failed to update the integrated object in multiple_datasets_object.")
          }
          
          # Close the modal dialog
          removeModal()
          
          # Mise à jour des éléments de l'interface utilisateur
          updateUIElements()
        }, error = function(e) {
          removeModal() # Close the modal dialog in case of an error
          showNotification(paste0("Error during integration: ", e$message), type = "error")
        })
      })
      
      # Add metadata field
      observeEvent(input$add_field, {
        rv_metadata$num_fields <- rv_metadata$num_fields + 1
      })
      
      # Create a responsive list to store metadata field names
      reactive_metadata_fields <- reactive({
        req(multiple_datasets_object())
        all_metadata_fields <- colnames(multiple_datasets_object()@meta.data)
        
        # Define a pattern to exclude fields that start with "RNA_snn"
        exclude_pattern <- "^RNA_snn"
        
        # Use grepl to filter out columns that match the pattern
        exclude_fields <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt", "integrated_snn_res.0.5")
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
              textInput(paste0("metadata_value_", dataset_name, "_", j), paste0("Value for ", dataset_name, " (", input[[paste0("dataset_name", which(datasets == dataset_name))]], ") ", j), value = "")
            }))
          )
        })
      })
      
      # Add metadata by the user
      observeEvent(input$add_metadata, {
        tryCatch({
          req(multiple_datasets_object()) # Assurez-vous que l'objet Seurat est chargé
          datasets <- unique(multiple_datasets_object()@meta.data$dataset)
          req(datasets) # Assurez-vous que les datasets sont présents
          
          # Créer une copie temporaire pour la modification
          seurat_temp <- multiple_datasets_object()
          
          for (j in 1:rv_metadata$num_fields) {
            metadata_field_name <- input[[paste0("metadata_name_", j)]]
            req(metadata_field_name) # Assurez-vous que le nom du champ de métadonnées n'est pas nul
            
            if (!metadata_field_name %in% colnames(seurat_temp@meta.data)) {
              seurat_temp@meta.data[[metadata_field_name]] <- NA
            }
            
            for (dataset in datasets) {
              metadata_field_value <- input[[paste0("metadata_value_", dataset, "_", j)]]
              req(metadata_field_value) # Assurez-vous que la valeur du champ de métadonnées n'est pas nulle
              
              rows <- which(seurat_temp@meta.data$dataset == dataset)
              seurat_temp@meta.data[rows, metadata_field_name] <- metadata_field_value
            }
          }
          
          # Mettre à jour l'objet Seurat intégré avec les nouvelles métadonnées
          multiple_datasets_object(seurat_temp)
          
          # Après mise à jour
          print(head(multiple_datasets_object()@meta.data))
        }, error = function(e) {
          showNotification(paste0("Error adding metadata: ", e$message), type = "error")
        })
      })
      
      # Identify and extract the first loaded Seurat object to obtain the RNA assay gene list
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


    # Observer forfind Clusters and display UMAP
    clustering_plot_merge <- reactiveVal()


    # Neighbors calculation
    observeEvent(input$runFindNeighbors, {
      tryCatch({
        # Show modal dialog
        showModal(modalDialog(
          title = "Please Wait",
          "Finding neighbors and running UMAP...",
          easyClose = FALSE,
          footer = NULL
        ))

        showNotification("Finding neighbors started...", type = "message")
        req(multiple_datasets_object())

        seurat_object_temp <- multiple_datasets_object()
        seurat_object_temp <- FindNeighbors(seurat_object_temp, dims = 1:input$dimension_2)
        seurat_object_temp <- RunUMAP(seurat_object_temp, dims = 1:input$dimension_2)
        multiple_datasets_object(seurat_object_temp)
        req(multiple_datasets_object()[["umap"]])

        clustering_plot_merge(DimPlot(multiple_datasets_object(), group.by="orig.ident") + ggtitle(NULL))

        showNotification("Finding neighbors and UMAP completed.", type = "message")

        # Close the modal dialog
        removeModal()
      }, error = function(e) {
        removeModal() # Close the modal dialog in case of an error
        showNotification(paste0("Error during find_neighbors: ", e$message), type = "error")
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
        
        seurat_integrated_temp <- FindClusters(multiple_datasets_object(), resolution = input$resolution_step1, algorithm = as.integer(input$algorithm_select))
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


    # Handler pour télécharger l'UMAP
    output$downloadUMAP_merge <- downloadHandler(
      filename = function() {
        paste("UMAP_plot", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
        tryCatch({
            req(clustering_plot_merge())
          ggsave(file, plot = clustering_plot_merge(), dpi = input$dpi_umap_merge, width = 10, height = 6)
        }, error = function(e) {
          showNotification(paste("Error saving UMAP plot: ", e$message), type = "error")
        })
      }
    )



############################## Calculate differential expressed genes for each cluster in the merged dataset ##############################


    gene_tables_merge <- reactiveValues()
    all_markers_merge <- reactiveVal()

    shinyjs::disable("download_DE_merged")

    # Function to calculate all differential markers once only
    clean_gene_names_for_html <- function(gene_names) {
      # Utiliser gsub pour remplacer tout ce qui suit un point et le point par rien
      cleaned_names <- gsub("\\..*$", "", gene_names) # Enlève le suffixe après le point
      cleaned_names <- gsub("\\.", "", cleaned_names) # Enlève aussi les points restants
      return(cleaned_names)
    }


    # Function to calculate all differential markers for the merged dataset
    calculate_merged_markers <- function() {
      tryCatch({
        req(multiple_datasets_object())
        markers <- FindAllMarkers(multiple_datasets_object(), min.pct = input$min_pct_all_multiple,
                                  logfc.threshold = input$logfc_threshold_all_multiple,)
        markers <- as.data.frame(markers)

        # Nettoyer les noms de gènes avant de générer le HTML
        cleaned_gene_names <- clean_gene_names_for_html(rownames(markers))
        markers$gene <- paste0('<span class="gene-name">', cleaned_gene_names, '</span>')
        markers$p_val <- format(markers$p_val, scientific = TRUE, digits = 3)
        markers$p_val_adj <- format(markers$p_val_adj, scientific = TRUE, digits = 3)
        markers$gene <- paste0('<a href="#" class="gene-name" data-gene="',cleaned_gene_names, '">', cleaned_gene_names, '</a>')
        all_markers_merge(markers)
      }, error = function(e) {
        showNotification(paste0("Error calculating markers:", e$message), type = "error")
        NULL
      })
    }

    update_gene_tables_display_merge <- function() {
      tryCatch({
        markers <- all_markers_merge()
        num_genes_to_display <- input$number_genes_merge
        for (cluster in unique(markers$cluster)) {
          gene_tables_merge[[paste0("table_", cluster)]] <- head(markers[markers$cluster == cluster, ], n = num_genes_to_display)
        }

        output$diff_genes_tables_merge <- renderUI({
          tagList(
            lapply(names(gene_tables_merge), function(name) {
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
        showNotification(paste0("Error updating gene table display:", e$message), type = "error")
      })
    }

    # Observer for run_DE_merged button
    observeEvent(input$run_DE_merged, {
      tryCatch({
        # Show modal dialog
        showModal(modalDialog(
          title = "Please Wait",
          "Calculating merged differential expression markers...",
          easyClose = FALSE,
          footer = NULL
        ))

        calculate_merged_markers()
        update_gene_tables_display_merge()
        shinyjs::enable("download_DE_merged")

        # Close the modal dialog
        removeModal()
      }, error = function(e) {
        removeModal() # Close the modal dialog in case of an error
        showNotification(paste0("Error running run_DE_merged:", e$message), type = "error")
      })
    })

    # Observer for updating the display only when the number of genes changes
    observeEvent(input$number_genes_merge, {
      tryCatch({
        if (!is.null(input$number_genes_merge) && !is.null(all_markers_merge())) {
          update_gene_tables_display_merge()
        }
      }, error = function(e) {
        showNotification(paste0("Error during gene update:", e$message), type = "error")
      })
    })


    # Observer for gene tables output
    observe({
      tryCatch({
        lapply(names(gene_tables_merge), function(name) {
          output[[name]] <- renderDT({
            datatable(gene_tables_merge[[name]], escape = FALSE)
          })
        })
      }, error = function(e) {
        showNotification(paste0("Error rendering gene tables:", e$message), type = "error")
      })
    })

    output$previous_tab_notification <- renderUI({
      tryCatch({
        if (!is.null(session$userData$previous_tab_notification_msg)) {
          list(
            tags$hr(),
            tags$strong(session$userData$previous_tab_notification_msg),
            tags$hr()
          )
        }
      }, error = function(e) {
        showNotification(paste0("Error when creating the notification:", e$message), type = "error")
      })
    })


    output$save_seurat_merge <- downloadHandler(
      filename = function() {
        return("seurat_object.rds")
      },
      content = function(file) {
        showModal(modalDialog(
          title = "Please Wait",
          "Preparing the seurat object for download...",
          easyClose = FALSE,
          footer = NULL
        ))
        tryCatch({
          saveRDS(multiple_datasets_object(), file)
          removeModal()
        }, error = function(e) {
          removeModal()
          showNotification(paste0("Error saving Seurat object: ", e$message), type = "error")
        })
      }
    )


    # Download markers as a CSV
    output$download_DE_merged <- downloadHandler(
      filename = function() {
        paste("DE_genes_merged_", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        tryCatch({
          markers <- all_markers_merge()
          write.csv(markers, file)
        }, error = function(e) {
          showNotification(paste0("Error downloading markers:", e$message), type = "error")
        })
      }
    )

############################## Visualize genes expressions ##############################

    # Fonction pour mettre à jour les choix du sélecteur 'group_by_select' en excluant certains champs et en ajoutant des noms personnalisés
    update_group_by_choices <- function() {
      req(multiple_datasets_object())
      
      metadata_fields <- colnames(multiple_datasets_object()@meta.data)
      
      # Exclure les champs spécifiques
      excluded_fields <- c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt")
      metadata_fields <- metadata_fields[!metadata_fields %in% excluded_fields]
      
      # Exclure les champs commençant par "integrated_snn"
      metadata_fields <- metadata_fields[!grepl("^integrated_snn", metadata_fields)]
      
      # Ajouter un nom personnalisé pour "ClusterIdents"
      combined_choices <- c(metadata_fields, "ClusterIdents")
      combined_display_names <- c(metadata_fields, "Cluster renamed")
      
      # Vérifier que les longueurs des noms et des valeurs correspondent
      if (length(combined_choices) != length(combined_display_names)) {
        stop("La longueur des noms et des valeurs ne correspond pas.")
      }
      
      names(combined_choices) <- combined_display_names
      
      print("Combined choices for group by select:")
      print(combined_choices)
      
      updateSelectInput(session, "group_by_select", choices = combined_choices, selected = "ClusterIdents")
    }
    
    # Call the function to update the group by choices
    observe({
      update_group_by_choices()
    })
    
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

    # Feature Plot
    observeEvent(input$runFeaturePlot, {
      tryCatch({
        req(input$gene_list_feature_merge)
        gene_list <- unique(trimws(strsplit(input$gene_list_feature_merge, ",")[[1]]))
        
        seurat_object <- multiple_datasets_object()
        req(seurat_object)
        
        selected_group_by <- input$group_by_select
        if (all(gene_list %in% rownames(LayerData(seurat_object, assay = "RNA", layer = "counts")))) {
          DefaultAssay(seurat_object) <- "RNA"
          
          min_cutoff <- ifelse(is.na(input$min_cutoff_feature_merge), NA, input$min_cutoff_feature_merge)
          max_cutoff <- ifelse(is.na(input$max_cutoff_feature_merge), NA, input$max_cutoff_feature_merge)
          
          plot <- FeaturePlot(
            seurat_object,
            features = gene_list,
            min.cutoff = min_cutoff,
            max.cutoff = max_cutoff
          )
          
          if (input$show_coexpression_merge && length(gene_list) == 2) {
            plot <- FeaturePlot(
              seurat_object,
              features = gene_list,
              blend = TRUE,
              blend.threshold = 1,
              order = TRUE,
              min.cutoff = min_cutoff,
              max.cutoff = max_cutoff
            )
          }
          
          if (input$add_noaxes_feature_merge) { plot <- plot + NoAxes() }
          if (input$add_nolegend_feature_merge) { plot <- plot + NoLegend() }
          
          plot <- plot + theme(
            text = element_text(size = input$text_size, face = if (input$bold_text) "bold" else "plain")
          )
          
          feature_plot_merge(plot)
        } else {
          showNotification("Requested genes are not present in the dataset.", type = "error")
        }
      }, error = function(e) {
        showNotification(paste0("Error during FeaturePlot generation: ", e$message), type = "error")
      })
    })
    
    output$FeaturePlot2 <- renderPlot({
      req(feature_plot_merge())
      feature_plot_merge()
    })
    

    # VlnPlot of the selected gene
    vln_plot_merge <- reactiveVal()
    
    # Violin Plot
    observeEvent(input$runVlnPlot, {
      tryCatch({
        req(input$gene_list_vln_merge)
        gene_list <- unique(trimws(strsplit(input$gene_list_vln_merge, ",")[[1]]))
        
        seurat_object <- multiple_datasets_object()
        req(seurat_object)
        
        selected_group_by <- input$group_by_select
        if (all(gene_list %in% rownames(LayerData(seurat_object, assay = "RNA", layer = "counts")))) {
          DefaultAssay(seurat_object) <- "RNA"
          
          plot <- VlnPlot(seurat_object, features = gene_list, group.by = selected_group_by, pt.size = ifelse(input$hide_vln_points_merge, 0, 1)) +
            NoLegend() +
            theme(
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
          
          vln_plot_merge(plot)
        } else {
          showNotification("Requested genes are not present in the dataset.", type = "error")
        }
      }, error = function(e) {
        showNotification(paste0("Error during Violin Plot generation: ", e$message), type = "error")
      })
    })
    
    output$VlnPlot2 <- renderPlot({
      req(vln_plot_merge())
      vln_plot_merge()
    })
    
    

    # DotPlot
    dot_plot_merge <- reactiveVal()

    # Dot Plot
    observeEvent(input$runDotPlot, {
      tryCatch({
        req(input$gene_list_dot_merge)
        gene_list <- unique(trimws(strsplit(input$gene_list_dot_merge, ",")[[1]]))
        
        seurat_object <- multiple_datasets_object()
        req(seurat_object)
        
        selected_group_by <- input$group_by_select
        if (selected_group_by %in% colnames(seurat_object@meta.data)) {
          if (all(gene_list %in% rownames(LayerData(seurat_object, assay = "RNA", layer = "counts")))) {
            DefaultAssay(seurat_object) <- "RNA"
            
            plot <- DotPlot(seurat_object, features = gene_list, group.by = selected_group_by) + RotatedAxis()
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
            
            dot_plot_merge(plot)
          } else {
            showNotification("Some of the requested genes are not present in the RNA assay of the dataset.", type = "error")
          }
        } else {
          showNotification("Selected group by option is not available in Seurat object metadata.", type = "error")
        }
      }, error = function(e) {
        showNotification(paste0("Error during DotPlot generation: ", e$message), type = "error")
      })
    })
    
    output$DotPlot2 <- renderPlot({ req(dot_plot_merge()); dot_plot_merge() })


    # Ridge Plot
   ridge_plot_merge <- reactiveVal()

   # Ridge Plot
   observeEvent(input$runRidgePlot, {
     tryCatch({
       req(input$gene_list_ridge_merge)
       gene_list <- unique(trimws(strsplit(input$gene_list_ridge_merge, ",")[[1]]))
       
       seurat_object <- multiple_datasets_object()
       req(seurat_object)
       
       selected_group_by <- input$group_by_select
       if (all(gene_list %in% rownames(LayerData(seurat_object, assay = "RNA", layer = "counts")))) {
         DefaultAssay(seurat_object) <- "RNA"
         
         plot <- RidgePlot(seurat_object, features = gene_list, group.by = selected_group_by)
         plot <- plot + theme(
           text = element_text(size = input$text_size, face = if (input$bold_text) "bold" else "plain")
         )
         
         ridge_plot_merge(plot)
       } else {
         showNotification("Requested genes are not present in the dataset.", type = "error")
       }
     }, error = function(e) {
       showNotification(paste0("Error during RidgePlot generation: ", e$message), type = "error")
     })
   })

    output$downloadVlnPlotMerge <- downloadHandler(
      filename = function() { paste("VlnPlot", Sys.Date(), ".png", sep = "") },
      content = function(file) {
        tryCatch({
          ggsave(file, plot = vln_plot_merge(), dpi = input$dpi_input_merge, width = 10, height = 6)
        }, error = function(e) {
          showNotification(paste("Error saving VlnPlot: ", e$message), type = "error")
        })
      }
    )

    output$downloadFeaturePlotMerge <- downloadHandler(
      filename = function() { paste("FeaturePlot", Sys.Date(), ".png", sep = "") },
      content = function(file) {
        tryCatch({
          ggsave(file, plot = feature_plot_merge(), dpi = input$dpi_input_merge, width = 10, height = 6)
        }, error = function(e) {
          showNotification(paste("Error saving FeaturePlot: ", e$message), type = "error")
        })
      }
    )

    output$downloadDotPlotMerge <- downloadHandler(
      filename = function() { paste("DotPlot", Sys.Date(), ".png", sep = "") },
      content = function(file) {
        tryCatch({
          ggsave(file, plot = dot_plot_merge(), dpi = input$dpi_input_merge, width = 10, height = 6)
        }, error = function(e) {
          showNotification(paste("Error saving DotPlot: ", e$message), type = "error")
        })
      }
    )

    output$downloadRidgePlotMerge <- downloadHandler(
      filename = function() { paste("RidgePlot", Sys.Date(), ".png", sep = "") },
      content = function(file) {
        tryCatch({
          ggsave(file, plot = ridge_plot_merge(), dpi = input$dpi_input_merge, width = 10, height = 6)
        }, error = function(e) {
          showNotification(paste("Error saving RidgePlot: ", e$message), type = "error")
        })
      }
    )


    number_of_nuclei_merge <- reactiveVal(NULL)


    observeEvent(input$analyze_btn_genes_expression_merge, {
      tryCatch({
        req(input$gene_list_genes_expression_merge, multiple_datasets_object())

        # Récupération de l'objet Seurat
        seurat_obj <- multiple_datasets_object()

        # Vérification de la présence des gènes dans l'objet Seurat
        gene_data <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
        available_genes <- rownames(gene_data)
        missing_genes <- setdiff(input$gene_list_genes_expression_merge, available_genes)
        if (length(missing_genes) > 0) {
          showNotification(paste("The following genes are not present in the dataset:", paste(missing_genes, collapse = ", ")), type = "warning")
        }
        selected_genes <- intersect(input$gene_list_genes_expression_merge, available_genes)
        if (length(selected_genes) == 0) {
          showNotification("No selected genes are present in the dataset.", type = "error")
          return()
        }

        # Vérifier que le champ 'dataset' existe dans les métadonnées
        if (!"dataset" %in% colnames(seurat_obj@meta.data)) {
          showNotification("The 'dataset' column is not found in the Seurat object's metadata.", type = "error")
          return()
        }

        # Vérifier que le seuil est défini
        if (is.null(input$logfc_threshold_genes_expression_merge) || input$logfc_threshold_genes_expression_merge == "") {
          showNotification("The logFC threshold is not specified.", type = "error")
          return()
        }

        # Extraction des données y compris les identités des clusters et des datasets
        data <- FetchData(seurat_obj, vars = c("ident", "dataset", selected_genes))

        # Vérification des identifiants de clusters et de datasets
        print("Available clusters in the dataset:")
        print(unique(data$ident))
        print("Available datasets in the dataset:")
        print(unique(data$dataset))

        # Calcul des métriques pour chaque gène
        expression_summary_list <- lapply(selected_genes, function(gene) {
          gene_data <- data[[gene]]
          print(paste("Processing gene:", gene))
          print(head(gene_data))  # Log pour vérifier les valeurs des données de comptage

          # Vérification du seuil
          print(paste("Threshold for expression:", input$logfc_threshold_genes_expression_merge))

          # Appliquer un seuil pour considérer un gène comme exprimé
          expressed_indices <- gene_data > input$logfc_threshold_genes_expression_merge
          print(paste("Gene data for", gene, ":", head(gene_data)))  # Imprimer les valeurs de gene_data
          print(paste("Expressed indices for gene", gene, ":", which(expressed_indices)))  # Log des indices exprimés
          print(paste("Expressed indices values for gene", gene, ":", gene_data[expressed_indices]))  # Imprimer les valeurs exprimées

          # Comparaison directe des valeurs
          comparison_results <- gene_data > input$logfc_threshold_genes_expression_merge
          print(paste("Comparison results for gene", gene, ":", head(comparison_results)))

          cells_expressed <- sum(expressed_indices, na.rm = TRUE)  # Pour gérer les valeurs NA
          total_cells <- length(gene_data)

          print(paste("Total cells for gene", gene, ":", total_cells))
          print(paste("Cells expressed for gene", gene, ":", cells_expressed))

          if (cells_expressed > 0) {
            cluster_info <- data$ident[expressed_indices]
            dataset_info <- data$dataset[expressed_indices]

            cluster_dataset_info <- paste(cluster_info, dataset_info, sep = "_")
            cluster_dataset_counts <- table(cluster_dataset_info)
            pct_per_cluster_dataset <- prop.table(cluster_dataset_counts) * 100

            df <- data.frame(
              Gene = gene,
              Cluster = sapply(strsplit(names(cluster_dataset_counts), "_"), `[`, 1),
              Dataset = sapply(strsplit(names(cluster_dataset_counts), "_"), `[`, 2),
              Cells_Expressed = as.numeric(cluster_dataset_counts),
              Percentage_of_Cells = as.numeric(pct_per_cluster_dataset),
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
          number_of_nuclei_merge(expression_df)  # Mise à jour de la variable réactive pour le téléchargement
        } else {
          expression_df <- data.frame()
          showNotification("No results meet the thresholds specified.", type = "warning")
        }

        # Affichage du tableau de résultats
        output$expression_summary_merge <- renderDataTable({
          datatable(expression_df, options = list(pageLength = 10, scrollX = TRUE))
        })

      }, error = function(e) {
        showNotification(paste("Error processing expression data: ", e$message), type = "error")
        print(paste("Error detail:", e$message))  # Log de détail de l'erreur pour le débogage
      })
    })

    # Téléchargement du tableau de comparaison des gènes
    output$download_genes_number_expression_merge <- downloadHandler(
      filename = function() {
        paste("number_of_gene_expression", Sys.Date(), ".csv", sep = "")
      },
      content = function(file) {
        tryCatch({
          req(number_of_nuclei_merge())
          write.csv(number_of_nuclei_merge(), file, row.names = FALSE)
        }, error = function(e) {
          showNotification(paste0("Error downloading number of gene expression comparison table: ", e$message), type = "error")
        })
      },
      contentType = "text/csv"
    )



############################## Heatmap and dual expression multi dataset ##############################

    # Tab 5: Heatmap and dual expression multi dataset

    heatmap_plot_multidataset <- reactiveVal()



    # Updated selectInput to choose how to group data
    observe({
      updateSelectInput(session, "dataset_select_heatmap", choices = reactive_metadata_fields())
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
      req(multiple_datasets_object(), input$gene_select_heatmap_multi)
      selected_genes <- input$gene_select_heatmap_multi
      seurat_object <- multiple_datasets_object()
      selected_group_by <- input$dataset_select_heatmap

      tryCatch({
        if (!input$select_all_clusters_merge) {
          if (nchar(input$text_clusters_merge) > 0) {
            specified_clusters <- unlist(strsplit(trimws(input$text_clusters_merge), ",\\s*"))
            valid_clusters <- specified_clusters %in% levels(Idents(seurat_object))
            if (sum(valid_clusters) == 0) {
              showNotification("None of the specified clusters is valid.", type = "error")
              return()
            }
            # Réordonner les niveaux du facteur des clusters selon l'ordre spécifié
            Idents(seurat_object) <- factor(Idents(seurat_object), levels = specified_clusters)

            seurat_object <- subset(seurat_object, idents = specified_clusters)
          } else {
            selected_clusters <- input$cluster_selector_merge
            Idents(seurat_object) <- factor(Idents(seurat_object), levels = selected_clusters)
            seurat_object <- subset(seurat_object, idents = selected_clusters)
          }
        }

        # Vérification finale des gènes valides dans l'objet Seurat filtré
        final_valid_genes <- selected_genes %in% rownames(LayerData(seurat_object, assay="RNA", layer='counts'))

        if(all(final_valid_genes)) {
          # Génération et affichage de la heatmap
          plot <- DoHeatmap(seurat_object, features = selected_genes, group.colors = c("lightgrey", "blue", "red"), group.by = selected_group_by) + NoLegend()
          heatmap_plot_multidataset(plot)
        } else {
          showNotification("Some of the selected genes are not found in the filtered dataset.", type = "error")
        }
      }, error = function(e) {
        showNotification("Genes are not found", type = "error")
      })
    })

    # Rendre la heatmap dans l'interface utilisateur
    output$heatmap_plot_multi <- renderPlot({
      req(heatmap_plot_multidataset())
      heatmap_plot_multidataset()
    })



    # Handler pour télécharger la heatmap
    output$download_heatmap_multi <- downloadHandler(
      filename = function() {
        paste0("heatmap_multi_", Sys.Date(), ".png")
      },
      content = function(file) {
        req(heatmap_plot_multidataset())
        tryCatch({
          ggsave(file, plot = heatmap_plot_multidataset(), width = 10, height = 8, dpi = input$dpi_heatmap_multi)
        }, error = function(e) {
          showNotification(paste("Erreur lors du téléchargement de la heatmap : ", e$message), type = "error")
        })
      }
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

    # Handler pour télécharger la heatmap
    output$download_scatter_multi <- downloadHandler(
      filename = function() {
        paste0("scatter_multi_", Sys.Date(), ".png")
      },
      content = function(file) {
        req(scatter_plot_multidataset())
        tryCatch({
          ggsave(file, plot = scatter_plot_multidataset(), width = 10, height = 8, dpi = input$dpi_scatter_multi)
        }, error = function(e) {
          showNotification(paste("Error donwloading scatter plot : ", e$message), type = "error")
        })
      }
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


    # Reactive variable for that tab
    cluster_colours_merge <- reactiveVal()

    observe({
      if (!is.null(multiple_datasets_object())) {
        updateSelectInput(session, "select_cluster_merge", choices = unique(Idents(multiple_datasets_object())))
        update_cluster_colours(multiple_datasets_object())
      }
    })

    # Function to update cluster colors
    update_cluster_colours <- function(seurat_object) {
      unique_idents <- unique(Idents(seurat_object))
      current_colours <- scales::hue_pal()(length(unique_idents))
      names(current_colours) <- sort(unique_idents)
      cluster_colours_merge(current_colours)
    }

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
      update_cluster_colours(updated_seurat)
      updateSelectInput(session, "select_cluster_merge", choices = unique(Idents(updated_seurat)))
      
      # Mise à jour de 'ClusterIdents' pour refléter les nouveaux noms
      updated_seurat$ClusterIdents <- Idents(updated_seurat)
      multiple_datasets_object(updated_seurat)
      
      # Mise à jour des choix du sélecteur 'group_by_select' pour inclure 'ClusterIdents'
      update_group_by_choices()
    })
    
    observeEvent(input$update_colour_merge_button, {
      req(input$select_color_merge, input$select_cluster_merge_color, cluster_colours_merge())
      current_colours <- cluster_colours_merge()
      current_colours[input$select_color_merge] <- input$select_cluster_merge_color
      cluster_colours_merge(current_colours)
    })


    # Finale UMAP
    output$umap_finale_merge <- renderPlotly({
      req(multiple_datasets_object(), cluster_colours_merge())
      tryCatch({
        plot_data <- DimPlot(multiple_datasets_object(), group.by = "ident", pt.size = input$pt_size_merge, label = TRUE, label.size = input$label_font_size_merge) +
          theme(axis.line = element_line(size = 0.5)) +
          scale_color_manual(values = cluster_colours_merge()) +
          NoLegend()
        interactive_plot <- ggplotly(plot_data, tooltip = "text")
        interactive_plot <- interactive_plot %>%
          layout(
            title = list(text = input$plot_title_merge, font = list(size = 24)),
            hovermode = "closest"
          )
        return(interactive_plot)
      }, error = function(e) {
        # Gestion d'erreur
        showNotification(paste("Erreur lors de la création de la UMAP finale : ", e$message), type = "error")
        return(NULL)
      })
    })



############################## Calculation of differentially expressed genes ##############################


    # Reactive variable for that tab
    diff_genes_compare <- reactiveVal() # Create a new reactive value to store markers_df
    diff_genes_compare_cluster <- reactiveVal()   # Define a new reactive value to store comparison tables
    diff_genes_compare_datasets <- reactiveVal()   # Define a new reactive value to store comparison tables
    filtered_umap_plot <- reactiveVal()

    shinyjs::disable("download_markers_single_cluster_merge")
    shinyjs::disable("download_markers_multiple_clusters_merge")
    shinyjs::disable("download_diff_dataset_cluster")

    # Drop-down menu to filter by Dataset
    output$dataset_filter_ui <- renderUI({
      req(multiple_datasets_object())
      selectInput("dataset_filter",
                  label = "Filter by Dataset",
                  choices = c(unique(multiple_datasets_object()@meta.data$dataset)),
                  selected = unique(multiple_datasets_object()@meta.data$dataset),
                  multiple = TRUE)
    })

    # Display of filtered UMAP graph
    output$filtered_umap_plot <- renderPlot({
      req(multiple_datasets_object(), input$dataset_filter)

      if (!"dataset" %in% colnames(multiple_datasets_object()@meta.data)) {
        stop("The dataset column was not found in the Seurat object's metadata.")
      }

      plot <- if ("Tous les Datasets" %in% input$dataset_filter || is.null(input$dataset_filter)) {
        DimPlot(multiple_datasets_object(), group.by = "ident", label = TRUE, label.size = 5) + NoLegend() + ggtitle(NULL)
      } else {
        valid_datasets <- input$dataset_filter %in% unique(multiple_datasets_object()@meta.data$dataset)
        subset_seurat <- subset(multiple_datasets_object(), subset = dataset %in% input$dataset_filter[valid_datasets])
        DimPlot(subset_seurat, group.by = "ident", label = TRUE, label.size = 5) + NoLegend() + ggtitle(NULL)
      }
      filtered_umap_plot(plot)  # Stockage du plot dans la variable réactive
      plot  # Affichage du plot
    })


    output$download_filtered_umap_plot <- downloadHandler(
      filename = function() {
        paste("UMAP_plot", Sys.Date(), ".png", sep = "")
      },
      content = function(file) {
        tryCatch({
          req(filtered_umap_plot())
          # Définir une valeur par défaut pour DPI si non spécifié dans l'interface utilisateur
          dpi_value <- if (!is.null(input$dpi_filtered_umap_plot)) input$dpi_filtered_umap_plot else 300
          ggsave(file, plot = filtered_umap_plot(), dpi = dpi_value, width = 10, height = 6)
        }, error = function(e) {
          showNotification(paste("Error in UMAP Plot: ", e$message), type = "error")
        })
      }
    )


    # Reactively update the cluster drop-down menu
    observe({
      req(multiple_datasets_object())
      cluster_choices <- unique(Idents(multiple_datasets_object()))
      updateSelectInput(session, "selected_cluster", choices = cluster_choices, selected = cluster_choices[1])
    })

    # Reactive function for markers
    observeEvent(input$calculate_DE, {
      tryCatch({
        # Show modal dialog
        showModal(modalDialog(
          title = "Please Wait",
          "Calculating differentially expressed genes...",
          easyClose = FALSE,
          footer = NULL
        ))

        showNotification("Calculating differentially expressed genes...", type = "message")

        subset_seurat <- multiple_datasets_object()

        # Try to calculate markers
        markers <- FindMarkers(subset_seurat,
                               ident.1 = input$selected_cluster,
                               min.pct = input$min_pct_merge,
                               logfc.threshold = input$logfc_threshold_merge, slot='data', assay='RNA')

        if (nrow(markers) > 0) {
          markers$p_val <- format(markers$p_val, scientific = TRUE, digits = 3)
          markers$p_val_adj <- format(markers$p_val_adj, scientific = TRUE, digits = 3)
          cleaned_gene_names <- clean_gene_names_for_html(rownames(markers))

          markers$gene <- paste0('<a href="#" class="gene-name" data-gene="', cleaned_gene_names, '">', cleaned_gene_names, '</a>')
          markers_df <- as.data.frame(markers)
          diff_genes_compare(markers_df)
          shinyjs::enable("download_markers_single_cluster_merge")
        } else {
          showNotification("No differentially expressed genes found for the selected cluster.", type = "message")
        }

        # Close the modal dialog
        removeModal()
      }, error = function(e) {
        removeModal() # Close the modal dialog in case of an error
        showNotification(paste0("Error when calculating differentially expressed genes: ", e$message), type = "error")
      })
    })



    output$DE_genes_table <- renderDataTable({
      tryCatch({
        datatable(diff_genes_compare(), escape = FALSE)
      }, error = function(e) {
        showNotification(paste0("Error when rendering the gene table: ", e$message), type = "error")
      })
    })


    # Download gene comparison table
    output$download_markers_single_cluster_merge <- downloadHandler(
      filename = function() {
        paste("genes-differents-", Sys.Date(), ".csv", sep="")      },
      content = function(file) {
        tryCatch({
          req(!is.null(diff_genes_compare))
          diff_genes_compare_DL <- diff_genes_compare()
          diff_genes_compare_DL_sorted <-diff_genes_compare_DL[order(diff_genes_compare_DL$p_val), ]
          write.csv(diff_genes_compare_DL_sorted, file)
        }, error = function(e) {
          showNotification(paste0("Error downloading gene comparison table: ", e$message), type = "error")
        })
      },
      contentType = "text/csv"
    )


    # Drop-down menus to select clusters for comparison
    output$cluster1_compare_ui <- renderUI({
      req(multiple_datasets_object())
      selectInput("cluster1_compare",
                  label = "Selec first cluster",
                  choices = unique(Idents(multiple_datasets_object())),
                  selected = NULL)
    })

    output$cluster2_compare_ui <- renderUI({
      req(multiple_datasets_object())
      selectInput("cluster2_compare",
                  label = "Select second cluster",
                  choices = unique(Idents(multiple_datasets_object())),
                  selected = NULL)
    })

    #  Cluster comparison
    observeEvent(input$compare_clusters_button, {
      tryCatch({
        # Show modal dialog
        showModal(modalDialog(
          title = "Please Wait",
          "Finding differentially expressed features...",
          easyClose = FALSE,
          footer = NULL
        ))

        showNotification("Finding differentially expressed features...", type = "message")
        req(input$cluster1_compare, input$cluster2_compare, multiple_datasets_object())
        if (input$cluster1_compare == input$cluster2_compare) {
          showNotification("Please select two different clusters for comparison!", type = "error")
          removeModal() # Close the modal dialog in case of same clusters
          return()
        }

        # Try to find markers
        temp_res <- FindMarkers(multiple_datasets_object(), ident.1 = input$cluster1_compare, ident.2 = input$cluster2_compare,
                                min.pct = input$min_pct_compare_merge,
                                logfc.threshold = input$logfc_threshold_compare_merge)

        temp_res$p_val <- format(temp_res$p_val, scientific = TRUE, digits = 3)
        temp_res$p_val_adj <- format(temp_res$p_val_adj, scientific = TRUE, digits = 3)
        temp_res$gene <- paste0('<a href="#" class="gene-name" data-gene="', rownames(temp_res), '">', rownames(temp_res), '</a>')
        shinyjs::enable("download_markers_multiple_clusters_merge")
        diff_genes_compare_cluster(temp_res)

        # Close the modal dialog
        removeModal()
      }, error = function(e) {
        removeModal() # Close the modal dialog in case of an error
        showNotification(paste0("Cluster comparison error: ", e$message), type = "error")
      })
    })


    output$diff_genes_table_compare <- renderDataTable({
      req(diff_genes_compare_cluster())
      tryCatch({
        datatable(diff_genes_compare_cluster(), escape = FALSE)
      }, error = function(e) {
        showNotification(paste0("Error when rendering the gene table: ", e$message), type = "error")
      })
    })

    # Download gene comparison table
    output$download_markers_multiple_clusters_merge <- downloadHandler(
      filename = function() {
        paste("diff-genes-comparison-", input$cluster1_compare, "-VS-", input$cluster2_compare, "-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        tryCatch({
          diff_genes_compare_cluster_DL <- diff_genes_compare_cluster()
          req(!is.null(diff_genes_compare_cluster_DL))
          diff_genes_compare_cluster_DL_sorted <- diff_genes_compare_cluster_DL[order(diff_genes_compare_cluster_DL$p_val), ]
          write.csv(diff_genes_compare_cluster_DL_sorted, file)
        }, error = function(e) {
          showNotification(paste0("Error downloading gene comparison table: ", e$message), type = "error")
        })
      },
      contentType = "text/csv"
    )



    output$dataset1_compare_ui <- renderUI({
      req(multiple_datasets_object())
      selectInput("dataset1_compare",
                  label = "Select the first dataset",
                  choices = unique(multiple_datasets_object()@meta.data$dataset),
                  selected = NULL)
    })

    output$dataset2_compare_ui <- renderUI({
      req(multiple_datasets_object())
      selectInput("dataset2_compare",
                  label = "Select the second dataset",
                  choices = unique(multiple_datasets_object()@meta.data$dataset),
                  selected = NULL)
    })

    output$cluster_compare_ui <- renderUI({
      req(multiple_datasets_object())
      selectInput("cluster_compare",
                  label = "Select a cluster",
                  choices = unique(Idents(multiple_datasets_object())),
                  selected = NULL)
    })
    # Observer for comparing datasets
    observeEvent(input$compare_datasets_button, {
      tryCatch({
        # Show modal dialog
        showModal(modalDialog(
          title = "Please Wait",
          "Comparing datasets...",
          easyClose = FALSE,
          footer = NULL
        ))

        req(input$dataset1_compare, input$dataset2_compare, multiple_datasets_object())

        if (input$dataset1_compare == input$dataset2_compare) {
          showNotification("Please select two different datasets for comparison!", type = "error")
          removeModal() # Close the modal dialog in case of same datasets
          return()
        }

        # Use subset.ident to filter directly in FindMarkers if all clusters are not selected
        temp_res <- FindMarkers(
          object = multiple_datasets_object(),
          ident.1 = input$dataset1_compare,
          ident.2 = input$dataset2_compare,
          group.by = "dataset",
          subset.ident = if (input$all_clusters) NULL else input$cluster_compare,
          min.pct = input$min_pct_compare_dataset_merge,
          logfc.threshold = input$logfc_threshold_datasets
        )

        temp_res$p_val <- format(temp_res$p_val, scientific = TRUE, digits = 3)
        temp_res$p_val_adj <- format(temp_res$p_val_adj, scientific = TRUE, digits = 3)
        temp_res$gene <- paste0('<a href="#" class="gene-name" data-gene="', rownames(temp_res), '">', rownames(temp_res), '</a>')
        diff_genes_compare_datasets(temp_res)

        print("Differentially expressed genes calculation completed.")

        if (nrow(temp_res) == 0) {
          showNotification("No differentially expressed genes found!", type = "error")
          removeModal() # Close the modal dialog in case no DEGs are found
          return()
        }

        # Update UI with results
        output$diff_dataset_cluster <- renderDataTable({
          datatable(diff_genes_compare_datasets(), escape = FALSE)
        })

        shinyjs::enable("download_diff_dataset_cluster")

        # Close the modal dialog
        removeModal()
      }, error = function(e) {
        removeModal() # Close the modal dialog in case of an error
        showNotification(paste0("Error when comparing datasets for a specific cluster: ", e$message), type = "error")
      })
    })


    # Download gene comparison table
    output$download_diff_dataset_cluster <- downloadHandler(
      filename = function() {
        paste("diff-genes-comparison-", input$cluster1_compare, "-VS-", input$cluster2_compare, "-", Sys.Date(), ".csv", sep="")
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
    
    
    
    
    
    # Function to generate the cluster table
    generate_cluster_table <- function() {
      req(multiple_datasets_object())
      seurat_object <- multiple_datasets_object()
      
      cluster_data <- table(Idents(seurat_object), seurat_object$orig.ident)
      cluster_df <- as.data.frame.matrix(cluster_data)
      
      cluster_df$Cluster <- rownames(cluster_df)
      cluster_df <- cluster_df[, c(ncol(cluster_df), 1:(ncol(cluster_df) - 1))]
      
      return(cluster_df)
    }
    
    observeEvent(input$generate_cluster_table, {
      output$cluster_table <- renderDataTable({
        generate_cluster_table()
      })
    })
    

############################## Subseting seurat object ##############################


    # Variable réactive pour cet onglet
    subset_seurat_merge <- reactiveVal(NULL)
    shinyjs::disable("download_subset_merge")

    # Mise à jour des choix d'identité de cellules pour le sous-ensemble
    observe({
      if (!is.null(multiple_datasets_object())) {
        updateSelectInput(session, "select_ident_subset_merge", choices = unique(Idents(multiple_datasets_object())))
      }
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
    output$download_subset_merge <- downloadHandler(
      filename = function() {
        paste("seurat_subset-", Sys.Date(), ".rds", sep="")
      },
      content = function(file) {
        showModal(modalDialog(
          title = "Please Wait",
          "Preparing the seurat object for download...",
          easyClose = FALSE,
          footer = NULL
        ))

        tryCatch({
          req(subset_seurat_merge())
          saveRDS(subset_seurat_merge(), file)
          removeModal()
        }, error = function(e) {
          removeModal()
          showNotification(paste0("Error while downloading Seurat subset: ", e$message), type = "error")
        })
      }
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

}
