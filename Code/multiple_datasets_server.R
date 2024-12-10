
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

      output$datasets_loaded <- reactive({
        !is.null(multiple_datasets_object())
      })
      outputOptions(output, 'datasets_loaded', suspendWhenHidden = FALSE)

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
        }, error = function(e) {
          showNotification(paste("An error occurred: ", e), type = "error")
        })
      })

      # Processing of loaded data with progress bar updates at the end of each step
      observe_file_input <- function(index) {
        # Ajouter isolate pour éviter les déclenchements multiples
        observeEvent(input[[paste0("merge", index)]], {
          # Vérifier si le fichier a déjà été traité
          if (!is.null(data_loaded[[paste0("loaded", index)]]) &&
              data_loaded[[paste0("loaded", index)]] == TRUE) {
            print(paste("Dataset", index, "already processed, skipping..."))
            return()
          }

          req(input[[paste0("merge", index)]])

          withProgress(message = paste("Processing dataset", index), value = 0, {
            tryCatch({
              print(paste("Starting data processing for dataset", index))
              dataset_type <- input[[paste0("dataset_type_merge", index)]]
              seurat_object <- NULL
              mt_pattern_merge <- ifelse(input$species_choice_merge == "mouse", "^mt-", "^MT-")

              # Get dataset name
              dataset_name <- input[[paste0("dataset_name", index)]]
              if (is.null(dataset_name) || dataset_name == "") {
                dataset_name <- paste0("Dataset_", index)
              }

              if (dataset_type == "seurat_object_merge") {
                print(paste("Loading pre-processed Seurat object for dataset", index))
                seurat_object <- readRDS(input[[paste0("merge", index)]]$datapath)

                if (!inherits(seurat_object, "Seurat")) {
                  stop("Loaded file is not a valid Seurat object")
                }

                DefaultAssay(seurat_object) <- "RNA"

              } else {
                # Raw data workflow
                print(paste("Processing raw data for dataset", index))

                # Create and clean unzip directory
                unzip_dir <- paste0("unzipped", index)
                if (dir.exists(unzip_dir)) {
                  unlink(unzip_dir, recursive = TRUE)
                }
                dir.create(unzip_dir)

                # Unzip and read data
                print(paste("Reading 10X data for dataset", index))
                unzip(input[[paste0("merge", index)]]$datapath, exdir = unzip_dir)
                data <- Read10X(unzip_dir)

                if (dataset_type == "snRNA_merge") {
                  print(paste("Creating Seurat object for snRNA data, dataset", index))
                  seurat_object <- CreateSeuratObject(
                    counts = data,
                    project = dataset_name,
                    min.cells = 3,
                    min.features = 200
                  )

                } else if (dataset_type == "multiome_merge") {
                  print(paste("Creating Seurat object for multiome data, dataset", index))
                  seurat_object <- CreateSeuratObject(
                    counts = data$`Gene Expression`,
                    project = dataset_name,
                    min.cells = 3,
                    min.features = 200
                  )
                }

                # Cleanup
                unlink(unzip_dir, recursive = TRUE)

                # Basic preprocessing for raw data
                print(paste("Processing raw data for dataset", index))
                seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = mt_pattern_merge)
                seurat_object <- NormalizeData(seurat_object, verbose = FALSE)
                seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 4000, verbose = FALSE)
                seurat_object <- ScaleData(seurat_object, features = rownames(seurat_object), verbose = FALSE)
                seurat_object <- RunPCA(seurat_object, features = VariableFeatures(seurat_object), verbose = FALSE)
              }

              # Add dataset info
              seurat_object$dataset <- dataset_name
              seurat_object$orig.ident <- dataset_name

              # Store the object
              seurat_objects[[paste0("seurat_object", index)]] <- seurat_object
              data_loaded[[paste0("loaded", index)]] <- TRUE

              print(paste("Dataset", index, "processed and stored successfully"))
              showNotification(paste("Dataset", index, "processed successfully!"), type = "message")

              # Enable metadata fields
              shinyjs::enable("add_field")
              shinyjs::enable("add_metadata")

            }, error = function(e) {
              print(paste("Error in dataset", index, ":", e$message))
              showNotification(paste("Error processing dataset", index, ":", e$message), type = "error")

              if (dir.exists(paste0("unzipped", index))) {
                unlink(paste0("unzipped", index), recursive = TRUE)
              }
            })
          })
        }, ignoreInit = TRUE, once = TRUE)  # Ajouter ces options pour éviter les déclenchements multiples
      }

      # Open the modal to load datasets
      observeEvent(input$open_file_input_modal, {
        showModal(modalDialog(
          title = "Load and Integrate Datasets",
          easyClose = FALSE,  # Prevent closing until integration is done
          footer = NULL,      # Footer controlled by the Integrate button
          tagList(
            numericInput('num_datasets', 'Number of datasets to upload:', value = 2, min = 1),
            radioButtons("species_choice_merge", label = "Select Species:", choices = list("Mouse" = "mouse", "Human" = "human"), selected = "mouse"),
            uiOutput("fileInputs"),  # Dynamic file input fields
            actionButton('integrate', 'Integrate', class = 'btn-primary', disabled = TRUE)  # Integrate button
          )
        ))

        # Generate file input fields dynamically based on the number of datasets
        observeEvent(input$num_datasets, {
          req(input$num_datasets > 0)
          output$fileInputs <- renderUI({
            lapply(1:input$num_datasets, function(i) {
              fluidRow(
                column(4, textInput(paste0('dataset_name', i), paste0('Dataset name ', i), value = paste0("Dataset ", i))),
                column(4, selectInput(paste0("dataset_type_merge", i), "Data type:",
                                      choices = list("snRNA-seq" = "snRNA_merge", "Multiome" = "multiome_merge", "Seurat Object" = "seurat_object_merge"))),
                column(4, fileInput(paste0('merge', i), paste0('Choose a Seurat or .gz file ', i), accept = c('.rds', '.zip')))
              )
            })
          })
        })

        # Observe file inputs to trigger preprocessing
        observeEvent(input$num_datasets, {
          lapply(1:input$num_datasets, function(i) {
            observe_file_input(i)  # This calls your existing observe_file_input function
          })
        })

        # Enable the Integrate button once all files are uploaded
        observe({
          req(input$num_datasets)

          # Pour les données brutes
          if (any(sapply(1:input$num_datasets, function(i) {
            type <- input[[paste0("dataset_type_merge", i)]]
            return(!is.null(type) && type != "seurat_object_merge")
          }))) {
            # Active le bouton que si tous les datasets sont traités
            num_loaded_datasets <- sum(sapply(1:input$num_datasets, function(i) {
              !is.null(data_loaded[[paste0("loaded", i)]]) && data_loaded[[paste0("loaded", i)]] == TRUE
            }))
            shinyjs::toggleState("integrate", condition = num_loaded_datasets == input$num_datasets)
          } else {
            # Pour les objets Seurat .rds, active le bouton dès que les fichiers sont chargés
            files_loaded <- all(sapply(1:input$num_datasets, function(i) {
              !is.null(input[[paste0("merge", i)]])
            }))
            shinyjs::toggleState("integrate", condition = files_loaded)
          }
        })
      })

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
          dims = 1:30,
          k.filter = 30,
          k.score = 20,
          k.anchor = 5,
          anchor.features = features
        )
        print("Integrating Data")
        seurat_integrated_temp <- IntegrateData(
          anchorset = anchors,
          dims = 1:30,
          features.to.integrate = features,
          k.weight = 50
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
          # Show modal dialog during integration (only once, at the integration step)
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
            print(paste("Added Seurat object", i, "to the integration list"))
          }
          print(paste("seurat_list length:", length(seurat_list)))

          if (length(seurat_list) < 2) {
            stop("Please upload at least two datasets for integration")
          }

          integrated_object <- integrate_data(seurat_list)
          multiple_datasets_object(integrated_object)  # Update the integrated object
          print("Finished data integration")

          if (!is.null(multiple_datasets_object())) {
            print("The integrated object is updated in multiple_datasets_object.")
            shinyjs::enable("add_field")
            shinyjs::enable("add_metadata")
          } else {
            print("Failed to update the integrated object in multiple_datasets_object.")
          }

          # Close the modal dialog after integration
          removeModal()

          # Update UI elements after integration
          updateUIElements()
        }, error = function(e) {
          removeModal()  # Close the modal dialog in case of an error
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
        paste("UMAP_plot", Sys.Date(), ".tiff", sep = "")
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

    # Call the update function when a new dataset is loaded or when the app starts
    observe({
      req(multiple_datasets_object())
      get_cluster_colors_merge(multiple_datasets_object())  # Initialize the cluster colors
    })


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

    # Observer unifié pour la gestion des métadonnées et clusters
    observe({
      req(multiple_datasets_object())
      seurat_object <- multiple_datasets_object()

      # 1. Définition des colonnes à ne pas montrer dans les sélecteurs
      excluded_patterns <- c(
        "^RNA_snn_res",
        "^pANN",
        "^DF",
        "^integrated",
        "^integrated_snn"
      )
      pattern <- paste(excluded_patterns, collapse = "|")

      # 2. Obtenir toutes les colonnes de métadonnées disponibles
      metadata_fields <- colnames(seurat_object@meta.data)
      # Filtrer seulement pour l'affichage
      display_fields <- metadata_fields[!grepl(pattern, metadata_fields)]

      # 3. Mise à jour des sélecteurs
      # 3.1 Group by select (garder toutes les options valides)
      updateSelectInput(session, "group_by_select",
                        choices = display_fields,
                        selected = input$group_by_select  # Garder la sélection actuelle
      )

      # 3.2 Metadata to compare (quand dataset est sélectionné)
      if(input$group_by_select == "dataset") {
        updateSelectInput(session, "metadata_to_compare",
                          choices = display_fields[display_fields != "dataset"],
                          selected = input$metadata_to_compare  # Garder la sélection actuelle
        )

        # 3.3 Mise à jour des clusters selon la colonne sélectionnée
        if(!is.null(input$metadata_to_compare)) {
          clusters <- unique(seurat_object@meta.data[[input$metadata_to_compare]])
          updateSelectInput(session, "cluster_order",
                            choices = clusters,
                            selected = clusters
          )
        }
      } else {
        # Si pas groupé par dataset, utiliser la colonne sélectionnée
        clusters <- unique(seurat_object@meta.data[[input$group_by_select]])
        updateSelectInput(session, "cluster_order",
                          choices = clusters,
                          selected = clusters
        )
      }
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
        req(input$gene_list_feature_merge, multiple_datasets_object())
        genes <- unique(trimws(strsplit(input$gene_list_feature_merge, ",")[[1]]))
        seurat_object <- multiple_datasets_object()
        req(seurat_object)
        DefaultAssay(seurat_object) <- input$viz_assay_merge
        
        # Préparer les paramètres de seuils
        min_cut <- if (!is.na(input$min_cutoff_feature_merge)) input$min_cutoff_feature_merge else NA
        max_cut <- if (!is.na(input$max_cutoff_feature_merge)) input$max_cutoff_feature_merge else NA
        
        if(input$group_by_select == "dataset") {
          # Créer un plot pour chaque dataset
          plot_list <- lapply(unique(seurat_object$dataset), function(ds) {
            # Subset pour le dataset actuel
            subset_obj <- subset(seurat_object, subset = dataset == ds)
            # Créer le FeaturePlot pour ce dataset avec seuils
            p <- FeaturePlot(
              subset_obj,
              features = genes,
              ncol = ifelse(length(genes) > 1, 2, 1),
              reduction = "umap",
              min.cutoff = min_cut,
              max.cutoff = max_cut
            ) +
              ggtitle(ds) +
              theme(plot.title = element_text(size = 14, face = "bold"))
            
            # Appliquer les options
            if(input$add_nolegend_feature_merge) p <- p + NoLegend()
            if(input$add_noaxes_feature_merge) p <- p + NoAxes()
            return(p)
          })
          # Combiner tous les plots
          combined_plot <- wrap_plots(plot_list, ncol = 1)
        } else {
          # Plot normal si pas groupé par dataset
          combined_plot <- FeaturePlot(
            seurat_object,
            features = genes,
            ncol = ifelse(length(genes) > 1, 2, 1),
            reduction = "umap",
            min.cutoff = min_cut,
            max.cutoff = max_cut
          )
          if(input$add_nolegend_feature_merge) combined_plot <- combined_plot + NoLegend()
          if(input$add_noaxes_feature_merge) combined_plot <- combined_plot + NoAxes()
        }
        feature_plot_merge(combined_plot)
        print(paste("Seuils appliqués - Min:", min_cut, "Max:", max_cut))
      }, error = function(e) {
        showNotification(paste("Error in FeaturePlot:", e$message), type = "error")
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
        req(input$gene_list_vln_merge, multiple_datasets_object())
        genes <- unique(trimws(strsplit(input$gene_list_vln_merge, ",")[[1]]))
        seurat_object <- multiple_datasets_object()
        req(seurat_object)
        DefaultAssay(seurat_object) <-  input$viz_assay_merge
        
        # Plot de base avec group.by
        plot <- VlnPlot(
          seurat_object,
          features = genes,
          group.by = input$group_by_select,
          pt.size = ifelse(input$hide_vln_points_merge, 0, 1)
        )
        
        # Options de base uniquement
        if (input$add_noaxes_vln_merge) plot <- plot + NoAxes()
        if (input$add_nolegend_vln_merge) plot <- plot + NoLegend()
        
        vln_plot_merge(plot)
      }, error = function(e) {
        showNotification(paste("Error in VlnPlot: ", e$message), type = "error")
      })
    })
    
    output$VlnPlot2 <- renderPlot({
      req(vln_plot_merge())
      vln_plot_merge()
    })
    
    
    

    # DotPlot reactive value
    dot_plot_merge <- reactiveVal()



    # Dot Plot Generation
    observeEvent(input$runDotPlot, {
      tryCatch({
        req(input$gene_list_dot_merge, multiple_datasets_object())
        genes <- unique(trimws(strsplit(input$gene_list_dot_merge, ",")[[1]]))
        seurat_object <- multiple_datasets_object()
        req(seurat_object)
        DefaultAssay(seurat_object) <- input$viz_assay_merge

        # Si groupé par dataset
        if(input$group_by_select == "dataset") {
          req(input$metadata_to_compare)

          # Créer colonne combinée dataset_category
          seurat_object$dataset_category <- paste(
            seurat_object$dataset,
            seurat_object@meta.data[[input$metadata_to_compare]],
            sep = "_"
          )

          # Créer combinaisons désirées
          if (!is.null(input$cluster_order) && length(input$cluster_order) > 0) {
            valid_combinations <- unlist(lapply(unique(seurat_object$dataset), function(ds) {
              paste(ds, input$cluster_order, sep = "_")
            }))

            # Filtrer les cellules
            cells_to_keep <- seurat_object$dataset_category %in% valid_combinations
            seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[cells_to_keep])

            # Réordonner niveaux
            seurat_object$dataset_category <- factor(
              seurat_object$dataset_category,
              levels = valid_combinations
            )
          }

          # Plot avec colonne combinée
          plot <- DotPlot(
            seurat_object,
            features = genes,
            group.by = "dataset_category"
          )

        } else {
          # Filtrer par clusters sélectionnés si pas groupé par dataset
          if (!is.null(input$cluster_order) && length(input$cluster_order) > 0) {
            cells_to_keep <- seurat_object@meta.data[[input$group_by_select]] %in% input$cluster_order
            seurat_object <- subset(seurat_object, cells = colnames(seurat_object)[cells_to_keep])

            # Réordonner les clusters
            seurat_object@meta.data[[input$group_by_select]] <- factor(
              seurat_object@meta.data[[input$group_by_select]],
              levels = input$cluster_order
            )
          }

          plot <- DotPlot(
            seurat_object,
            features = genes,
            group.by = input$group_by_select
          )
        }

        # Rotation des axes
        if (input$invert_axes) {
          plot <- plot + coord_flip()
        } else {
          plot <- plot + RotatedAxis()
        }

        # Options supplémentaires
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

    output$downloadVlnPlotMerge <- downloadHandler(
      filename = function() { paste("VlnPlot", Sys.Date(), ".tiff", sep = "") },
      content = function(file) {
        tryCatch({
          ggsave(file, plot = vln_plot_merge(), dpi = input$dpi_input_merge, width = 10, height = 6)
        }, error = function(e) {
          showNotification(paste("Error saving VlnPlot: ", e$message), type = "error")
        })
      }
    )

    output$downloadFeaturePlotMerge <- downloadHandler(
      filename = function() { paste("FeaturePlot", Sys.Date(), ".tiff", sep = "") },
      content = function(file) {
        tryCatch({
          ggsave(file, plot = feature_plot_merge(), dpi = input$dpi_input_merge, width = 10, height = 6)
        }, error = function(e) {
          showNotification(paste("Error saving FeaturePlot: ", e$message), type = "error")
        })
      }
    )

    output$downloadDotPlotMerge <- downloadHandler(
      filename = function() { paste("DotPlot", Sys.Date(), ".tiff", sep = "") },
      content = function(file) {
        tryCatch({
          ggsave(file, plot = dot_plot_merge(), dpi = input$dpi_input_merge, width = 10, height = 6)
        }, error = function(e) {
          showNotification(paste("Error saving DotPlot: ", e$message), type = "error")
        })
      }
    )

    output$downloadRidgePlotMerge <- downloadHandler(
      filename = function() { paste("RidgePlot", Sys.Date(), ".tiff", sep = "") },
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
        gene_data <- GetAssayData(seurat_obj, assay = input$viz_assay_merge, slot = "counts")
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

      tryCatch({
        # Gestion des clusters
        if (!input$select_all_clusters_merge) {
          if (nchar(input$text_clusters_merge) > 0) {
            specified_clusters <- unlist(strsplit(trimws(input$text_clusters_merge), ",\\s*"))
            valid_clusters <- specified_clusters %in% levels(Idents(seurat_object))
            if (sum(valid_clusters) == 0) {
              showNotification("None of the specified clusters is valid.", type = "error")
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

        # Sélection des gènes selon l'option choisie
        if(input$use_top10_genes_merge) {
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
          selected_genes <- input$gene_select_heatmap_multi
        }

        # Vérification finale des gènes
        final_valid_genes <- selected_genes %in% rownames(LayerData(seurat_object, assay="RNA", layer='counts'))
        if(any(final_valid_genes)) {
          # Génération et affichage de la heatmap
          plot <- DoHeatmap(seurat_object,
                            features = selected_genes[final_valid_genes],
                            group.colors = c("lightgrey", "blue", "red"),
                            group.by = selected_group_by) + NoLegend()
          heatmap_plot_multidataset(plot)
        } else {
          showNotification("No valid genes found in the dataset.", type = "error")
        }

      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error")
      })
      removeModal()
    })


    # Rendre la heatmap dans l'interface utilisateur
    output$heatmap_plot_multi <- renderPlot({
      req(heatmap_plot_multidataset())
      heatmap_plot_multidataset()
    })



    # Handler pour télécharger la heatmap
    output$download_heatmap_multi <- downloadHandler(
      filename = function() {
        paste0("heatmap_multi_", Sys.Date(), ".tiff")
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
        paste0("scatter_multi_", Sys.Date(), ".tiff")
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

    # Sauvegarder les couleurs des clusters dans l'objet Seurat après mise à jour
    observeEvent(input$update_colour_merge_button, {
      message("Update the color of the selected cluster for multiple datasets")

      # Charger les couleurs actuelles à partir de l'objet Seurat
      updated_seurat <- multiple_datasets_object()

      # Récupérer les couleurs actuelles ou initialiser les couleurs
      cluster_colors <- get_cluster_colors_merge(updated_seurat)

      # Mettre à jour la couleur du cluster sélectionné
      if (!is.null(input$select_color_merge) && input$select_color_merge %in% names(cluster_colors)) {
        cluster_colors[input$select_color_merge] <- input$select_cluster_merge_color
        message(paste("Updating color for cluster", input$select_color_merge, "to", input$select_cluster_merge_color))
      } else {
        showNotification("Selected cluster is not valid.", type = "error")
        return()
      }

      # Sauvegarder les nouvelles couleurs dans @misc
      updated_seurat@misc$cluster_colors <- cluster_colors
      multiple_datasets_object(updated_seurat)  # Mettre à jour l'objet Seurat avec les nouvelles couleurs

      # Vérifier la mise à jour des couleurs
      message("Current cluster colors:")
      print(updated_seurat@misc$cluster_colors)

      showNotification("Cluster colors saved in Seurat object.", type = "message")
    })


    # UMAP pour l'affichage final avec couleurs personnalisées
    output$umap_finale_merge <- renderPlotly({
      req(multiple_datasets_object())
      message("Generating the final UMAP for multiple datasets")

      updated_seurat <- multiple_datasets_object()

      # Récupérer les couleurs depuis @misc
      cluster_colors <- get_cluster_colors_merge(updated_seurat)

      # Générer le plot UMAP avec les couleurs mises à jour
      plot_data <- DimPlot(
        updated_seurat,
        group.by = "ident",
        pt.size = input$pt_size_merge,
        label = TRUE,
        label.size = input$label_font_size_merge
      ) +
        scale_color_manual(values = cluster_colors) +
        NoLegend() +
        theme(axis.line = element_line(size = 0.5))

      # Convertir en plot interactif avec plotly
      interactive_plot <- ggplotly(plot_data, tooltip = "text") %>%
        layout(
          title = list(text = input$plot_title_merge, font = list(size = 24)),
          hovermode = "closest"
        )

      message("UMAP final generated with updated colors for multiple datasets.")
      return(interactive_plot)
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

      # Vérifier si la colonne dataset existe dans les métadonnées
      if (!"dataset" %in% colnames(multiple_datasets_object()@meta.data)) {
        stop("The dataset column was not found in the Seurat object's metadata.")
      }

      # Récupérer les couleurs des clusters depuis l'objet Seurat
      cluster_colors <- get_cluster_colors_merge(multiple_datasets_object())

      # Utiliser la valeur de la case à cocher pour afficher ou masquer les noms des clusters
      show_labels <- input$show_labels_merge

      # Créer le graphique UMAP basé sur les datasets sélectionnés
      plot <- if ("Tous les Datasets" %in% input$dataset_filter || is.null(input$dataset_filter)) {
        # Afficher tous les clusters
        DimPlot(multiple_datasets_object(), group.by = "ident", label = show_labels, label.size = 5) +
          scale_color_manual(values = cluster_colors) +  # Appliquer les couleurs personnalisées
          NoAxes() +             # Supprimer les axes
          NoLegend() +           # Supprimer la légende
          ggtitle(NULL)          # Supprimer le titre
      } else {
        # Filtrer les datasets sélectionnés
        valid_datasets <- input$dataset_filter %in% unique(multiple_datasets_object()@meta.data$dataset)
        subset_seurat <- subset(multiple_datasets_object(), subset = dataset %in% input$dataset_filter[valid_datasets])

        # Afficher les clusters filtrés
        DimPlot(subset_seurat, group.by = "ident", label = show_labels, label.size = 5) +
          scale_color_manual(values = cluster_colors) +  # Appliquer les couleurs personnalisées
          NoAxes() +             # Supprimer les axes
          NoLegend() +           # Supprimer la légende
          ggtitle(NULL)          # Supprimer le titre
      }

      filtered_umap_plot(plot)  # Stocker le graphique dans la variable réactive (si nécessaire)
      plot  # Afficher le graphique
    })

    output$download_filtered_umap_plot <- downloadHandler(
      filename = function() {
        paste("UMAP_plot", Sys.Date(), ".tiff", sep = "")
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


    # UI pour sélectionner les clusters
    output$cluster1_compare_ui <- renderUI({
      req(multiple_datasets_object())
      selectInput("cluster1_compare",
                  label = "Select clusters for group 1",
                  choices = unique(Idents(multiple_datasets_object())),
                  selected = NULL,
                  multiple = TRUE)  # Permettre la sélection multiple
    })

    output$cluster2_compare_ui <- renderUI({
      req(multiple_datasets_object())
      # Exclure les clusters déjà sélectionnés dans group 1
      remaining_clusters <- setdiff(unique(Idents(multiple_datasets_object())), input$cluster1_compare)
      selectInput("cluster2_compare",
                  label = "Select clusters for group 2",
                  choices = remaining_clusters,
                  selected = NULL,
                  multiple = TRUE)  # Permettre la sélection multiple
    })

    # Comparaison des clusters
    observeEvent(input$compare_clusters_button, {
      tryCatch({
        showModal(modalDialog(
          title = "Please Wait",
          "Finding differentially expressed features...",
          easyClose = FALSE,
          footer = NULL
        ))

        req(input$cluster1_compare, input$cluster2_compare, multiple_datasets_object())

        # Vérifier qu'il n'y a pas de chevauchement entre les groupes
        if(any(input$cluster1_compare %in% input$cluster2_compare)) {
          showNotification("Groups must not have overlapping clusters!", type = "error")
          removeModal()
          return()
        }

        seurat_obj <- multiple_datasets_object()

        # Créer une nouvelle identité pour les groupes de comparaison
        new_idents <- as.character(Idents(seurat_obj))
        new_idents[new_idents %in% input$cluster1_compare] <- "group1"
        new_idents[new_idents %in% input$cluster2_compare] <- "group2"
        Idents(seurat_obj) <- new_idents

        # Trouver les marqueurs
        temp_res <- FindMarkers(seurat_obj,
                                ident.1 = "group1",
                                ident.2 = "group2",
                                min.pct = input$min_pct_compare_merge,
                                logfc.threshold = input$logfc_threshold_compare_merge)

        # Formater les résultats
        temp_res$p_val <- format(temp_res$p_val, scientific = TRUE, digits = 3)
        temp_res$p_val_adj <- format(temp_res$p_val_adj, scientific = TRUE, digits = 3)
        temp_res$gene <- paste0('<a href="#" class="gene-name" data-gene="', rownames(temp_res), '">', rownames(temp_res), '</a>')

        # Ajouter l'information sur les groupes comparés
        temp_res$comparison <- sprintf("Group1(%s) vs Group2(%s)",
                                       paste(input$cluster1_compare, collapse=","),
                                       paste(input$cluster2_compare, collapse=","))

        shinyjs::enable("download_markers_multiple_clusters_merge")
        diff_genes_compare_cluster(temp_res)
        removeModal()

      }, error = function(e) {
        removeModal()
        showNotification(paste("Cluster comparison error:", e$message), type = "error")
      })
    })

    # Observer pour mettre à jour les choix du deuxième groupe
    observeEvent(input$cluster1_compare, {
      req(multiple_datasets_object())
      remaining_clusters <- setdiff(unique(Idents(multiple_datasets_object())), input$cluster1_compare)
      updateSelectInput(session, "cluster2_compare",
                        choices = remaining_clusters)
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


    ################### Comparison of multiples clusters betwin multiples datasets

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

    # Observer pour la comparaison
    observeEvent(input$compare_datasets_button, {
      tryCatch({
        showModal(modalDialog(
          title = "Please Wait",
          "Comparing datasets and clusters...",
          easyClose = FALSE,
          footer = NULL
        ))

        req(input$dataset1_compare, input$dataset2_compare, multiple_datasets_object())

        # Vérifier qu'il n'y a pas de chevauchement entre les groupes de datasets
        if(any(input$dataset1_compare %in% input$dataset2_compare)) {
          showNotification("Groups must not have overlapping datasets!", type = "error")
          removeModal()
          return()
        }

        seurat_obj <- multiple_datasets_object()

        # Créer une nouvelle colonne pour identifier les combinaisons dataset-cluster
        seurat_obj$dataset_cluster <- paste(seurat_obj@meta.data$dataset,
                                            Idents(seurat_obj),
                                            sep = "_")

        # Créer les identifiants des groupes
        group1_combinations <- expand.grid(
          dataset = input$dataset1_compare,
          cluster = if(input$all_clusters) unique(Idents(seurat_obj)) else input$cluster_compare,
          stringsAsFactors = FALSE
        )
        group2_combinations <- expand.grid(
          dataset = input$dataset2_compare,
          cluster = if(input$all_clusters) unique(Idents(seurat_obj)) else input$cluster_compare,
          stringsAsFactors = FALSE
        )

        # Créer les identifiants de groupe
        group1_ids <- apply(group1_combinations, 1, function(x) paste(x[1], x[2], sep = "_"))
        group2_ids <- apply(group2_combinations, 1, function(x) paste(x[1], x[2], sep = "_"))

        # Trouver les marqueurs
        temp_res <- FindMarkers(
          object = seurat_obj,
          ident.1 = group1_ids,
          ident.2 = group2_ids,
          group.by = "dataset_cluster",
          min.pct = input$min_pct_compare_dataset_merge,
          logfc.threshold = input$logfc_threshold_datasets
        )

        # Formater les résultats
        temp_res$p_val <- format(temp_res$p_val, scientific = TRUE, digits = 3)
        temp_res$p_val_adj <- format(temp_res$p_val_adj, scientific = TRUE, digits = 3)
        temp_res$gene <- paste0('<a href="#" class="gene-name" data-gene="', rownames(temp_res), '">', rownames(temp_res), '</a>')

        # Ajouter les informations de comparaison
        temp_res$comparison <- sprintf(
          "Group1(Datasets: %s, Clusters: %s) vs Group2(Datasets: %s, Clusters: %s)",
          paste(input$dataset1_compare, collapse=","),
          paste(if(input$all_clusters) "all" else input$cluster_compare, collapse=","),
          paste(input$dataset2_compare, collapse=","),
          paste(if(input$all_clusters) "all" else input$cluster_compare, collapse=",")
        )

        diff_genes_compare_datasets(temp_res)

        if (nrow(temp_res) == 0) {
          showNotification("No differentially expressed genes found!", type = "error")
          removeModal()
          return()
        }

        output$diff_dataset_cluster <- renderDataTable({
          datatable(diff_genes_compare_datasets(), escape = FALSE)
        })

        shinyjs::enable("download_diff_dataset_cluster")
        removeModal()

      }, error = function(e) {
        removeModal()
        showNotification(paste("Error in comparison:", e$message), type = "error")
      })
    })

    # Observer pour mettre à jour les choix du deuxième groupe de datasets
    observeEvent(input$dataset1_compare, {
      req(multiple_datasets_object())
      remaining_datasets <- setdiff(unique(multiple_datasets_object()@meta.data$dataset),
                                    input$dataset1_compare)
      updateSelectInput(session, "dataset2_compare",
                        choices = remaining_datasets)
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




    # Function to generate the cluster table with counts and proportions
    generate_cluster_table <- function() {
      req(multiple_datasets_object())
      seurat_object <- multiple_datasets_object()

      # Get absolute counts
      cluster_data <- table(Idents(seurat_object), seurat_object$orig.ident)
      cluster_df <- as.data.frame.matrix(cluster_data)

      # Calculate total nuclei per dataset
      dataset_totals <- colSums(cluster_data)

      # Calculate proportions
      proportion_df <- sweep(cluster_data, 2, dataset_totals, "/") * 100
      proportion_df <- as.data.frame.matrix(proportion_df)

      # Combine counts and proportions
      result_df <- data.frame(Cluster = rownames(cluster_df))

      # Add counts and proportions for each dataset
      for(dataset in colnames(cluster_df)) {
        result_df[[paste0(dataset, "_count")]] <- cluster_df[[dataset]]
        result_df[[paste0(dataset, "_percent")]] <- round(proportion_df[[dataset]], 2)
      }

      return(result_df)
    }

    # Update the observer
    observeEvent(input$generate_cluster_table, {
      output$cluster_table <- renderDataTable({
        df <- generate_cluster_table()
        # Format the table to show counts and percentages
        datatable(df,
                  options = list(
                    pageLength = 25,
                    scrollX = TRUE
                  ),
                  caption = htmltools::tags$caption(
                    style = 'caption-side: bottom; text-align: center;',
                    'Count and percentage of nuclei per cluster and dataset'
                  )
        ) %>%
          formatStyle(
            columns = grep("_percent$", colnames(df)),
            suffix = "%"
          )
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
