############################## MultiNicheNet Server mise à jour ##############################

multinichenet_server <- function(input, output, session) {

  ############################## NicheNet Server ##############################
  # Reactive values for storing NicheNet data
  data_nichenet <- reactiveValues(
    seuratObj = NULL,
    lr_network = NULL,
    ligand_target_matrix = NULL,
    weighted_networks = NULL,
    nichenet_output = NULL
  )

  # Modal to load data for NicheNet
  observeEvent(input$open_modal_nichenet, {
    showModal(modalDialog(
      title = "Load Dataset and Ligand-Receptor Data for NicheNet",
      easyClose = FALSE,
      footer = tagList(
        modalButton("Cancel"),
        actionButton("start_load_nichenet", "Load Data")
      ),
      fluidRow(
        column(6, fileInput("seurat_file_nichenet", "Choose Seurat Object (.rds)", accept = ".rds")),
        column(6, fileInput("lr_network_file_nichenet", "Upload Ligand-Receptor Network (.rds)", accept = ".rds")),
        column(6, fileInput("ligand_target_matrix_file_nichenet", "Upload Ligand-Target Matrix (.rds)", accept = ".rds")),
        column(6, fileInput("weighted_networks_file_nichenet", "Upload Weighted Networks (.rds)", accept = ".rds"))
      )
    ))
  })



      # Load data and check columns
      observeEvent(input$start_load_nichenet, {
        showModal(modalDialog(
          title = "Loading data...",
          "Please wait while the data is being loaded.",
          footer = NULL,
          easyClose = FALSE
        ))

        tryCatch({
          # Load Seurat object
          data_nichenet$seuratObj <- readRDS(input$seurat_file_nichenet$datapath)

          # Load other necessary files
          data_nichenet$lr_network <- readRDS(input$lr_network_file_nichenet$datapath)
          data_nichenet$ligand_target_matrix <- readRDS(input$ligand_target_matrix_file_nichenet$datapath)
          data_nichenet$weighted_networks <- readRDS(input$weighted_networks_file_nichenet$datapath)

          # Check the column names in lr_network
          print("Checking columns in lr_network:")
          print(colnames(data_nichenet$lr_network))

          # Only rename columns if 'from' and 'to' are present, otherwise assume 'ligand' and 'receptor' are correct
          if ("from" %in% colnames(data_nichenet$lr_network) & "to" %in% colnames(data_nichenet$lr_network)) {
            data_nichenet$lr_network <- data_nichenet$lr_network %>%
              dplyr::rename(ligand = from, receptor = to)
          }

          # Ensure ligand and receptor names are valid R names
          data_nichenet$lr_network <- data_nichenet$lr_network %>%
            distinct(ligand, receptor) %>%
            mutate(ligand = make.names(ligand), receptor = make.names(receptor))

          # Apply make.names() to ligand_target_matrix
          showNotification("Checking and renaming ligand-target matrix...", type = "message")
          colnames(data_nichenet$ligand_target_matrix) <- colnames(data_nichenet$ligand_target_matrix) %>% make.names()
          rownames(data_nichenet$ligand_target_matrix) <- rownames(data_nichenet$ligand_target_matrix) %>% make.names()

          # Filter ligands in lr_network based on ligand_target_matrix
          valid_ligands <- colnames(data_nichenet$ligand_target_matrix)
          data_nichenet$lr_network <- data_nichenet$lr_network %>%
            dplyr::filter(ligand %in% valid_ligands)

          # Check if there are any valid ligands left
          if (nrow(data_nichenet$lr_network) == 0) {
            stop("No valid ligands in lr_network after filtering with ligand_target_matrix")
          }
          # If the columns 'from' and 'to' are missing, rename 'ligand' and 'receptor'
          if (!"from" %in% colnames(data_nichenet$lr_network)) {
            data_nichenet$lr_network <- data_nichenet$lr_network %>%
              dplyr::rename(from = ligand, to = receptor)
          }

          # Convert Seurat object and adjust gene names for mouse
          showNotification("Updating Seurat object...", type = "message")
          data_nichenet$seuratObj <- UpdateSeuratObject(data_nichenet$seuratObj)
          data_nichenet$seuratObj <- alias_to_symbol_seurat(data_nichenet$seuratObj, "mouse")
          print("Seurat object updated and gene names converted for mouse.")

          # Retrieve metadata columns for selection inputs
          columns_metadata <- colnames(data_nichenet$seuratObj@meta.data)

          # Update selectInput for cell identity column (senders and receivers)
          updateSelectInput(session, "cell_identity_column", choices = columns_metadata)

          # Update selectInput for condition column
          updateSelectInput(session, "condition_column", choices = columns_metadata)

          removeModal()  # Remove loading modal after data is loaded
          showNotification("Data loaded successfully!", type = "message")

        }, error = function(e) {
          removeModal()  # Remove modal if an error occurs
          # Show full error message
          showNotification(paste("Error loading NicheNet data:", e$message), type = "error")
        })
      })


  # Dynamically update select inputs for sender and receiver based on selected cell identity column
      # Dynamically update select inputs for sender and receiver based on selected cell identity column
      observeEvent(input$cell_identity_column, {
        req(input$cell_identity_column)

        # Ensure the Seurat object and metadata are loaded
        req(data_nichenet$seuratObj)

        # Get unique cell identities from the selected column
        cell_types <- unique(data_nichenet$seuratObj@meta.data[[input$cell_identity_column]])

        # Update sender and receiver select inputs to allow multiple selections
        updateSelectizeInput(session, "sender_select", choices = cell_types, server = TRUE)
        updateSelectizeInput(session, "receiver_select", choices = cell_types, server = TRUE)
      })


  # Dynamically update select inputs for condition_oi and condition_reference based on selected condition column
  observeEvent(input$condition_column, {
    req(input$condition_column)

    # Ensure the Seurat object and metadata are loaded
    req(data_nichenet$seuratObj)

    # Get unique conditions from the selected column
    conditions <- unique(data_nichenet$seuratObj@meta.data[[input$condition_column]])

    # Update select inputs for condition_oi and condition_reference
    updateSelectInput(session, "condition_oi_select", choices = conditions)
    updateSelectInput(session, "condition_reference_select", choices = conditions)
  })

  # Run NicheNet analysis when all parameters are defined
  observeEvent(input$run_nichenet, {
    req(input$sender_select,
        input$receiver_select,
        input$condition_oi_select,
        input$condition_reference_select,
        input$expression_pct,
        input$top_n_ligands_nichenet)

    # Modal pour indiquer que l'analyse est en cours
    showModal(modalDialog(
      title = "Running NicheNet analysis...",
      "Please wait while the analysis is being performed.",
      footer = NULL,
      easyClose = FALSE
    ))

    tryCatch({
      # Set the RNA assay as the default assay for the Seurat object
      DefaultAssay(data_nichenet$seuratObj) <- "RNA"

      # Logs for tracking execution
      showNotification("Starting NicheNet analysis...", type = "message")

      # Run the NicheNet analysis with the provided parameters
      data_nichenet$nichenet_output <- nichenet_seuratobj_aggregate(
        seurat_obj = data_nichenet$seuratObj,
        sender = input$sender_select,
        receiver = input$receiver_select,
        condition_colname = input$condition_column,
        condition_oi = input$condition_oi_select,
        condition_reference = input$condition_reference_select,
        expression_pct = input$expression_pct,
        top_n_ligands = input$top_n_ligands_nichenet,  # Nouveau paramètre
        ligand_target_matrix = data_nichenet$ligand_target_matrix,
        lr_network = data_nichenet$lr_network,
        weighted_networks = data_nichenet$weighted_networks
      )

      # Logs during execution
      showNotification("NicheNet analysis in progress...", type = "message")

      removeModal()  # Remove modal after analysis completion
      showNotification("NicheNet analysis completed!", type = "message")

    }, error = function(e) {
      removeModal()  # Remove modal if an error occurs
      showNotification(paste("Error during NicheNet analysis:", e$message), type = "error")
    })
  })

  # Display NicheNet results
  observeEvent(input$view_nichenet_results, {
    req(data_nichenet$nichenet_output)

    output$nichenet_output_display <- renderUI({
      output_type <- input$nichenet_output_type

      # Créer une liste pour stocker les outputs
      output_list <- list()

      # Afficher le contenu selon le type
      output_list[["content"]] <- if (output_type == "ligand_activities") {
        renderTable({
          data_nichenet$nichenet_output$ligand_activities
        })
      } else if (output_type == "top_ligands") {
        renderText({
          paste(data_nichenet$nichenet_output$top_ligands, collapse = ", ")
        })
      } else if (output_type == "top_targets") {
        renderText({
          paste(data_nichenet$nichenet_output$top_targets, collapse = ", ")
        })
      } else if (output_type %in% c("ligand_expression_dotplot",
                                    "ligand_differential_expression_heatmap",
                                    "ligand_target_heatmap",
                                    "ligand_receptor_heatmap")) {
        renderPlot({
          data_nichenet$nichenet_output[[output_type]]
        })
      } else if (output_type %in% c("ligand_target_matrix", "ligand_target_df")) {
        renderTable({
          data_nichenet$nichenet_output[[output_type]]
        })
      }

      # Retourner la liste d'outputs
      do.call(tagList, output_list)
    })

    # Gestionnaire de téléchargement simplifié
    output$download_nichenet_result <- downloadHandler(
      filename = function() {
        output_type <- input$nichenet_output_type
        date_str <- format(Sys.time(), "%Y%m%d_%H%M%S")

        # Extension selon le type de sortie
        ext <- if (output_type %in% c("ligand_expression_dotplot",
                                      "ligand_differential_expression_heatmap",
                                      "ligand_target_heatmap",
                                      "ligand_receptor_heatmap")) {
          "tiff"
        } else if (output_type %in% c("ligand_activities",
                                      "ligand_target_matrix",
                                      "ligand_target_df")) {
          "csv"
        } else {
          "txt"
        }

        paste0("nichenet_", output_type, "_", date_str, ".", ext)
      },
      content = function(file) {
        output_type <- input$nichenet_output_type
        result <- data_nichenet$nichenet_output[[output_type]]

        tryCatch({
          if (output_type %in% c("ligand_expression_dotplot",
                                 "ligand_differential_expression_heatmap",
                                 "ligand_target_heatmap",
                                 "ligand_receptor_heatmap")) {
            # Export en TIFF pour les plots
            width <- 8
            height <- 6
            dpi <- as.numeric(input$plot_dpi)

            tiff(file,
                 width = width,
                 height = height,
                 units = 'in',
                 res = dpi,
                 compression = "lzw")
            print(result)
            dev.off()

          } else if (output_type %in% c("ligand_activities",
                                        "ligand_target_matrix",
                                        "ligand_target_df")) {
            # Export en CSV pour les tableaux
            write.csv(result, file, row.names = FALSE)

          } else {
            # Export en texte pour les listes simples
            writeLines(paste(result, collapse = ", "), file)
          }
        }, error = function(e) {
          if (dev.cur() > 1) dev.off()
          cat("Error in export:", e$message, "\n")
        })
      }
    )
  })

  ############################## Circos Plot Server ##############################
  # Step 1: Assign ligands to sender cells
  observeEvent(input$assign_ligands, {
    req(data_nichenet$seuratObj, input$cell_identity_column, data_nichenet$nichenet_output)

    tryCatch({
      # Extract the chosen cell type column from user input
      celltype_col <- input$cell_identity_column

      # Convert the values in the selected cell type column to character for proper comparisons
      data_nichenet$seuratObj@meta.data[[celltype_col]] <- as.character(data_nichenet$seuratObj@meta.data[[celltype_col]])

      # Assign ligands to sender cells based on expression and function settings
      ligand_type_indication_df <- assign_ligands_to_celltype(
        seuratObj = data_nichenet$seuratObj,
        ligands = data_nichenet$nichenet_output$top_ligands,
        celltype_col = celltype_col,
        func.agg = mean,
        func.assign = function(x) mean(x) + 0.5 * sd(x)
      )

      # Ensure columns are correctly named for downstream functions
      colnames(ligand_type_indication_df) <- c("ligand_type", "ligand")

      # Log the ligand assignment summary to help with debugging
      cat("Log: Ligand assignment results (first few rows):\n")
      print(head(ligand_type_indication_df))
      cat("Log: Ligand assignment summary:\n")
      print(table(ligand_type_indication_df$ligand_type))

      # Store the resulting DataFrame for further steps
      data_nichenet$ligand_type_indication_df <- ligand_type_indication_df
      showNotification("Ligands assigned to sender cells successfully.", type = "message")
    }, error = function(e) {
      showNotification(paste("Error during ligand assignment:", e$message), type = "error")
    })
  })

  # Step 2: Define ligand-target links of interest
  observeEvent(input$define_links, {
    req(data_nichenet$ligand_type_indication_df, data_nichenet$nichenet_output)

    tryCatch({
      cat("Log: Starting to define ligand-target links...\n")

      # Retrieve ligand-type data and ensure columns are as expected
      ligand_type_indication_df <- data_nichenet$ligand_type_indication_df
      cat("Log: Ligand-type indication dataframe (first few rows):\n")
      print(head(ligand_type_indication_df))

      # Retrieve ligand-target links data frame from the NicheNet output
      active_ligand_target_links_df <- data_nichenet$nichenet_output$ligand_target_df
      active_ligand_target_links_df$target_type <- input$target_type

      cat("Log: Original active_ligand_target_links_df (first few rows):\n")
      print(head(active_ligand_target_links_df))

      # Join ligand_type_indication_df with active_ligand_target_links_df on the 'ligand' column
      circos_links <- active_ligand_target_links_df %>%
        inner_join(ligand_type_indication_df, by = "ligand")

      cat("Log: Circos links after join (first few rows):\n")
      print(head(circos_links))

      # Verify that the resulting dataframe has the expected columns
      if (!all(c("ligand", "target", "weight", "target_type", "ligand_type") %in% colnames(circos_links))) {
        stop("The resulting dataframe is missing required columns.")
      }

      # Store the resulting circos links for visualization
      data_nichenet$circos_links <- circos_links
      showNotification("Ligand-target links defined successfully.", type = "message")
    }, error = function(e) {
      showNotification(paste("Error defining ligand-target links:", e$message), type = "error")
      cat("Error defining ligand-target links:", e$message, "\n")
    })
  })




  # Step 3: Generate dynamic inputs for customizing ligand and target colors
  # Générer les inputs de couleur pour les donneurs
  output$donor_color_inputs <- renderUI({
    req(input$sender_select)

    # Créer un color picker pour chaque donneur sélectionné
    tagList(
      h4("Donor Cell Colors"),
      lapply(input$sender_select, function(donor) {
        colourInput(
          inputId = paste0("donor_color_", make.names(donor)),
          label = paste("Color for", donor),
          value = sample(rainbow(length(input$sender_select)), 1)  # Couleur aléatoire par défaut
        )
      })
    )
  })

  # Générer l'input de couleur pour le receveur
  output$receiver_color_input <- renderUI({
    req(input$receiver_select)

    tagList(
      h4("Receiver Cell Color"),
      colourInput(
        inputId = "receiver_color",
        label = paste("Color for receiver:", input$receiver_select),
        value = "#999999"
      )
    )
  })


  # Création d'une valeur réactive pour vis_circos_obj
  vis_circos_obj <- reactiveVal(NULL)

  # Observateur pour la préparation des données du Circos Plot
  observeEvent(input$draw_circos_plot, {
    req(data_nichenet$circos_links)

    tryCatch({
      circos.clear()

      # Obtenir les couleurs des ligands à partir des entrées utilisateur
      ligand_colors <- sapply(unique(data_nichenet$circos_links$ligand_type), function(type) {
        input[[paste0("ligand_color_", type)]]
      })
      names(ligand_colors) <- unique(data_nichenet$circos_links$ligand_type)

      # Obtenir la couleur de la cible à partir de l'entrée utilisateur
      target_colors <- sapply(unique(data_nichenet$circos_links$target_type), function(type) {
        input[[paste0("target_color_", type)]]
      })
      names(target_colors) <- unique(data_nichenet$circos_links$target_type)

      cat("Log: Ligand colors:\n")
      print(ligand_colors)
      cat("Log: Target colors:\n")
      print(target_colors)

      # Préparer l'objet de visualisation circos
      circos_obj <- prepare_circos_visualization(
        data_nichenet$circos_links,
        ligand_colors = ligand_colors,
        target_colors = target_colors
      )

      # Stocker l'objet dans la valeur réactive
      vis_circos_obj(circos_obj)

      cat("Log: vis_circos_obj structure:\n")
      str(vis_circos_obj())

      showNotification("Circos Plot data prepared successfully!", type = "message")
    }, error = function(e) {
      showNotification(paste("Error preparing Circos Plot data:", e$message), type = "error")
      cat("Error preparing Circos Plot data:", e$message, "\n")
    })
  })
  make_circos_plot <- function(vis_circos_obj,
                               transparency = FALSE,
                               top_n_ligands = 30,
                               top_n_receptors = 30,
                               link_size = 1,
                               text_size = 0.6,
                               text_position = "Outside",
                               sort_mode = "Descending",
                               gap_degree = 1,
                               donor_colors = NULL,
                               receiver_color = NULL) {
    tryCatch({
      req(vis_circos_obj$links_circle)
      req(donor_colors)
      req(receiver_color)

      circos.clear()

      # Nettoyage des données initiales
      links_df <- vis_circos_obj$links_circle %>%
        dplyr::mutate(
          ligand = trimws(ligand),
          target = as.character(target),
          weight = as.numeric(weight)
        )

      # Joindre les informations de type
      links_df <- links_df %>%
        dplyr::left_join(
          data_nichenet$ligand_type_indication_df %>%
            dplyr::mutate(ligand = trimws(ligand)),
          by = "ligand"
        )

      # Calculer les poids totaux par ligand et type
      ligand_weights <- links_df %>%
        dplyr::group_by(ligand, ligand_type) %>%
        dplyr::summarise(total_weight = sum(weight, na.rm = TRUE), .groups = 'drop') %>%
        dplyr::arrange(ligand_type, desc(total_weight))

      # Sélectionner les top ligands de manière équilibrée
      n_types <- length(unique(ligand_weights$ligand_type))
      ligands_per_type <- ceiling(top_n_ligands / n_types)

      top_ligands <- ligand_weights %>%
        dplyr::group_by(ligand_type) %>%
        dplyr::slice_max(order_by = total_weight, n = ligands_per_type) %>%
        dplyr::ungroup() %>%
        dplyr::slice_max(order_by = total_weight, n = top_n_ligands) %>%
        dplyr::pull(ligand)

      # Sélectionner les top récepteurs
      target_weights <- links_df %>%
        dplyr::group_by(target) %>%
        dplyr::summarise(total_weight = sum(weight, na.rm = TRUE)) %>%
        dplyr::arrange(desc(total_weight))

      top_targets <- target_weights$target[1:min(top_n_receptors, nrow(target_weights))]

      # Filtrer les données
      links_df <- links_df %>%
        dplyr::filter(ligand %in% top_ligands & target %in% top_targets)

      # Normaliser les poids
      links_df <- links_df %>%
        dplyr::mutate(
          weight = (weight - min(weight)) / (max(weight) - min(weight)),
          weight = weight * 0.9 + 0.1
        )

      # Information sur les ligands et leurs types
      ligand_info <- links_df %>%
        dplyr::select(ligand, ligand_type) %>%
        dplyr::distinct()

      # Assigner les couleurs aux ligands
      ligand_colors <- sapply(top_ligands, function(l) {
        l_type <- ligand_info$ligand_type[ligand_info$ligand == l]
        if(length(l_type) > 0 && !is.na(l_type) && l_type %in% names(donor_colors)) {
          donor_colors[l_type]
        } else {
          "#808080"  # Couleur par défaut
        }
      })
      names(ligand_colors) <- top_ligands

      # Couleurs pour les récepteurs
      target_colors <- rep(receiver_color, length(top_targets))
      names(target_colors) <- top_targets

      # Combiner les couleurs
      color_data_circos <- c(ligand_colors, target_colors)

      # Créer le vecteur de groupes
      groups <- c(
        setNames(ligand_info$ligand_type[match(names(ligand_colors), ligand_info$ligand)],
                 names(ligand_colors)),
        setNames(rep("Receiver", length(target_colors)),
                 names(target_colors))
      )

      # Configuration du circos
      circos.par(gap.degree = as.numeric(gap_degree),
                 start.degree = 90)

      # Création du diagramme
      chordDiagram(
        links_df %>% dplyr::select(ligand, target, weight),
        grid.col = color_data_circos,
        transparency = ifelse(transparency, 0.7, 0.3),
        directional = 1,
        direction.type = c("diffHeight", "arrows"),
        link.arr.type = "big.arrow",
        link.sort = TRUE,
        link.decreasing = (sort_mode == "Descending"),
        link.lwd = link_size * 2,
        annotationTrack = "grid",
        preAllocateTracks = list(track.height = max(strwidth(unlist(dimnames(links_df))))),
        annotationTrackHeight = c(0.05, 0.1),
        group = groups
      )

      # Ajout des étiquettes
      if (text_position %in% c("Inside", "Both", "Outside")) {
        circos.track(track.index = 1, panel.fun = function(x, y) {
          sector.index <- get.cell.meta.data("sector.index")
          xlim <- get.cell.meta.data("xlim")
          ylim <- get.cell.meta.data("ylim")

          if (text_position == "Inside" || text_position == "Both") {
            circos.text(mean(xlim), ylim[1], sector.index,
                        facing = "clockwise", niceFacing = TRUE,
                        adj = c(0, 0.5), cex = text_size,
                        col = ifelse(sector.index %in% names(ligand_colors), "black", "darkgray"))
          }
          if (text_position == "Outside" || text_position == "Both") {
            circos.text(mean(xlim), ylim[2], sector.index,
                        facing = "clockwise", niceFacing = TRUE,
                        adj = c(0, 0.5), cex = text_size,
                        col = ifelse(sector.index %in% names(ligand_colors), "black", "darkgray"))
          }
        }, bg.border = NA)
      }

      receiver_type <- unique(links_df$target_type)[1]

      # Ajout de la légende avec le vrai nom du receveur
      legend_colors <- c(donor_colors, setNames(receiver_color, receiver_type))
      legend("bottomright",
             legend = names(legend_colors),
             fill = unname(legend_colors),
             border = NA,
             bty = "n",
             cex = 0.8,
             title = "Cell Types")

    }, error = function(e) {
      cat("Error in make_circos_plot:", e$message, "\n")
      plot.new()
      text(0.5, 0.5, paste("Error generating Circos Plot:", e$message), cex = 1.2)
    })
  }
  # Rendu du plot
  output$circos_plot_output <- renderPlot({
    req(input$draw_circos_plot)
    req(vis_circos_obj())

    par(mar = c(2, 2, 2, 2))
    layout(matrix(1, 1, 1))

    # Récupérer les couleurs des donneurs
    donor_colors <- sapply(input$sender_select, function(donor) {
      input[[paste0("donor_color_", make.names(donor))]]
    })
    names(donor_colors) <- input$sender_select

    tryCatch({
      make_circos_plot(
        vis_circos_obj(),
        transparency = as.logical(input$transparency_circos),
        top_n_ligands = length(unique(vis_circos_obj()$links_circle$ligand)),
        top_n_receptors = as.integer(input$top_n_receptors),
        text_size = as.numeric(input$text_size %||% 0.6),
        text_position = input$text_position %||% "Outside",
        sort_mode = input$sort_mode %||% "Descending",
        gap_degree = as.numeric(input$gap_degree %||% 1),
        link_size = as.numeric(input$link_size %||% 1),
        donor_colors = donor_colors,
        receiver_color = input$receiver_color %||% "#999999"
      )
    }, error = function(e) {
      cat("Error generating Circos Plot:", e$message, "\n")
      plot.new()
      text(0.5, 0.5, paste("Error generating Circos Plot:", e$message), cex = 1.2)
    })
  }, res = 96,
  height = function() {
    w <- session$clientData$output_circos_plot_output_width
    return(w)
  },
  width = function() {
    return(session$clientData$output_circos_plot_output_width)
  })

  # Handler de téléchargement
  output$download_circo_plot <- downloadHandler(
    filename = function() {
      paste("circos_plot_", format(Sys.time(), "%Y%m%d"), ".tiff", sep = "")
    },
    content = function(file) {
      tryCatch({
        tiff(file,
             width = 10,
             height = 8,
             units = 'in',
             res = input$plot_dpi,
             compression = "lzw")

        par(mar = c(2, 2, 2, 2))

        # Récupérer les couleurs des donneurs
        donor_colors <- sapply(input$sender_select, function(donor) {
          input[[paste0("donor_color_", make.names(donor))]]
        })
        names(donor_colors) <- input$sender_select

        make_circos_plot(
          vis_circos_obj(),
          transparency = as.logical(input$transparency_circos),
          top_n_ligands = length(unique(vis_circos_obj()$links_circle$ligand)),
          top_n_receptors = as.integer(input$top_n_receptors),
          text_size = as.numeric(input$text_size %||% 0.6),
          text_position = input$text_position %||% "Outside",
          sort_mode = input$sort_mode %||% "Descending",
          gap_degree = as.numeric(input$gap_degree %||% 1),
          link_size = as.numeric(input$link_size %||% 1),
          donor_colors = donor_colors,
          receiver_color = input$receiver_color %||% "#999999"
        )

        dev.off()
      }, error = function(e) {
        if (dev.cur() > 1) dev.off()
        cat("Error saving plot:", e$message, "\n")
      })
    }
  )


##################### Multinichenet #################################



  # Reactive values pour stocker le dataset et autres données
  data <- reactiveValues(
    sce = NULL,
    metadata = NULL,
    lr_network = NULL,
    ligand_target_matrix = NULL,
    contrast_tbl = NULL,
    contrasts_oi = NULL,
    senders_oi = NULL,
    receivers_oi = NULL,
    frq_list = NULL,
    abundance_info = NULL,
    genes_oi = NULL,
    sample_id = NULL,
    group_id = NULL,
    selected_groups = NULL,  # Stocker les groupes sélectionnés
    celltype_id = NULL
  )

  metadata_columns <- reactiveVal(NULL)  # Stocke les colonnes des métadonnées

  group_levels <- reactiveVal(NULL)  # Stocke les niveaux de groupes pour les contrastes
  contrast_formula <- reactiveVal(NULL)  # Stocke la formule de contraste générée

  ############################### BOÎTE MODALE BLOQUANTE ##################################

  show_loading_modal <- function(message) {
    showModal(modalDialog(
      title = "Processing...",
      paste(message, "Please wait."),
      easyClose = FALSE,
      footer = NULL
    ))
  }

  hide_loading_modal <- function() {
    removeModal()
  }

  ############################### ÉTAPE 1: CHARGEMENT DES DONNÉES ##################################


  # Ouvrir une boîte modale pour sélectionner le fichier et les données de ligands-récepteurs
  observeEvent(input$open_modal, {
    showModal(modalDialog(
      title = "Load Dataset and Ligand-Receptor Data",
      easyClose = FALSE,
      footer = tagList(
        modalButton("Cancel"),
        actionButton("start_load", "Load Data")
      ),
      fluidRow(
        column(6, fileInput("sce_file", "Choose Seurat or SingleCellExperiment Object (.rds)", accept = ".rds")),
        column(12, h4("Upload Ligand-Receptor Network and Ligand-Target Matrix Files")),
        column(6, fileInput("lr_network_file", "Upload Ligand-Receptor Network (.rds)", accept = ".rds")),
        column(6, fileInput("ligand_target_matrix_file", "Upload Ligand-Target Matrix (.rds)", accept = ".rds"))
      )
    ))
  })
  # Fonction pour rendre les niveaux des colonnes de colData valides
  sanitize_colData <- function(sce) {
    col_data <- colData(sce)

    for (col_name in colnames(col_data)) {
      # Appliquer make.names aux niveaux des colonnes qui sont des facteurs
      if (is.factor(col_data[[col_name]])) {
        levels(col_data[[col_name]]) <- make.names(levels(col_data[[col_name]]), unique = TRUE)
      }

      # Si ce n'est pas un facteur mais une chaîne de caractères, convertir en facteur et appliquer make.names
      if (is.character(col_data[[col_name]])) {
        col_data[[col_name]] <- factor(col_data[[col_name]])
        levels(col_data[[col_name]]) <- make.names(levels(col_data[[col_name]]), unique = TRUE)
      }
    }

    colData(sce) <- col_data  # Mettre à jour colData avec les noms nettoyés
    return(sce)
  }

  # Observer pour démarrer le chargement des données
  observeEvent(input$start_load, {
    # Valider que l'utilisateur a bien sélectionné un fichier
    if (is.null(input$sce_file)) {
      showNotification("Please select a dataset!", type = "error")
      return(NULL)
    }

    # Valider que l'utilisateur a bien sélectionné les fichiers Ligand-Receptor et Ligand-Target Matrix
    if (is.null(input$lr_network_file) || is.null(input$ligand_target_matrix_file)) {
      showNotification("Please upload the ligand-receptor network and ligand-target matrix!", type = "error")
      return(NULL)
    }

    # Afficher une boîte modale pendant le chargement
    show_loading_modal("Loading Dataset")

    tryCatch({
      # Étape 1 : Charger l'objet Seurat ou SingleCellExperiment
      sce_data <- readRDS(input$sce_file$datapath)

      # Vérifier si c'est un objet Seurat et le convertir en SingleCellExperiment
      if ("Seurat" %in% class(sce_data)) {
        print("Converting Seurat object to SingleCellExperiment...")
        sce_data <- Seurat::as.SingleCellExperiment(sce_data, assay = "RNA")

        # Ajouter conversion des noms de gènes
        sce_data <- alias_to_symbol_SCE(sce_data, "mouse") %>% makenames_SCE()

        data$metadata <- colData(sce_data)  # Extraire les métadonnées après conversion
        showNotification("Seurat object converted to SingleCellExperiment successfully!", type = "message")
      } else if ("SingleCellExperiment" %in% class(sce_data)) {
        # Ajouter conversion des noms de gènes
        sce_data <- alias_to_symbol_SCE(sce_data, "mouse") %>% makenames_SCE()
        data$metadata <- colData(sce_data)  # Extraire les métadonnées directement
      } else {
        stop("Unknown dataset format. Please provide a Seurat or SingleCellExperiment object.")
      }

      # Appliquer make.names à toutes les colonnes de colData qui sont des facteurs ou des caractères
      sce_data <- sanitize_colData(sce_data)

      # Stocker l'objet SingleCellExperiment
      data$sce <- sce_data

      # Stocker les colonnes des métadonnées
      metadata_columns(colnames(data$metadata))
      showNotification("Dataset loaded successfully!", type = "message")

      # Étape 2: Charger le réseau ligand-récepteur et la matrice ligand-cible
      data$lr_network <- readRDS(input$lr_network_file$datapath)
      data$ligand_target_matrix <- readRDS(input$ligand_target_matrix_file$datapath)

      # Nettoyage des noms de colonnes et de lignes pour ligand_target_matrix
      colnames(data$ligand_target_matrix) <- make.names(colnames(data$ligand_target_matrix))
      rownames(data$ligand_target_matrix) <- make.names(rownames(data$ligand_target_matrix))

      showNotification("Ligand-receptor network and ligand-target matrix loaded successfully!", type = "message")


    }, error = function(e) {
      # En cas d'erreur, afficher une notification et imprimer l'erreur complète
      showNotification(paste("Error during data loading:", e$message), type = "error")
      print(paste("Full error:", capture.output(str(e))))  # Afficher l'erreur complète pour le débogage

    }, finally = {
      # Cacher la boîte modale après le chargement
      hide_loading_modal()
    })
  })

  # Fonction pour associer les symboles MGI aux identifiants Ensembl dans l'objet SCE
  update_gene_symbols <- function(sce, conversion_table) {
    # Vérifier que le fichier de conversion contient bien les colonnes requises
    if (!all(c("Gene.stable.ID", "Marker.Symbol") %in% colnames(conversion_table))) {
      stop("Les colonnes 'Gene.stable.ID' ou 'Marker.Symbol' sont manquantes dans la table de conversion.")
    }

    # Associer les nouveaux symboles de gènes (Marker.Symbol) aux identifiants Ensembl dans l'objet sce
    gene_symbols <- conversion_table$Marker.Symbol[match(rownames(sce), conversion_table$Gene.stable.ID)]

    # Remplacer les noms de gènes dans l'objet SCE
    rownames(sce) <- ifelse(!is.na(gene_symbols) & gene_symbols != "", gene_symbols, rownames(sce))

    # Retourner l'objet SCE mis à jour
    return(sce)
  }


  ############################### GROUPES: SÉLECTION ET FORMULE DE CONTRASTE ##################################

  # Mettre à jour les colonnes de métadonnées après le chargement
  observe({
    columns <- metadata_columns()
    if (!is.null(columns)) {
      updateSelectInput(session, "sample_id", choices = columns)
      updateSelectInput(session, "group_id", choices = columns)
      updateSelectInput(session, "celltype_id", choices = columns)
    }
  })

  # Observer pour sauvegarder les paramètres de métadonnées après que l'utilisateur clique sur le bouton "Set Metadata Parameters"
  observeEvent(input$define_params, {
    req(input$sample_id, input$group_id, input$celltype_id)  # Vérification des paramètres définis

    # Sauvegarder dans des variables réactives
    data$sample_id <- input$sample_id
    data$group_id <- input$group_id
    data$celltype_id <- input$celltype_id

    # Récupérer les données
    sce <- data$sce
    col_data_names <- colnames(SummarizedExperiment::colData(sce))

    # Vérification si les colonnes existent
    if (!(input$group_id %in% col_data_names)) {
      showNotification("Error: Selected group_id not found in colData.", type = "error")
      return(NULL)
    }
    if (!(input$celltype_id %in% col_data_names)) {
      showNotification("Error: Selected celltype_id not found in colData.", type = "error")
      return(NULL)
    }
    if (!(input$sample_id %in% col_data_names)) {
      showNotification("Error: Selected sample_id not found in colData.", type = "error")
      return(NULL)
    }

    # Si tout est bon, afficher une notification
    showNotification("Metadata parameters set and valid.", type = "message")

    # Assurez-vous que les colonnes sont bien des facteurs
    SummarizedExperiment::colData(sce)[, input$group_id] <- factor(SummarizedExperiment::colData(sce)[, input$group_id])
    SummarizedExperiment::colData(sce)[, input$celltype_id] <- factor(SummarizedExperiment::colData(sce)[, input$celltype_id])
    SummarizedExperiment::colData(sce)[, input$sample_id] <- factor(SummarizedExperiment::colData(sce)[, input$sample_id])

    # Assurez-vous que les informations sont sauvegardées
    data$sce <- sce
  })



  # Mettre à jour les niveaux de groupe lorsqu'une colonne group_id est sélectionnée
  observeEvent(input$group_id, {
    req(data$sce)  # Vérifier que les données sont chargées
    sce_data <- data$sce
    selected_group_column <- input$group_id

    tryCatch({
      # Vérifier que la colonne group_id est bien dans colData
      if (!selected_group_column %in% colnames(colData(sce_data))) {
        stop(paste("The column", selected_group_column, "does not exist in colData."))
      }

      # Extraire les niveaux du groupe et afficher les premières lignes de colData pour vérifier son contenu
      print(head(colData(sce_data)))  # Debugging: afficher les premières lignes de colData

      groups <- unique(colData(sce_data)[[selected_group_column]])  # Extraire les niveaux du groupe

      # Mettre à jour les niveaux de groupe
      group_levels(groups)

    }, error = function(e) {
      showNotification(paste("Error extracting group levels:", e$message), type = "error")
      print(paste("Full error:", capture.output(str(e))))  # Debugging: afficher l'erreur complète
    })
  })

  # UI dynamique pour sélectionner les groupes (2 ou 3 groupes à comparer)
  output$group_select_ui <- renderUI({
    groups <- group_levels()

    if (!is.null(groups) && length(groups) > 0) {
      selectInput("group_select", "Select Groups for Contrast", choices = groups, multiple = TRUE)  # Ici, multiple = TRUE est correct
    } else {
      p("No groups available for selection. Please select a valid group ID column.")
    }
  })


  ############################### STEP 3: DEFINE CONTRASTS ##################################
  generate_contrast_formula <- function(selected_groups, sce, group_id) {
    # Get exact group names from data
    actual_groups <- levels(factor(SummarizedExperiment::colData(sce)[[group_id]]))

    print("Available groups in data (from levels):")
    print(actual_groups)
    print("Selected groups:")
    print(selected_groups)

    result <- list()

    # For 2 groups
    if (length(selected_groups) == 2) {
      g1 <- selected_groups[1]
      g2 <- selected_groups[2]

      result$contrasts_oi <- sprintf("'%s-%s','%s-%s'", g1, g2, g2, g1)
      result$contrast_tbl <- tibble(
        contrast = c(paste0(g1, "-", g2), paste0(g2, "-", g1)),
        group = c(g1, g2)
      )
    }
    # For 3 groups
    else if (length(selected_groups) == 3) {
      g1 <- selected_groups[1]
      g2 <- selected_groups[2]
      g3 <- selected_groups[3]

      # Format for 3 groups with division by 2
      contrast1 <- sprintf("%s-(%s+%s)/2", g1, g2, g3)
      contrast2 <- sprintf("%s-(%s+%s)/2", g2, g1, g3)
      contrast3 <- sprintf("%s-(%s+%s)/2", g3, g1, g2)

      result$contrasts_oi <- sprintf("'%s','%s','%s'", contrast1, contrast2, contrast3)
      result$contrast_tbl <- tibble(
        contrast = c(contrast1, contrast2, contrast3),
        group = c(g1, g2, g3)
      )
    } else {
      showNotification("Please select 2 or 3 groups.", type = "error")
      return(NULL)
    }

    print("Final contrast table:")
    print(result$contrast_tbl)
    print("Final contrasts_oi:")
    print(result$contrasts_oi)

    return(result)
  }

  # Save contrast observer
  observeEvent(input$save_contrast, {
    req(data$sce, input$group_id)

    # Get exact group names from data
    all_groups <- levels(factor(SummarizedExperiment::colData(data$sce)[[input$group_id]]))
    selected_groups <- input$group_select

    if (is.null(selected_groups) || length(selected_groups) < 2) {
      showNotification("Please select 2 or 3 groups.", type = "error")
      return(NULL)
    }

    # Match selected groups to actual levels
    selected_groups <- all_groups[match(selected_groups, all_groups)]

    tryCatch({
      contrast_info <- generate_contrast_formula(selected_groups, data$sce, input$group_id)
      if (is.null(contrast_info)) return(NULL)

      # Store results
      data$contrasts_oi <- contrast_info$contrasts_oi
      data$contrast_tbl <- contrast_info$contrast_tbl

      # Display formula
      output$generated_contrast_formula <- renderText({
        contrast_info$contrasts_oi
      })

      # Set other parameters
      data$senders_oi <- unique(SummarizedExperiment::colData(data$sce)[[input$celltype_id]])
      data$receivers_oi <- unique(SummarizedExperiment::colData(data$sce)[[input$celltype_id]])

      showNotification("Contrast and parameters saved successfully!", type = "message")

    }, error = function(e) {
      showNotification(paste("Error saving contrast:", e$message), type = "error")
    })
  })
####################### ANALYSIS #####################################

  ######## Cell-type filtering: determine which cell types are sufficiently present #############

  observeEvent(input$calculate_abundance_info, {
    req(data$sce)

    tryCatch({
      print("Calculating Abundance Info...")

      # Call get_abundance_info function
      abundance_info <- get_abundance_info(
        sce = data$sce,
        sample_id = input$sample_id,
        group_id = input$group_id,
        celltype_id = input$celltype_id,
        min_cells = input$min_cells,
        senders_oi = data$senders_oi,
        receivers_oi = data$receivers_oi,
        batches = NA
      )

      # Store abundance info in reactive values
      data$abundance_info <- abundance_info

      print("Abundance Info calculated.")
      showNotification("Abundance Info calculated and saved successfully!", type = "message")

    }, error = function(e) {
      showNotification(paste("Error calculating Abundance Info:", e$message), type = "error")
      print(paste("Full error:", capture.output(str(e))))
    })
  })


######## Gene filtering: determine which genes are sufficiently expressed in each present cell type #############


   observeEvent(input$calculate_frq_list, {
    req(data$sce)  # Vérifier que les données sont chargées

    tryCatch({
      print("Calculating Fraction Expression List...")

      # Appel à la fonction get_frac_exprs pour calculer frq_list
      frq_list <- get_frac_exprs(
        sce = data$sce,
        sample_id = input$sample_id,
        celltype_id = input$celltype_id,
        group_id = input$group_id,
        batches = NA,
        min_cells = input$min_cells,
        fraction_cutoff = input$fraction_cutoff,
        min_sample_prop = input$min_sample_prop
      )

      # Stockage de frq_list dans une variable réactive
      data$frq_list <- frq_list

      print("Fraction Expression List calculated.")
      print(frq_list)

      # Extraction des gènes d'intérêt en utilisant dplyr::filter
      genes_oi <- frq_list$expressed_df %>%
        dplyr::filter(expressed == TRUE) %>%
        dplyr::pull(gene) %>%
        unique()

      # Mettre à jour l'objet sce avec les gènes d'intérêt
      sce <- data$sce[genes_oi, ]
      data$genes_oi <- genes_oi
      data$sce <- sce

      showNotification("Fraction Expression List calculated and saved successfully!", type = "message")

    }, error = function(e) {
      showNotification(paste("Error calculating Fraction Expression List:", e$message), type = "error")
      print(paste("Full error:", capture.output(str(e))))  # Log complet
    })
  })





  # Observer for abundance and expression information analysis
  observeEvent(input$run_abundance_expression_info, {
    req(data$sce, data$frq_list, data$abundance_info, data$lr_network)

    tryCatch({
      print("Starting abundance and expression info analysis...")

      # Process abundance and expression info with suppressMessages
      abundance_expression_info <- suppressMessages(
        process_abundance_expression_info(
          sce = data$sce,
          sample_id = input$sample_id,
          group_id = input$group_id,
          celltype_id = input$celltype_id,
          min_cells = input$min_cells,
          senders_oi = data$senders_oi,
          receivers_oi = data$receivers_oi,
          lr_network = data$lr_network,
          frq_list = data$frq_list,
          abundance_info = data$abundance_info,
          batches = NA
        )
      )

      # Store the results
      data$abundance_expression_info <- abundance_expression_info

      print("Abundance and Expression Info analysis completed.")
      print("Results structure:")
      print(str(abundance_expression_info))

      # Verify that data was stored
      if (!is.null(data$abundance_expression_info)) {
        showNotification("Abundance expression info stored successfully!", type = "message")
      } else {
        stop("Failed to store abundance expression info")
      }

    }, error = function(e) {
      showNotification(paste("Error during analysis:", e$message), type = "error")
      print("Full error trace:")
      print(capture.output(str(e)))
    })
  })

  # Add validation check
  observe({
    if (!is.null(data$abundance_expression_info)) {
      print("Validating abundance_expression_info:")
      print("Components present:")
      print(names(data$abundance_expression_info))
    }
  })

  # Observer for Differential Expression analysis and preparing data for ligand analysis
  observeEvent(input$run_DE_analysis, {
    req(data$sce, data$lr_network)

    tryCatch({
      print("Running DE analysis...")

      # Run DE analysis
      DE_info <- get_DE_info(
        sce = data$sce,
        sample_id = input$sample_id,
        group_id = input$group_id,
        celltype_id = input$celltype_id,
        contrasts_oi = data$contrasts_oi,
        min_cells = 1,
        expressed_df = data$frq_list$expressed_df,
        batches = NA,
        covariates = NA
      )

      # Store complete DE results
      data$de_output <- DE_info

      print("Processing DE results...")

      # Get celltype_de based on empirical p-values setting
      if(input$empirical_pval_multinichenet) {
        print("Calculating empirical p-values...")
        DE_info_emp <- get_empirical_pvals(DE_info$celltype_de$de_output_tidy)
        celltype_de <- DE_info_emp$de_output_tidy_emp %>%
          dplyr::select(-p_val, -p_adj) %>%
          dplyr::rename(
            p_val = p_emp,
            p_adj = p_adj_emp
          )
      } else {
        print("Using standard p-values...")
        celltype_de <- DE_info$celltype_de$de_output_tidy
      }

      # Validate celltype_de
      if(is.null(celltype_de)) {
        stop("celltype_de is NULL after processing")
      }

      print("Dimensions of celltype_de:")
      print(dim(celltype_de))
      print("Columns in celltype_de:")
      print(colnames(celltype_de))

      # Store celltype_de
      data$celltype_de <- celltype_de

      # Combine sender and receiver DE information
      print("Combining sender and receiver DE information...")
      sender_receiver_de <- combine_sender_receiver_de(
        sender_de = celltype_de,
        receiver_de = celltype_de,
        senders_oi = data$senders_oi,
        receivers_oi = data$receivers_oi,
        lr_network = data$lr_network
      )

      # Store sender_receiver_de
      data$sender_receiver_de <- sender_receiver_de

      print("DE analysis complete. Summary:")
      print(paste("Number of celltype DE results:", nrow(celltype_de)))
      print(paste("Number of sender-receiver pairs:", nrow(sender_receiver_de)))

      showNotification("DE analysis and sender-receiver pairing completed successfully!", type = "message")

    }, error = function(e) {
      showNotification(paste("Error during DE analysis:", e$message), type = "error")
      print("Full error trace:")
      print(capture.output(str(e)))
    })
  })

  # Validation observer for DE results
  observe({
    if(!is.null(data$celltype_de)) {
      print("Validating DE results:")
      print("celltype_de structure:")
      str(data$celltype_de)
      print("sender_receiver_de structure:")
      str(data$sender_receiver_de)
    }
  })




  observeEvent(input$calculate_abundance_and_frq, {
    req(data$sce)

    showModal(modalDialog(
      title = "Calculating Analysis Parameters",
      "Please wait while calculations are running...",
      footer = NULL,
      easyClose = FALSE
    ))

    tryCatch({
      # Step 1: Calculate Abundance Info
      print("Calculating Abundance Info...")
      abundance_info <- get_abundance_info(
        sce = data$sce,
        sample_id = input$sample_id,
        group_id = input$group_id,
        celltype_id = input$celltype_id,
        min_cells = input$min_cells,
        senders_oi = data$senders_oi,
        receivers_oi = data$receivers_oi,
        batches = NA
      )

      # Store abundance info
      data$abundance_info <- abundance_info
      print("Abundance Info calculated")

      # Step 2: Calculate Fraction Expression List
      print("Calculating Fraction Expression List...")
      frq_list <- get_frac_exprs(
        sce = data$sce,
        sample_id = input$sample_id,
        celltype_id = input$celltype_id,
        group_id = input$group_id,
        batches = NA,
        min_cells = input$min_cells,
        fraction_cutoff = input$fraction_cutoff,
        min_sample_prop = input$min_sample_prop
      )

      # Store frq_list
      data$frq_list <- frq_list

      # Extract genes of interest
      genes_oi <- frq_list$expressed_df %>%
        dplyr::filter(expressed == TRUE) %>%
        dplyr::pull(gene) %>%
        unique()

      # Update sce with genes of interest
      sce <- data$sce[genes_oi, ]
      data$genes_oi <- genes_oi
      data$sce <- sce

      print("Fraction Expression List calculated")
      print(paste("Number of genes of interest:", length(genes_oi)))

      showNotification("Abundance and Expression calculations completed successfully!", type = "message")

    }, error = function(e) {
      showNotification(paste("Error during calculations:", e$message), type = "error")
      print("Full error trace:")
      print(capture.output(str(e)))
    }, finally = {
      removeModal()
    })
  })

  # Dans le server, créer un nouvel observer qui combine les deux analyses
  observeEvent(input$run_expression_and_DE, {
    req(data$sce, data$frq_list, data$abundance_info, data$lr_network)

    showModal(modalDialog(
      title = "Running Analysis",
      "Please wait while the analyses are running...",
      footer = NULL,
      easyClose = FALSE
    ))

    print("Starting combined analysis...")

    tryCatch({
      # Abundance & Expression Info analysis
      print("Starting Abundance and Expression Info analysis...")
      abundance_expression_info <- suppressMessages(
        process_abundance_expression_info(
          sce = data$sce,
          sample_id = input$sample_id,
          group_id = input$group_id,
          celltype_id = input$celltype_id,
          min_cells = input$min_cells,
          senders_oi = data$senders_oi,
          receivers_oi = data$receivers_oi,
          lr_network = data$lr_network,
          frq_list = data$frq_list,
          abundance_info = data$abundance_info,
          batches = NA
        )
      )

      data$abundance_expression_info <- abundance_expression_info
      print("Abundance and Expression Info analysis completed")

      # DE Analysis
      print("Starting DE analysis...")
      de_output <- get_DE_info(
        sce = data$sce,
        sample_id = input$sample_id,
        group_id = input$group_id,
        celltype_id = input$celltype_id,
        contrasts_oi = data$contrasts_oi,
        min_cells = 1,
        expressed_df = data$frq_list$expressed_df,
        batches = NA,
        covariates = NA
      )

      data$de_output <- de_output
      print("DE output structure:")
      print(str(de_output))

      # Process DE results
      print("Processing DE results...")
      # On vérifie si celltype_de est une liste
      if (!is.null(de_output$celltype_de) && is.list(de_output$celltype_de) &&
          !is.null(de_output$celltype_de$de_output_tidy)) {

        celltype_de <- de_output$celltype_de$de_output_tidy
        # Convertir en data frame si ce n'est pas déjà le cas
        if (!is.data.frame(celltype_de)) {
          celltype_de <- as.data.frame(celltype_de)
        }

        # Ajouter direction_regulation
        celltype_de <- celltype_de %>%
          dplyr::mutate(
            logFC = as.numeric(logFC),
            direction_regulation = ifelse(logFC > 0, "up", "down")
          )

        # Assurer que direction_regulation est un facteur avec les bons niveaux
        celltype_de$direction_regulation <- factor(
          celltype_de$direction_regulation,
          levels = c("up", "down")
        )

        data$celltype_de <- celltype_de
        print("DE results processed successfully")
        print(paste("Number of DE results:", nrow(celltype_de)))

      } else {
        stop("Invalid structure in DE output")
      }

      # Update plots
      if (!is.null(de_output$hist_pvals)) {
        output$de_pvals_hist <- renderPlot({
          de_output$hist_pvals
        }, res = 96)

        output$download_de_pvals_hist <- downloadHandler(
          filename = function() {
            paste("de_pvals_histogram_", format(Sys.time(), "%Y%m%d"), ".tiff", sep = "")
          },
          content = function(file) {
            tryCatch({
              tiff(file, width = 10, height = 8, units = 'in', res = 300, compression = "lzw")
              print(de_output$hist_pvals)
              dev.off()
            }, error = function(e) {
              if (dev.cur() > 1) dev.off()
              showNotification(paste("Error saving plot:", e$message), type = "error")
            })
          }
        )
      }

      showNotification("Both analyses completed successfully!", type = "message")

    }, error = function(e) {
      showNotification(paste("Error during analysis:", e$message), type = "error")
      print("Full error trace:")
      print(capture.output(str(e)))
    }, finally = {
      removeModal()
    })
  })


  # Render DE p-values histogram
  output$de_pvals_hist <- renderPlot({
    req(data$de_output)
    data$de_output$hist_pvals
  }, res = 96)

  # Download handler for DE p-values histogram
  output$download_de_pvals_hist <- downloadHandler(
    filename = function() {
      paste("de_pvals_histogram_", format(Sys.time(), "%Y%m%d"), ".tiff", sep = "")
    },
    content = function(file) {
      tryCatch({
        # Save the plot as TIFF
        tiff(file,
             width = 10,
             height = 8,
             units = 'in',
             res = 300,
             compression = "lzw")

        # Print the plot
        print(data$de_output$hist_pvals)

        # Close the device
        dev.off()

      }, error = function(e) {
        if (dev.cur() > 1) dev.off()
        showNotification(paste("Error saving plot:", e$message), type = "error")
      })
    }
  )




  ################# Ligand receptor analysis ################


  # Function to process geneset data and calculate the ratio of geneset vs background genes
  process_geneset_data <- function(contrast, celltype_de, logFC_threshold, p_val_adj, p_val_threshold) {
    tryCatch({
      # Process each cell type to get background and geneset information
      celltype_de$cluster_id %>% unique() %>%
        lapply(function(cluster) {
          # Get background genes for current cell type
          background_genes <- celltype_de %>%
            dplyr::filter(cluster_id == cluster) %>%
            pull(gene) %>% unique()
          n_background <- length(background_genes)

          # Get differentially expressed genes based on thresholds
          de_data <- celltype_de %>%
            dplyr::filter(
              cluster_id == cluster,
              contrast == !!contrast,
              if(p_val_adj) p_adj <= p_val_threshold else p_val <= p_val_threshold
            )

          # Split into up and down regulated genes
          genes_up <- de_data %>%
            dplyr::filter(logFC >= logFC_threshold) %>%
            pull(gene) %>% unique()

          genes_down <- de_data %>%
            dplyr::filter(logFC <= -logFC_threshold) %>%
            pull(gene) %>% unique()

          # Calculate proportions and check if within recommended range (1/200 to 1/10)
          prop_up <- length(genes_up) / n_background
          prop_down <- length(genes_down) / n_background

          in_range_up <- prop_up >= 1/200 && prop_up <= 1/10
          in_range_down <- prop_down >= 1/200 && prop_down <= 1/10

          # Return results as tibble
          tibble(
            cluster_id = cluster,
            n_background = n_background,
            n_geneset_up = length(genes_up),
            n_geneset_down = length(genes_down),
            prop_geneset_up = prop_up,
            prop_geneset_down = prop_down,
            in_range_up = in_range_up,
            in_range_down = in_range_down,
            contrast = contrast,
            logFC_threshold = logFC_threshold,
            p_val_threshold = p_val_threshold,
            adjusted = p_val_adj
          )
        }) %>% bind_rows()
    }, error = function(e) {
      showNotification(paste("Error in geneset processing:", e$message), type = "error")
      NULL
    })
  }

  # Observer for combined geneset and ligand activity analysis
  observeEvent(input$run_combined_analysis, {
    req(data$sce, data$lr_network, data$ligand_target_matrix, data$abundance_expression_info)

    showModal(modalDialog(
      title = "Analysis in Progress",
      "Please wait while the analyses are running...",
      footer = NULL,
      easyClose = FALSE
    ))

    tryCatch({
      # Step 1: Geneset Assessment
      print("Starting geneset assessment...")
      geneset_assessment <- data$contrast_tbl$contrast %>%
        lapply(
          process_geneset_data,
          celltype_de = data$celltype_de,
          logFC_threshold = input$logFC_threshold_geneset_multinichenet,
          p_val_adj = FALSE,
          p_val_threshold = input$p_val_threshold_geneset_multinichenet
        ) %>% bind_rows()

      data$geneset_assessment <- geneset_assessment

      # Step 2: Create sender_receiver_de
      print("Creating sender-receiver DE information...")
      sender_receiver_de <- combine_sender_receiver_de(
        sender_de = data$celltype_de,
        receiver_de = data$celltype_de,
        senders_oi = data$senders_oi,
        receivers_oi = data$receivers_oi,
        lr_network = data$lr_network
      )

      # Store sender_receiver_de
      data$sender_receiver_de <- sender_receiver_de

      # Step 3: Ligand Activity Analysis
      print("Starting ligand activity analysis...")
      receivers_oi <- intersect(
        data$receivers_oi,
        data$celltype_de$cluster_id %>% unique()
      )

      ligand_activities_targets_DEgenes <- suppressMessages(suppressWarnings(
        get_ligand_activities_targets_DEgenes(
          receiver_de = data$celltype_de,
          receivers_oi = receivers_oi,
          ligand_target_matrix = data$ligand_target_matrix,
          logFC_threshold = input$logFC_threshold_geneset_multinichenet,
          p_val_threshold = input$p_val_threshold_geneset_multinichenet,
          p_val_adj = input$p_val_adj_geneset_multinichenet,
          top_n_target = as.numeric(input$top_n_target_multinichenet),
          verbose = TRUE,
          n.cores = as.numeric(input$n_cores_multinichenet)
        )
      ))

      data$ligand_activities_targets_DEgenes <- ligand_activities_targets_DEgenes
      print("Ligand activity analysis completed")

      # Step 4: Prioritization Analysis
      print("Starting prioritization analysis...")

      # Prepare metadata for prioritization
      metadata_combined <- SummarizedExperiment::colData(data$sce) %>%
        tibble::as_tibble()

      grouping_tbl <- metadata_combined[,c(input$sample_id, input$group_id)] %>%
        tibble::as_tibble() %>%
        distinct()
      colnames(grouping_tbl) <- c("sample", "group")
      data$grouping_tbl <- grouping_tbl
      # Create sender-receiver table
      sender_receiver_tbl <- data$sender_receiver_de %>%
        distinct(sender, receiver)

      # Generate prioritization tables
      prioritization_tables <- suppressMessages(
        generate_prioritization_tables(
          sender_receiver_info = data$abundance_expression_info$sender_receiver_info,
          sender_receiver_de = data$sender_receiver_de,
          ligand_activities_targets_DEgenes = ligand_activities_targets_DEgenes,
          contrast_tbl = data$contrast_tbl,
          sender_receiver_tbl = sender_receiver_tbl,
          grouping_tbl = grouping_tbl,
          scenario = "regular",
          fraction_cutoff = input$fraction_cutoff_multinichenet,
          abundance_data_receiver = data$abundance_expression_info$abundance_data_receiver,
          abundance_data_sender = data$abundance_expression_info$abundance_data_sender,
          ligand_activity_down = input$ligand_activity_down_multinichenet
        )
      )
      data$prioritization_tables <- prioritization_tables
      print("Prioritization analysis completed")

      # Step 5: Calculate LR-Target Correlations
      print("Calculating ligand-receptor-target correlations...")
      lr_target_prior_cor <- lr_target_prior_cor_inference(
        receivers_oi = prioritization_tables$group_prioritization_tbl$receiver %>% unique(),
        abundance_expression_info = data$abundance_expression_info,
        celltype_de = data$celltype_de,
        grouping_tbl = grouping_tbl,
        prioritization_tables = prioritization_tables,
        ligand_target_matrix = data$ligand_target_matrix,
        logFC_threshold = input$logFC_threshold_geneset_multinichenet,
        p_val_threshold = input$p_val_threshold_geneset_multinichenet,
        p_val_adj = input$p_val_adj_geneset_multinichenet
      )

      # Update displays
      output$ligand_activities_table_multinichenet <- renderTable({
        ligand_activities_targets_DEgenes$ligand_activities %>% head(20)
      })

      output$prioritization_table_multinichenet <- renderTable({
        prioritization_tables$group_prioritization_tbl %>% head(20)
      })

      data$multinichenet_output <- list(
        celltype_info = data$abundance_expression_info$celltype_info,
        celltype_de = data$celltype_de,
        sender_receiver_info = data$abundance_expression_info$sender_receiver_info,
        sender_receiver_de = data$sender_receiver_de,
        ligand_activities_targets_DEgenes = data$ligand_activities_targets_DEgenes,
        prioritization_tables = data$prioritization_tables,
        grouping_tbl = grouping_tbl,
        lr_target_prior_cor = lr_target_prior_cor
      )

      print("Stored results in multinichenet_output:")
      print(names(data$multinichenet_output))

      showNotification("Complete analysis pipeline finished successfully!", type = "message")

    }, error = function(e) {
      showNotification(paste("Error during analysis:", e$message), type = "error")
      print("Full error trace:")
      print(capture.output(str(e)))
    }, finally = {
      removeModal()
    })
  })

  # Handler pour les ligand activities
  output$download_ligand_activities <- downloadHandler(
    filename = function() {
      paste("ligand_activities_", format(Sys.time(), "%Y%m%d"), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(data$ligand_activities_targets_DEgenes$ligand_activities, file, row.names = FALSE)
    }
  )

  # Handler pour la table de prioritization
  output$download_prioritization <- downloadHandler(
    filename = function() {
      paste("prioritization_table_", format(Sys.time(), "%Y%m%d"), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(data$prioritization_tables$group_prioritization_tbl, file, row.names = FALSE)
    }
  )

  # Handler pour télécharger l'objet RDS complet
  output$download_all_results <- downloadHandler(
    filename = function() {
      paste("multinichenet_output_", format(Sys.time(), "%Y%m%d"), ".rds", sep = "")
    },
    content = function(file) {
      req(data$multinichenet_output)

      showModal(modalDialog(
        title = "Saving Analysis Results",
        "Please wait while saving the results...",
        footer = NULL,
        easyClose = FALSE
      ))

      tryCatch({
        # D'abord créer une copie de l'objet
        multinichenet_output <- data$multinichenet_output

        # Utiliser make_lite_output pour compresser
        print("Compressing output...")
        multinichenet_output <- make_lite_output(multinichenet_output)

        # Sauvegarder l'objet compressé
        print("Saving to file...")
        saveRDS(multinichenet_output, file = file)

        showNotification("Results saved successfully!", type = "message")

      }, error = function(e) {
        showNotification(paste("Error saving results:", e$message), type = "error")
        print("Error details:")
        print(str(e))
      }, finally = {
        removeModal()
      })
    }
  )

  ################### Circos Plot Functions ###################



  create_circos_plots <- function(prioritization_tables,
                                  top_n = 50,
                                  rank_per_group = FALSE,
                                  color_palette = "Spectral",
                                  track_height = 0.05,  # Reduced default height
                                  track_margin = c(0.01, 0.01),
                                  gap_degree = 2) {
    # Get top interactions
    prioritized_tbl_oi_all <- get_top_n_lr_pairs(
      prioritization_tables,
      top_n = top_n,
      rank_per_group = rank_per_group
    )

    # Create prioritized table
    prioritized_tbl_oi <- prioritization_tables$group_prioritization_tbl %>%
      dplyr::filter(id %in% prioritized_tbl_oi_all$id) %>%
      distinct(id, sender, receiver, ligand, receptor, group) %>%
      left_join(prioritized_tbl_oi_all)

    # Handle NA values
    prioritized_tbl_oi$prioritization_score[is.na(prioritized_tbl_oi$prioritization_score)] <- 0

    # Get unique senders and receivers
    senders_receivers <- union(
      prioritized_tbl_oi$sender %>% unique(),
      prioritized_tbl_oi$receiver %>% unique()
    ) %>% sort()

    # Generate colors, handling the case where number of colors needed exceeds palette size
    n_colors_needed <- length(senders_receivers)
    n_colors_available <- min(n_colors_needed, 11)  # RColorBrewer palettes have max 11 colors

    colors_sender <- RColorBrewer::brewer.pal(
      n = n_colors_available,
      name = color_palette
    )

    # If we need more colors, interpolate
    if(n_colors_needed > n_colors_available) {
      colors_sender <- colorRampPalette(colors_sender)(n_colors_needed)
    }

    colors_sender <- setNames(colors_sender, senders_receivers)
    colors_receiver <- colors_sender

    # Set circos parameters
    circos.clear()
    circos.par(
      track.height = track_height,
      track.margin = track_margin,
      gap.degree = gap_degree,
      cell.padding = c(0, 0, 0, 0)  # Minimize cell padding
    )

    # Create circos plots
    circos_list <- make_circos_group_comparison(
      prioritized_tbl_oi,
      colors_sender,
      colors_receiver
    )

    return(circos_list)
  }
  # Dans l'observer
  observeEvent(input$generate_circos_multinichenet, {
    req(data$multinichenet_output$prioritization_tables)

    showModal(modalDialog(
      title = "Generating Circos Plots",
      "Please wait while the plots are being generated...",
      footer = NULL,
      easyClose = FALSE
    ))

    tryCatch({
      # Create circos plots with user parameters
      data$circos_plots <- create_circos_plots(
        prioritization_tables = data$multinichenet_output$prioritization_tables,
        top_n = input$top_n_interactions_multinichenet,
        color_palette = input$color_palette_multinichenet
      )

      # Update group choices for download
      updateSelectInput(session, "selected_group",
                        choices = names(data$circos_plots))

      # Display plots
      output$circos_plots_ui_multinichenet <- renderUI({
        box(
          title = "Circos Plots",
          width = 12,
          lapply(names(data$circos_plots), function(group) {
            plotname <- paste0("circos_plot_", make.names(group))
            div(
              h4(paste("Group:", group)),
              plotOutput(plotname),
              hr()
            )
          })
        )
      })

      # Create plot outputs
      for(group in names(data$circos_plots)) {
        local({
          group_local <- group
          output[[paste0("circos_plot_", make.names(group_local))]] <- renderPlot({
            par(mar = c(1, 1, 1, 1))
            replayPlot(data$circos_plots[[group_local]])
          })
        })
      }

      showNotification("Circos plots generated successfully!", type = "message")

    }, error = function(e) {
      showNotification(paste("Error generating Circos plots:", e$message), type = "error")
    }, finally = {
      removeModal()
    })
  })


  # Handler pour télécharger le circos plot
  output$download_circos_plot <- downloadHandler(
    filename = function() {
      paste0("circos_plot_", input$selected_group, "_", format(Sys.time(), "%Y%m%d"), ".tiff")
    },
    content = function(file) {
      req(data$circos_plots, input$selected_group)

      tryCatch({
        # Créer le fichier TIFF
        tiff(file,
             width = 10,
             height = 10,
             units = "in",
             res = input$circos_plot_dpi,
             compression = "lzw")

        # Recréer le plot
        replayPlot(data$circos_plots[[input$selected_group]])

        # Fermer le device
        dev.off()

        showNotification("Plot downloaded successfully!", type = "message")

      }, error = function(e) {
        showNotification(paste("Error saving plot:", e$message), type = "error")
        if (dev.cur() > 1) dev.off()  # S'assurer que le device est fermé en cas d'erreur
      })
    }
  )

  # Update observer for select inputs
  observe({
    req(data$multinichenet_output)

    # Get unique values
    groups <- unique(data$multinichenet_output$prioritization_tables$group_prioritization_tbl$group)
    senders <- unique(data$multinichenet_output$prioritization_tables$group_prioritization_tbl$sender)
    receivers <- unique(data$multinichenet_output$prioritization_tables$group_prioritization_tbl$receiver)

    # Add "All" option at the beginning of the lists
    senders_choices <- c("All Senders" = "all", setNames(senders, senders))
    receivers_choices <- c("All Receivers" = "all", setNames(receivers, receivers))

    # Update select inputs
    updateSelectInput(session, "selected_group_viz", choices = groups)
    updateSelectInput(session, "specific_sender",
                      choices = senders_choices,
                      selected = "all")
    updateSelectInput(session, "specific_receiver",
                      choices = receivers_choices,
                      selected = "all")
  })

  # Update plot generation event
  observeEvent(input$generate_lr_plot, {
    req(data$multinichenet_output, input$selected_group_viz)

    showModal(modalDialog(
      title = "Generating Plot",
      "Please wait...",
      footer = NULL,
      easyClose = FALSE
    ))

    tryCatch({
      # Initialize parameters for get_top_n_lr_pairs
      params <- list(
        prioritization_tables = data$multinichenet_output$prioritization_tables,
        top_n = input$top_n_interactions,
        groups_oi = input$selected_group_viz
      )

      # Add sender/receiver parameters only if specific ones are selected
      if (!is.null(input$specific_sender) && input$specific_sender != "all") {
        params$senders_oi <- input$specific_sender
      }
      if (!is.null(input$specific_receiver) && input$specific_receiver != "all") {
        params$receivers_oi <- input$specific_receiver
      }

      # Get prioritized interactions
      prioritized_interactions <- do.call(get_top_n_lr_pairs, params)

      # Generate plot
      data$current_lr_plot <- make_sample_lr_prod_activity_plots_Omnipath(
        data$multinichenet_output$prioritization_tables,
        prioritized_interactions %>% dplyr::inner_join(data$lr_network)
      )

      # Display plot
      output$lr_activity_plot <- renderPlot({
        data$current_lr_plot
      })

      showNotification("Plot generated successfully!", type = "message")

    }, error = function(e) {
      showNotification(paste("Error generating plot:", e$message), type = "error")
      print(paste("Full error:", capture.output(str(e))))
    }, finally = {
      removeModal()
    })
  })

  # Download handlers
  output$download_lr_plot <- downloadHandler(
    filename = function() {
      paste0("lr_activity_plot_", format(Sys.time(), "%Y%m%d"), ".tiff")
    },
    content = function(file) {
      tryCatch({
        req(data$current_lr_plot)  # Ensure the plot exists
        ggsave(file, plot = data$current_lr_plot, device = "tiff", width = 10, height = 8, dpi = input$lr_plot_dpi)
      }, error = function(e) {
        showNotification(paste("Download error:", e$message), type = "error")
      })
    }
  )

  # Observer pour mettre à jour les choix de sender
  observe({
    req(data$multinichenet_output$prioritization_tables)

    # Get unique senders
    senders <- unique(data$multinichenet_output$prioritization_tables$group_prioritization_tbl$sender)
    updateSelectInput(session, "sender_to_color", choices = senders)
  })

  # Create dynamic color selectors
  output$color_selectors <- renderUI({
    req(data$multinichenet_output)

    senders <- unique(data$multinichenet_output$prioritization_tables$group_prioritization_tbl$sender)

    lapply(senders, function(sender) {
      colourInput(
        inputId = paste0("color_", make.names(sender)),
        label = paste("Color for", sender),
        value = "pink"
      )
    })
  })

  # Network Generation
  observeEvent(input$generate_network, {
    req(data$multinichenet_output, input$sender_to_color)

    showModal(modalDialog(
      title = "Generating Network",
      "Please wait...",
      footer = NULL,
      easyClose = FALSE
    ))

    tryCatch({
      prioritized_tbl_oi_all <- get_top_n_lr_pairs(
        data$multinichenet_output$prioritization_tables,
        input$top_n_network,
        rank_per_group = FALSE
      )

      # Différencier le traitement selon l'utilisation des corrélations
      if(input$use_correlations && !is.null(data$multinichenet_output$lr_target_prior_cor)) {
        # Filtrer par corrélation
        lr_target_prior_cor_filtered <- data$multinichenet_output$prioritization_tables$group_prioritization_tbl$group %>%
          unique() %>%
          lapply(function(group_oi) {
            lr_target_prior_cor_filtered <- data$multinichenet_output$lr_target_prior_cor %>%
              dplyr::inner_join(
                data$multinichenet_output$ligand_activities_targets_DEgenes$ligand_activities %>%
                  dplyr::distinct(ligand, target, direction_regulation, contrast)
              ) %>%
              dplyr::inner_join(data$contrast_tbl) %>%
              dplyr::filter(group == group_oi)

            # Filter up-regulated
            lr_target_prior_cor_filtered_up <- lr_target_prior_cor_filtered %>%
              dplyr::filter(direction_regulation == "up") %>%
              dplyr::filter((rank_of_target < input$top_n_target) &
                              (pearson > input$correlation_threshold))

            # Filter down-regulated
            lr_target_prior_cor_filtered_down <- lr_target_prior_cor_filtered %>%
              dplyr::filter(direction_regulation == "down") %>%
              dplyr::filter((rank_of_target < input$top_n_target) &
                              (pearson < -input$correlation_threshold))

            # Combine results
            dplyr::bind_rows(
              lr_target_prior_cor_filtered_up,
              lr_target_prior_cor_filtered_down
            )
          }) %>% dplyr::bind_rows()

        lr_target_df <- lr_target_prior_cor_filtered %>%
          dplyr::distinct(group, sender, receiver, ligand, receptor, id, target, direction_regulation)

      } else {
        # Utiliser la méthode standard sans corrélation
        lr_target_prior <- prioritized_tbl_oi_all %>%
          dplyr::inner_join(
            data$multinichenet_output$ligand_activities_targets_DEgenes$ligand_activities %>%
              dplyr::distinct(ligand, target, direction_regulation, contrast) %>%
              dplyr::inner_join(data$contrast_tbl) %>%
              dplyr::ungroup()
          )

        lr_target_df <- lr_target_prior %>%
          dplyr::distinct(group, sender, receiver, ligand, receptor, id, target, direction_regulation)
      }

      # Generate network
      network <- infer_intercellular_regulatory_network(lr_target_df, prioritized_tbl_oi_all)

      # Préparer les couleurs
      unique_senders <- unique(network$nodes$celltype)
      colors_sender <- setNames(
        RColorBrewer::brewer.pal(n = min(length(unique_senders), 8), name = "Set2"),
        unique_senders
      )
      colors_sender[input$sender_to_color] <- input$sender_color

      # Visualize network
      network_graph <- visualize_network(network, colors_sender)
      data$current_network <- network_graph$plot

      # Display network
      output$regulatory_network_plot <- renderPlot({
        data$current_network
      })

      showNotification("Network generated successfully!", type = "message")

    }, error = function(e) {
      showNotification(paste("Error generating network:", e$message), type = "error")
      print("Full error:")
      print(str(e))
    }, finally = {
      removeModal()
    })
  })


  output$download_network <- downloadHandler(
    filename = function() {
      paste0("regulatory_network_", format(Sys.time(), "%Y%m%d"), ".tiff")
    },
    content = function(file) {
      tryCatch({
        req(data$current_network)  # Ensure the network exists
        ggsave(file, plot = data$current_network, device = "tiff", width = 12, height = 10, dpi = input$network_plot_dpi)
      }, error = function(e) {
        showNotification(paste("Download error:", e$message), type = "error")
      })
    }
  )



}
