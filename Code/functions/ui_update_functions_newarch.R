# ui_update_functions_newarch.R

updateAssayChoices <- function(session, seurat_object, assay_input_id = "viz_assay") {
  # Update assay choices based on available assays in Seurat object
  # Args:
  #   session: Shiny session object
  #   seurat_object: Seurat object
  #   assay_input_id: ID of the selectInput to update
  
  if (is.null(seurat_object)) {
    return(NULL)
  }
  
  tryCatch({
    # Get available assays
    available_assays <- names(seurat_object@assays)
    
    # Debug message
    message(paste("Available assays:", paste(available_assays, collapse = ", ")))
    
    # Update selectInput
    updateSelectInput(session, assay_input_id,
                      choices = available_assays,
                      selected = if("RNA" %in% available_assays) "RNA" else available_assays[1])
    
  }, error = function(e) {
    message(paste("Error updating assay choices:", e$message))
  })
}

updateGeneChoices <- function(session, seurat_object, selected_assay, gene_input_ids = NULL) {
  # Update gene choices based on selected assay
  # Args:
  #   session: Shiny session object
  #   seurat_object: Seurat object
  #   selected_assay: Currently selected assay
  #   gene_input_ids: Vector of input IDs to update (if NULL, updates default ones)
  
  if (is.null(seurat_object) || is.null(selected_assay)) {
    return(NULL)
  }
  
  # Default gene input IDs if not specified
  if (is.null(gene_input_ids)) {
    gene_input_ids <- c("gene_select", "gene_select_heatmap", "gene_select_genes_analysis")
  }
  
  tryCatch({
    # Get gene list from selected assay
    gene_list <- sort(rownames(LayerData(seurat_object,
                                         assay = selected_assay,
                                         layer = 'counts')))
    
    # Update all specified gene inputs
    for (input_id in gene_input_ids) {
      updatePickerInput(session, input_id, choices = gene_list)
    }
    
    message(paste("Updated", length(gene_input_ids), "gene input(s) with", length(gene_list), "genes"))
    
  }, error = function(e) {
    message(paste("Error updating gene choices:", e$message))
  })
}

updateClusterChoices <- function(session, seurat_object, cluster_input_ids = NULL) {
  # Update cluster choices based on Seurat object identities
  # Args:
  #   session: Shiny session object
  #   seurat_object: Seurat object
  #   cluster_input_ids: List with input IDs and their types
  
  if (is.null(seurat_object)) {
    return(NULL)
  }
  
  # Default cluster input IDs if not specified
  if (is.null(cluster_input_ids)) {
    cluster_input_ids <- list(
      select = c("select_cluster", "ident_1"),
      checkbox = c("ident_2")
    )
  }
  
  tryCatch({
    # Get cluster choices
    cluster_choices <- unique(Idents(seurat_object))
    
    message(paste("Available clusters:", paste(cluster_choices, collapse = ", ")))
    
    # Update selectInput types
    if ("select" %in% names(cluster_input_ids)) {
      for (input_id in cluster_input_ids$select) {
        updateSelectInput(session, input_id, choices = cluster_choices)
      }
    }
    
    # Update checkboxGroupInput types
    if ("checkbox" %in% names(cluster_input_ids)) {
      for (input_id in cluster_input_ids$checkbox) {
        updateCheckboxGroupInput(session, input_id, choices = cluster_choices, selected = "")
      }
    }
    
  }, error = function(e) {
    message(paste("Error updating cluster choices:", e$message))
  })
}

createGeneUpdateFunction <- function(input) {
  # Returns a function that can access current input values
  function(session, new_gene, gene_text_input_ids = NULL) {
    if (is.null(new_gene) || new_gene == "") {
      return(NULL)
    }
    
    # Default gene text input IDs if not specified
    if (is.null(gene_text_input_ids)) {
      gene_text_input_ids <- list(
        gene_list_vln = input$gene_list_vln,
        gene_list_feature = input$gene_list_feature,
        gene_list_dotplot = input$gene_list_dotplot,
        gene_list_ridge_plot = input$gene_list_ridge_plot
      )
    }
    
    # Helper function to update genes
    update_genes <- function(existing_genes, new_gene) {
      if (is.null(existing_genes) || existing_genes == "") {
        return(new_gene)
      } else {
        genes <- strsplit(existing_genes, ",")[[1]]
        genes <- unique(c(trimws(genes), new_gene))
        return(paste(genes, collapse = ", "))
      }
    }
    
    # Update all specified text inputs
    for (input_id in names(gene_text_input_ids)) {
      current_value <- gene_text_input_ids[[input_id]]
      updateTextInput(session, input_id, value = update_genes(current_value, new_gene))
    }
  }
}


# Dans ui_update_functions_newarch.R, ajouter :
updateClusterTextInputs <- function(session, seurat_object, text_input_ids) {
  if (is.null(seurat_object)) return(NULL)
  tryCatch({
    cluster_list <- getClusters(seurat_object)
    cluster_text <- paste(cluster_list, collapse = ",")
    for (input_id in text_input_ids) {
      updateTextInput(session, input_id, value = cluster_text)
    }
    message(paste("Updated", length(text_input_ids), "cluster text input(s)"))
  }, error = function(e) {
    message(paste("Error updating cluster text inputs:", e$message))
  })
}
