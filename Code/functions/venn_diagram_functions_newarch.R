# venn_diagram_functions_newarch.R

safeIntersect <- function(x, y) {
  # Safe intersection function
  if (length(x) > 0 && length(y) > 0) {
    return(intersect(x, y))
  } else {
    return(character(0))
  }
}

safeSetdiff <- function(x, y) {
  # Safe set difference function
  if (length(x) > 0) {
    if (length(y) > 0) {
      return(setdiff(x, y))
    } else {
      return(x)
    }
  } else {
    return(character(0))
  }
}

calculateOverlaps <- function(gene_lists, list_names) {
  # Function to calculate overlaps between sets
  if (length(gene_lists) < 2) {
    return(NULL)
  }
  
  # All individual sets
  all_combinations <- list()
  for (i in 1:length(gene_lists)) {
    name <- list_names[i]
    all_combinations[[name]] <- gene_lists[[i]]
  }
  
  # Intersections between sets
  if (length(gene_lists) >= 2) {
    for (i in 2:length(gene_lists)) {
      # All combinations of i elements
      combos <- combn(1:length(gene_lists), i, simplify = FALSE)
      for (combo in combos) {
        # Start with first set
        intersection <- gene_lists[[combo[1]]]
        name_parts <- list_names[combo[1]]
        
        # Add intersections with other sets
        for (j in 2:length(combo)) {
          intersection <- safeIntersect(intersection, gene_lists[[combo[j]]])
          name_parts <- c(name_parts, list_names[combo[j]])
        }
        
        combo_name <- paste(name_parts, collapse = " & ")
        all_combinations[[combo_name]] <- intersection
      }
    }
  }
  
  return(all_combinations)
}

cleanGeneNamesForHtml <- function(genes) {
  # Clean gene names for HTML display
  cleaned <- gsub("[^a-zA-Z0-9_.-]", "_", genes)
  return(cleaned)
}

storeDETable <- function(table_storage, table_data, table_name, description, type, parameters) {
  # Store differential expression table for Venn diagram analysis
  if (is.null(table_data) || nrow(table_data) == 0) {
    showNotification("Cannot store an empty table.", type = "warning")
    return(list(success = FALSE, storage = table_storage))
  }
  
  # Create entry for storage
  table_entry <- list(
    name = table_name,
    description = description,
    data = table_data,
    type = type,
    gene_list = rownames(table_data),
    timestamp = Sys.time(),
    parameters = parameters
  )
  
  # Add to storage
  current_tables <- table_storage
  current_tables[[table_name]] <- table_entry
  
  showNotification(paste0("Stored ", nrow(table_data), " genes in '", description, "'"), type = "message")
  return(list(success = TRUE, storage = current_tables))
}

filterGeneList <- function(table_storage, table_name, significant_only, log_fc_threshold, p_val_threshold, use_adjusted_p, direction) {
  # Filter gene list based on user parameters
  if (table_name == "" || is.null(table_storage[[table_name]])) {
    return(character(0))
  }
  
  # Get the gene table
  gene_table <- table_storage[[table_name]]$data
  
  # Convert scientifically formatted p-values back to numeric
  gene_table$p_val <- as.numeric(gene_table$p_val)
  gene_table$p_val_adj <- as.numeric(gene_table$p_val_adj)
  
  # Select genes based on significance
  if (significant_only) {
    p_val_col <- if (use_adjusted_p) "p_val_adj" else "p_val"
    gene_table <- gene_table[gene_table[[p_val_col]] < p_val_threshold, ]
  }
  
  # Select genes based on fold change direction
  if (direction == "up") {
    gene_table <- gene_table[gene_table$avg_log2FC > log_fc_threshold, ]
  } else if (direction == "down") {
    gene_table <- gene_table[gene_table$avg_log2FC < -log_fc_threshold, ]
  } else { # "both"
    gene_table <- gene_table[abs(gene_table$avg_log2FC) > log_fc_threshold, ]
  }
  
  # Return gene names
  if (nrow(gene_table) > 0) {
    return(rownames(gene_table))
  } else {
    return(character(0))
  }
}

generateVennDiagram <- function(gene_lists, list_names, colors = c("#FF6B6B", "#4ECDC4", "#45B7D1")) {
  # Generate Venn diagram object
  n_sets <- length(list_names)
  
  if (n_sets == 2) {
    # For 2-set Venn diagram
    venn_obj <- venn.diagram(
      x = list(
        A = gene_lists[[1]],
        B = gene_lists[[2]]
      ),
      category.names = list_names,
      filename = NULL,
      output = TRUE,
      lwd = 2,
      lty = 'blank',
      fill = colors[1:2],
      alpha = 0.5,
      cex = 2,
      fontfamily = "sans",
      cat.cex = 1.5,
      cat.fontfamily = "sans",
      cat.default.pos = "outer",
      cat.dist = 0.05,
      cat.pos = c(-20, 20),
      euler.d = FALSE,
      scaled = FALSE
    )
  } else if (n_sets == 3) {
    # For 3-set Venn diagram
    venn_obj <- venn.diagram(
      x = list(
        A = gene_lists[[1]],
        B = gene_lists[[2]],
        C = gene_lists[[3]]
      ),
      category.names = list_names,
      filename = NULL,
      output = TRUE,
      lwd = 2,
      lty = 'blank',
      fill = colors[1:3],
      alpha = 0.5,
      cex = 1.5,
      fontfamily = "sans",
      cat.cex = 1.2,
      cat.fontfamily = "sans",
      cat.default.pos = "outer",
      euler.d = FALSE,
      scaled = FALSE
    )
  } else {
    stop("Venn diagrams support only 2 or 3 sets")
  }
  
  return(venn_obj)
}

updateVennSelectInputs <- function(session, table_storage, input_ids) {
  # Update Venn diagram select inputs with available tables
  if (length(table_storage) > 0) {
    table_choices <- c("None" = "", names(table_storage))
    for (input_id in input_ids) {
      updateSelectInput(session, input_id, choices = table_choices)
    }
  }
}

processVennGeneration <- function(table_storage, selected_tables, filter_params, colors) {
  # Process complete Venn diagram generation
  tryCatch({
    # Validate inputs
    selected_tables <- selected_tables[selected_tables != ""]
    if (length(selected_tables) < 2) {
      return(list(success = FALSE, message = "Please select at least 2 gene lists for comparison"))
    }
    
    # Filter gene lists based on user parameters
    gene_lists <- list()
    list_names <- character(length(selected_tables))
    
    for (i in 1:length(selected_tables)) {
      table_name <- selected_tables[i]
      gene_lists[[i]] <- filterGeneList(
        table_storage, 
        table_name,
        filter_params$significant_only[[i]],
        filter_params$log_fc_threshold[[i]],
        filter_params$p_val_threshold,
        filter_params$use_adjusted_p,
        filter_params$direction
      )
      list_names[i] <- table_storage[[table_name]]$description
    }
    
    # Calculate all overlaps
    overlaps <- calculateOverlaps(gene_lists, list_names)
    
    # Generate Venn diagram object
    venn_obj <- generateVennDiagram(gene_lists, list_names, colors)
    
    return(list(
      success = TRUE,
      venn_plot = venn_obj,
      overlaps = overlaps,
      list_names = list_names
    ))
    
  }, error = function(e) {
    return(list(success = FALSE, message = paste("Error generating Venn diagram:", e$message)))
  })
}