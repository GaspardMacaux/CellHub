# cells_genes_expressions_newarch.R
# Modular functions for gene expression, co-expression, and cluster composition analysis
# Works with both single and integrated datasets

# Function to create visual percentage bar for expression
create_expression_percentage_bar <- function(percentage) {
  # Color gradient based on expression percentage
  if (percentage > 75) {
    color <- "#28a745"  # Green for high expression
  } else if (percentage > 50) {
    color <- "#17a2b8"  # Cyan for moderate-high
  } else if (percentage > 25) {
    color <- "#ffc107"  # Yellow for moderate
  } else if (percentage > 10) {
    color <- "#fd7e14"  # Orange for low
  } else {
    color <- "#dc3545"  # Red for very low
  }
  
  bar_html <- paste0(
    '<div style="width: 100%; height: 20px; border: 1px solid #ddd; ',
    'border-radius: 4px; overflow: hidden; background-color: #f8f9fa;">',
    '<div style="width: ', percentage, '%; height: 100%; ',
    'background: linear-gradient(90deg, ', color, ' 0%, ', 
    adjustcolor(color, alpha.f = 0.7), ' 100%); ',
    'display: flex; align-items: center; justify-content: center; ',
    'transition: width 0.3s ease;" ',
    'title="', percentage, '% of cells express this gene">',
    '<span style="font-size: 11px; color: white; font-weight: bold; ',
    'text-shadow: 1px 1px 2px rgba(0,0,0,0.5);">',
    if (percentage > 15) paste0(percentage, '%') else '',
    '</span>',
    '</div></div>'
  )
  
  return(bar_html)
}

# Function to create size visualization bar for cluster composition
create_size_bar <- function(percentage) {
  # Color gradient based on size
  if (percentage > 20) {
    color <- "#2E86AB"  # Dark blue for large clusters
  } else if (percentage > 10) {
    color <- "#A23B72"  # Purple for medium clusters
  } else if (percentage > 5) {
    color <- "#F18F01"  # Orange for small clusters
  } else {
    color <- "#C73E1D"  # Red for very small clusters
  }
  
  bar_html <- paste0(
    '<div style="width: 100%; height: 20px; border: 1px solid #ccc; ',
    'display: flex; align-items: center; background-color: #f5f5f5;">',
    '<div style="width: ', min(percentage * 4, 100), '%; height: 18px; ',
    'background-color: ', color, '; margin: 1px; ',
    'display: flex; align-items: center; justify-content: center;" ',
    'title="', percentage, '% of total cells">',
    if (percentage > 5) paste0('<span style="font-size: 10px; color: white; font-weight: bold;">', 
                               percentage, '%</span>') else '',
    '</div></div>'
  )
  
  return(bar_html)
}

# Helper function to create co-expression visualization bar
create_coexpression_bar <- function(pct_both, pct_only1, pct_only2, pct_neither) {
  bar_html <- paste0(
    '<div style="width: 150px; height: 20px; display: flex; border: 1px solid #ddd; border-radius: 3px; overflow: hidden;">',
    '<div style="width: ', pct_neither, '%; background: #e0e0e0;" title="Neither: ', pct_neither, '%"></div>',
    '<div style="width: ', pct_only2, '%; background: #87CEEB;" title="Only Gene 2: ', pct_only2, '%"></div>',
    '<div style="width: ', pct_only1, '%; background: #FFA07A;" title="Only Gene 1: ', pct_only1, '%"></div>',
    '<div style="width: ', pct_both, '%; background: #32CD32;" title="Both: ', pct_both, '%"></div>',
    '</div>'
  )
  return(bar_html)
}

# Main function to analyze gene expression with enhanced metrics
analyze_gene_expression <- function(seurat_obj, 
                                    selected_genes, 
                                    assay_name = "RNA",
                                    expression_threshold = 0.1,
                                    is_integrated = FALSE) {
  
  tryCatch({
    # Set assay
    DefaultAssay(seurat_obj) <- assay_name
    
    # Get normalized data
    gene_data <- GetAssayData(seurat_obj, assay = assay_name, slot = "data")
    
    # Get available genes
    available_genes <- rownames(gene_data)
    selected_genes <- intersect(selected_genes, available_genes)
    
    if (length(selected_genes) == 0) {
      stop("No selected genes are present in the dataset")
    }
    
    # Get metadata
    cluster_info <- data.frame(
      cell = colnames(seurat_obj),
      cluster = as.character(Idents(seurat_obj)),
      stringsAsFactors = FALSE
    )
    
    # For integrated datasets, add dataset info
    if (is_integrated && "dataset" %in% colnames(seurat_obj@meta.data)) {
      cluster_info$dataset <- seurat_obj@meta.data$dataset
    }
    
    # Count total cells
    if (is_integrated) {
      # Count by cluster and dataset
      cluster_counts_total <- cluster_info %>%
        group_by(cluster, dataset) %>%
        summarise(n = n(), .groups = 'drop')
    } else {
      # Count by cluster only
      cluster_counts_total <- table(cluster_info$cluster)
    }
    
    # Process each gene
    expression_summary_list <- lapply(selected_genes, function(gene) {
      gene_expression <- gene_data[gene, ]
      expressed_indices <- gene_expression > expression_threshold
      
      if (sum(expressed_indices) == 0) {
        return(NULL)
      }
      
      # Get cells expressing the gene
      expressing_cells <- cluster_info[expressed_indices, ]
      
      if (is_integrated) {
        # Group by cluster and dataset for integrated data
        expression_summary <- expressing_cells %>%
          group_by(cluster, dataset) %>%
          summarise(
            Cells_Expressed = n(),
            .groups = 'drop'
          )
        
        # Calculate mean expression for expressing cells
        mean_expr_by_group <- sapply(1:nrow(expression_summary), function(i) {
          cells_in_group <- which(
            cluster_info$cluster == expression_summary$cluster[i] & 
              cluster_info$dataset == expression_summary$dataset[i] &
              expressed_indices
          )
          mean(gene_expression[cells_in_group])
        })
        expression_summary$Mean_Expression <- round(mean_expr_by_group, 3)
        
        # Merge with total counts
        expression_summary <- merge(
          expression_summary,
          cluster_counts_total,
          by = c("cluster", "dataset"),
          all.x = TRUE
        )
        names(expression_summary)[which(names(expression_summary) == "n")] <- "Total_Cells_in_Cluster"
        
        # Calculate dataset-specific totals
        dataset_totals <- cluster_info %>%
          group_by(dataset) %>%
          summarise(Total_Cells_Dataset = n(), .groups = 'drop')
        
        expression_summary <- merge(expression_summary, dataset_totals, by = "dataset")
        
        # Calculate percentages
        expression_summary$Percentage_of_Cluster <- round(
          (expression_summary$Cells_Expressed / expression_summary$Total_Cells_in_Cluster) * 100, 2
        )
        expression_summary$Percentage_of_Dataset <- round(
          (expression_summary$Cells_Expressed / expression_summary$Total_Cells_Dataset) * 100, 2
        )
        
      } else {
        # Standard processing for non-integrated data
        cluster_counts_expressed <- table(expressing_cells$cluster)
        
        # Calculate mean expression by cluster for expressing cells
        mean_expression_by_cluster <- aggregate(
          gene_expression[expressed_indices],
          by = list(Cluster = expressing_cells$cluster),
          FUN = mean
        )
        
        expression_summary <- data.frame(
          cluster = names(cluster_counts_expressed),
          Cells_Expressed = as.numeric(cluster_counts_expressed),
          Total_Cells_in_Cluster = as.numeric(cluster_counts_total[names(cluster_counts_expressed)]),
          Mean_Expression = round(mean_expression_by_cluster$x[match(
            names(cluster_counts_expressed), 
            mean_expression_by_cluster$Cluster
          )], 3),
          stringsAsFactors = FALSE
        )
        
        expression_summary$Percentage_of_Cluster <- round(
          (expression_summary$Cells_Expressed / expression_summary$Total_Cells_in_Cluster) * 100, 2
        )
      }
      
      # Add gene name
      expression_summary$Gene <- gene
      
      # Add visual elements
      expression_summary$Expression_Bar <- sapply(
        expression_summary$Percentage_of_Cluster, 
        create_expression_percentage_bar
      )
      
      # Create cell count summary
      expression_summary$Cell_Summary <- paste0(
        expression_summary$Cells_Expressed, "/", 
        expression_summary$Total_Cells_in_Cluster, " cells"
      )
      
      return(expression_summary)
    })
    
    # Combine results
    valid_summaries <- Filter(Negate(is.null), expression_summary_list)
    
    if (length(valid_summaries) == 0) {
      stop("No cells express the selected genes above the threshold")
    }
    
    expression_df <- do.call(rbind, valid_summaries)
    
    # Reorder columns based on dataset type
    if (is_integrated) {
      col_order <- c("Gene", "cluster", "dataset", "Cell_Summary", 
                     "Cells_Expressed", "Total_Cells_in_Cluster", 
                     "Percentage_of_Cluster", "Percentage_of_Dataset",
                     "Expression_Bar", "Mean_Expression")
    } else {
      col_order <- c("Gene", "cluster", "Cell_Summary",
                     "Cells_Expressed", "Total_Cells_in_Cluster", 
                     "Percentage_of_Cluster", "Expression_Bar", 
                     "Mean_Expression")
    }
    
    expression_df <- expression_df[, col_order]
    
    # Sort by gene and percentage
    expression_df <- expression_df[order(expression_df$Gene, -expression_df$Percentage_of_Cluster), ]
    
    return(list(
      data = expression_df,
      threshold = expression_threshold,
      is_integrated = is_integrated
    ))
    
  }, error = function(e) {
    stop(paste("Error in analyze_gene_expression:", e$message))
  })
}

# Enhanced function to analyze gene co-expression with better visualizations
analyze_gene_coexpression <- function(seurat_obj, 
                                      genes, 
                                      assay_name = "RNA",
                                      expression_threshold = 0,
                                      is_integrated = FALSE) {
  
  tryCatch({
    # Ensure genes is a character vector
    if (is.character(genes) && length(genes) == 1) {
      genes <- trimws(unlist(strsplit(genes, ",")))
    }
    genes <- genes[genes != ""]
    
    if (length(genes) < 2) {
      stop("Please provide at least 2 genes for co-expression analysis")
    }
    
    # Set assay
    DefaultAssay(seurat_obj) <- assay_name
    
    # Get expression data
    expr_data <- GetAssayData(seurat_obj, slot = "data", assay = assay_name)
    
    # Validate genes
    available_genes <- rownames(expr_data)
    missing_genes <- setdiff(genes, available_genes)
    
    if (length(missing_genes) > 0) {
      warning(paste("Missing genes:", paste(missing_genes, collapse = ", ")))
      genes <- intersect(genes, available_genes)
    }
    
    if (length(genes) < 2) {
      stop("Less than 2 valid genes found in dataset")
    }
    
    # Extract expression for selected genes
    gene_expr <- as.matrix(expr_data[genes, , drop = FALSE])
    
    # Create binary expression matrix
    gene_binary <- gene_expr > expression_threshold
    
    # Get metadata
    metadata <- data.frame(
      cell = colnames(seurat_obj),
      cluster = as.character(Idents(seurat_obj)),
      stringsAsFactors = FALSE
    )
    
    # Add dataset info for integrated data
    if (is_integrated && "dataset" %in% colnames(seurat_obj@meta.data)) {
      metadata$dataset <- seurat_obj@meta.data$dataset
    }
    
    # Function to analyze co-expression patterns
    analyze_group <- function(group_cells, group_name, dataset_name = NULL) {
      group_binary <- gene_binary[, group_cells, drop = FALSE]
      
      results <- data.frame(
        Group = group_name,
        N_Cells_Total = length(group_cells),
        stringsAsFactors = FALSE
      )
      
      if (!is.null(dataset_name)) {
        results$Dataset <- dataset_name
      }
      
      # Add individual gene counts with visual bars
      for (gene in genes) {
        col_name <- paste0("N_", gene, "_Positive")
        results[[col_name]] <- sum(group_binary[gene, ])
        
        # Add percentage
        pct_col <- paste0("Pct_", gene, "_Positive")
        pct_value <- round(100 * results[[col_name]] / results$N_Cells_Total, 1)
        results[[pct_col]] <- pct_value
        
        # Add visual bar for each gene
        bar_col <- paste0(gene, "_Bar")
        results[[bar_col]] <- create_expression_percentage_bar(pct_value)
      }
      
      # Calculate co-expression patterns
      if (length(genes) == 2) {
        # For 2 genes
        results$N_Both_Positive <- sum(group_binary[genes[1], ] & group_binary[genes[2], ])
        results$N_Only_Gene1 <- sum(group_binary[genes[1], ] & !group_binary[genes[2], ])
        results$N_Only_Gene2 <- sum(!group_binary[genes[1], ] & group_binary[genes[2], ])
        results$N_Neither <- sum(!group_binary[genes[1], ] & !group_binary[genes[2], ])
        
        # Add percentages
        results$Pct_Both <- round(100 * results$N_Both_Positive / results$N_Cells_Total, 1)
        results$Pct_Only_Gene1 <- round(100 * results$N_Only_Gene1 / results$N_Cells_Total, 1)
        results$Pct_Only_Gene2 <- round(100 * results$N_Only_Gene2 / results$N_Cells_Total, 1)
        results$Pct_Neither <- round(100 * results$N_Neither / results$N_Cells_Total, 1)
        
        # Add summary text
        results$Coexpression_Summary <- paste0(
          "Both: ", results$N_Both_Positive, " (", results$Pct_Both, "%) | ",
          "Only ", genes[1], ": ", results$N_Only_Gene1, " (", results$Pct_Only_Gene1, "%) | ",
          "Only ", genes[2], ": ", results$N_Only_Gene2, " (", results$Pct_Only_Gene2, "%)"
        )
        
      } else {
        # For multiple genes
        n_genes_expressed <- colSums(group_binary)
        
        results$N_All_Positive <- sum(n_genes_expressed == length(genes))
        results$N_Any_Positive <- sum(n_genes_expressed > 0)
        results$N_None_Positive <- sum(n_genes_expressed == 0)
        
        results$Pct_All <- round(100 * results$N_All_Positive / results$N_Cells_Total, 1)
        results$Pct_Any <- round(100 * results$N_Any_Positive / results$N_Cells_Total, 1)
        results$Pct_None <- round(100 * results$N_None_Positive / results$N_Cells_Total, 1)
        
        # Distribution of number of genes expressed
        for (i in 0:length(genes)) {
          col_name <- paste0("N_Expressing_", i, "_Genes")
          results[[col_name]] <- sum(n_genes_expressed == i)
        }
        
        # Add summary
        results$Coexpression_Summary <- paste0(
          "All genes: ", results$N_All_Positive, " (", results$Pct_All, "%) | ",
          "Any gene: ", results$N_Any_Positive, " (", results$Pct_Any, "%)"
        )
      }
      
      return(results)
    }
    
    # Analyze by groups
    if (is_integrated) {
      # Analyze by cluster and dataset
      all_results <- metadata %>%
        group_by(cluster, dataset) %>%
        group_map(~ analyze_group(.x$cell, .y$cluster, .y$dataset)) %>%
        bind_rows()
      
      # Add overall by dataset
      dataset_overall <- metadata %>%
        group_by(dataset) %>%
        group_map(~ analyze_group(.x$cell, paste("Overall", .y$dataset), .y$dataset)) %>%
        bind_rows()
      
      all_results <- rbind(all_results, dataset_overall)
      
    } else {
      # Analyze by cluster only
      all_clusters <- unique(metadata$cluster)
      cluster_results <- lapply(all_clusters, function(cl) {
        cells <- metadata$cell[metadata$cluster == cl]
        analyze_group(cells, cl)
      })
      all_results <- do.call(rbind, cluster_results)
      
      # Add overall
      overall_results <- analyze_group(metadata$cell, "Overall")
      all_results <- rbind(overall_results, all_results)
    }
    
    # Add visual elements for key percentages
    if (length(genes) == 2) {
      all_results$Coexpression_Visual <- mapply(
        function(both, only1, only2, neither) {
          create_coexpression_bar(both, only1, only2, neither)
        },
        all_results$Pct_Both,
        all_results$Pct_Only_Gene1,
        all_results$Pct_Only_Gene2,
        all_results$Pct_Neither
      )
    }
    
    return(list(
      data = all_results,
      genes_analyzed = genes,
      threshold_used = expression_threshold,
      is_integrated = is_integrated
    ))
    
  }, error = function(e) {
    stop(paste("Error in analyze_gene_coexpression:", e$message))
  })
}

# Function to create cluster composition table with support for integrated datasets
create_cluster_composition_table <- function(seurat_obj, is_integrated = FALSE) {
  # Get cluster information
  cluster_data <- data.frame(
    cluster = as.character(Idents(seurat_obj)),
    nFeature = seurat_obj$nFeature_RNA,
    nCount = seurat_obj$nCount_RNA,
    percent_mt = seurat_obj$percent.mt,
    stringsAsFactors = FALSE
  )
  
  # Add dataset info for integrated data
  if (is_integrated && "dataset" %in% colnames(seurat_obj@meta.data)) {
    cluster_data$dataset <- seurat_obj@meta.data$dataset
  }
  
  # Calculate statistics
  if (is_integrated) {
    # Group by cluster and dataset
    cluster_stats <- cluster_data %>%
      group_by(cluster, dataset) %>%
      summarise(
        Cell_Count = n(),
        Mean_Genes = round(mean(nFeature, na.rm = TRUE), 1),
        Median_Genes = round(median(nFeature, na.rm = TRUE), 1),
        Mean_UMI = round(mean(nCount, na.rm = TRUE), 1),
        Median_UMI = round(median(nCount, na.rm = TRUE), 1),
        Mean_MT_Percent = round(mean(percent_mt, na.rm = TRUE), 2),
        .groups = 'drop'
      )
    
    # Calculate percentages within dataset
    dataset_totals <- cluster_data %>%
      group_by(dataset) %>%
      summarise(Total_in_Dataset = n(), .groups = 'drop')
    
    cluster_stats <- merge(cluster_stats, dataset_totals, by = "dataset")
    cluster_stats$Percent_of_Dataset <- round(
      (cluster_stats$Cell_Count / cluster_stats$Total_in_Dataset) * 100, 1
    )
    
  } else {
    # Standard single dataset processing
    cluster_stats <- cluster_data %>%
      group_by(cluster) %>%
      summarise(
        Cell_Count = n(),
        Mean_Genes = round(mean(nFeature, na.rm = TRUE), 1),
        Median_Genes = round(median(nFeature, na.rm = TRUE), 1),
        Mean_UMI = round(mean(nCount, na.rm = TRUE), 1),
        Median_UMI = round(median(nCount, na.rm = TRUE), 1),
        Mean_MT_Percent = round(mean(percent_mt, na.rm = TRUE), 2),
        .groups = 'drop'
      )
  }
  
  # Calculate total cells for percentage
  total_cells <- nrow(cluster_data)
  cluster_stats$Percent_of_Total <- round((cluster_stats$Cell_Count / total_cells) * 100, 1)
  
  # Create visual representation of cluster sizes
  cluster_stats$Size_Bar <- sapply(cluster_stats$Percent_of_Total, create_size_bar)
  
  # Create quality summary
  cluster_stats$Quality_Summary <- paste0(
    "Genes: ", cluster_stats$Mean_Genes, " | ",
    "UMI: ", cluster_stats$Mean_UMI, " | ",
    "MT: ", cluster_stats$Mean_MT_Percent, "%"
  )
  
  # Reorder columns based on dataset type
  if (is_integrated) {
    cluster_stats <- cluster_stats[, c("cluster", "dataset", "Cell_Count", 
                                       "Percent_of_Total", "Percent_of_Dataset", "Size_Bar", 
                                       "Mean_Genes", "Median_Genes", "Mean_UMI", "Median_UMI", 
                                       "Mean_MT_Percent", "Quality_Summary")]
  } else {
    cluster_stats <- cluster_stats[, c("cluster", "Cell_Count", "Percent_of_Total", "Size_Bar", 
                                       "Mean_Genes", "Median_Genes", "Mean_UMI", "Median_UMI", 
                                       "Mean_MT_Percent", "Quality_Summary")]
  }
  
  return(cluster_stats)
}

# Enhanced render function for expression analysis table
render_expression_table <- function(expression_data, table_id) {
  
  is_integrated <- expression_data$is_integrated
  df <- expression_data$data
  
  # Define column visibility and names
  if (is_integrated) {
    hidden_cols <- c("Cells_Expressed", "Total_Cells_in_Cluster")
    col_names <- c('Gene', 'Cluster', 'Dataset', 'Cells Expressing', 
                   'Cells Expressed', 'Total Cells', '% in Cluster', 
                   '% in Dataset', 'Expression Visual', 'Mean Expression')
  } else {
    hidden_cols <- c("Cells_Expressed", "Total_Cells_in_Cluster")
    col_names <- c('Gene', 'Cluster', 'Cells Expressing',
                   'Cells Expressed', 'Total Cells', '% Expressing', 
                   'Expression Visual', 'Mean Expression')
  }
  
  # Create column indices for hiding
  all_cols <- names(df)
  hidden_indices <- which(all_cols %in% hidden_cols) - 1
  
  dt <- datatable(
    df,
    options = list(
      pageLength = 15,
      scrollX = TRUE,
      columnDefs = list(
        list(visible = FALSE, targets = hidden_indices),
        list(className = 'dt-center', targets = '_all'),
        list(width = '120px', targets = which(all_cols == "Expression_Bar") - 1),
        list(width = '100px', targets = which(all_cols == "Cell_Summary") - 1)
      ),
      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel'),
      order = if (is_integrated) {
        list(list(0, 'asc'), list(6, 'desc'))
      } else {
        list(list(0, 'asc'), list(5, 'desc'))
      }
    ),
    rownames = FALSE,
    caption = paste("Gene Expression Analysis (Threshold:", expression_data$threshold, ")"),
    escape = FALSE,
    colnames = col_names
  ) %>%
    formatStyle(
      'Mean_Expression',
      backgroundColor = styleInterval(
        c(0.5, 1, 2, 5), 
        c('#f0f0f0', '#d4edda', '#fff3cd', '#f8d7da', '#d1ecf1')
      )
    )
  
  # Add styling for dataset column if integrated
  if (is_integrated && "dataset" %in% names(df)) {
    dt <- dt %>%
      formatStyle(
        'dataset',
        backgroundColor = styleEqual(
          unique(df$dataset),
          rainbow(length(unique(df$dataset)), alpha = 0.3)
        )
      )
  }
  
  return(dt)
}

# Enhanced render function for co-expression table
render_coexpression_table <- function(coexpression_data, table_id) {
  
  df <- coexpression_data$data
  is_integrated <- coexpression_data$is_integrated
  genes <- coexpression_data$genes_analyzed
  
  # Define which columns to hide (the N_ columns)
  hidden_cols <- names(df)[grepl("^N_", names(df))]
  hidden_indices <- which(names(df) %in% hidden_cols) - 1
  
  dt <- datatable(
    df,
    options = list(
      pageLength = 20,
      scrollX = TRUE,
      columnDefs = list(
        list(visible = FALSE, targets = hidden_indices),
        list(className = 'dt-center', targets = '_all'),
        list(className = 'dt-left', targets = which(names(df) == "Coexpression_Summary") - 1)
      ),
      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel')
    ),
    rownames = FALSE,
    caption = paste("Co-expression Analysis -", paste(genes, collapse = ", ")),
    escape = FALSE
  ) %>%
    formatStyle(
      'Group',
      target = 'row',
      backgroundColor = styleEqual(
        grep("Overall", df$Group, value = TRUE),
        rep('#ffffcc', length(grep("Overall", df$Group)))
      ),
      fontWeight = styleEqual(
        grep("Overall", df$Group, value = TRUE),
        rep('bold', length(grep("Overall", df$Group)))
      )
    )
  
  # Add percentage bars for all Pct_ columns
  pct_cols <- names(df)[grepl("^Pct_", names(df))]
  for (col in pct_cols) {
    if (col %in% names(df)) {
      dt <- dt %>%
        formatStyle(
          col,
          background = styleColorBar(c(0, 100), '#e8f4fd'),
          backgroundSize = '98% 88%',
          backgroundRepeat = 'no-repeat',
          backgroundPosition = 'center'
        )
    }
  }
  
  if (is_integrated && "Dataset" %in% names(df)) {
    dt <- dt %>%
      formatStyle(
        'Dataset',
        backgroundColor = styleEqual(
          unique(df$Dataset),
          rainbow(length(unique(df$Dataset)), alpha = 0.3)
        )
      )
  }
  
  return(dt)
}

# Render function for cluster composition table
render_cluster_composition_table <- function(cluster_data, is_integrated = FALSE) {
  
  # Define columns based on dataset type
  if (is_integrated) {
    col_names <- c('Cluster', 'Dataset', 'Cell Count', '% Total', '% Dataset', 'Size Visual', 
                   'Mean Genes', 'Median Genes', 'Mean UMI', 'Median UMI', 
                   'Mean MT%', 'Quality Summary')
    center_targets <- c(0, 1, 2, 3, 4, 6, 7, 8, 9, 10)
    left_targets <- c(5, 11)
  } else {
    col_names <- c('Cluster', 'Cell Count', '% Total', 'Size Visual', 
                   'Mean Genes', 'Median Genes', 'Mean UMI', 'Median UMI', 
                   'Mean MT%', 'Quality Summary')
    center_targets <- c(0, 1, 2, 4, 5, 6, 7, 8)
    left_targets <- c(3, 9)
  }
  
  dt <- datatable(
    cluster_data,
    options = list(
      pageLength = 15,
      scrollX = TRUE,
      columnDefs = list(
        list(className = 'dt-center', targets = center_targets),
        list(className = 'dt-left', targets = left_targets),
        list(width = '100px', targets = which(names(cluster_data) == "Size_Bar") - 1),
        list(width = '200px', targets = which(names(cluster_data) == "Quality_Summary") - 1)
      ),
      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel')
    ),
    rownames = FALSE,
    caption = "Cluster Composition and Quality Metrics",
    escape = FALSE,
    colnames = col_names
  ) %>%
    formatStyle(
      'Cell_Count',
      background = styleColorBar(range(cluster_data$Cell_Count), '#e8f4fd'),
      backgroundSize = '80% 90%',
      backgroundRepeat = 'no-repeat',
      backgroundPosition = 'center'
    ) %>%
    formatStyle(
      'Mean_MT_Percent',
      backgroundColor = styleInterval(c(5, 10, 20), c('#90EE90', '#FFFF99', '#FFB347', '#FF6B6B'))
    )
  
  # Add dataset coloring for integrated data
  if (is_integrated && "dataset" %in% names(cluster_data)) {
    dt <- dt %>%
      formatStyle(
        'dataset',
        backgroundColor = styleEqual(
          unique(cluster_data$dataset),
          rainbow(length(unique(cluster_data$dataset)), alpha = 0.3)
        )
      )
  }
  
  return(dt)
}

# Function to create co-expression plot
create_coexpression_plot <- function(results_data, genes) {
  
  # Remove "Overall" for plotting
  plot_data <- results_data[!grepl("Overall", results_data$Group), ]
  
  if (length(genes) == 2) {
    # For 2 genes, create stacked bar plot
    plot_data_long <- plot_data %>%
      dplyr::select(Group, N_Both_Positive, N_Only_Gene1, N_Only_Gene2, N_Neither) %>%
      tidyr::pivot_longer(cols = -Group, names_to = "Pattern", values_to = "Count")
    
    # Create proper labels after pivoting
    plot_data_long <- plot_data_long %>%
      dplyr::mutate(
        Pattern_Label = dplyr::case_when(
          Pattern == "N_Neither" ~ "Neither",
          Pattern == "N_Only_Gene1" ~ paste("Only", genes[1]),
          Pattern == "N_Only_Gene2" ~ paste("Only", genes[2]),
          Pattern == "N_Both_Positive" ~ "Both",
          TRUE ~ Pattern
        ),
        Pattern_Label = factor(Pattern_Label, 
                               levels = c("Neither", 
                                          paste("Only", genes[2]), 
                                          paste("Only", genes[1]), 
                                          "Both"))
      )
    
    # Create the plot with proper color mapping
    p <- ggplot(plot_data_long, aes(x = Group, y = Count, fill = Pattern_Label)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(
        values = setNames(
          c("gray80", "lightblue", "lightcoral", "darkgreen"),
          c("Neither", paste("Only", genes[2]), paste("Only", genes[1]), "Both")
        ),
        name = "Expression Pattern"
      ) +
      labs(title = paste("Co-expression of", genes[1], "and", genes[2]),
           y = "Number of cells",
           x = "Cluster") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
  } else {
    # For multiple genes, show distribution of number of genes expressed
    cols_to_plot <- grep("N_Expressing_\\d+_Genes", names(plot_data), value = TRUE)
    
    plot_data_long <- plot_data %>%
      dplyr::select(Group, dplyr::all_of(cols_to_plot)) %>%
      tidyr::pivot_longer(cols = -Group, names_to = "N_Genes", values_to = "Count") %>%
      dplyr::mutate(N_Genes = as.numeric(gsub("N_Expressing_(\\d+)_Genes", "\\1", N_Genes)))
    
    p <- ggplot(plot_data_long, aes(x = Group, y = Count, fill = factor(N_Genes))) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_viridis_d(name = "# Genes\nExpressed") +
      labs(title = paste("Expression of", length(genes), "genes"),
           subtitle = paste(genes, collapse = ", "),
           y = "Number of cells",
           x = "Cluster") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  return(p)
}


# Function to create co-expression summary statistics
create_coexpression_summary_stats <- function(data, genes) {
  if (length(genes) == 2) {
    # For 2 genes
    overall_data <- data[grepl("Overall", data$Group), ]
    summary_stats <- data.frame(Dataset = gsub("Overall ", "", overall_data$Group), Total_Cells = overall_data$N_Cells_Total, Both_Genes = paste0(overall_data$N_Both_Positive, " (", overall_data$Pct_Both, "%)"), Only_Gene1 = paste0(overall_data$N_Only_Gene1, " (", overall_data$Pct_Only_Gene1, "%)"), Only_Gene2 = paste0(overall_data$N_Only_Gene2, " (", overall_data$Pct_Only_Gene2, "%)"), Neither = paste0(overall_data$N_Neither, " (", overall_data$Pct_Neither, "%)"))
    names(summary_stats)[3:6] <- c("Both Genes", paste("Only", genes[1]), paste("Only", genes[2]), "Neither")
  } else {
    # For multiple genes
    overall_data <- data[grepl("Overall", data$Group), ]
    summary_stats <- data.frame(Dataset = gsub("Overall ", "", overall_data$Group), Total_Cells = overall_data$N_Cells_Total, All_Genes = paste0(overall_data$N_All_Positive, " (", overall_data$Pct_All, "%)"), Any_Gene = paste0(overall_data$N_Any_Positive, " (", overall_data$Pct_Any, "%)"), No_Genes = paste0(overall_data$N_None_Positive, " (", overall_data$Pct_None, "%)"))
  }
  return(summary_stats)
}


# Wrapper functions for backward compatibility
create_single_cluster_composition_table <- function(seurat_obj) {
  create_cluster_composition_table(seurat_obj, is_integrated = FALSE)
}

analyze_gene_coexpression_single <- function(seurat_obj, gene_text, expression_threshold = 0) {
  analyze_gene_coexpression(
    seurat_obj = seurat_obj,
    genes = gene_text,
    assay_name = DefaultAssay(seurat_obj),
    expression_threshold = expression_threshold,
    is_integrated = FALSE
  )
}

create_coexpression_plot_single <- function(results_data, genes) {
  create_coexpression_plot(results_data, genes)
}

# Helper function for formatted percentage display
format_percentage_compact <- function(n, total, color = "blue") {
  pct <- round(100 * n / total, 1)
  sprintf(
    '<div style="display: inline-block; position: relative; width: 60px;">
    <div style="position: absolute; width: %.1f%%; height: 20px; background-color: %s; opacity: 0.3;"></div>
    <div style="position: relative; text-align: center; line-height: 20px; font-weight: bold;">%.1f%%</div>
  </div>',
    pct, color, pct
  )
}