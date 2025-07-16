
############################## User Interface ##############################


ui <- dashboardPage(

  ############################## HTML et CSS ##############################

  dashboardHeader(title = tags$img(src = "header.png",   height = "53px", width = "218px")),
  dashboardSidebar(
    tags$head(
      tags$style(HTML("
                        /* Change header color */
                        .main-header .navbar {
                          background-color: rgb(122, 207, 176) !important;
                        }
                        /* Aligner le contenu de la sidebar avec le contenu principal */
                        .sidebar .sidebar-menu {
                          margin-top: 5px;
                        }

                        /* Augmenter la largeur de la bulle d'information */
                        .popover {
                          width: 300px;
                        }

                        /* Style pour ajouter l'image au header */
                        .main-header {
                          background-image: url('www/header.png');
                          background-repeat: no-repeat;
                          background-position: center center;
                        }

                        .info-box {
                          width: 100%;
                        }
/* Style pour le plot carré réactif */
                .reactive-square-plot {
                    width: 100%;
                    height: auto;
                    aspect-ratio: 1 / 1;
                }
                      ")),
      tags$script(src = "my_script.js")
    ),
    useShinyjs(),

    ############################## Sidebar Menu ##############################
    sidebarMenu(
      menuItem("Introduction", tabName = "introduction", icon = icon("mug-hot")),
      menuItem(
        "Single Dataset Analysis", icon = icon("virus-covid"),
        menuItem("Load dataset", tabName = "load_dataset", icon = icon("download")),
        menuItem("Data cleanup & Variable features", tabName = "qc", icon = icon("bug")),
        menuItem("Dimensional reduction", tabName = "dimensional_reduction", icon = icon("home")),
        menuItem("Clustering", tabName = "clustering", icon = icon("brain")),
        menuItem("Doublet Detection", tabName = "doublet_detection",icon = icon("dna")),
        menuItem("Plot Gene expressions", tabName = "plots_genes_expressions", icon = icon("flask")),
        menuItem("Heat maps & Dual expression", tabName = "heatmaps_dualexpression", icon = icon("thermometer-half")),
        menuItem("Assigning cell types identity", tabName = "assigning_cell_type_identity", icon = icon("id-card")),
        menuItem("Clusters comparison", tabName = "cluster_comparison", icon = icon("earth-americas")),
        menuItem("Subset", tabName = "subset", icon = icon("project-diagram"))
      ),
      menuItem(
        "Multiple Datasets Analysis", icon = icon("viruses"),
        menuItem("Load datasets", tabName = "load_datasets_merge", icon = icon("download")),
        menuItem("Clustering", tabName = "clustering_merge", icon = icon("brain")),
        menuItem("Plot Gene expressions", tabName = "plots_genes_expressions_merge", icon = icon("flask")),
        menuItem("Heatmaps & Dual expression", tabName = "heatmaps_dualexpression_merge", icon = icon("thermometer-half")),
        menuItem("Assigning cell types identity", tabName = "assigning_cell_type_identity_merge", icon = icon("id-card")),
        menuItem("Clusters comparison", tabName = "cluster_comparison_merge", icon = icon("earth-americas")),
        menuItem("Subset", tabName = "subset_merge", icon = icon("project-diagram"))
      ),
      menuItem(
        "Cell Chat", icon = icon("envelope"),
        menuItem("Load dataset", tabName = "load_data_cellchat", icon = icon("download")),
        menuItem("Ligand-Receptor", tabName = "ligand_receptor_cellchat", icon = icon("link")),
        menuItem("Circle Plot", tabName = "circle_plot_cellchat", icon = icon("circle")),
        menuItem("Circle Plot Global", tabName = "circle_plot_global_cellchat", icon = icon("circle"))

      ),
      
      menuItem(
        "Transcriptomic Spatial", tabName = "transcriptomic_spatial", icon = icon("map"),
              menuItem("Load Spatial dataset", tabName = "load_spatial_dataset", icon = icon("download")),
               menuItem("Normalisation Step", tabName = "spatial_normalisation", icon = icon("chart-line")),
               menuItem("Clustering", tabName = "spatial_clustering", icon = icon("circle")),
               menuItem("Interactive Spatial Visualization", tabName = "spatial_interactive_visualization", icon = icon("magnifying-glass")),
              menuItem("Gene expression visualisation", tabName = "spatial_gene_expression_visualisation", icon = icon("database")),
               menuItem("Cluster Annotation", tabName = "spatial_cluster_annotation", icon = icon("tags")),
              menuItem(" Spatial Tissue Viewer", tabName = "spatial_tissue_viewer", icon = icon("magnifying-glass")),
               menuItem("Marker Analysis", tabName = "spatial_marker_analysis", icon = icon("dna"))
     
         ),
      menuItem(
        "Monocle", icon = icon("magnifying-glass"),
        menuItem("Load dataset", tabName = "trajectory", icon = icon("download")),
        menuItem("Trajectory Analysis", tabName = "genes_pseudotime", icon = icon("route"))
        #menuItem("Gene modules", tabName = "gene_modules", icon = icon("route"))

      ),
      menuItem("Acknowledgement & Licence", tabName = "acknowledgement", icon = icon("heart"))
    )
  ),

  ############################## Single/Introduction ##############################
  dashboardBody(
    tabItems(
      tabItem(tabName = "introduction",
              div(class = "well", style = "background: linear-gradient(rgba(0,0,0,0.7), rgba(0,0,0,0.7)), url('muscle.png');
                                   background-size: cover;
                                   background-position: center;
                                   color: white;
                                   padding: 40px;
                                   min-height: 800px;",
              div(class = "container-fluid", style = "padding: 0;",
               div(class = "well", style = "margin-bottom: 30px; background-color: #f8f9fa;",
                  h1("Welcome to Cell-Hub, the single-Cell and single-nucleus analysis tool", style = "color: #2c3e50; text-align: center; margin-bottom: 25px;"),
                  div(class = "row",
                      div(class = "col-md-8 col-md-offset-2",
                          h3("About This Application", style = "color: #34495e; margin-bottom: 20px;"),
                          div(style = "font-size: 16px; line-height: 1.6; color: #555;",
                              p(
                                tags$strong("Single-Cell and Single-Nucleus RNA Sequencing Analysis"),
                                "This application is designed for analyzing scRNA-seq and snRNA-seq data produced using 10x Genomics technology."
                              ),

                              h4("Key Features:", style = "color: #2980b9; margin-top: 20px;"),
                              tags$ul(
                                tags$li("Analysis of gene expression at single-cell resolution"),
                                tags$li("Built on powerful Seurat and Monocle libraries"),
                                tags$li("Complete suite of quality control and analysis tools"),
                                tags$li("Trajectory analysis for cell fate determination"),
                                tags$li("User-friendly interface requiring no coding experience")
                              ),

                              h4("Technologies Used:", style = "color: #2980b9; margin-top: 20px;"),
                              p("The application leverages four main libraries:"),
                              tags$ul(
                                tags$li(tags$strong("Seurat:"), "A comprehensive tool for scRNA-seq analysis, including quality control, visualization, and clustering."),
                                tags$li(tags$strong("Monocle:"), "Advanced trajectory analysis for understanding cell fate decisions and developmental processes."),
                                tags$li(tags$strong("NicheNet & MultiNicheNet:"), "Specialized tool for ligand-receptor analysis, predicting which ligands from sender cells influence target gene expression in receiver cells.")
                              )
                          )
                      )
                  ))
              ),
                  div(class = "well", style = "margin-top: 20px; background-color: white;",
                  h3("Analysis Pipeline", style = "color: #2c3e50; text-align: center; margin-bottom: 20px;"),
                  div(style = "text-align: center; max-width: 100%; overflow-x: auto;",
                      tags$img(src = 'pipeline.png',
                               style = "max-width: 110%; height: auto; border-radius: 5px; box-shadow: 0 2px 4px rgba(0,0,0,0.1);")
                  ),

              )
      ))
      ,
      tabItem(tabName = "load_dataset",
              fluidPage(
                div(class = "container-fluid",
                    div(class = "row",
                        div(class = "col-12 text-center",
                            div(style = "background-color: #f8f9fa; padding: 20px; margin-bottom: 20px; border-radius: 5px;",
                                h2("Single Dataset Analysis", style = "color: #2c3e50; font-weight: 500; margin-bottom: 10px;"),
                                p("Load and analyze your single-cell RNA sequencing data", style = "color: #7f8c8d; font-size: 16px;")
                            )
                        )
                    )
                ),
                div(class = "well",
                    h4("Data Loading Options", style = "color: #666;"),
                    p("This section allows you to load and analyze a single dataset. You have two options:"),
                    h5("Option 1: Load Raw 10X Data", style = "font-weight: bold; color: #444;"),
                    tags$ul(
                      tags$li("Create a ZIP file containing the three 10X Genomics output files:"),
                      tags$ul(
                        tags$li(tags$code("barcodes.tsv.gz"), "- Contains cell barcodes"),
                        tags$li(tags$code("features.tsv.gz"), "- Contains gene information"),
                        tags$li(tags$code("matrix.mtx.gz"), "- Contains the expression matrix")
                      ),
                      tags$li("All files must be at the root of the ZIP file, not in a subfolder")
                    ),
                    h5("Option 2: Load Processed Data", style = "font-weight: bold; color: #444;"),
                    tags$ul(
                      tags$li("If you have already processed your data with Seurat, you can load the saved RDS file directly"),
                      tags$li("This option allows you to resume analysis from a previously saved point")
                    ),
                    div(style = "border-top: 1px solid #ddd; margin-top: 20px; padding-top: 10px; color: #666; font-size: 0.9em;",
                        p("Sources:",
                          tags$ul(
                            tags$li(
                              tags$a(href="https://satijalab.org/seurat/articles/get_started_v5_new",
                                     "Seurat V5 Introduction Vignette",
                                     target="_blank")
                            ),
                            tags$li(
                              "Hao et al., Dictionary learning for integrative, multimodal and scalable single-cell analysis. ",
                              tags$i("Nature Biotechnology"),
                              " 42, 293–304 (2024). ",
                              tags$a(href="https://doi.org/10.1038/s41587-023-02100-3",
                                     "DOI: 10.1038/s41587-023-02100-3",
                                     target="_blank")
                            )
                          )
                        )
                    )
                ),
                box(title = "Dataset Parameters", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                    fluidRow(
                      column(6,
                             div(class = "well",
                                 selectInput("species_choice",  "Select species:",choices = c("Mouse" = "mouse",  "Human" = "human", "Rat" = "rat"),selected = "mouse")                             )
                      ),
                      column(6,
                             div(class = "well",
                                 radioButtons("dataset_type", "Choose Dataset Type", choices = list("snRNA-seq" = "snRNA", "Multiome" = "multiome"))
                             )
                      )
                    )
                ),
                box(title = "Load Dataset", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                    fluidRow(
                      column(6,
                             h4("Load Raw 10X Data"),
                             fileInput('file', 'Select ZIP file containing 10X data', accept = c('.zip','.h5')),
                             p(class = "text-muted", "Upload a ZIP file containing barcodes.tsv.gz, matrix.mtx.gz, and features.tsv.gz")
                      ),
                      column(6,
                             h4("Load Processed Data"),
                             fileInput("load_seurat_file", "Select Seurat Object (.rds)", accept = ".rds"),
                             p(class = "text-muted", "Upload a pre-processed Seurat object saved as RDS file")
                      )
                    )
                )
                )
              )
      ,
      ############################## Single/QC metrics and normalization  ##############################
      tabItem(
        tabName = "qc",
        box(title = "Quality Control Overview",
            status = "info",solidHeader = TRUE, collapsible = TRUE, width = 14,
            div(class = "well",
                h4("Quality Control Steps:", style = "color: #2c3e50;"),
                p("Follow these steps to ensure high-quality data:"),
                tags$ol(
                  tags$li(strong("Feature Selection:"),
                          "Set the range of unique genes detected per cell. Generally, cells with very few genes may be empty droplets or poor quality cells,
          while those with too many genes might be doublets."),
                  tags$li(strong("Mitochondrial Filtering:"),
                          "Set the maximum percentage of mitochondrial genes. High mitochondrial content often indicates cell stress or death."),
                  tags$li(strong("Visual Inspection:"),
                          "Use the violin plots and scatter plots to evaluate the distribution of these metrics across your cells."),
                  tags$li(strong("Apply Filters:"),
                          "Click 'Apply QC Filters' to remove cells that don't meet the criteria."),
                  tags$li(strong("Normalization:"),
                          "Finally, normalize the filtered data to account for technical variations between cells.")
                ),
                div(style = "margin-top: 15px;",
                    tags$em("Note: The default parameters are general guidelines. Optimal values may vary depending on your specific experimental context.")
                )
            )
        ),
        sidebarLayout(
          sidebarPanel(
            infoBoxOutput("nuclei_count"),
            div(style = "display: inline-block; width: 80%;",
                sliderInput("nFeature_range", "Unique genes detected in each cell", min = 0, max = 5000, value = c(200, 2500))
            ),
            div(style = "display: inline-block; width: 18%;",actionButton("exclamation1",label = icon("exclamation-triangle"), `data-toggle` = "popover", `data-html` = "true",`data-content` = "The number of unique genes detected in each cell. <ul><li>Low-quality cells or empty droplets will often have very few genes</li><li>Cell doublets or multiplets may exhibit an aberrantly high gene count.</li><li>Be careful, in most cases it's better to stay between 200 and 2500</li></ul>")
            ) ,
            div(style = "display: inline-block; width: 80%;",sliderInput("percent.mt_max", "Maximum value for  mitochondrial genome:", value = 5, min = 0, max = 100)
            ),
            div(style = "display: inline-block; width: 18%;",actionButton("exclamation3", label = icon("exclamation-triangle"), `data-toggle` = "popover", `data-content` = "It is normal to have mitochondrial DNA contamination for single cell, but be careful with this value for single nuclei")
            ),
            div(style = "display: inline-block; width: 80%;",sliderInput(inputId = "scale_factor", label = "Scale factor", value = 10000, min = 500, max = 50000, step = 500)
            ),
            div(style = "display: inline-block; width: 18%;",actionButton("exclamation3", label = icon("exclamation-triangle"), `data-toggle` = "popover", `data-content` = "Feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result")
            ),
            actionButton("QCmetrics", "QC metrics plot"),tags$br(),
            br(),
            actionButton("show_plots", "Feature Scatter Plots"),tags$br(),
            br(),
            actionButton("apply_qc", "Apply QC Filters"),tags$br(),
            br(),
            numericInput("num_var_features", "Number of Variable Features:", value = 2000, min = 1),
            br(),
            actionButton("normalize_data", "Normalize Data",class = "btn-primary")
          ),
          mainPanel(
            div(style = "height: 500px;", plotOutput("vlnplot")),
            fluidRow(
              column(6, div(style = "height: 400px;", plotOutput("scatter_plot1"))),
              column(6, div(style = "height: 400px;", plotOutput("scatter_plot2")))
            ),
            fluidRow(
              column(8, div(style = "height: 400px;", plotOutput("variable_feature_plot")))
            )
          )
        )),

      ############################## Single/Scaling, PCA and elbow plot ##############################
      tabItem(tabName = "dimensional_reduction",
              box(title = "Dimensional Reduction Analysis", status = "info", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  div(
                    h4("Understanding Data Scaling and PCA", style = "color: #2c3e50;"),
                    p("This step reduces the complexity of your data while preserving important biological signals:"),
                    tags$ol(
                      tags$li(strong("Data Scaling:"), "First, we normalize the expression of each gene to give equal weight to all genes, regardless of their expression level. This prevents highly-expressed genes from dominating the analysis."),
                      tags$li(strong("Principal Component Analysis (PCA):"), "PCA identifies the main sources of variation in your data, transforming thousands of gene expression measurements into a smaller set of meaningful components."),
                      tags$li(strong("Elbow Plot Analysis:"), "The elbow plot helps determine how many principal components to use in downstream analysis by showing the percentage of variance explained by each component.")
                    ),
                    div(style = "margin-top: 15px;",
                        tags$em("Note: The number of principal components chosen will affect all subsequent analyses.")
                    )
                  )
              ),
              box(title = "Scale Data and Run PCA", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  div(style = "margin-bottom: 20px;",
                      actionButton("scale_button", "Scale Data & Run PCA", class = "btn-primary", style = "width: 200px;")
                  ),
                  verbatimTextOutput("pca_results"),
                  plotOutput("loading_plot")
              ),
              box(title = "Determine Optimal Components", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  div(style = "margin-bottom: 20px;",
                      actionButton("run_elbow", "Show Elbow Plot", class = "btn-primary", style = "width: 200px;")
                  ),
                  plotOutput("elbow_plot")
              )
      ),

      ############################## Single/Neighbors calculation and clustering ##############################
      tabItem(tabName = "clustering",
              box(title = "Cell Clustering Analysis", status = "info", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  div(
                    h4("Understanding Cell Clustering", style = "color: #2c3e50;"),
                    p("This step groups cells based on their gene expression similarities to identify distinct cell populations:"),
                    tags$ol(
                      tags$li(strong("Neighbor Finding:"),
                              "First, we identify cells with similar gene expression profiles. The number of dimensions used affects how these similarities are calculated."),
                      tags$li(strong("Clustering Resolution:"),
                              "Then, we group similar cells together. The resolution parameter controls how finely the cells are grouped - higher values create more clusters."),
                      tags$li(strong("Algorithm Selection:"),
                              "Different clustering algorithms are available, each with its own approach to grouping cells.")
                    ),
                    div(style = "margin-top: 15px;",
                        tags$em("Note: Finding the right balance between too many and too few clusters is key for meaningful biological interpretation.")
                    )
                  )
              ),
              box(title = "Cluster cells", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  column(width = 3,
                         numericInput("dimension_1", "Number of dimensions:", value = 10, min = 1),
                         actionButton("run_neighbors", "Find neighbors", class = "btn-primary"),
                         checkboxInput("remove_axes", "Remove Axes", FALSE)
                  ),

                  column(width = 3,
                         numericInput("resolution_step1", "Resolution:", min = 0.1, max = 2, step = 0.1, value = 0.5),
                         actionButton("run_clustering", "Find clusters", class = "btn-primary"),
                         checkboxInput("remove_legend", "Remove Legend", FALSE)
                  ),
                  column(width = 3,
                         selectInput("algorithm_select", "Select Algorithm:", choices = list("Original Louvain" = 1, "Louvain with Multilevel Refinement" = 2, "SLM Algorithm" = 3)),
                         selectInput("umap_export_format", "File Format:",choices = c("TIFF" = "tiff", "PNG" = "png","PDF" = "pdf", "JPEG" = "jpeg", "SVG" = "svg"),selected = "tiff")

                  ),
                  column(width = 3,
                         numericInput("dpi_umap", "Resolution (DPI):", value = 300, min = 72, max = 1200, step = 72),
                         br(),   downloadButton("downloadUMAP", "Download UMAP")
                  )
              ),
              box(title = "Clustering Results", status = "primary", solidHeader = TRUE, width = 14,
                  plotOutput("clustering_plot")
              )
      ),

      ############################## Single/DoubletFinder ##############################
      tabItem(tabName = "doublet_detection",
              box(title = "Doublet Detection Overview", status = "info", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  div(
                    h4("About Doublet Detection", style = "color: #2c3e50;"),
                    p("Doublets are artifacts where two cells are captured together, leading to mixed expression profiles.
       DoubletFinder helps identify and remove these artifacts to ensure data quality."),
                    tags$ol(
                      tags$li(strong("Expected Rate:"), "Set based on your experimental protocol, typically 5-10% for 10x Genomics data."),
                      tags$li(strong("Parameter Settings:"), "pN and pK values control the algorithm's sensitivity and specificity."),
                      tags$li(strong("Visualization:"), "Results are shown on a UMAP plot where doublets are highlighted.")
                    )
                  )
              ),
              box(title = "Doublet Detection Parameters", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  column(4,
                         numericInput("doublet_rate", "Expected Doublet Rate (%):", value = 7.5, min = 0, max = 100, step = 0.5),
                         numericInput("pc_use", "Number of PCs:", value = 10, min = 1, max = 50)
                  ),
                  column(4,
                         numericInput("pN_value", "pN value:", value = 0.25, min = 0, max = 1, step = 0.05),
                         numericInput("pK_value", "pK value:", value = 0.09, min = 0, max = 1, step = 0.01)
                  ),
                  column(4,
                         div(style = "margin-top: 25px",
                             actionButton("run_doubletfinder", "Detect Doublets", class = "btn-primary", style = "margin-right: 10px;"),
                             actionButton("remove_doublets", "Remove Doublets", class = "btn-warning")
                         )
                  )
              ),
              fluidRow(
                column(6,
                       box(title = "Doublet UMAP", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                           plotOutput("doublet_umap")
                       )
                ),
                column(6,
                       box(title = "Statistics", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                           DTOutput("doublet_stats"),
                           div(style = "margin-top: 10px",
                               downloadButton("download_doublet_results", "Download Results", class = "btn-info")
                           )
                       )
                )
              )
      ),


      ############################## Single/Visualization of expressed genes ##############################
      tabItem(tabName = "plots_genes_expressions",
              box(title = "Gene Expression Visualization", status = "info", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  div(
                    h4("Multiple Ways to Visualize Gene Expression", style = "color: #2c3e50;"),
                    p("Choose from different visualization methods to understand gene expression patterns:"),
                    tags$ol(
                      tags$li(strong("Feature Plot:"), "Shows gene expression on UMAP projection. Useful for seeing spatial expression patterns."),
                      tags$li(strong("Violin Plot:"), "Displays expression distribution across clusters. Good for comparing expression levels."),
                      tags$li(strong("Dot Plot:"), "Shows both expression level and percentage of expressing cells. Perfect for comparing multiple genes."),
                      tags$li(strong("Ridge Plot:"), "Alternative to violin plots, better for visualizing many clusters.")
                    )
                  )
              ),
              box(title = "Gene Selection", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  fluidRow(
                    column(width = 4,
                           pickerInput("gene_select", "Select Genes:", choices = NULL, multiple = TRUE,
                                       options = list(`actions-box` = TRUE, `live-search` = TRUE)),
                                       numericInput("title_text_size", "Title Size:", value = 14, min = 8, max = 32)
                    ),
                    column(4,
                           selectInput("viz_assay", "Select Assay:",
                                       choices = NULL,  # On laisse vide, sera mis à jour dynamiquement
                                       selected = "RNA"),
                                       numericInput("axis_text_size", "Axis Text Size:", value = 12, min = 6, max = 30),
                           numericInput("axis_line_width", "Axis Line Width:", value = 1, min = 0.5, max = 3, step = 0.1)
                    ),
                    column(width = 4,
                           numericInput("dpi_plot", "Images resolution for download:", value = 300, min = 72, step = 72),
                           selectInput("plot_format", "Select Image Format:",
                                       choices =c("PNG" = "png", "JPEG" = "jpeg", "TIFF" = "tiff","SVG"="svg","PDF" = "pdf"), selected = "png"),
                             downloadButton("save_seurat_object_2", "Download Seurat Object")

                    )
                  )
              ),
              box(title = "Feature Plot", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  textInput("gene_list_feature", "Selected genes for FeaturePlot:", value = ""),
                  fluidRow(
                    column(width = 3,
                           actionButton("show_feature", "Display Plot", class = "btn-primary"),
                           checkboxInput("add_noaxes_feature", "Remove Axes", FALSE)
                    ),
                    column(width = 3,
                           selectInput("min_cutoff", "Minimum Cutoff:",
                                       choices = c("None" = NA, "q1" = "q1", "q10" = "q10", "q20" = "q20", "q30" = "q30",
                                                   "q40" = "q40", "q50" = "q50", "q60" = "q60", "q70" = "q70", "q80" = "q80",
                                                   "q90" = "q90", "q99" = "q99"), selected = NA),
                           checkboxInput("add_nolegend_feature", "Remove Legend", FALSE)
                    ),
                    column(width = 3,
                           selectInput("max_cutoff", "Maximum Cutoff:",
                                       choices = c("None" = NA, "q1" = "q1", "q10" = "q10", "q20" = "q20", "q30" = "q30",
                                                   "q40" = "q40", "q50" = "q50", "q60" = "q60", "q70" = "q70", "q80" = "q80",
                                                   "q90" = "q90", "q99" = "q99"), selected = NA),
                           checkboxInput("show_coexpression", "Show co-expression", value = FALSE)
                    ),
                    column(width = 3,
                           downloadButton("downloadFeaturePlot", "Download Plot")
                    )
                  ),
                  plotOutput("feature_plot")
              ),
              box(title = "Violin Plot",
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  width = 14,
                  textInput("gene_list_vln", "Selected genes:", value = ""),
                  selectInput("cluster_order_vln",
                              "Order of Clusters",
                              choices = NULL,
                              multiple = TRUE),
                  fluidRow(
                    column(width = 3, actionButton("show_vln", "Display Plot", class = "btn-primary")
                    ),
                    column(width = 3,  checkboxInput("hide_vln_points", "Hide Points", FALSE)
                    ),
                    column(width = 2,checkboxInput("add_noaxes_vln", "Remove Axes", FALSE)
                    ),
                    column(width = 2, checkboxInput("add_nolegend_vln", "Remove Legend", FALSE)
                    ),
                    column(width = 2, downloadButton("downloadVlnPlot", "Download Plot")
                    )
                  ),
                  plotOutput("vln_plot")
              ),
              box(title = "Dot Plot",
                  status = "primary",
                  solidHeader = TRUE,
                  collapsible = TRUE,
                  width = 14,
                  fluidRow(
                    column(width = 4,
                           textInput("gene_list_dotplot", "Selected genes:", value = ""),
                           actionButton("show_dot", "Display Plot", class = "btn-primary"),
                           checkboxInput("add_noaxes_dot", "Remove Axes", FALSE)
                    ),
                    column(width = 4,
                           selectInput("cluster_order_dot", "Order of Clusters", choices = NULL, multiple = TRUE),  # Changé ici
                           checkboxInput("invert_axes", "Invert Axes", value = FALSE),
                           checkboxInput("add_nolegend_dot", "Remove Legend", FALSE)
                    ),
                    column(width = 4,
                           downloadButton("downloadDotPlot", "Download Plot")
                    )
                  ),
                  plotOutput("dot_plot")
              ),
              box(title = "Ridge Plot", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  textInput("gene_list_ridge_plot", "Selected genes:", value = ""),
                  fluidRow(
                    column(width = 3,
                           actionButton("show_ridge", "Display Plot", class = "btn-primary")
                    ),
                    column(width = 3,
                           checkboxInput("add_noaxes_ridge", "Remove Axes", FALSE)
                    ),
                    column(width = 3,
                           checkboxInput("add_nolegend_ridge", "Remove Legend", FALSE)
                    ),
                    column(width = 3,
                           downloadButton("download_ridge_plot", "Download Plot")
                    )
                  ),
                  plotOutput("ridge_plot")
              ),
              box(title = "Cell Expression Analysis", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  fluidRow(
                    column(width = 4,
                           pickerInput("gene_select_genes_analysis", "Select Genes:", choices = NULL,
                                       multiple = TRUE, options = list(`actions-box` = TRUE, `live-search` = TRUE))
                    ),
                    column(width = 4,
                           numericInput("expression_threshold", "Minimum Expression Threshold:", value = 1, min = 0, max = 10),
                           actionButton("analyze_btn", "Analyze Expression", class = "btn-primary")
                    ),
                    column(width = 4,
                           downloadButton("download_genes_number_expresion", "Download Results")
                    )
                  ),
                  dataTableOutput("expression_summary")
              )
      ),

      ############################## Single/Heatmaps and scatter plots ##############################
      tabItem(tabName = "heatmaps_dualexpression",
              box(title = "Expression Visualization Tools", status = "info", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  div(
                    h4("Advanced Expression Analysis", style = "color: #2c3e50;"),
                    p("Use these tools to explore relationships between genes and clusters:"),
                    tags$ol(
                      tags$li(strong("Heatmap:"), "Visualize expression patterns across multiple genes and clusters simultaneously.
         Useful for identifying gene signatures and cluster-specific patterns."),
                      tags$li(strong("Feature Scatter:"), "Compare expression levels between two genes.
         Helps identify correlations and relationships between genes of interest.")
                    ),
                    div(style = "margin-top: 15px;",
                        tags$em("Note: You can specify clusters of interest or use automated selections for both visualizations.")
                    )
                  )
              ),
              box(title = "Heatmap Analysis", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  fluidRow(
                    column(width = 4,
                           pickerInput("gene_select_heatmap", "Select Genes:",
                                       choices = NULL, multiple = TRUE,
                                       options = list(`actions-box` = TRUE, `live-search` = TRUE)),
                           checkboxInput("use_top10_genes", "Use top 10 genes per cluster", FALSE),
                           actionButton("generateCustomHeatmap", "Generate Heatmap", class = "btn-primary")
                    ),
                    column(width = 4,
                           textInput("text_clusters", "Specify Clusters (comma-separated):", ""),
                           checkboxInput("select_all_clusters", "Select All Clusters", TRUE)
                    ),
                    column(width = 4,
                           numericInput("dpi_heatmap_single", "Resolution:", value = 300, min = 72, step = 72),
                           downloadButton("download_heatmap_single", "Download Heatmap")
                    )
                  ),
                  plotOutput("heatmap_single")
              ),
              box(title = "Gene Correlation Analysis", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  fluidRow(
                    column(width = 4,
                           pickerInput("feature1_select", "First Gene:",
                                       choices = NULL, multiple = FALSE,
                                       options = list(`actions-box` = TRUE, `live-search` = TRUE)),
                           pickerInput("feature2_select", "Second Gene:",
                                       choices = NULL, multiple = FALSE,
                                       options = list(`actions-box` = TRUE, `live-search` = TRUE))
                    ),
                    column(width = 4,
                           textInput("scatter_text_clusters", "Specify Clusters:", ""),
                           checkboxInput("scatter_select_all_clusters", "Select All Clusters", TRUE),
                           actionButton("generateFeatureScatter", "Generate Plot", class = "btn-primary")
                    ),
                    column(width = 4,
                           numericInput("dpi_scatter", "Resolution:", value = 300, min = 72, step = 72),
                           downloadButton("download_scatter_plot", "Download Plot")
                    )
                  ),
                  plotOutput("scatterPlot")
              )
      )  ,
      ############################## Single/Final UMAP ##############################
      tabItem(tabName = "assigning_cell_type_identity",
              box(title = "Cell Type Assignment", status = "info", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  div(
                    h4("Annotating Cell Clusters", style = "color: #2c3e50;"),
                    p("Assign biological identities to your clusters based on marker gene expression:"),
                    tags$ol(
                      tags$li(strong("Rename Clusters:"), "Use known cell type markers to give meaningful names to your clusters."),
                      tags$li(strong("Customize Visualization:"), "Adjust the appearance of your UMAP to highlight important features."),
                      tags$li(strong("Color Scheme:"), "Assign specific colors to clusters for consistent visualization across plots.")
                    ),
                    div(style = "margin-top: 15px;",
                        tags$em("Note: Consistent naming and coloring helps in communicating your findings effectively.")
                    )
                  )
              ),
              box(title = "Cluster Identity Assignment", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  fluidRow(
                    column(width = 4,
                           div(class = "well",
                               h4("Rename Clusters", style = "color: #2c3e50; margin-bottom: 15px;"),
                               selectInput("select_cluster", "Select cluster:", choices = NULL),
                               textInput("rename_single_cluster", "New name:"),
                               actionButton("rename_single_cluster_button", "Apply Name", class = "btn-primary")
                           )
                    ),
                    column(width = 4,
                           div(class = "well",
                               h4("Plot Settings", style = "color: #2c3e50; margin-bottom: 15px;"),
                               textInput("plot_title", "Plot title:", value = "UMAP Final"),
                               numericInput("label_font_size", "Label size:", value = 5, min = 1, max = 20, step = 0.5),
                               numericInput("pt_size", "Point size:", value = 0.3, min = 0.1, max = 3, step = 0.1)
                           )
                    ),
                    column(width = 4,
                           div(class = "well",
                               h4("Color Settings", style = "color: #2c3e50; margin-bottom: 15px;"),
                               selectInput("cluster_select", "Select cluster:", choices = NULL),
                               colourInput("cluster_colour", "Choose color:", value = "red"),
                               actionButton("update_colour", "Update Color", class = "btn-info"),
                               br(),
                               br(),
                               downloadButton("save_seurat_object_1", "Download Seurat Object")
                           )
                    )
                  )
              ),
              box(title = "Interactive UMAP", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  plotlyOutput("umap_finale")
              ),
              box(title = "Alternative Visualizations", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  fluidRow(
                    column(width = 3,
                           selectInput("plot_type_select", "Plot type:",
                                       choices = c("FeaturePlot", "VlnPlot", "DotPlot", "RidgePlot"))
                    )
                  ),
                  plotOutput("selected_plot_display")
              )
      ),

      ############################## Single/Find markers for a specific cluster ##############################
      tabItem(tabName = "cluster_comparison",
              box(title = "Cluster Biomarker Analysis", status = "info", solidHeader = TRUE, collapsible = TRUE,  width = 14,
                  div(
                    h4("Understanding Cluster Biomarkers"),
                    p("This section helps identify genes that characterize specific clusters:"),
                    tags$ol(
                      tags$li(strong("Global Comparison:"), "Compare one cluster against all others to find unique markers."),
                      tags$li(strong("Pairwise Comparison:"), "Compare two specific clusters to find distinguishing genes.")
                    )
                  )
              ),
              box(title = "Cluster Overview", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  plotOutput("umap_plot"),
                  fluidRow(
                    column(width = 6,
                           numericInput("dpi_umap_cluster", "Plot resolution:", value = 300, min = 72, step = 72)
                    ),
                    column(width = 6,
                           selectInput("umap_cluster_format", "File Format:",
                                       choices = c("TIFF" = "tiff", "PNG" = "png", "PDF" = "pdf",
                                                   "JPEG" = "jpeg", "SVG" = "svg"),
                                       selected = "tiff")
                    )
                  ),
                  fluidRow(
                    column(width = 3,
                           downloadButton("downloadUMAPCluster", "Download Plot")
                    ),
                    column(width = 3,
                           checkboxInput("show_labels", "Show Cluster Labels", value = TRUE)
                    ),
                    column(width = 3,
                           downloadButton("save_seurat_object_3", "Download Seurat Object")
                    ),
                    column(width = 3) # Colonne vide pour l'équilibre
                  )
              )
              ,
              box(title = "Global Cluster Comparison", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  fluidRow(
                    column(width = 4,
                           selectInput("cluster", "Target Cluster:", choices = NULL),
                           sliderInput("number_genes_10", "Number of Genes:", min = 10, max = 1000, value = 10, step = 50)
                    ),
                    column(width = 4,
                           numericInput("logfc_threshold_single", "Log2FC Threshold:", value = 0.1, min = 0, max = 5, step = 0.1),
                           actionButton("find_markers", "Find Markers", class = "btn-primary")
                    ),
                    column(width = 4,
                           numericInput("min_pct_single", "Min Expression %:", value = 0.01, min = 0, max = 1, step = 0.01),
                           downloadButton("download_markers_single_cluster", "Download Results")
                    )
                  ),
                  uiOutput("gene_tables_10")
              ),
              box(title = "Pairwise Cluster Comparison", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  fluidRow(
                    column(width = 4,
                           selectInput("cluster1", "Group 1 Clusters:", choices = NULL, multiple = TRUE),
                           selectInput("cluster2", "Group 2 Clusters:", choices = NULL, multiple = TRUE),
                           sliderInput("n_diff_markers", "Number of Genes:", min = 10, max = 1000, value = 10, step = 50)
                    ),
                    column(width = 4,
                           numericInput("logfc_threshold_comparison", "Log2FC Threshold:", value = 0.1, min = 0, max = 5, step = 0.1),
                           actionButton("compare_markers", "Compare Groups", class = "btn-primary")
                    ),
                    column(width = 4,
                           numericInput("min_pct_comparison", "Min Expression %:", value = 0.01, min = 0, max = 1, step = 0.01),
                           downloadButton("download_markers_multiple_clusters", "Download Results")
                    )
                  ),
                  uiOutput("gene_tables_new")
              ),
              box(
                title = "Venn Diagram Comparison", 
                status = "primary", 
                solidHeader = TRUE, 
                collapsible = TRUE, 
                width = 14,
                fluidRow(
                  column(4,
                         selectInput("venn_table_1_single", "Select first gene list:", 
                                     choices = c("None" = ""), 
                                     selected = ""),
                         checkboxInput("significant_only_venn_1_single", "Include only significant genes", value = TRUE),
                         numericInput("log_fc_threshold_venn_1_single", "Log2FC threshold:", value = 0.25, min = 0, step = 0.05)
                  ),
                  column(4,
                         selectInput("venn_table_2_single", "Select second gene list:", 
                                     choices = c("None" = ""), 
                                     selected = ""),
                         checkboxInput("significant_only_venn_2_single", "Include only significant genes", value = TRUE),
                         numericInput("log_fc_threshold_venn_2_single", "Log2FC threshold:", value = 0.25, min = 0, step = 0.05)
                  ),
                  column(4,
                         selectInput("venn_table_3_single", "Select third gene list (optional):", 
                                     choices = c("None" = ""), 
                                     selected = ""),
                         checkboxInput("significant_only_venn_3_single", "Include only significant genes", value = TRUE),
                         numericInput("log_fc_threshold_venn_3_single", "Log2FC threshold:", value = 0.25, min = 0, step = 0.05)
                  )
                ),
                fluidRow(
                  column(4,
                         selectInput("venn_direction_single", "Gene selection criteria:", 
                                     choices = c("Up-regulated (avg_log2FC > threshold)" = "up",
                                                 "Down-regulated (avg_log2FC < -threshold)" = "down",
                                                 "Both directions (|avg_log2FC| > threshold)" = "both"),
                                     selected = "up")
                  ),
                  column(4,
                         numericInput("p_val_threshold_venn_single", "p-value threshold:", 
                                      value = 0.05, min = 0, max = 1, step = 0.01)
                  ),
                  column(4,
                         checkboxInput("use_adjusted_p_venn_single", "Use adjusted p-values", value = TRUE),
                         actionButton("generate_venn_btn_single", "Generate Venn Diagram", 
                                      class = "btn-primary btn-lg", 
                                      style = "margin-top: 23px;")
                  )
                ),
                fluidRow(
                  column(4,
                         colourInput("venn_color_1_single", "Color for first set:", value = "#56B4E9")
                  ),
                  column(4,
                         colourInput("venn_color_2_single", "Color for second set:", value = "#E69F00")
                  ),
                  column(4,
                         colourInput("venn_color_3_single", "Color for third set (if used):", value = "#009E73")
                  )
                ),
                plotOutput("venn_plot_single", height = "400px"),
                fluidRow(
                  column(6,
                         selectInput("venn_diagram_format_single", "Export format:", 
                                     choices = c("PNG" = "png", "TIFF" = "tiff", "PDF" = "pdf", "JPEG" = "jpeg"),
                                     selected = "png"),
                         numericInput("venn_diagram_dpi_single", "DPI:", value = 300, min = 72, step = 72)
                  ),
                  column(6,
                         downloadButton("download_venn_diagram_single", "Download Venn Diagram"),
                         downloadButton("download_venn_gene_lists_single", "Download Gene Lists")
                  )
                ),
                fluidRow(
                  column(12,
                         selectInput("selected_gene_set_single", "View gene list:", choices = NULL),
                         DTOutput("venn_gene_table_single")
                  )
                )
              )
              ,
              box(title = "Cluster Composition Analysis", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  fluidRow(
                    column(6,
                           actionButton("generate_cluster_composition_single", "Generate Cluster Composition", class = "btn-primary")
                    ),
                    column(6,
                           downloadButton("download_cluster_composition_single", "Download Table", class = "btn-success")
                    )
                  ),
                  br(),
                  DTOutput("cluster_composition_single")
              ),
              box(
                title = "Gene Co-expression Analysis", 
                status = "primary", 
                solidHeader = TRUE, 
                collapsible = TRUE, 
                width = 14,
                
                fluidRow(
                  column(
                    width = 4,
                    textAreaInput("gene_text_coexpression_single", "Enter genes (comma-separated):", value = "", placeholder = "e.g., Pax7, Myod1", height = "80px", width = "100%"),
                    helpText("Enter exactly 2 genes for co-expression analysis")
                  ),
                  
                  column(
                    width = 4,
                    numericInput("gene1_threshold_single", "Gene 1 expression threshold:", value = 0, min = 0, max = 10, step = 0.1),
                    numericInput("gene2_threshold_single", "Gene 2 expression threshold:", value = 0, min = 0, max = 10, step = 0.1),
                    textInput("coexpr_text_clusters", "Specific Clusters (optional):", ""),
                    checkboxInput("coexpr_select_all_clusters", "Analyze All Clusters", TRUE)
                  ),
                  
                  column(
                    width = 4,
                    selectInput("coexpr_plot_type", "Visualization:", choices = c("Stacked Bar" = "bar", "Heatmap" = "heat"), selected = "bar"),
                    br(),
                    actionButton("analyze_coexpression_single", "Analyze Co-expression", class = "btn-primary", style = "width: 100%; margin-bottom: 10px;"),
                    fluidRow(
                      column(6, downloadButton("download_coexpression_table_single", "Download Table", style = "width: 100%;")),
                      column(6, downloadButton("download_coexpression_plot_single", "Download Plot", style = "width: 100%;"))
                    )
                  )
                ),
                
                tabsetPanel(
                  tabPanel("Results Table", DTOutput("gene_coexpression_table_single")),
                  tabPanel("Visualization", plotOutput("gene_coexpression_plot_single", height = "500px"))
                )
              )),


      ############################## Single/Subsetting ##############################
      tabItem(tabName = "subset",
              box(title = "Data Subsetting Tools", status = "info", solidHeader = TRUE, collapsible = TRUE,  width = 14,
                  div(
                    h4("Create Focused Datasets", style = "color: #2c3e50;"),
                    p("Select specific cell populations either by cluster identity or gene expression patterns. All previous analyses are preserved in the subset."),
                    tags$ol(
                      tags$li(strong("Cluster-based:"), "Select specific clusters to retain"),
                      tags$li(strong("Expression-based:"), "Filter cells based on expression of key genes"),
                      tags$li(strong("Results:"), "Compare original and subset UMAPs")
                    )
                  )
              ),

              box(title = "Original Dataset", status = "primary", solidHeader = TRUE, width = 14,
                  plotOutput("global_umap")
              ),
              fluidRow(
                column(width = 6,
                       box(title = "Subset by Clusters", status = "primary", solidHeader = TRUE, width = NULL,
                           selectInput("select_ident_subset", "Select clusters:", choices = NULL, multiple = TRUE),
                           actionButton("apply_subset", "Create Subset", class = "btn-primary")
                       )
                ),
                column(width = 6,
                       box(title = "Subset by Expression", status = "primary", solidHeader = TRUE, width = NULL,
                           numericInput("expression_threshold", "Expression threshold:", value = 0.1),
                           textInput("gene_list_subset", "Genes (comma-separated):", value = ""),
                           numericInput("num_genes_to_express", "Minimum expressed genes:", value = 1, min = 1),
                           actionButton("apply_gene_subset", "Create Subset", class = "btn-primary")
                       )
                )
              ),
              box(title = "Subset Preview", status = "primary", solidHeader = TRUE, width = 14,
                  plotOutput("subset_umap"),
                  div(style = "margin-top: 15px; text-align: center;",
                      downloadButton("download_subset_seurat", "Save as RDS", class = "btn-success")
                  )
              )
      ),
      ############################## Multiple/Loading Data ##############################
      tabItem(
        tabName = "load_datasets_merge",
        fluidRow(
          box(
            title = "Dataset Loading Options",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            fluidRow(
              # Option gauche: Charger plusieurs datasets
              column(6,
                     div(
                       class = "well",
                       style = "height: 350px; padding: 15px; border-right: 2px solid #ddd;",
                       h4("Load Multiple Raw Datasets", style = "color: #666; text-align: center;"),
                       p("Use this option to load and integrate multiple raw or processed datasets."),
                       p("Requirements:"),
                       tags$ul(
                         tags$li("Raw data: Compress 10X files into separate ZIP files"),
                         tags$li("Processed data: Load multiple Seurat objects (.rds files)"),
                         tags$li("You can mix both types in the same integration")
                       ),
                       div(style = "text-align: center; margin-top: 20px;",
                           actionButton("open_file_input_modal", "Load Multiple Datasets",
                                        class = "btn-primary btn-lg",
                                        icon = icon("upload"))
                       )
                     )
              ),

              # Option droite: Charger un objet pré-intégré
              column(6,
                     div(
                       class = "well",
                       style = "height: 350px; padding: 15px;",
                       h4("Load Pre-integrated Object", style = "color: #666; text-align: center;"),
                       p("Use this option if you already have a Seurat object with integrated datasets."),
                       p("Requirements:"),
                       tags$ul(
                         tags$li("A single .rds file containing an integrated Seurat object"),
                         tags$li("The object should contain multiple datasets already integrated")
                       ),
                       div(style = "text-align: center; margin-top: 20px;",
                           fileInput("load_seurat_file_merge",
                                     label = NULL,
                                     accept = ".rds",
                                     buttonLabel = "Browse for Seurat Object...",
                                     placeholder = "No file selected")
                       )
                     )
              )
            )
          )
        ),
        fluidRow(
          conditionalPanel(
            condition = "output.datasets_loaded",
            box(
              title = "Integration and Metadata Management",
              status = "primary",
              solidHeader = TRUE,
              collapsible = TRUE,
              width = 12,
              fluidRow(
                column(6,
                       h4("Add Metadata to Datasets"),
                       actionButton('add_field', 'Add New Metadata Field',
                                    icon = icon("plus"),
                                    class = "btn-info")
                ),
                column(6,
                       uiOutput("metadata_inputs"),
                       br(),
                       actionButton('add_metadata', 'Save Metadata',
                                    icon = icon("save"),
                                    class = "btn-success")
                )
              )
            )
          )
        )
      ),
      ############################## Multiple/Scaling and PCA reduction ##############################
      tabItem(
        tabName = "clustering_merge",
        fluidRow(
          box(title = "Scaling & PCA", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
              fluidRow(
                column(6,
                       actionButton("runScalePCA", "Run Scaling, PCA and Elbow Plot",class = "btn-primary")
                ),
                column(6,
                       actionButton("runHarmony", "Run Harmony Integration"),
                       actionButton("exclamation_harmony", label = icon("exclamation-triangle"),
                                    `data-toggle` = "popover",
                                    `data-html` = "true",
                                  `data-content` = "Harmony is an efficient algorithm for integrating multiple datasets. It corrects batch effects while preserving biological variation.")
                )
              ),
              plotOutput("elbow_plot2")
          ),
          box(title = "Harmony Settings", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
              fluidRow(
                column(4,
                       selectInput("harmony_vars", "Variables to integrate:", choices = NULL, multiple = TRUE),
                       numericInput("harmony_dims", "Number of dimensions:", value = 15, min = 1)
                )
              )
          ),
          box(title = "Cluster cells", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12, height = 650,
              column(width = 3,
                     numericInput("dimension_2", "Number of dimensions:", value = 15, min = 1),
                     actionButton("runFindNeighbors", "Find Neighbors and run UMAP",class = "btn-primary"),
                     checkboxInput("remove_axes_umap_merge", "Remove Axes", FALSE)
              ),
              column(width = 3,
                     numericInput("resolution_step2", "Resolution for clustering:", min = 0.01, step = 0.1, value = 0.5),
                     actionButton("runFindClusters", "Find clusters",class = "btn-primary"),
                     checkboxInput("remove_legend_umap_merge", "Remove Legend", FALSE)
              ),
              column(width = 3,
                     selectInput("algorithm_select", "Select Algorithm:",  choices = list("Original Louvain" = 1,"Louvain with Multilevel Refinement" = 2,"SLM Algorithm" = 3)),

                    selectInput("umap_merge_format", "File Format:", choices = c("TIFF" = "tiff",  "PNG" = "png",  "PDF" = "pdf",  "JPEG" = "jpeg",   "SVG" = "svg"), selected = "tiff"),
              ),
              column(width = 3,
                     numericInput("dpi_umap_merge", "Resolution (DPI):", value = 300, min = 72, max = 1200, step = 72),
                     br(),
                     downloadButton("downloadUMAP_merge", "Download UMAP")

              ),
              plotOutput("UMAPPlot_cluster_merge")
          )
        )),
      ############################## Multiple/Visualize genes expressions ##############################
      tabItem(tabName = "plots_genes_expressions_merge",
              fluidRow(
                box(title = "Select genes for plots", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
                    fluidRow(
                      column(4, selectInput("group_by_select", "Group By:", choices = c("dataset" = "dataset", "Seurat Clusters" = "seurat_clusters", "Cell Types" = "cell_type"))),
                      column(4, selectInput("viz_assay_merge", "Select Assay:", choices = c("RNA", "integrated"), selected = "RNA")),
                      column(4, numericInput("dpi_input_merge", "Images resolution for download:", value = 300, min = 72, step = 72))
                    ),
                    fluidRow(
                      column(4,
                             conditionalPanel(condition = "input.group_by_select == 'dataset'",
                                              fluidRow(
                                                column(12, selectInput("metadata_to_compare", "Optional: Split by:", choices = c("None" = ""), selected = "")),
                                                column(12, selectInput("comparison_mode", "Comparison Mode:", choices = c("Split by clusters" = "split", "Compare selected clusters" = "subset")))
                                              )
                             )
                      ),
                      column(4,
                             selectInput("plot_format_merge", "Select Image Format:",
                                         choices = c("PNG" = "png", "JPEG" = "jpeg", "TIFF" = "tiff","SVG"="svg","PDF" = "pdf"),  selected = "png"
                             )
                      ),
                      column(4,
                             pickerInput("geneInput_merge", "Select Genes:", choices = NULL, multiple = TRUE, options = list(`actions-box` = TRUE, `live-search` = TRUE)),
                             br(),
                             downloadButton("save_seurat_merge_2", "Download Seurat Object")
                      )
                    )
                )
              ),
              fluidRow(
                box(title = "Visualize with a Feature plot", status = "primary", collapsed = TRUE, solidHeader = TRUE, collapsible = TRUE, width = 12,
                    textInput("gene_list_feature_merge", "Selected genes for FeaturePlot:", value = ""),
                    fluidRow(
                      column(6,
                             fluidRow(
                               column(6, actionButton("runFeaturePlot", "Run Feature Plot",class = "btn-primary")),
                               column(6, checkboxInput("add_nolegend_feature_merge", "Remove Legend", value = FALSE))
                             ),
                             fluidRow(
                               column(6, checkboxInput("add_noaxes_feature_merge", "Remove Axes", value = FALSE)),
                               column(6, checkboxInput("show_coexpression_merge", "Show Co-expression", value = FALSE))
                             ),
                             fluidRow(
                               column(12, downloadButton("downloadFeaturePlotMerge", "Download FeaturePlot"))
                             )
                      ),
                      column(6,
                             fluidRow(
                               column(12, selectInput("min_cutoff_feature_merge", "Minimum Cutoff:", choices = c("None" = NA, "q1" = "q1", "q10" = "q10", "q20" = "q20", "q30" = "q30", "q40" = "q40", "q50" = "q50", "q60" = "q60", "q70" = "q70", "q80" = "q80", "q90" = "q90", "q99" = "q99"), selected = NA)),
                               column(12, selectInput("max_cutoff_feature_merge", "Maximum Cutoff:", choices = c("None" = NA, "q1" = "q1", "q10" = "q10", "q20" = "q20", "q30" = "q30", "q40" = "q40", "q50" = "q50", "q60" = "q60", "q70" = "q70", "q80" = "q80", "q90" = "q90", "q99" = "q99"), selected = NA))
                             )
                      )
                    ),
                    div(style = "overflow: auto;", plotOutput("FeaturePlot2"))
                )
              ),
              fluidRow(
                box(title = "Visualize with a Violin plot", status = "primary", collapsed = TRUE, solidHeader = TRUE, collapsible = TRUE, width = 12,
                    textInput("gene_list_vln_merge", "Selected genes for VlnPlot:", value = ""),
                    fluidRow(
                      column(4,
                             actionButton("runVlnPlot", "Run Vln Plot",class = "btn-primary"),
                             checkboxInput("add_nolegend_vln_merge", "Remove Legend", value = FALSE),
                             checkboxInput("add_noaxes_vln_merge", "Remove Axes", value = FALSE)
                      ),
                      column(4,
                             selectInput("cluster_order_vln_merge", "Select clusters to show:", choices = NULL,multiple = TRUE
                             ),
                             checkboxInput("hide_vln_points_merge", "Hide points", value = FALSE)
                      ),
                      column(4,
                             downloadButton("downloadVlnPlotMerge", "Download VlnPlot")
                      )
                    ),
                    div(style = "overflow: auto;", plotOutput("VlnPlot2"))
                )
              ),
              fluidRow(
                box(title = "Visualize with a Dot plot", status = "primary", collapsed = TRUE, solidHeader = TRUE, collapsible = TRUE, width = 12,
                    textInput("gene_list_dot_merge", "Selected genes for DotPlot:", value = ""),
                    fluidRow(
                      column(4,
                             actionButton("runDotPlot", "Generate DotPlot",class = "btn-primary"),
                             checkboxInput("add_nolegend_dot_merge", "Remove Legend", value = FALSE),
                             checkboxInput("add_noaxes_dot_merge", "Remove Axes", value = FALSE)
                      ),
                      column(4,
                             selectInput("cluster_order_dotplot_merge", "Select clusters to show:", choices = NULL, multiple = TRUE),
                             checkboxInput("invert_axes", "Invert Axes", value = FALSE)
                      ),
                      column(4, downloadButton("downloadDotPlotMerge", "Download DotPlot"))
                    ),
                    div(style = "overflow: auto;", plotOutput("DotPlot2"))
                )
              ),
              fluidRow(
                box(title = "Visualize with a Ridge plot", status = "primary", collapsed = TRUE, solidHeader = TRUE, collapsible = TRUE, width = 12,
                    textInput("gene_list_ridge_merge", "Selected genes for Ridge Plot:", value = ""),
                    fluidRow(
                      column(4,
                             actionButton("runRidgePlot", "Run Ridge Plot",class = "btn-primary"),
                             checkboxInput("add_nolegend_ridge_merge", "Remove Legend", value = FALSE),
                             checkboxInput("add_noaxes_ridge_merge", "Remove Axes", value = FALSE)
                      ),
                      column(4, downloadButton("downloadRidgePlotMerge", "Download Ridge Plot"))
                    ),
                    div(style = "overflow: auto;", plotOutput("Ridge_plot_merge"))
                )
              ),
              fluidRow(
                box(title = "Visualize gene expression", status = "primary", collapsed = TRUE, solidHeader = TRUE, collapsible = TRUE, width = 12,
                    fluidRow(
                      column(4, textInput("gene_list_genes_expression_merge", "Selected genes:", value = "")),
                      column(4, numericInput("logfc_threshold_genes_expression_merge", "Expression Threshold:", value = 1, min = 0, max = 10)),
                      column(4,
                             actionButton("analyze_btn_genes_expression_merge", "Analyze Expression",class = "btn-primary"),
                             downloadButton("download_genes_number_expression_merge", "Download Table")
                      )
                    ),
                    div(style = "overflow: auto;", dataTableOutput("expression_summary_merge"))
                )
              )
      ),

      ############################## Multiple/Heatmap and dual expression multi datasets ##############################
      tabItem(
        tabName = "heatmaps_dualexpression_merge",
        h2("Heatmaps for Multiple Datasets"),
        fluidRow(
          box(
            title = "Heatmap Settings",
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            width = 12,
            height = 650,
            fluidRow(
              column(
                width = 12,
                h4("Gene Selection", style = "margin-top:0; color:#666;"),
                fluidRow(
                  column(
                    width = 8,
                    pickerInput("gene_select_heatmap_multi", "Select Genes:",
                                choices = NULL,
                                multiple = TRUE,
                                options = list(`actions-box` = TRUE, `live-search` = TRUE))
                  ),
                  column(
                    width = 4,
                    checkboxInput("use_top10_genes_merge", "Use top 10 genes per cluster", FALSE)
                  )
                ),
                hr(style = "margin-top:10px; margin-bottom:15px;")
              )
            ),
            fluidRow(
              column(
                width = 6,
                h4("Grouping & Clusters", style = "margin-top:0; color:#666;"),
                selectInput("dataset_select_heatmap", "Group By:",
                            choices = c("dataset", "cluster")),
                selectInput("assay_select_heatmap", "Data Assay:",
                            choices = c("RNA", "integrated"),
                            selected = "RNA"),
                checkboxInput("select_all_clusters_merge", "Select All Clusters", TRUE),
                conditionalPanel(
                  condition = "!input.select_all_clusters_merge",
                  textInput("text_clusters_merge", "Specify Clusters (comma-separated):", "")
                )
              ),
              column(
                width = 6,
                h4("Display & Export", style = "margin-top:0; color:#666;"),
                numericInput("dpi_heatmap_multi", "Resolution (DPI):", value = 300, min = 72, step = 72),
                selectInput("heatmap_format", "Export Format:",
                            choices = c("PNG" = "png", "PDF" = "pdf", "TIFF" = "tiff", "JPEG" = "jpeg"),
                            selected = "png"),
                div(
                  style = "margin-top: 20px;",
                  actionButton("generateHeatmapMulti", "Generate Heatmap", class = "btn-primary"),
                  downloadButton("download_heatmap_multi", "Download Heatmap")
                )
              )
            ),
            plotOutput("heatmap_plot_multi")
          )
        )
        ,
        fluidRow(
          box(title = "Dual expression", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,height = 650,
              column(width = 3,
                     pickerInput("gene_select_scatter_multi1", "Select Genes:", choices = NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE)),
                     pickerInput("gene_select_scatter_multi2", "Select Genes:", choices = NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE)),
              ),
              column(width = 3,
                     textInput("text_clusters_merge_scatter", "Specify Clusters (comma-separated):", ""),
                     checkboxInput("select_all_clusters_merge_scatter", "Select All Clusters", TRUE)),


              column(width = 3,
                     selectInput("dataset_select_scatter", "Group By:", choices = c("dataset", "cluster")),
                     numericInput("dpi_scatter_multi", "Resolution for Download:", value = 300, min = 72, step = 72),
              ),
              column(width = 3,
                     br(),
                     downloadButton("download_scatter_multi", "Download Heatmap"),
                     br(),
                     br(),
                     br(),
                     actionButton("generateScatterMulti", "Generate Dual expression",class = "btn-primary"),
              ),
              plotOutput("scatter_plot_multi")
          ))) ,
      ############################## Multiple/Assigning cell identity merge ##############################
      tabItem(tabName = "assigning_cell_type_identity_merge",
              # Introduction Box
              box(title = "Cell Type Assignment for Merged Data", status = "info", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  div(
                    h4("Annotating Cell Clusters in Merged Datasets", style = "color: #2c3e50;"),
                    p("Assign biological identities to your clusters across integrated datasets:"),
                    tags$ol(
                      tags$li(strong("Rename Clusters:"), "Provide consistent naming across all datasets."),
                      tags$li(strong("Visualization Settings:"), "Customize UMAP appearance for merged data."),
                      tags$li(strong("Color Coordination:"), "Maintain consistent color schemes across datasets.")
                    ),
                    div(style = "margin-top: 15px;",
                        tags$em("Note: Consistent annotation across datasets is crucial for comparative analysis.")
                    )
                  )
              ),
              box(title = "Merged Cluster Identity Assignment", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  fluidRow(
                    column(width = 4,
                           div(class = "well",
                               h4("Rename Clusters", style = "color: #2c3e50; margin-bottom: 15px;"),
                               selectInput("select_cluster_merge", "Select cluster:", choices = NULL),
                               textInput("rename_single_cluster_merge", "New name:"),
                               actionButton("rename_single_cluster_merge_button", "Apply Name", class = "btn-primary")
                           )
                    ),
                    column(width = 4,
                           div(class = "well",
                               h4("Plot Settings", style = "color: #2c3e50; margin-bottom: 15px;"),
                               textInput("plot_title_merge", "Plot title:", value = "UMAP Final Merge"),
                               numericInput("label_font_size_merge", "Label size:", value = 5, min = 1, max = 20, step = 0.5),
                               numericInput("pt_size_merge", "Point size:", value = 1, min = 0.01, max = 5, step = 0.1),
                               checkboxInput("show_cluster_labels", "Show cluster labels", value = TRUE)
                           )
                    ),
                    column(width = 4,
                           div(class = "well",
                               h4("Color Settings", style = "color: #2c3e50; margin-bottom: 15px;"),
                               selectInput("select_color_merge", "Select cluster:", choices = NULL),
                               colourInput("select_cluster_merge_color", "Choose color:", value = "red"),
                               actionButton("update_colour_merge_button", "Update Color", class = "btn-info"),
                               br(),
                               br(),
                               downloadButton("save_seurat_merge_3", "Save Seurat Object"),

                           )
                    )
                  )
              ),
              box(title = "Interactive UMAP - Merged Data", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  plotlyOutput("umap_finale_merge")
              ),
              box(title = "Alternative Visualizations", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  fluidRow(
                    column(width = 3,
                           selectInput("plot_type_select_merge", "Plot type:",
                                       choices = c("FeaturePlot", "VlnPlot", "DotPlot", "RidgePlot"))
                    )
                  ),
                  plotOutput("selected_plot_display_merge")
              )
      ),

      ############################## Multiple/Calculation of differentially expressed genes ##############################
      tabItem(
        tabName = "cluster_comparison_merge",
        box( title = "UMAP Visualization",   status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
     fluidRow(
       plotOutput("filtered_umap_plot"),

            column(width = 4,
                   uiOutput("dataset_filter_ui"),
                   checkboxInput("bold_labels_merge", "Bold labels", value = FALSE),
                   checkboxInput("show_labels_merge", "Show labels", value = TRUE),
            ),
            column(width = 4,
           numericInput("filtered_umap_plot_dpi", "Resolution (DPI):", value = 300, min = 72, step = 72),
          downloadButton("save_seurat_merge_4", "Save Seurat Object")
            ),
            column(width = 4,
                   selectInput("filtered_umap_format", "File Format:",
                               choices = c("TIFF" = "tiff", "PNG" = "png", "PDF" = "pdf", "JPEG" = "jpeg", "SVG" = "svg"),
                               selected = "tiff"),
                   downloadButton("download_filtered_umap_plot", "Download UMAP plot")
            )
          )
        )
        ,
        box(  title = "Compares one cluster with all others", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
              column(width = 4,
                     selectInput("selected_cluster", label = "Select a cluster for comparison", choices = NULL),
                    actionButton("calculate_DE", "Start analysis",class = "btn-primary")
              ),
              column(width = 4,
                     numericInput("logfc_threshold_merge", label = "Log2 Fold Change threshold:", value = 0.25),
                     downloadButton('download_markers_single_cluster_merge', 'Download differentially expressed genes')
              ),
              column(width = 4,
                     numericInput("min_pct_merge", "Percentage threshold:", value = 0.01, min = 0, max = 1, step = 0.01)
              ),
              DTOutput('DE_genes_table')
        ),
        box(  title = "Compares one cluster with one other cluster", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
              column(width = 4,
                     uiOutput("cluster1_compare_ui"),
                     numericInput("min_pct_compare_merge", "Percentage threshold:", value = 0.01, min = 0, max = 1, step = 0.01),
              ),
              column(width = 4,
                     uiOutput("cluster2_compare_ui"),

                     actionButton("compare_clusters_button", "Start analysis between those Clusters",class = "btn-primary")
              ),
              column(width = 4,
                     numericInput("logfc_threshold_compare_merge", label = "Log2 Fold Change threshold:", value = 0.25),

                     downloadButton('download_markers_multiple_clusters_merge', 'Download differentially expressed genes')
              ),
              DTOutput("diff_genes_table_compare")
        ),
        box(title = "Compares a cluster between two datasets ", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
            column(width = 4,
                   uiOutput("dataset1_compare_ui"),
                   numericInput("logfc_threshold_datasets", label = "Log2 Fold Change threshold:", value = 0.25),
                   actionButton("compare_datasets_button", "start analysis  between those datasets",class = "btn-primary")
            ),
            column(width = 4,
                   uiOutput("dataset2_compare_ui"),
                   numericInput("min_pct_compare_dataset_merge", "Percentage threshold:", value = 0.01, min = 0, max = 1, step = 0.01),
                   checkboxInput("all_clusters", "Compare all clusters", value = FALSE)
            ),
            column(width = 4,
                   uiOutput("cluster_compare_ui"),
                   downloadButton("download_diff_dataset_cluster", "Download differentially expressed genes")
            ),
            DTOutput("diff_dataset_cluster")
        ),
     box(title = "Cluster composition", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
         fluidRow(
           column(6,
                  actionButton("generate_cluster_table_multiple", "Generate Cluster Table", class = "btn-primary")
           ),
           column(6,
                  downloadButton("download_cluster_composition_multiple", "Download Table", class = "btn-success")
           )
         ),
         br(),
         DTOutput("cluster_table_multiple")
     ),
     
     box( title = "Venn Diagram Comparison",  status = "primary",   solidHeader = TRUE,  collapsible = TRUE,   width = 14,
       fluidRow(
         column(4,
                selectInput("venn_table_1", "Select first gene list:", 
                            choices = c("None" = ""), 
                            selected = ""),
                checkboxInput("significant_only_venn_1", "Include only significant genes", value = TRUE),
                numericInput("log_fc_threshold_venn_1", "Log2FC threshold:", value = 0.25, min = 0, step = 0.05)
         ),
         column(4,
                selectInput("venn_table_2", "Select second gene list:", 
                            choices = c("None" = ""), 
                            selected = ""),
                checkboxInput("significant_only_venn_2", "Include only significant genes", value = TRUE),
                numericInput("log_fc_threshold_venn_2", "Log2FC threshold:", value = 0.25, min = 0, step = 0.05)
         ),
         column(4,
                selectInput("venn_table_3", "Select third gene list (optional):", 
                            choices = c("None" = ""), 
                            selected = ""),
                checkboxInput("significant_only_venn_3", "Include only significant genes", value = TRUE),
                numericInput("log_fc_threshold_venn_3", "Log2FC threshold:", value = 0.25, min = 0, step = 0.05)
         )
       ),
       fluidRow(
         column(4,
                selectInput("venn_direction", "Gene selection criteria:", 
                            choices = c("Up-regulated (avg_log2FC > threshold)" = "up",
                                        "Down-regulated (avg_log2FC < -threshold)" = "down",
                                        "Both directions (|avg_log2FC| > threshold)" = "both"),
                            selected = "up")
         ),
         column(4,
                numericInput("p_val_threshold_venn", "p-value threshold:", 
                             value = 0.05, min = 0, max = 1, step = 0.01)
         ),
         column(4,
                checkboxInput("use_adjusted_p_venn", "Use adjusted p-values", value = TRUE),
                actionButton("generate_venn_btn", "Generate Venn Diagram", 
                             class = "btn-primary btn-lg", 
                             style = "margin-top: 23px;")
         )
       ),
       fluidRow(
         column(4,
                colourInput("venn_color_1", "Color for first set:", value = "#56B4E9")
         ),
         column(4,
                colourInput("venn_color_2", "Color for second set:", value = "#E69F00")
         ),
         column(4,
                colourInput("venn_color_3", "Color for third set (if used):", value = "#009E73")
         )
       ),
       plotOutput("venn_plot", height = "400px"),
       fluidRow(
         column(6,
                selectInput("venn_diagram_format", "Export format:", 
                            choices = c("PNG" = "png", "TIFF" = "tiff", "PDF" = "pdf", "JPEG" = "jpeg"),
                            selected = "png"),
                numericInput("venn_diagram_dpi", "DPI:", value = 300, min = 72, step = 72)
         ),
         column(6,
                downloadButton("download_venn_diagram", "Download Venn Diagram"),
                downloadButton("download_venn_gene_lists", "Download Gene Lists")
         )
       ),
       fluidRow(
         column(12,
                selectInput("selected_gene_set", "View gene list:", choices = NULL),
                DTOutput("venn_gene_table")
         )
       )
     ),
     box(
       title = "Gene Co-expression Analysis", 
       status = "primary", 
       solidHeader = TRUE, 
       collapsible = TRUE, 
       width = 14,
       
       fluidRow(
         column(
           width = 4,
           textAreaInput(
             "gene_text_coexpression_multiple",
             "Enter genes (comma-separated):",
             value = "",
             placeholder = "e.g., Pax7, Myod1, Myog",
             height = "80px",
             width = "100%"
           ),
           numericInput(
             "coexpr_threshold_multiple",
             "Expression threshold:",
             value = 0,
             min = 0,
             max = 10,
             step = 0.1
           )
         ),
         
         column(
           width = 4,
           selectInput(
             "coexpr_group_by_multiple",
             "Primary grouping:",
             choices = c(
               "By Dataset" = "dataset",
               "By Cluster" = "cluster"
             ),
             selected = "dataset"
           ),
           
           conditionalPanel(
             condition = "input.coexpr_group_by_multiple == 'cluster'",
             selectInput(
               "coexpr_secondary_split_multiple",
               "Show dataset breakdown:",
               choices = c(
                 "No" = "none",
                 "Yes" = "dataset"
               ),
               selected = "none"
             )
           ),
           
           conditionalPanel(
             condition = "input.coexpr_group_by_multiple == 'dataset'",
             selectInput(
               "coexpr_secondary_split_multiple_alt",
               "Show cluster breakdown:",
               choices = c(
                 "No" = "none",
                 "Yes" = "cluster"
               ),
               selected = "none"
             )
           )
         ),
         
         column(
           width = 4,
           br(),
           actionButton(
             "analyze_coexpression_multiple",
             "Analyze Co-expression",
             class = "btn-primary",
             style = "width: 100%; margin-bottom: 10px;"
           ),
           fluidRow(
             column(6,
                    downloadButton(
                      "download_coexpression_table_multiple",
                      "Download Table",
                      style = "width: 100%;"
                    )
             ),
             column(6,
                    downloadButton(
                      "download_coexpression_plot_multiple",
                      "Download Plot",
                      style = "width: 100%;"
                    )
             )
           )
         )
       ),
       
       hr(),
       
       tabsetPanel(
         tabPanel(
           "Results Table",
           DTOutput("gene_coexpression_table_multiple")
         ),
         tabPanel(
           "Visualization",
           plotOutput("gene_coexpression_plot_multiple", height = "700px")
         ),
         tabPanel(
           "Summary",
           verbatimTextOutput("coexpression_summary_multiple")
         )
       )
     )
     ),
      ############################## Multiple subset##############################
      
     tabItem(tabName = "subset_merge",
             box(width = 14, plotOutput("global_umap_merge")),
             fluidRow(
               column(width = 4,
                      box(title = "Subset by clusters", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = NULL,
                          selectInput("select_ident_subset_merge", "Select the clusters to include:", choices = NULL, multiple = TRUE),
                          actionButton("apply_subset_merge", "Apply cluster based subset", class = "btn-primary")
                      )
               ),
               column(width = 4,
                      box(title = "Subset by genes expressions", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = NULL,
                          numericInput("expression_threshold_merge", "Expression threshold", value = 1, min = 0),
                          textInput("gene_list_merge", "Enter genes separated by comma (,)", value = ""),
                          numericInput("num_genes_to_express_merge", "Number of genes to be expressed from the list:", value = 1, min = 1),
                          actionButton("apply_gene_subset_merge", "Apply Gene Subsetting", class = "btn-primary")
                      )
               ),
               column(width = 4,
                      box(title = "Subset by metadata", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = NULL,
                          selectInput("metadata_column_subset", "Select metadata column:", choices = NULL),
                          uiOutput("metadata_values_ui"),
                          actionButton("apply_metadata_subset", "Apply Metadata Subset", class = "btn-primary")
                      )
               )
             ),
             box(width = 14, plotOutput("subset_umap_merge")),
             box(title = "Save Subset", status = "success", solidHeader = TRUE, width = 14, align = "center",
                 downloadButton("download_subset_merge", "Save subset as .RDS", class = "btn-lg")
             )
     ),
      ############################## Cell Chat load data ##############################
      tabItem(
        tabName = "load_data_cellchat",
        fluidRow(
          box(title = "About CellChat Analysis", status = "info", solidHeader = TRUE, collapsible = TRUE, width = 12,
              p("CellChat is a tool for analyzing cell-cell communication through ligand-receptor interactions in single-cell transcriptomics data."),
              h4("Key Features:"),
              tags$ul(
                tags$li(strong("Data Input:"), "Uses Seurat objects for cell-cell communication inference"),
                tags$li(strong("Database:"), "Curated ligand-receptor interactions from multiple sources"),
                tags$li(strong("Analysis:"), "Identifies significant interactions and communication patterns"),
                tags$li(strong("Visualization:"), "Multiple ways to visualize cell-cell communication networks")
              ),
              tags$ul(
                tags$li(icon("lightbulb"), " Tip: Ensure your Seurat object has clear cell type annotations"),
                tags$li(icon("book"), " Publication: ", tags$a(href="https://www.nature.com/articles/s41467-021-21246-9", "Jin et al., Nature Communications, 2021", target="_blank"))
              ),
              tags$small(
                "Documentation: ",
                tags$a(href="https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html",
                       "CellChat Tutorial",
                       target="_blank")
              )
          )
        ),
        fluidRow(
          box(title = "Load Data", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 6,
              fileInput("seurat_file_cellchat", "Upload Seurat RDS file")
          ),
          box(title = "GaspouDB Settings", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 6,
              selectInput("species_cellchat", "Species", choices = c("Mouse" = "mouse", "Human" = "human"), selected = "mouse"),
              actionButton("load_db_cellchat", "Load Database"),
              p(strong("GaspouDB Components:")),
              tags$ul(
                tags$li("CellChatDB base interactions"),
                tags$li("MultiNicheDB database"),
                tags$li("CellPhoneDB database"),
                tags$li("CellTalkDB database")
              ),
          ),
          box(title = "Processing Logs", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
              verbatimTextOutput("loading_logs")
          ),
        )
      ),
      ############################## Ligand receptor analysis##############################
      tabItem(tabName = "ligand_receptor_cellchat",
              fluidRow(
                box(title = "Step 1: Create CellChat Objects", status = "primary", solidHeader = TRUE, width = 12,
                    column(6,
                           selectInput("group_by_cellchat", "Group cells by (e.g. cell type)", choices = NULL),
                           checkboxInput("use_subset_cellchat", "Create subset by condition", value = FALSE),
                           conditionalPanel(
                             condition = "input.use_subset_cellchat == true",
                             selectInput("subset_column_cellchat", "Select condition to compare", choices = NULL)
                           ),
                           actionButton("create_cellchat", "Create Objects", class = "btn-primary")
                    )
                )
              ),
              fluidRow(
                box(title = "Step 2: Analyze CellChat Objects", status = "primary", solidHeader = TRUE, width = 12,
                    selectInput("select_objects_cellchat", "Select objects to analyze", choices = NULL, multiple = TRUE),
                    actionButton("analyze_cellchat", "Run Analysis", class = "btn-primary")
                )
              ),
              fluidRow(
                box(title = "Bubble Plot Parameters", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
                    fluidRow(
                      column(3,
                             selectInput("subset_choose_cellchat", "Select CellChat object", choices = NULL),
                             checkboxInput("flip_axes", "Flip axes (swap X and Y)", value = FALSE),
                             checkboxInput("exclude_intra", "Exclude intra-group communication", value = TRUE)
                      ),
                      column(3,
                             selectizeInput("sources_use_cellchat", "Source cell types:", choices = NULL, multiple = TRUE,
                                            options = list(maxItems = NULL, plugins = list('remove_button'),
                                                           placeholder = 'Select source cell types')
                             ),
                             selectizeInput("targets_use_cellchat", "Target cell types:", choices = NULL, multiple = TRUE,
                                            options = list(maxItems = NULL, plugins = list('remove_button'),
                                                           placeholder = 'Select target cell types')
                             ),
                      ),
                      column(3,
                             numericInput("threshold_cellchat", "Threshold:", value = 0.05, min = 0, max = 1, step = 0.01),
                             numericInput("plot_width", "Plot width (px):", value = 900, min = 400, max = 3000, step = 100),
                             numericInput("plot_height", "Plot height (px):", value = 1500, min = 400, max = 3000, step = 100)
                      ),
                      column(3,
                             selectInput("bubble_plot_format", "Format:",
                                         choices = c("PNG" = "png", "JPEG" = "jpeg", "TIFF" = "tiff","SVG"="svg","PDF" = "pdf"),
                                         selected = "png"
                             ), br(),
                             downloadButton("download_bubble_plot", "Download Plot"),
                             downloadButton("download_info_table", "Download Info Table"),
                             br(), br(), br(),
                             actionButton("generate_plot_cellchat", "Generate Plot", class = "btn-primary")
                      )
                    ),
                    # Main display area with plot and info side by side
                    fluidRow(
                      column(9,
                             # Main plot area
                             plotOutput("bubble_plot_cellchat")
                      ),
                      column(3,
                             # Two info panels stacked
                             h4("Interaction Type Summary"),
                             tableOutput("interaction_type_summary"),
                             hr(),
                             h4("All Interactions"),
                             div(
                               tableOutput("interaction_info_table"),
                               style = "max-height: 500px; overflow-y: auto;" # Make table scrollable
                             )
                      )
                    )
                )
              )),
     # UI pour Tab 3 - Chord Plot
     tabItem(
       tabName = "circle_plot_cellchat",
       fluidRow(
         box(
           title = "Chord Plot - Ligand-Receptor Interactions", 
           status = "primary", 
           solidHeader = TRUE, 
           collapsible = TRUE, 
           width = 12,
           
           # Info box
           div(class = "well",
               h4("About Chord Plots", style = "color: #2c3e50;"),
               p("Chord diagrams visualize communication networks between cell types. The width of each chord represents the strength of communication through specific ligand-receptor pairs."),
               p("Tips: Adjust plot dimensions for better visibility when many interactions are present.")
           ),
           
           # Parameters
           fluidRow(
             column(3,
                    h5("Data Selection"),
                    selectInput("cellchat_obj_chord", "CellChat object:", choices = NULL),
                    selectizeInput("sender_groups", "Sender cell types:", 
                                   choices = NULL, 
                                   multiple = TRUE,
                                   options = list(placeholder = "Select senders")),
                    selectizeInput("receiver_groups", "Receiver cell types:", 
                                   choices = NULL, 
                                   multiple = TRUE,
                                   options = list(placeholder = "Select receivers"))
             ),
             column(3,
                    h5("Plot Parameters"),
                    numericInput("prob_threshold_chord", "Min probability:", 
                                 value = 0.05, min = 0, max = 1, step = 0.01),
                    sliderInput("chord_label_distance", "Label distance:", 
                                value = 0.1, min = -0.5, max = 1, step = 0.05),
                    sliderInput("chord_label_size", "Label size:", 
                                value = 0.8, min = 0.4, max = 2, step = 0.1)
             ),
             column(3,
                    h5("Display Settings"),
                    sliderInput("chord_plot_dims", "Plot size (px):", 
                                value = 800, min = 400, max = 2000, step = 100),
                    checkboxInput("use_custom_colors", "Custom colors", value = FALSE),
                    conditionalPanel(
                      condition = "input.use_custom_colors == true",
                      uiOutput("color_inputs_chord")
                    )
             ),
             column(3,
                    h5("Export Options"),
                    selectInput("chord_plot_format", "Format:",
                                choices = c("PNG" = "png", "PDF" = "pdf", 
                                            "SVG" = "svg", "TIFF" = "tiff"),
                                selected = "png"),
                    numericInput("chord_export_width", "Export width (px):", 
                                 value = 1200, min = 600, max = 3000),
                    br(),
                    actionButton("generate_chord_plot", "Generate Plot", 
                                 class = "btn-primary btn-block"),
                    br(),
                    downloadButton("download_chord_plot", "Download Plot",
                                   class = "btn-success btn-block")
             )
           ),
           
           # Dynamic plot area
           hr(),
           uiOutput("chord_plot_ui")
         )
       )
     ),
     
     # UI pour Tab 4 - Circle Plots
     tabItem(
       tabName = "circle_plot_global_cellchat",
       
       # Global plot section
       fluidRow(
         box(
           title = "Network Overview - Global Communication", 
           status = "primary", 
           solidHeader = TRUE, 
           collapsible = TRUE, 
           width = 12,
           
           div(class = "well",
               h4("Global Network Visualization", style = "color: #2c3e50;"),
               p("Shows overall communication patterns between all cell types. Circle size represents cell type abundance, edge thickness shows communication strength.")
           ),
           
           fluidRow(
             column(3,
                    selectInput("cellchat_obj_global", "CellChat object:", choices = NULL),
                    selectInput("plot_type_global", "Metric:",
                                choices = c("Number of interactions" = "count",
                                            "Interaction strength" = "weight"),
                                selected = "weight")
             ),
             column(3,
                    sliderInput("global_vertex_size", "Cell type size:", 
                                value = 15, min = 5, max = 30, step = 1),
                    sliderInput("global_edge_width", "Max edge width:", 
                                value = 10, min = 1, max = 20, step = 1)
             ),
             column(3,
                    sliderInput("global_label_size", "Label size:", 
                                value = 0.8, min = 0.4, max = 2, step = 0.1),
                    sliderInput("global_margin", "Plot margin:", 
                                value = 0.1, min = -0.2, max = 0.5, step = 0.05)
             ),
             column(3,
                    numericInput("global_plot_size", "Display size (px):", 
                                 value = 600, min = 400, max = 1000),
                    br(),
                    actionButton("generate_global_plot", "Generate Plot", 
                                 class = "btn-primary btn-block"),
                    br(),
                    downloadButton("download_global_plot", "Download Plot",
                                   class = "btn-success btn-block")
             )
           ),
           
           hr(),
           uiOutput("global_plot_ui")
         )
       ),
       
       # Cell type specific plots
       fluidRow(
         box(
           title = "Cell Type Focus - Outgoing & Incoming Signals", 
           status = "primary", 
           solidHeader = TRUE, 
           collapsible = TRUE, 
           width = 12,
           
           div(class = "well",
               h4("Cell Type Specific Networks", style = "color: #2c3e50;"),
               p("Visualize communication patterns for specific cell types. Each plot shows signals sent from (or received by) the selected cell type.")
           ),
           
           fluidRow(
             column(4,
                    selectInput("cellchat_obj_specific", "CellChat object:", choices = NULL),
                    selectizeInput("cell_types_to_show", "Select cell types:",
                                   choices = NULL,
                                   multiple = TRUE,
                                   options = list(
                                     placeholder = 'Choose cell types to analyze',
                                     plugins = list('remove_button')
                                   )),
                    radioButtons("signal_direction", "Signal direction:",
                                 choices = c("Outgoing (as sender)" = "outgoing",
                                             "Incoming (as receiver)" = "incoming"),
                                 selected = "outgoing")
             ),
             column(4,
                    sliderInput("specific_label_size", "Label size:", 
                                value = 0.8, min = 0.4, max = 2, step = 0.1),
                    sliderInput("specific_margin", "Plot margin:", 
                                value = 0.15, min = 0, max = 0.5, step = 0.05),
                    numericInput("n_cols_specific", "Columns:", 
                                 value = 3, min = 1, max = 6)
             ),
             column(4,
                    numericInput("specific_plot_height", "Plot height (px):", 
                                 value = 250, min = 150, max = 500),
                    br(),
                    actionButton("generate_specific_plots", "Generate Plots", 
                                 class = "btn-primary btn-block"),
                    br(),
                    downloadButton("download_specific_plots", "Download Plots",
                                   class = "btn-success btn-block")
             )
           ),
           
           hr(),
           uiOutput("specific_plots_ui")
         )
       )
     ),
  
  ############################## Spatial/Load Dataset ##############################
  tabItem(
    tabName = "load_spatial_dataset",
    div(class = "container-fluid",
        div(class = "row",
            div(class = "col-12 text-center",
                div(style = "background-color: #f8f9fa; padding: 20px; margin-bottom: 20px; border-radius: 5px;",
                    h2("Spatial Transcriptomics Analysis", style = "color: #2c3e50; font-weight: 500; margin-bottom: 10px;"),
                    p("Load and analyze your spatial transcriptomics data", style = "color: #7f8c8d; font-size: 16px;")
                )
            )
        )
    ),
    box(title = "Spatial Data Overview", status = "info", solidHeader = TRUE, collapsible = TRUE, width = 14,
        div(class = "well",
            h4("About Spatial Transcriptomics", style = "color: #2c3e50;"),
            p("Spatial transcriptomics combines gene expression analysis with spatial location information, 
            allowing you to understand how gene expression varies across tissue sections."),
            h5("Supported Technologies:", style = "font-weight: bold; color: #444;"),
            tags$ul(
              tags$li(strong("10x Visium:"), "High-resolution spatial gene expression profiling"),
              tags$li(strong("Slide-seq:"), "Single-cell resolution spatial transcriptomics"),
              tags$li(strong("STARmap:"), "In situ sequencing-based spatial transcriptomics"),
              tags$li(strong("seqFISH:"), "Sequential fluorescence in situ hybridization")
            ),
            h5("Data Requirements:", style = "font-weight: bold; color: #444;"),
            tags$ul(
              tags$li("Expression matrix (genes x spots/cells)"),
              tags$li("Spatial coordinates for each spot/cell"),
              tags$li("Optional: Histological image of the tissue section"),
              tags$li("Optional: Spot/cell metadata")
            )
        )
    ),
    box(title = "Supported Data Formats", status = "info", solidHeader = TRUE, collapsible = TRUE, width = 14,
        div(class = "well",
            h5("10x Visium Data Structure:", style = "font-weight: bold; color: #444;"),
            tags$pre(
              "your_spatial_data.zip/
├── filtered_feature_bc_matrix.h5     # OR filtered_feature_bc_matrix/ folder
└── spatial/
    ├── scalefactors_json.json
    ├── tissue_positions.parquet      # New format
    ├── tissue_positions_list.csv     # Old format (also supported)
    ├── tissue_hires_image.png        # Optional
    └── tissue_lowres_image.png       # Optional"
            ),
            p(class = "text-muted", 
              "The application automatically detects whether your data uses the new Parquet format or the older CSV format for tissue positions.")
        )
    ),
    box(title = "Dataset Parameters", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
        fluidRow(
          column(3,
                 div(class = "well",
                     selectInput("spatial_species_choice", "Select Species:",
                                 choices = c("Mouse" = "mouse", "Human" = "human", "Rat" = "rat"),
                                 selected = "mouse")
                 )
          ),
          column(3,
                 div(class = "well",
                     selectInput("spatial_technology", "Spatial Technology:",
                                 choices = c("10x Visium" = "visium",
                                             "10x Visium HD" = "visium_hd",
                                             "Slide-seq" = "slideseq",
                                             "STARmap" = "starmap",
                                             "Custom" = "custom"),
                                 selected = "visium")
                 )
          ),
          column(3,
                 div(class = "well",
                     selectInput("spatial_resolution", "Data Resolution:",
                                 choices = c("Auto-detect" = "auto",
                                             "2µm (very high)" = "2um", 
                                             "8µm (high)" = "8um",
                                             "16µm (standard)" = "16um"),
                                 selected = "auto"),  # CHANGE THIS: était probablement "0.5"
                     helpText("Select resolution or use auto-detect")
                 )
          ),
          column(3,
                 div(class = "well",
                     checkboxInput("has_image", "Include tissue image", value = TRUE),
                     conditionalPanel(
                       condition = "input.has_image == true",
                       selectInput("image_format", "Image format:",
                                   choices = c("PNG" = "png", "JPEG" = "jpg", "TIFF" = "tif"),
                                   selected = "png")
                     )
                 )
          )
        )
    ),
    box(title = "Load Spatial Dataset", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
        fluidRow(
          column(6,
                 h4("Load 10x Visium Data"),
                 fileInput('spatial_file', 'Select ZIP file containing spatial data', 
                           accept = c('.zip')),
                 p(class = "text-muted", 
                   "Upload a ZIP file containing: filtered_feature_bc_matrix.h5, spatial/scalefactors_json.json, 
                 spatial/tissue_positions_list.csv, and optionally tissue images")
          ),
          column(6,
                 h4("Load Processed Spatial Object"),
                 fileInput("load_spatial_seurat_file", "Select Spatial Seurat Object (.rds)", 
                           accept = ".rds"),
                 p(class = "text-muted", 
                   "Upload a pre-processed spatial Seurat object saved as RDS file")
               
          )
        )
    ),
    box(title = "Dataset Information", status = "success", solidHeader = TRUE, collapsible = TRUE, width = 14,
        fluidRow(
          column(4,
                 infoBoxOutput("spatial_spots_count", width = 12)
          ),
          column(4,
                 infoBoxOutput("spatial_genes_count", width = 12)
          ),
          column(4,
                 infoBoxOutput("spatial_samples_count", width = 12)
          )
        ),
        DTOutput("spatial_dataset_summary")
    )
  ),
  
  ############################## Spatial/Normalization & QC ##############################
  tabItem(
    tabName = "spatial_normalisation",
    box(title = "Spatial Quality Control Overview", status = "info", solidHeader = TRUE, collapsible = TRUE, width = 14,
        div(class = "well",
            h4("Spatial Data Quality Assessment", style = "color: #2c3e50;"),
            p("Quality control for spatial data involves both standard scRNA-seq metrics and spatial-specific considerations:"),
            tags$ol(
              tags$li(strong("Standard QC Metrics:"), "Number of genes, UMI counts, mitochondrial content"),
              tags$li(strong("Spatial Coverage:"), "Distribution of spots across the tissue section"),
              tags$li(strong("Tissue Detection:"), "Identification of spots over tissue vs. background"),
              tags$li(strong("Spatial Continuity:"), "Assessment of expression continuity across neighboring spots"),
              tags$li(strong("Data Sketching:"), "Optional subsampling for large datasets to improve processing speed")
            )
        )
    ),
    sidebarLayout(
      sidebarPanel(
        h4("QC Parameters"),
        infoBoxOutput("spatial_spots_info"),
        div(style = "display: inline-block; width: 80%;",
            sliderInput("spatial_nFeature_range", "Genes per spot", min = 0, max = 8000, value = c(50, 6000))
        ),
        div(style = "display: inline-block; width: 18%;",
            actionButton("spatial_genes_info", label = icon("exclamation-triangle"), 
                         `data-toggle` = "popover", `data-html` = "true",
                         `data-content` = "Number of genes detected per spot. Spatial spots typically detect more genes than single cells.")
        ),
        div(style = "display: inline-block; width: 80%;",
            sliderInput("spatial_nCount_range", "UMI counts per spot", min = 0, max = 50000, value = c(50, 40000))
        ),
        div(style = "display: inline-block; width: 18%;",
            actionButton("spatial_umi_info", label = icon("exclamation-triangle"), 
                         `data-toggle` = "popover", 
                         `data-content` = "Total UMI counts per spot. Higher than single cells due to multiple cells per spot.")
        ),
        div(style = "display: inline-block; width: 80%;",
            sliderInput("spatial_mt_max", "Maximum mitochondrial %", value = 20, min = 0, max = 100)
        ),
        div(style = "display: inline-block; width: 18%;",
            actionButton("spatial_mt_info", label = icon("exclamation-triangle"), 
                         `data-toggle` = "popover", 
                         `data-content` = "Mitochondrial gene percentage. May be higher in spatial data.")
        ),
        checkboxInput("filter_tissue_spots", "Keep only spots over tissue", value = TRUE),
        actionButton("spatial_qc_plots", "Generate QC Plots", class = "btn-primary"),
        br(), br(),
        actionButton("apply_spatial_qc", "Apply QC Filters", class = "btn-warning"),
        br(), br(),
        numericInput("spatial_var_features", "Number of Variable Features:", value = 3000, min = 1000, max = 8000),
        actionButton("normalize_spatial_data", "Normalize Spatial Data", class = "btn-success"),
        
        # Data Sketching Section
        conditionalPanel(
          hr(),
          h4("Data Sketching", style = "color: #666; margin-top: 20px;"),
          p("Sketching creates a representative subset for faster processing.", 
            style = "font-size: 0.9em; color: #666; margin-bottom: 15px;"),
          
          # Show recommendations
          div(style = "background-color: #f8f9fa; padding: 10px; border-radius: 5px; margin-bottom: 15px;",
              verbatimTextOutput("sketching_recommendations", placeholder = TRUE)
          ),
          
          fluidRow(
            column(12, 
                   numericInput("sketch_ncells", "Number of cells to sketch:", 
                                value = 50000, min = 1000, step = 1000)
            )
          ),
          fluidRow(
            column(12, 
                   selectInput("sketch_method", "Sketching method:", 
                               choices = c("Uniform" = "Uniform", "LeverageScore" = "LeverageScore"), 
                               selected = "Uniform")
            )
          ),
          
          # Dynamic sketch size recommendation
          div(style = "margin-bottom: 15px;",
              textOutput("sketch_size_recommendation")
          ),
          
          actionButton("apply_sketching", "Apply Sketching", 
                       class = "btn-info", style = "margin-bottom: 15px; width: 100%;")
        )
      
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("Spatial Plots",
                   fluidRow(
                     column(6, plotOutput("spatial_tissue_plot")),
                     column(6, plotOutput("spatial_qc_violin"))
                   ),
                   fluidRow(
                     column(12, plotOutput("spatial_feature_scatter"))
                   )
          ),
          tabPanel("Standard QC",
                   div(style = "height: 500px;", plotOutput("spatial_vlnplot")),
                   fluidRow(
                     column(6, plotOutput("spatial_scatter1")),
                     column(6, plotOutput("spatial_scatter2"))
                   )
          ),
          tabPanel("Variable Features",
                   plotOutput("spatial_variable_features")
          )
        )
      )
    )
  ),
  
  
  ############################## Spatial/Clustering ##############################
  tabItem(
    tabName = "spatial_clustering",
    box(title = "Spatial Clustering Overview", status = "info", solidHeader = TRUE, collapsible = TRUE, width = 14,
        div(
          h4("Spatial-Aware Clustering", style = "color: #2c3e50;"),
          p("Clustering spatial transcriptomics data can incorporate both expression similarity and spatial proximity:"),
          tags$ol(
            tags$li(strong("Expression-based Clustering:"), "Traditional clustering based on gene expression similarity"),
            tags$li(strong("Spatial Clustering:"), "Methods that consider spatial relationships between spots"),
            tags$li(strong("Domain Detection:"), "Identification of spatially coherent regions with similar expression")
          )
        )
    ),
    
    # Processing Status Box
    box(title = "Processing Status", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
        fluidRow(
          column(6,
                 infoBoxOutput("spatial_processing_status")
          ),
          column(6,
                 infoBoxOutput("spatial_assay_status")
          )
        ),

    ),
    # Add this box to the spatial gene expression UI
    box(title = "Assay Selection & Management", status = "warning", solidHeader = TRUE, collapsible = TRUE, width = 14,
        fluidRow(
          column(4,
                 selectInput("spatial_default_assay", "Select Default Assay:", 
                             choices = c("Spatial"), selected = "Spatial"),
                 actionButton("set_default_assay", "Set as Default", class = "btn-warning")
          ),
          column(4,
                 infoBoxOutput("spatial_current_assay_box", width = 12)
          ),
          column(4,
                 verbatimTextOutput("current_assay_info", placeholder = TRUE)
          )
        )
    ),
    # Dimensional Reduction Box
    box(title = "Dimensional Reduction", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
        fluidRow(
          column(3,
                 selectInput("spatial_assay_select", "Analysis Assay:", choices = c("Spatial"), selected = "Spatial"),
                 actionButton("spatial_scale_pca", "Scale Data & Run PCA", class = "btn-primary", style = "width: 100%; margin-bottom: 10px;")
          ),
          column(3,
                 numericInput("spatial_pca_dims", "Number of PCA dimensions:", value = 30, min = 10, max = 50),
                 numericInput("spatial_umap_dims", "Dimensions for UMAP:", value = 30, min = 10, max = 50)
          ),
          column(3,
                 selectInput("spatial_reduction_method", "Reduction method:",
                             choices = c("PCA" = "pca", "ICA" = "ica", "LSI" = "lsi"),
                             selected = "pca")
          ),
          column(3,
                 h5("Current Status:", style = "margin-top: 0;"),
                 verbatimTextOutput("spatial_pca_status", placeholder = TRUE),
                 verbatimTextOutput("spatial_umap_status", placeholder = TRUE)
          )
        ),
        plotOutput("spatial_elbow_plot", height = "400px")
    ),
    # Clustering Analysis Box
    box(title = "Clustering Analysis", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
        fluidRow(
          column(3,
                 numericInput("spatial_cluster_dims", "Clustering dimensions:", value = 5, min = 1, max = 50),
                 actionButton("spatial_run_neighbors_umap", "Run Neighbors + UMAP", 
                              class = "btn-primary", style = "width: 100%; margin-bottom: 10px;"),
                 actionButton("spatial_find_clusters", "Find Clusters", class = "btn-success", style = "width: 100%;")
          ),
          column(3,
                 numericInput("spatial_resolution", "Clustering resolution:", 
                              value = 0.5, min = 0.1, max = 150, step = 0.05),
                 selectInput("spatial_cluster_algorithm", "Algorithm:",
                             choices = list("Louvain" = 1, "Louvain (multilevel)" = 2, "SLM" = 3),
                             selected = 1)
          ),
          column(3,
                 checkboxInput("spatial_plot_labels", "Show cluster labels", value = TRUE),
                 numericInput("spatial_plot_dpi", "Plot resolution (DPI):", value = 300, min = 72, step = 72),
          ),
          column(3,
                 h5("Clustering Info:", style = "margin-top: 0;"),
                 verbatimTextOutput("spatial_neighbors_status", placeholder = TRUE),
                 verbatimTextOutput("spatial_clusters_status", placeholder = TRUE),
                 br(),
                 downloadButton("save_spatial_object", "Save spatial object")
          )
        )
    ),
    
    # Visualization Results
    fluidRow(
      column(6,
             box(title = "UMAP Clustering", status = "primary", solidHeader = TRUE, width = 12, height = "600px",
                 plotOutput("spatial_umap_clusters", height = "520px")
             )
      ),
      column(6,
             box(title = "Spatial Clustering", status = "primary", solidHeader = TRUE, width = 12, height = "600px",
                 plotOutput("spatial_tissue_clusters", height = "520px")
             )
      )
    ),
    
  ),
  
  
  
  ############################## interactive visualization ##############################
  tabItem(
    tabName = "spatial_interactive_visualization",
    box(title = "Interactive Spatial Visualization", status = "success", solidHeader = TRUE, collapsible = TRUE, width = 14,
        fluidRow(
          column(3,
                 selectInput("spatial_interactive_feature", "Display feature:",
                             choices = c("Clusters" = "seurat_clusters", "UMI Count" = "nCount_Spatial", 
                                         "Gene Count" = "nFeature_Spatial", "Mitochondrial %" = "percent.mt"),
                             selected = "seurat_clusters")
          ),
          column(3,
                 numericInput("spatial_interactive_pt_size", "Point size:", value = 0.8, min = 0.1, max = 3, step = 0.1),
                 numericInput("spatial_interactive_alpha", "Transparency:", value = 0.7, min = 0.1, max = 1, step = 0.1)
          ),
          column(3,
                 numericInput("spatial_interactive_max_points", "Max points:", 
                              value = 10000, min = 1000, max = 50000, step = 1000),
                 div(style = "font-size: 0.8em; color: #666;",
                     "Reduce for better performance")
          ),
          column(3,
                 actionButton("refresh_interactive_spatial", "Refresh Plot", 
                              class = "btn-primary", style = "width: 100%; margin-bottom: 10px;"),
                 downloadButton("download_interactive_spatial", "Download Plot", 
                                class = "btn-success", style = "width: 100%;")
          )
        ),
        fluidRow(
          column(12,
                 verbatimTextOutput("spatial_plot_info", placeholder = TRUE)
          )
        ),
        div(style = "height: 700px; width: 100%; overflow: auto; border: 1px solid #ddd; border-radius: 5px;",
            plotlyOutput("interactive_spatial_plot", height = "650px")
        )
    )
  ),
  
  
  
  ############################## Spatial/Gene Expression Visualization ##############################
  tabItem(
    tabName = "spatial_gene_expression_visualisation",
    box(title = "Spatial Gene Expression", status = "info", solidHeader = TRUE, collapsible = TRUE, width = 14,
        div(
          h4("Visualizing Gene Expression in Space", style = "color: #2c3e50;"),
          p("Spatial transcriptomics allows visualization of gene expression patterns directly on tissue sections:"),
          tags$ul(
            tags$li(strong("Spatial Feature Plots:"), "Show gene expression overlaid on tissue coordinates"),
            tags$li(strong("Violin Plots:"), "Display expression distribution across clusters"),
            tags$li(strong("Dot Plots:"), "Show both expression level and percentage of expressing cells"),
            tags$li(strong("Co-expression Analysis:"), "Identify genes with similar spatial patterns")
          )
        )
    ),
    box(title = "Gene Selection & Parameters", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
        fluidRow(
          column(4,
                 h6("Gene Selection (Global)", style = "color: #666; margin-bottom: 10px;"),
                 pickerInput("spatial_genes_select", "Select Genes:", choices = NULL, multiple = TRUE,
                             options = list(`actions-box` = TRUE, `live-search` = TRUE)),
                 div(style = "font-size: 0.8em; color: #666;",
                     "Selected genes will auto-fill in plot sections below")
          ),
          column(4,
                 h6("Global Parameters", style = "color: #666; margin-bottom: 10px;"),
                 selectInput("spatial_viz_assay", "Assay:", choices = c("Spatial"), selected = "Spatial"),
                 numericInput("spatial_export_dpi", "Export DPI:", value = 300, min = 72, step = 72)
          ),
          column(4,
                 h6("Dataset Information", style = "color: #666; margin-bottom: 10px;"),
                 verbatimTextOutput("spatial_assay_info", placeholder = TRUE)
          )
        )
    ),
   
    box(title = "2. Violin Plot by Clusters", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
        fluidRow(
          column(12,
                 textInput("spatial_violin_genes_text", "Genes for Violin Plot (comma-separated):", 
                           placeholder = "e.g., Myh4, Myh7, Tnnt3",
                           width = "100%")
          )
        ),
        
        fluidRow(
          column(4,
                 actionButton("show_spatial_violin", "Generate Violin Plot", class = "btn-success"),
                 checkboxInput("spatial_violin_split", "Split violin plot", value = FALSE),
                 checkboxInput("spatial_violin_points", "Show points", value = TRUE)
          ),
          column(4,
                 selectInput("spatial_cluster_order", "Cluster order:", choices = NULL, multiple = TRUE),
                 numericInput("spatial_violin_pt_size", "Point size:", value = 0.1, min = 0, max = 2, step = 0.1)
          ),
          column(4,
                 downloadButton("download_spatial_violin", "Download Violin Plot"),
                 br(), br(),
                 checkboxInput("spatial_violin_log", "Log scale", value = FALSE)
          )
        ),
        plotOutput("spatial_violin_plot", height = "500px")
    ),
        box(title = "Dot Plot by Clusters", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
        fluidRow(
          column(12,
                 textInput("spatial_dotplot_genes_text", "Genes for Dot Plot (comma-separated):", 
                           placeholder = "e.g., Myh4, Myh7, Tnnt3",
                           width = "100%")
          )
        ),
        
        fluidRow(
          column(4,
                 actionButton("show_spatial_dotplot", "Generate Dot Plot", class = "btn-info"),
                 checkboxInput("spatial_dot_scale", "Scale dot size", value = TRUE),
                 selectInput("spatial_dot_scale_by", "Scale by:", choices = c("radius", "size"), selected = "radius")
          ),
          column(4,
                 numericInput("spatial_dot_min", "Min dot size:", value = 0, min = 0, max = 10),
                 numericInput("spatial_dot_max", "Max dot size:", value = 6, min = 1, max = 20),
                 checkboxInput("spatial_dot_cluster_idents", "Cluster on y-axis", value = FALSE)
          ),
          column(4,
                 downloadButton("download_spatial_dotplot", "Download Dot Plot"),
                 br(), br(),
                 selectInput("spatial_dot_colors", "Color palette:", 
                             choices = c("RdYlBu", "viridis", "plasma", "Blues", "Reds"), 
                             selected = "RdYlBu")
          )
        ),
        plotOutput("spatial_dotplot", height = "500px")
    ),
    box(title = "Spatial Feature Plots", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
        fluidRow(
          column(4,
                 textInput("spatial_genes_text", "Genes for this plot:", 
                           placeholder = "Auto-filled from selection above or type manually",
                           width = "100%"),
                 sliderInput("spatial_image_size", "Image size:", 
                             min = 0.5, max = 2, value = 1, step = 0.1,
                             ticks = TRUE)
          ),
          column(4,
                 numericInput("spatial_pt_size", "Spot size:", value = 1.5, min = 0.1, max = 5, step = 0.1),
                 checkboxInput("spatial_combine_plots", "Combine multiple genes", value = TRUE),
                 sliderInput("spatial_brightness", "Brightness:", min = 0.5, max = 3, value = 1, step = 0.1)
          ),
          column(4,
                 actionButton("show_spatial_features", "Generate Plots", class = "btn-primary", style = "width: 100%;"),
                 downloadButton("download_spatial_features", "Download", class = "btn-success", style = "width: 100%; margin-top: 10px;")
          )
        ),
        plotOutput("spatial_feature_plots", height = "600px")
    )),
  

  ############################## Spatial/Interactive Tissue Viewer ##############################
  tabItem(
    tabName = "spatial_tissue_viewer",
    
    # Header
    box(title = "Spatial Tissue Viewer", status = "primary", solidHeader = TRUE, width = 12,
        p("Compare high-resolution histology images with transcriptomic data side by side")
    ),
    
    # Controls row
    fluidRow(
      # Left controls for H&E
      column(4,
             box(title = "H&E Image Controls", status = "info", solidHeader = TRUE, width = 12,
                 fileInput("load_hne_image", "Load H&E Image:", 
                           accept = c(".png", ".jpg", ".jpeg", ".tiff", ".tif",".mrxs")),
                 hr(),
                 h5("H&E Zoom Controls:"),
                 actionButton("hne_zoom_in", "Zoom In", icon = icon("plus"), 
                              class = "btn-primary", style = "width: 45%; margin-right: 5px;"),
                 actionButton("hne_zoom_out", "Zoom Out", icon = icon("minus"), 
                              class = "btn-warning", style = "width: 45%;"),
                 br(), br(),
                 actionButton("hne_reset_zoom", "Reset View", icon = icon("home"), 
                              class = "btn-success", style = "width: 92%;"),
                 br(), br(),
                 sliderInput("hne_zoom_factor", "Zoom Factor:", 
                             min = 1, max = 5, value = 2, step = 0.5)
             )
      ),
      
      # Center info
      column(4,
             box(title = "Viewer Info", status = "warning", solidHeader = TRUE, width = 12,
                 h5("Navigation Instructions:"),
                 tags$ul(
                   tags$li("Click on image: Center and zoom"),
                   tags$li("Use buttons for controlled zoom"),
                   tags$li("Drag sliders to pan when zoomed")
                 ),
                 hr(),
                 verbatimTextOutput("viewer_status_info")
             )
      ),
      
      # Right controls for Transcriptomics
      column(4,
             box(title = "Transcriptomics Controls", status = "success", solidHeader = TRUE, width = 12,
                 selectInput("transcriptomics_feature", "Display Feature:",
                             choices = c("Clusters" = "seurat_clusters", 
                                         "Gene Count" = "nFeature_Spatial"),
                             selected = "seurat_clusters"),
                 hr(),
                 h5("Transcriptomics Zoom:"),
                 actionButton("trans_zoom_in", "Zoom In", icon = icon("plus"), 
                              class = "btn-primary", style = "width: 45%; margin-right: 5px;"),
                 actionButton("trans_zoom_out", "Zoom Out", icon = icon("minus"), 
                              class = "btn-warning", style = "width: 45%;"),
                 br(), br(),
                 actionButton("trans_reset_zoom", "Reset View", icon = icon("home"), 
                              class = "btn-success", style = "width: 92%;"),
                 br(), br(),
                 sliderInput("trans_zoom_factor", "Zoom Factor:", 
                             min = 1, max = 5, value = 2, step = 0.5)
             )
      )
    ),
    
    # Main plots row
    fluidRow(
      # H&E Plot
      column(6,
             box(title = "H&E Histology", status = "info", solidHeader = TRUE, width = 12,
                 plotOutput("hne_plot", 
                            click = "hne_plot_click",
                            brush = brushOpts(id = "hne_plot_brush", resetOnNew = TRUE),
                            height = "600px")
             )
      ),
      
      # Transcriptomics Plot
      column(6,
             box(title = "Spatial Transcriptomics", status = "success", solidHeader = TRUE, width = 12,
                 plotOutput("transcriptomics_plot", 
                            click = "trans_plot_click",
                            brush = brushOpts(id = "trans_plot_brush", resetOnNew = TRUE),
                            height = "600px")
             )
      )
    ),
    
    # Export row
    fluidRow(
      column(12,
             box(title = "Export Options", status = "primary", width = 12, collapsible = TRUE,
                 fluidRow(
                   column(3,
                          selectInput("export_format_viewer", "Format:",
                                      choices = c("PNG" = "png", "JPEG" = "jpeg", 
                                                  "TIFF" = "tiff", "PDF" = "pdf"),
                                      selected = "png")
                   ),
                   column(3,
                          numericInput("export_dpi_viewer", "DPI:", 
                                       value = 300, min = 72, max = 600, step = 72)
                   ),
                   column(3,
                          downloadButton("download_hne_view", "Download H&E View", 
                                         class = "btn-info", style = "width: 100%;")
                   ),
                   column(3,
                          downloadButton("download_trans_view", "Download Transcriptomics", 
                                         class = "btn-success", style = "width: 100%;")
                   )
                 )
             )
      )
    )
  ),
  ############################## Spatial/Marker Analysis ##############################
  tabItem(
    tabName = "spatial_marker_analysis",
    box(title = "Spatial Cluster Marker Analysis", status = "info", solidHeader = TRUE, collapsible = TRUE, width = 14,
        div(
          h4("Find Marker Genes for Clusters", style = "color: #2c3e50;"),
          p("Identify genes that specifically characterize spatial clusters:"),
          tags$ul(
            tags$li(strong("One vs All:"), "Compare one cluster against all other clusters"),
            tags$li(strong("One vs Selected:"), "Compare one cluster against specific clusters"),
            tags$li(strong("Statistical Analysis:"), "Find differentially expressed genes using Wilcoxon test"),
            tags$li(strong("Export Results:"), "Download marker genes as CSV file")
          )
        )
    ),
    
    fluidRow(
      # Left panel for controls
      column(4,
             box(title = "Analysis Parameters", status = "primary", solidHeader = TRUE, width = 12,
                 h4("Comparison Type", style = "color: #666;"),
                 radioButtons("spatial_comparison_type", "Choose comparison:",
                              choices = c("One vs All" = "one_vs_all",
                                          "One vs Selected" = "one_vs_selected"),
                              selected = "one_vs_all"),
                 
                 hr(),
                 
                 h4("Select Clusters", style = "color: #666;"),
                 selectInput("spatial_marker_cluster", "Target cluster (ident.1):", 
                             choices = NULL),
                 
                 conditionalPanel(
                   condition = "input.spatial_comparison_type == 'one_vs_selected'",
                   selectInput("spatial_comparison_clusters", "Compare against (ident.2):", 
                               choices = NULL,
                               multiple = TRUE),
                   p("Select one or more clusters to compare against", 
                     style = "font-size: 0.9em; color: #666;")
                 ),
                 
                 hr(),
                 
                 h4("Parameters", style = "color: #666;"),
                 numericInput("spatial_marker_logfc", "Min Log2 Fold Change:", 
                              value = 0.25, min = 0, max = 2, step = 0.1),
                 numericInput("spatial_marker_min_pct", "Min % cells expressing:", 
                              value = 0.1, min = 0, max = 1, step = 0.05),
                 checkboxInput("spatial_marker_only_pos", "Only positive markers", value = TRUE),
                 
                 hr(),
                 
                 actionButton("run_spatial_markers", "Find Markers", 
                              class = "btn-primary", style = "width: 100%;"),
                 br(), br(),
                 downloadButton("download_spatial_markers", "Download Results", 
                                style = "width: 100%;")
             )
      ),
      
      # Right panel for results
      column(8,
             box(title = "Marker Gene Results", status = "primary", solidHeader = TRUE, width = 12,
                 h4(textOutput("marker_results_title"), style = "color: #666;"),
                 DTOutput("spatial_marker_results")
             )
      )
    )
  ),
  
  
      ############################## Trajectory/Monocle Conversion and trajectory ##############################
  tabItem(
    tabName = "trajectory",
    fluidRow(
      # Information box - width réduit pour moins dominer la page
      box(title = "About Trajectory Analysis with Monocle", status = "info", solidHeader = TRUE, collapsible = TRUE, width = 12,
          p("Monocle is a tool for analyzing single-cell expression data that helps reconstruct the developmental/differentiation
                trajectories of cells. It orders cells based on their transcriptional similarity to create a 'pseudotime' trajectory."),

          h4("Key Features:"),
          tags$ul(
            tags$li(strong("Convert to Monocle:"), "Transforms your Seurat object into a Monocle cell_data_set object"),
            tags$li(strong("Construct Graph:"), "Creates a trajectory graph based on transcriptional similarities"),
            tags$li(strong("Root Selection:"), "Identifies the starting point of the biological process"),
            tags$li(strong("Pseudotime Analysis:"), "Orders cells along a continuous trajectory representing their progress through the biological process")
          ),

          # Conseils simplifiés
          tags$ul(
            tags$li(icon("lightbulb"), " Tip: Choose a root cell from a population you believe represents the starting state of your process."),
            tags$li(icon("exclamation-triangle"), " Note: The number of clusters (k) affects trajectory construction - start with a moderate value.")
          ),

          # Source en plus petit
          tags$small(
            "Source: ",
            tags$a(href="http://cole-trapnell-lab.github.io/monocle-release/",
                   "Monocle Documentation",
                   target="_blank"),
            " - Trapnell Lab"
          )
      )),

    fluidRow(
      # Box 1: Data Loading
      box(title = "Load Data", status = "primary", solidHeader = TRUE, width = 12,
          fluidRow(
            column(6,
                   fileInput('load_seurat_file_monocle', 'Choose Seurat file (.rds)', accept = ".rds", width = "100%")
            ),
            column(6,
                   fluidRow(
                     column(6,
                            infoBoxOutput("monocle_cells_info", width = 12)
                     ),
                     column(6,
                            infoBoxOutput("monocle_clusters_info", width = 12)
                     )
                   )
            )
          )
      ),
      
      # Box 2: Trajectory Analysis
      box(title = "Trajectory Analysis", status = "success", solidHeader = TRUE, width = 12,
          fluidRow(
            column(4,
                   h5("Step 1: Convert", style = "color: #666;"),
                   actionButton("convertToMonocle", "Convert to Monocle", 
                                class = "btn-primary", icon = icon("exchange-alt"), 
                                style = "width: 100%;")
            ),
            column(4,
                   h5("Step 2: Build Trajectory", style = "color: #666;"),
                   actionButton("constructGraph", "Construct Graph", 
                                class = "btn-success", icon = icon("project-diagram"), 
                                style = "width: 100%;")
            ),
            column(4,
                   h5("Step 3: Set Root", style = "color: #666;"),
                   selectInput("root_cluster_select", label = NULL, choices = NULL, width = "100%"),
                   actionButton("set_root_cell", "Calculate Pseudotime", 
                                class = "btn-warning", icon = icon("flag"), 
                                style = "width: 100%; margin-top: 5px;")
            )
          )
      
    )
    ,



      box(title = "Trajectory Plot", status = "primary", solidHeader = TRUE, width = 6,
          plotOutput("trajectoryPlot", height = "500px")
      ),
      box(title = "Pseudotime Distribution", status = "primary", solidHeader = TRUE, width = 6,
          plotOutput("pseudotimePlot", height = "500px")
      ),

      box(title = "Download Options", status = "primary", solidHeader = TRUE, width = 12,
          column(3,
                 selectInput("trajectory_download_format", "File Format:",
                             choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg", "TIFF" = "tiff"),
                             selected = "png")
          ),
          column(3,
                 numericInput("trajectory_download_dpi", "Resolution (DPI):",
                              value = 300, min = 72, max = 1200, step = 10)
          ),
          column(3,
                 downloadButton("download_trajectory_umap", "Download Trajectory")
          ),
          column(3,
                 downloadButton("download_pseudotime_umap", "Download Pseudotime")
          )
      ),
    )
  )
  ,
      ############################## Trajectory/Differential expressed genes along the trajectory ##############################
  tabItem(
    tabName = "genes_pseudotime",
    box(title = "About Pseudotime Gene Expression Analysis", status = "info", solidHeader = TRUE, collapsible = TRUE, width = 14,
        p("This analysis identifies genes that change their expression patterns along the trajectory, revealing the molecular progression of cells through biological processes."),
        h4("Analysis Steps:"),
        tags$ul(
          tags$li(strong("Differential Gene Test:"), "Identifies genes that significantly change expression over pseudotime"),
          tags$li(strong("Expression Visualization:"), "Plots expression patterns of selected genes along the trajectory"),
          tags$li(strong("Result Export:"), "Download differential expression results for further analysis")
        ),
        div(style="background-color: #f8f9fa; padding: 10px; border-radius: 5px; margin-top: 10px;",
            tags$ul(
              tags$li(icon("lightbulb"), " Tip: Focus on genes with low q-values for the most reliable results"),
              tags$li(icon("info-circle"), " The q-value cutoff controls the false discovery rate in your analysis"),
              tags$li(icon("chart-line"), " Visualize multiple genes together to identify co-regulated patterns")
            )
        ),
        div(style = "border-top: 1px solid #ddd; margin-top: 20px; padding-top: 10px; color: #666; font-size: 0.9em;",
            p("For more information, see: ",
              tags$a(href="http://cole-trapnell-lab.github.io/monocle-release/docs/#differentialgetest-details-and-options",
                     "Monocle Documentation - Differential Expression Testing", target="_blank")
            )
        )
    ),
    fluidRow(
      box(title = "Differential Expression Analysis", status = "primary", solidHeader = TRUE, width = 12,
          column(4,
                 actionButton("run_diff_gene_pseudotime", "Run Differential Gene Test", class = "btn-primary"),
                 div(style="margin-top: 5px;", "Identifies genes that change significantly over pseudotime")
          ),
          column(4,
                 downloadButton("download_pseudotime_diff_genes", "Download Results (CSV)", class = "btn-info"),
                 div(style="margin-top: 5px;", "Export complete analysis results")
          ),
          column(4,
                 numericInput("sig_genes_cutoff", "q-value cutoff:", value = 0.05, min = 0, max = 1, step = 0.01),
                 div(style="margin-top: 5px;", "Adjust significance threshold for differential expression")
          ),
          DTOutput("diffGeneTable")
      )
    ),
    
    fluidRow(
      box(title = "Gene Expression Visualization", status = "primary", solidHeader = TRUE, width = 12,
          fluidRow(
            column(4,
                   pickerInput("gene_selection", 
                               "Select Genes to Visualize:", 
                               choices = NULL, 
                               multiple = TRUE, 
                               options = list(
                                 `actions-box` = TRUE, 
                                 `live-search` = TRUE,
                                 `selected-text-format` = "count > 3",
                                 title = "Search and select genes"
                               ))
            ),
            column(4,
                   actionButton("generate_gene_path", "Generate Plot",  class = "btn-primary", icon = icon("chart-line"), style = "width: 100%; margin-top: 25px;"),
                   br(),
                   downloadButton("download_gene_path_plot",  "Download",class = "btn-primary", icon = icon("download"),style = "width: 100%; margin-top: 25px;")
            ),
            column(4,
             
                   
                   selectInput("plot_format_monocle",  "Export Format:",  choices = c("PNG" = "png", "PDF" = "pdf", "TIFF" = "tiff"), selected = "png"),
                   numericInput("plot_dpi_monocle",  "Resolution (DPI):", value = 300, min = 72, max = 600, step = 72)
                    )
          ),
          hr(style = "border-top: 1px solid #ccc; margin: 15px 0;"),
          fluidRow(
            column(12,
                   div(style = "position: relative;",
                       div(style = "position: absolute; right: 0; top: -40px; z-index: 10;",
                       ),
                       plotOutput("gene_path_plot_output", height = "600px")
                   )
            )
          )
      )
    )
  ),

                  #####################Gene modules#############################
  tabItem(
    tabName = "gene_modules",
    fluidRow(
      box(title = "Gene Module Analysis", status = "info", solidHeader = TRUE, collapsible = TRUE, width = 12,
          p("This analysis groups genes with similar expression patterns along the trajectory into modules."),
          p(strong("Prerequisites:"), " Run the differential gene test first to identify genes that change along pseudotime."),
          tags$ul(
            tags$li("Modules help identify co-regulated genes"),
            tags$li("Each module represents a distinct expression pattern"),
            tags$li("Useful for understanding biological processes along the trajectory")
          )
      )
    ),
    fluidRow(
      box(title = "Find Gene Modules", status = "primary", solidHeader = TRUE, width = 6,
          h4("Module Detection Parameters"),
          p("Resolution controls the number of modules: lower values = fewer modules"),
          numericInput("resolution_value", "Resolution:", value = 0.001, min = 0.0001, max = 1, step = 0.0001),
          br(),
          actionButton("find_gene_modules", "Find Gene Modules", class = "btn-primary", icon = icon("search")),
          br(), br(),
          downloadButton("download_module_genes", "Download Module Assignments", class = "btn-info")
      ),
      box(title = "Module Summary", status = "primary", solidHeader = TRUE, width = 6,
          DTOutput("module_summary_table")
      )
    ),
    fluidRow(
      box(title = "Module Heatmap", status = "primary", solidHeader = TRUE, width = 6,
          actionButton("generate_module_heatmap", "Generate Module Heatmap", class = "btn-success"),
          br(), br(),
          plotOutput("module_heatmap", height = "500px"),
          br(),
          downloadButton("download_module_heatmap", "Download Heatmap")
      ),
      box(title = "Module Visualization", status = "primary", solidHeader = TRUE, width = 6,
          selectInput("selected_modules", "Select Modules to Visualize:",
                      choices = NULL, multiple = TRUE),
          numericInput("genes_per_module", "Top genes per module:", value = 3, min = 1, max = 10),
          actionButton("visualize_modules", "Show Module Genes", class = "btn-success"),
          br(), br(),
          plotOutput("module_visualization", height = "500px"),
          br(),
          downloadButton("download_module_visualization", "Download Plot")
      )
    ),
    fluidRow(
      box(title = "Export Settings", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
          column(4,
                 selectInput("module_download_format", "File Format:",
                             choices = c("PNG" = "png", "PDF" = "pdf", "TIFF" = "tiff"),
                             selected = "png")
          ),
          column(4,
                 numericInput("module_plot_dpi", "Resolution (DPI):", value = 300, min = 72, max = 600)
          ),
          column(4,
                 numericInput("module_plot_width", "Width (inches):", value = 10, min = 5, max = 20),
                 numericInput("module_plot_height", "Height (inches):", value = 8, min = 5, max = 20)
          )
      )
    )
  ),
    
    
      ############################## Acknowlegment and Licence ##############################
      tabItem(tabName = "acknowledgement",
              div(class = "container-fluid", style = "padding: 0;",
                  div(class = "well", style = "background: linear-gradient(rgba(0,0,0,0.7), rgba(0,0,0,0.7)), url('muscle.png');
                                   background-size: cover;
                                   background-position: center;
                                   color: white;
                                   padding: 40px;
                                   min-height: 800px;",
                      h1("Acknowledgements & License", class = "text-center", style = "font-size: 36px; margin-bottom: 40px;"),
                      div(class = "row",
                          div(class = "col-md-10 col-md-offset-1",
                              div(class = "section", style = "margin-bottom: 30px;",
                                  h3("Development", style = "color: #7ACFB0;"),
                                  p("This application was developed by Gaspard Macaux and is the property of the Neuromuscular Development, Genetics and Physiopathology laboratory directed by Dr. Pascal Maire.")
                              ),
                              div(class = "section", style = "margin-bottom: 30px;",
                                  h3("Contributors", style = "color: #7ACFB0;"),
                                  p("Special thanks to:"),
                                  tags$ul(
                                    tags$li("Edgar Jauliac, Léa Delivry and Hugues Escoffier for their expertise in transcriptomic analysis"),
                                    tags$li("Maxime Di Gallo for testing the application")
                                  )
                              ),
                              div(class = "section", style = "margin-bottom: 30px;",
                                  h3("Technologies", style = "color: #7ACFB0;"),
                                  tags$ul(
                                    tags$li("Seurat - Comprehensive single-cell analysis toolkit"),
                                    tags$li("Shiny - Web application framework"),
                                    tags$li("R Studio - Development environment"),
                                    tags$li("Cell Chat - Ligand-Receptor analysis"),
                                    tags$li("Monocle - Trajectory analysis")
                                  )
                              ),
                              div(class = "section", style = "margin-top: 40px;",
                                  h3("License", style = "color: #7ACFB0;"),
                                  p("This application is licensed under the GPL3."),
                                  tags$a(href = "https://www.gnu.org/licenses/gpl-3.0.html", "Learn more about GPL3", style = "color: #7ACFB0;")
                              )
                          )
                      )
                  )
              )
      ) # Fin du dernier tabItem (acknowledgement)
  ), # Fermeture de tabItems
# Script pour les popovers
tags$script(HTML('$(function () { $("[data-toggle=\'popover\']").popover(); });'))
) # Fermeture de dashboardBody
)
