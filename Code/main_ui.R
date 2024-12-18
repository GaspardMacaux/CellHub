
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
    # Ajouter des onglets à la barre latérale
    sidebarMenu(
      menuItem("Introduction", tabName = "introduction", icon = icon("mug-hot")),
      menuItem(
        "Single Dataset Analysis", icon = icon("virus-covid"),
        menuItem("Load dataset", tabName = "intro", icon = icon("file-arrow-down")),
        menuItem("Data cleanup & Variable features", tabName = "qc", icon = icon("bug")),
        menuItem("Dimensional reduction", tabName = "perform_linear_dimensional_reduction", icon = icon("home")),
        menuItem("Clustering", tabName = "cluster_the_cells", icon = icon("brain")),
        menuItem("Doublet Detection", tabName = "doubletfinder",icon = icon("dna")),
        menuItem("Differentially expressed genes", tabName = "gene_cluster", icon = icon("react")),
        menuItem("Genes expressions", tabName = "visualizing_marker_expression", icon = icon("flask")),
        menuItem("Heat maps & Dual expression", tabName = "heat_maps", icon = icon("thermometer-half")),
        menuItem("Assigning cell type identity", tabName = "assigning_cell_type_identity", icon = icon("id-card")),
        menuItem("Cluster biomarkers", tabName = "cluster_biomarkers", icon = icon("earth-americas")),
        menuItem("Subset", tabName = "subset", icon = icon("project-diagram"))
      ),
      menuItem(
        "Multiple Datasets Analysis", icon = icon("viruses"),
        menuItem("Load datasets", tabName = "merge_dataset", icon = icon("file-arrow-down")),
        menuItem("Clustering", tabName = "plot_merge", icon = icon("brain")),
        menuItem("Differentially expressed genes", tabName = "DE_merged_dataset", icon = icon("react")),
        menuItem("Genes expressions", tabName = "visualization_merge", icon = icon("flask")),
        menuItem("Heatmaps", tabName = "heatmaps_multi", icon = icon("thermometer-half")),

        menuItem("Assigning cell type identity", tabName = "assigning_cell_type_identity_merge", icon = icon("id-card")),
        menuItem("Cluster biomarkers", tabName = "combined_analysis", icon = icon("earth-americas")),
        menuItem("Subset", tabName = "subset_merge", icon = icon("project-diagram"))
      ),

      menuItem(
        "Trajectory Analysis", icon = icon("project-diagram"),
        menuItem("Trajectory Analysis", tabName = "trajectory", icon = icon("clock")),
        menuItem("Genes pseudotime", tabName = "genes_pseudotime", icon = icon("hourglass"))
      ),

      menuItem(
        "Ligand-Receptor Analysis", icon = icon("dna"),
        menuItem("NicheNet", tabName = "nichenet_load_and_define", icon = icon("upload")),
        menuItem("Results Viewer", tabName = "nichenet_run_and_view", icon = icon("chart-line")),
        menuItem("Circos Plot", tabName = "circos_plot", icon = icon("circle")),

        menuItem("MultiNicheNet", tabName = "multinichenet_load_and_define", icon = icon("database")),
        menuItem("MultiNicheNet Results", tabName = "multinichenet_run_and_view", icon = icon("chart-line")),
        menuItem("Ligand Analysis", tabName = "genes_set_ligan_inferance", icon = icon("circle-nodes")),
        menuItem("Activity Plots", tabName = "lr_activity_plots", icon = icon("chart-area")),
        menuItem("Network View", tabName = "regulatory_network", icon = icon("diagram-project"))
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
                  h1("Welcome to Single-Cell Analysis Tool", style = "color: #2c3e50; text-align: center; margin-bottom: 25px;"),
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
      tabItem(tabName = "intro",
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
                # Section explicative
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

                # Options box
                box(title = "Dataset Parameters", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                    fluidRow(
                      column(6,
                             div(class = "well",
                                 radioButtons("species_choice", label = "Select Species:", choices = list("Mouse" = "mouse", "Human" = "human"), selected = "mouse")
                             )
                      ),
                      column(6,
                             div(class = "well",
                                 radioButtons("dataset_type", "Choose Dataset Type", choices = list("snRNA-seq" = "snRNA", "Multiome" = "multiome"))
                             )
                      )
                    )
                ),

                # Loading box
                box(title = "Load Dataset", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                    fluidRow(
                      column(6,
                             h4("Load Raw 10X Data"),
                             fileInput('file', 'Select ZIP file containing 10X data', accept = c('.zip')),
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
      ),
      ############################## Single/QC metrics and normalization  ##############################

      tabItem(
        tabName = "qc",
        # Introduction Box
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
            actionButton("normalize_data", "Normalize Data")
          ),
          mainPanel(
            div(style = "height: 500px;", plotOutput("vlnplot")),  # Augmenter la hauteur pour vlnplot
            fluidRow(
              column(6, div(style = "height: 400px;", plotOutput("scatter_plot1"))),  # Hauteur pour scatter_plot1
              column(6, div(style = "height: 400px;", plotOutput("scatter_plot2")))   # Hauteur pour scatter_plot2
            ),
            fluidRow(
              column(8, div(style = "height: 400px;", plotOutput("variable_feature_plot")))  # Hauteur pour variable_feature_plot
            )
          )
        )),

      ############################## Single/Scaling, PCA and elbow plot ##############################
      tabItem(tabName = "perform_linear_dimensional_reduction",
              # Introduction Box
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

              # Scaling and PCA Box
              box(title = "Scale Data and Run PCA", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  div(style = "margin-bottom: 20px;",
                      actionButton("scale_button", "Scale Data & Run PCA", class = "btn-primary", style = "width: 200px;")
                  ),
                  verbatimTextOutput("pca_results"),
                  plotOutput("loading_plot")
              ),

              # Elbow Plot Box
              box(title = "Determine Optimal Components", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  div(style = "margin-bottom: 20px;",
                      actionButton("run_elbow", "Show Elbow Plot", class = "btn-primary", style = "width: 200px;")
                  ),
                  plotOutput("elbow_plot")
              )
      ),

      ############################## Single/Neighbors calculation and clustering ##############################
      tabItem(tabName = "cluster_the_cells",
              # Introduction Box
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
              # Main clustering box
              box(title = "Cluster cells", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  column(width = 4,
                         numericInput("dimension_1", "Number of dimensions:", value = 10, min = 1),
                         actionButton("run_neighbors", "Find neighbors", class = "btn-primary"),
                         checkboxInput("remove_axes", "Remove Axes", FALSE)
                  ),

                  column(width = 4,
                         numericInput("resolution_step1", "Resolution:", min = 0.1, max = 2, step = 0.1, value = 0.5),
                         actionButton("run_clustering", "Find clusters", class = "btn-primary"),
                         checkboxInput("remove_legend", "Remove Legend", FALSE)
                  ),

                  column(width = 4,
                         selectInput("algorithm_select", "Select Algorithm:",
                                     choices = list("Original Louvain" = 1, "Louvain with Multilevel Refinement" = 2, "SLM Algorithm" = 3)),
                         numericInput("dpi_umap", "Image resolution:", value = 300, min = 72, step = 72),
                         downloadButton("downloadUMAP", "Download UMAP")
                  )
              ),
              # Plot output
              box(title = "Clustering Results", status = "primary", solidHeader = TRUE, width = 14,
                  plotOutput("clustering_plot")
              )
      ),

      ############################## Single/DoubletFinder ##############################
      tabItem(tabName = "doubletfinder",
              # Introduction Box
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

              # Parameters Box
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

              # Results Boxes
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

      ############################## Single/Calculate differential expressed genes for each cluster ##############################
      tabItem(tabName = "gene_cluster",
              # Introduction Box
              box(title = "Differential Expression Analysis", status = "info", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  div(
                    h4("Understanding Differential Expression", style = "color: #2c3e50;"),
                    p("This analysis identifies genes that are uniquely expressed in each cluster, helping to determine cluster identity and function:"),
                    tags$ol(
                      tags$li(strong("Fold Change Threshold:"),
                              "Sets the minimum difference in expression between clusters. Higher values indicate stronger cluster-specific expression."),
                      tags$li(strong("Percentage Threshold:"),
                              "Defines the minimum percentage of cells that must express the gene in either cluster."),
                      tags$li(strong("Number of Genes:"),
                              "Controls how many top differential genes to return per cluster.")
                    ),
                    div(style = "margin-top: 15px;",
                        tags$em("Note: Adjusting these parameters helps balance sensitivity and specificity in identifying marker genes.")
                    )
                  )
              ),

              # Parameters Box
              box(title = "Differentially Expressed Genes", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  column(width = 4,
                         numericInput("logfc_threshold_all_single", "Log2 Fold Change threshold:", value = 0.25),
                         numericInput("min_pct_all_single", "Percentage threshold:", value = 0.01, min = 0, max = 1, step = 0.01)
                  ),
                  column(width = 4,
                         sliderInput("number_genes", "Number of genes:", min = 10, max = 1000, value = 10, step = 50),
                         actionButton("run_DE", "Find Differential Genes", class = "btn-primary")
                  ),
                  column(width = 4,
                         div(style = "margin-top: 25px",
                             downloadButton("download_DE", "Download Results", class = "btn-info"),
                             div(style = "margin-top: 15px",
                                 downloadButton("save_seurat", "Save Seurat Object", class = "btn-success")
                             )
                         )
                  )
              ),

              # Results Box
              box(title = "Results Table", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  uiOutput("diff_genes_tables")
              )
      ),

      ############################## Single/Visualization of expressed genes ##############################
      tabItem(tabName = "visualizing_marker_expression",
              # Introduction Box
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

              # Gene Selection Box
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
                           numericInput("dpi_plot", "Images resolution for download:", value = 300, min = 72, step = 72)
                    )
                  )
              ),

              # Feature Plot Box
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

              # Dans l'UI, ajoutons le checkbox manquant
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
                    column(width = 3,
                           actionButton("show_vln", "Display Plot", class = "btn-primary")
                    ),
                    column(width = 3,
                           checkboxInput("hide_vln_points", "Hide Points", FALSE)
                    ),
                    column(width = 2,
                           checkboxInput("add_noaxes_vln", "Remove Axes", FALSE)
                    ),
                    column(width = 2,
                           checkboxInput("add_nolegend_vln", "Remove Legend", FALSE)  # Ajouté ici
                    ),
                    column(width = 2,
                           downloadButton("downloadVlnPlot", "Download Plot")
                    )
                  ),
                  plotOutput("vln_plot")
              ),
              # Dot Plot Box
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

              # Ridge Plot Box
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

              # Cell Number Box
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
      tabItem(tabName = "heat_maps",
              # Introduction Box
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

              # Heatmap Box
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

              # Feature Scatter Box
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
              # Introduction Box
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

              # Settings Box
              box(title = "Cluster Identity Assignment", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  fluidRow(
                    # Cluster Renaming
                    column(width = 4,
                           div(class = "well",
                               h4("Rename Clusters", style = "color: #2c3e50; margin-bottom: 15px;"),
                               selectInput("select_cluster", "Select cluster:", choices = NULL),
                               textInput("rename_single_cluster", "New name:"),
                               actionButton("rename_single_cluster_button", "Apply Name", class = "btn-primary")
                           )
                    ),

                    # Plot Settings
                    column(width = 4,
                           div(class = "well",
                               h4("Plot Settings", style = "color: #2c3e50; margin-bottom: 15px;"),
                               textInput("plot_title", "Plot title:", value = "UMAP Final"),
                               numericInput("label_font_size", "Label size:", value = 5, min = 1, max = 20, step = 0.5),
                               numericInput("pt_size", "Point size:", value = 0.3, min = 0.1, max = 3, step = 0.1)
                           )
                    ),

                    # Color Settings
                    column(width = 4,
                           div(class = "well",
                               h4("Color Settings", style = "color: #2c3e50; margin-bottom: 15px;"),
                               selectInput("cluster_select", "Select cluster:", choices = NULL),
                               colourInput("cluster_colour", "Choose color:", value = "red"),
                               actionButton("update_colour", "Update Color", class = "btn-info")
                           )
                    )
                  )
              ),

              # UMAP Plot
              box(title = "Interactive UMAP", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  plotlyOutput("umap_finale")
              ),

              # Additional Visualizations
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
      tabItem(tabName = "cluster_biomarkers",
              # Introduction Box
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

              # UMAP Display
              box(title = "Cluster Overview", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                  plotOutput("umap_plot"),
                  fluidRow(
                    column(width = 4,
                           numericInput("dpi_umap_cluster", "Plot resolution:", value = 300, min = 72, step = 72)
                    ),
                    column(width = 4,
                           downloadButton("downloadUMAPCluster", "Download Plot")
                    ),
                    column(width = 4,
                           checkboxInput("show_labels", "Show Cluster Labels", value = TRUE)
                    )
                  )
              ),

              # Global Comparison Box
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

              # Pairwise Comparison Box
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
              )
      ),
      ############################## Single/Subsetting ##############################
      tabItem(tabName = "subset",
              # Introduction Box
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

              # Original UMAP
              box(title = "Original Dataset", status = "primary", solidHeader = TRUE, width = 14,
                  plotOutput("global_umap")
              ),

              # Subsetting Options
              fluidRow(
                # Cluster-based subsetting
                column(width = 6,
                       box(title = "Subset by Clusters", status = "primary", solidHeader = TRUE, width = NULL,
                           selectInput("select_ident_subset", "Select clusters:", choices = NULL, multiple = TRUE),
                           actionButton("apply_subset", "Create Subset", class = "btn-primary")
                       )
                ),

                # Gene-based subsetting
                column(width = 6,
                       box(title = "Subset by Expression", status = "primary", solidHeader = TRUE, width = NULL,
                           numericInput("expression_threshold", "Expression threshold:", value = 0.1),
                           textInput("gene_list_subset", "Genes (comma-separated):", value = ""),
                           numericInput("num_genes_to_express", "Minimum expressed genes:", value = 1, min = 1),
                           actionButton("apply_gene_subset", "Create Subset", class = "btn-primary")
                       )
                )
              ),

              # Results & Download
              box(title = "Subset Preview", status = "primary", solidHeader = TRUE, width = 14,
                  plotOutput("subset_umap"),
                  div(style = "margin-top: 15px; text-align: center;",
                      downloadButton("download_subset_seurat", "Save as RDS", class = "btn-success")
                  )
              )
      ),

      ############################## Multiple/Loading Data ##############################
      tabItem(
        tabName = "merge_dataset",
        # Box pour le chargement des données - toujours visible
        fluidRow(
          box(
            title = "Dataset Loading Options",
            status = "primary",
            solidHeader = TRUE,
            width = 12,
            tabsetPanel(
              # Option 1: Load raw datasets
              tabPanel(
                "Load Multiple Raw Datasets",
                br(),
                div(
                  class = "well",
                  h4("Load Multiple 10X Datasets for Integration", style = "color: #666;"),
                  p("Use this option if you want to load and integrate multiple raw datasets."),
                  p("Requirements:"),
                  tags$ul(
                    tags$li("Raw data option: Compress the three 10X files (barcodes.tsv.gz, matrix.mtx.gz, and features.tsv.gz) for each dataset into separate ZIP files"),
                    tags$li("Processed data option: You can also load multiple processed Seurat objects (.rds files) that will be integrated together"),
                    tags$li("You can mix both types: some datasets as ZIP files and others as RDS files")
                  ),
                  br(),
                  actionButton("open_file_input_modal", "Load Datasets",
                               class = "btn-primary",
                               icon = icon("upload"))
                )
              ),
              # Option 2: Load pre-integrated object
              tabPanel(
                "Load Pre-integrated Object",
                br(),
                div(
                  class = "well",
                  h4("Load Pre-integrated Seurat Object", style = "color: #666;"),
                  p("Use this option if you already have a Seurat object that contains multiple integrated datasets."),
                  p("Requirements:"),
                  tags$ul(
                    tags$li("A single .rds file containing an integrated Seurat object"),
                    tags$li("The object should contain multiple datasets already integrated")
                  ),
                  br(),
                  fileInput("load_seurat_file_merge",
                            label = NULL,
                            accept = ".rds",
                            buttonLabel = "Browse...",
                            placeholder = "No file selected")
                )
              )
            )
          )
        ),

        # Box pour les métadonnées - apparaît quand des données sont chargées
        fluidRow(
          conditionalPanel(
            condition = "output.datasets_loaded",
            box(
              title = "Integration and Metadata Management",
              status = "primary",
              solidHeader = TRUE,
              collapsible = TRUE,
              width = 14,
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
      )
      ,


      ############################## Multiple/Scaling and PCA reduction ##############################

      tabItem(
        tabName = "plot_merge",
        fluidRow(

          box(title = "Scaling & PCA", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
              fluidRow(
                column(6,
                       actionButton("runScalePCA", "Run Scaling, PCA and Elbow Plot")
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
              column(width = 4,
                     numericInput("dimension_2", "Number of dimensions:", value = 15, min = 1),
                     actionButton("runFindNeighbors", "Find Neighbors and run UMAP"),
                     checkboxInput("remove_axes_umap_merge", "Remove Axes", FALSE)

              ),
              column(width = 4,
                     numericInput("resolution_step2", "Resolution for clustering:", min = 0.01, step = 0.1, value = 0.5),
                     actionButton("runFindClusters", "Find clusters"),
                     checkboxInput("remove_legend_umap_merge", "Remove Legend", FALSE)
              ),
              column(width = 4,
                     selectInput("algorithm_select", "Select Algorithm:", choices = list("Original Louvain" = 1, "Louvain with Multilevel Refinement" = 2, "SLM Algorithm" = 3)),
                     numericInput("dpi_umap_merge", "DPI for UMAP Download:", value = 300, min = 72, max = 1200, step = 72),
                     downloadButton("downloadUMAP_merge", "Download UMAP")

              ),
              plotOutput("UMAPPlot_cluster_merge")

          )

        )),

      ############################## Multiple/Calculate differential expressed genes for each cluster in the merged dataset ##############################

      tabItem(
        tabName = "DE_merged_dataset",
        box(  title = "Differentially expressed genes", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
              column(width = 4,

                     numericInput("logfc_threshold_all_multiple", label = "Log2FC threshold", value = 0.25),
                     numericInput("min_pct_all_multiple", "Pct threshold ", value = 0.01, min = 0, max = 1, step = 0.01)
              ),
              column(width = 4,
                     numericInput("number_genes_merge", "Number of genes to display:", min = 1, max = 2000, value = 10, step = 10),
                     br(),

                     actionButton("run_DE_merged", "Differentially expressed genes")
              ),
              column(width = 4,

                     br(),
                     downloadButton("save_seurat_merge", "Save Seurat Object"),
                     br(),
                     br(),
                     downloadButton('download_DE_merged', 'Download differentially expressed genes')
              )
        )
        ,


        uiOutput("diff_genes_tables_merge"),
        uiOutput("previous_tab_notification")
      ),

      ############################## Multiple/Visualize genes expressions ##############################

      tabItem(
        tabName = "visualization_merge",
        fluidRow(
          # Options générales
          box(title = "Select genes for plots", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
              fluidRow(
                # Première ligne
                column(4, selectInput("group_by_select", "Group By:", choices = c("dataset", "cluster"))),
                column(4, selectInput("viz_assay_merge", "Select Assay:", choices = c("RNA", "integrated"), selected = "RNA")),
                column(4, numericInput("dpi_input_merge", "Images resolution for download:", value = 300, min = 72, step = 72))
              ),
              fluidRow(
                # Deuxième ligne
                column(6, 
                       conditionalPanel(
                         condition = "input.group_by_select == 'dataset'",
                         selectInput("metadata_to_compare", "Compare by:", choices = NULL)
                       )
                ),
                column(6, pickerInput("geneInput_merge", "Select Genes:", choices = NULL, multiple = TRUE, options = list(`actions-box` = TRUE, `live-search` = TRUE)))
              )
          ),

          # Feature plot
          box(title = "Visualize with a Feature plot", status = "primary", collapsed = TRUE, solidHeader = TRUE, collapsible = TRUE, width = 12,
              textInput("gene_list_feature_merge", "Selected genes for FeaturePlot:", value = ""),
              fluidRow(
                column(6,
                       fluidRow(
                         column(6, actionButton("runFeaturePlot", "Run Feature Plot")),
                         column(6, checkboxInput("add_nolegend_feature_merge", "Remove Legend", FALSE))
                       ),
                       fluidRow(
                         column(6, checkboxInput("add_noaxes_feature_merge", "Remove Axes", FALSE)),
                         column(6, checkboxInput("show_coexpression_merge", "Show Co-expression", value = FALSE))
                       ),
                       fluidRow(
                         column(12, downloadButton("downloadFeaturePlotMerge", "Download FeaturePlot"))
                       )
                ),
                column(6,
                       fluidRow(
                         column(12,
                                selectInput("min_cutoff_feature_merge", "Minimum Cutoff:",
                                            choices = c("None" = NA, "q1" = "q1", "q10" = "q10", "q20" = "q20", "q30" = "q30", "q40" = "q40", "q50" = "q50", "q60" = "q60", "q70" = "q70", "q80" = "q80", "q90" = "q90", "q99" = "q99"), 
                                            selected = NA
                                )
                         )
                       ),
                       fluidRow(
                         column(12,
                                selectInput("max_cutoff_feature_merge", "Maximum Cutoff:",
                                            choices = c("None" = NA, "q1" = "q1", "q10" = "q10", "q20" = "q20", "q30" = "q30", "q40" = "q40", "q50" = "q50", "q60" = "q60", "q70" = "q70", "q80" = "q80", "q90" = "q90", "q99" = "q99"), 
                                            selected = NA
                                )
                         )
                       )
                )
              ),
              plotOutput("FeaturePlot2")
          ),
          # Violin plot
          box(title = "Visualize with a Violin plot", status = "primary",collapsed = TRUE,  solidHeader = TRUE, collapsible = TRUE, width = 12,
              textInput("gene_list_vln_merge", "Selected genes for VlnPlot:", value = ""),
              fluidRow(
                column(4,
                       actionButton("runVlnPlot", "Run Vln Plot"),
                       checkboxInput("add_nolegend_vln_merge", "Remove Legend", FALSE),
                       checkboxInput("add_noaxes_vln_merge", "Remove Axes", FALSE)
                ),
                column(4,
                       checkboxInput("hide_vln_points_merge", "Hide points", FALSE),
                       downloadButton("downloadVlnPlotMerge", "Download VlnPlot")
                )
              ),
              plotOutput("VlnPlot2")
          ),

          # Dot plot
          box(title = "Visualize with a Dot plot", status = "primary", collapsed = TRUE, solidHeader = TRUE, collapsible = TRUE, width = 12,
              textInput("gene_list_dot_merge", "Selected genes for DotPlot:", value = ""),
              fluidRow(
                column(4,
                       actionButton("runDotPlot", "Generate DotPlot"),
                       checkboxInput("add_nolegend_dot_merge", "Remove Legend", FALSE),
                       checkboxInput("add_noaxes_dot_merge", "Remove Axes", FALSE)
                ),
                column(4,
                       selectInput("cluster_order", "Select clusters to show:", choices = NULL, multiple = TRUE),
                       checkboxInput("invert_axes", "Invert Axes", FALSE)
                ),
                column(4,
                       downloadButton("downloadDotPlotMerge", "Download DotPlot")
                )
              ),
              plotOutput("DotPlot2")
          ),

          # Ridge plot
          box(title = "Visualize with a Ridge plot",status = "primary",collapsed = TRUE,solidHeader = TRUE, collapsible = TRUE,width = 12,
          textInput("gene_list_ridge_merge", "Selected genes for Ridge Plot:", value = ""),
              fluidRow(
                column(4,
                       actionButton("runRidgePlot", "Run Ridge Plot"),
                       checkboxInput("add_nolegend_ridge_merge", "Remove Legend", FALSE),
                       checkboxInput("add_noaxes_ridge_merge", "Remove Axes", FALSE)
                ),
                column(4,
                       downloadButton("downloadRidgePlotMerge", "Download Ridge Plot")
                )
              ),
              plotOutput("Ridge_plot_merge")
          ),
          box(title = "Visualize gene expression", status = "primary", collapsed = TRUE, solidHeader = TRUE, collapsible = TRUE, width = 12,
              fluidRow(
                column(4,
                       textInput("gene_list_genes_expression_merge", "Selected genes:", value = "")),

                column(4,
                       numericInput("logfc_threshold_genes_expression_merge", "Expression Threshold:", value = 1, min = 0, max = 10)),
                column(4,
                       actionButton("analyze_btn_genes_expression_merge", "Analyze Expression"),
                       downloadButton("download_genes_number_expression_merge", "Download Table")
                )),
              dataTableOutput("expression_summary_merge")
          )


        )),

      ############################## Multiple/Heatmap and dual expression multi datasets ##############################

      tabItem(
        tabName = "heatmaps_multi",
        h2("Heatmaps for Multiple Datasets"),
        fluidRow(
          box(title = "Heatmap Settings", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,height = 650,
              column(width = 4,
                     pickerInput("gene_select_heatmap_multi", "Select Genes:", choices = NULL, multiple = TRUE, options = list(`actions-box` = TRUE, `live-search` = TRUE)),
                     selectInput("dataset_select_heatmap", "Group By:", choices = c("dataset", "cluster")),
                     checkboxInput("use_top10_genes_merge", "Use top 10 genes per cluster", FALSE)                              ),
              column(width = 4,
                     actionButton("generateHeatmapMulti", "Generate Heatmap"),
                     textInput("text_clusters_merge", "Specify Clusters (comma-separated):", ""),
                     checkboxInput("select_all_clusters_merge", "Select All Clusters", TRUE)
              ),
              column(width = 4,
                     numericInput("dpi_heatmap_multi", "Resolution for Download:", value = 300, min = 72, step = 72),
                     downloadButton("download_heatmap_multi", "Download Heatmap")
              ),
              plotOutput("heatmap_plot_multi")
          )
        ),
        fluidRow(
          box(title = "Dual expression", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,height = 650,
              column(width = 4,
                     pickerInput("gene_select_scatter_multi1", "Select Genes:", choices = NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE)),
                     pickerInput("gene_select_scatter_multi2", "Select Genes:", choices = NULL, options = list(`actions-box` = TRUE, `live-search` = TRUE)),
              ),
              column(width = 4,
                     actionButton("generateScatterMulti", "Generate Dual expression"),
                     textInput("text_clusters_merge_scatter", "Specify Clusters (comma-separated):", ""),
                     checkboxInput("select_all_clusters_merge_scatter", "Select All Clusters", TRUE)),


              column(width = 4,
                     selectInput("dataset_select_scatter", "Group By:", choices = c("dataset", "cluster")),
                     numericInput("dpi_scatter_multi", "Resolution for Download:", value = 300, min = 72, step = 72),

                     downloadButton("download_scatter_multi", "Download Heatmap")
              ),
              plotOutput("scatter_plot_multi")

          )))
      ,

      ############################## Multiple/Assigning cell identity merge ##############################


      tabItem(
        tabName = "assigning_cell_type_identity_merge",
        box(  title = "Assign cell type identity", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
              column(width = 4,
                     selectInput("select_cluster_merge", "Select Cluster:", choices = NULL),
                     textInput("rename_single_cluster_merge", "New name for selected cluster:"),
                     actionButton("rename_single_cluster_merge_button", "Rename the selected cluster")
              ),
              column(width = 4,
                     textInput("plot_title_merge", "Plot title:", value = "UMAP Finale Merge"),
                     numericInput("label_font_size_merge", "Label font size", value = 5, min = 1, max = 20, step = 0.5),
                     numericInput("pt_size_merge", "Points size:", value = 0.1, min = 0.01, max = 5, step = 0.1),
                     checkboxInput("show_cluster_labels", "Show cluster labels", value = TRUE)
              ),
              column(width = 4,
                     selectInput("select_color_merge", "Select the cluster to be modified:", choices = NULL),
                     colourInput("select_cluster_merge_color", "Choose a new color for the cluster:", value = "red"),
                     actionButton("update_colour_merge_button", "Update cluster color")
              )
        ),
        plotlyOutput("umap_finale_merge"),
        selectInput("plot_type_select_merge", "Select Plot Type:",  choices = c("FeaturePlot", "VlnPlot", "DotPlot", "RidgePlot")),
        plotOutput("selected_plot_display_merge")

      ),

      ############################## Multiple/Calculation of differentially expressed genes ##############################

      tabItem(
        tabName = "combined_analysis",
        plotOutput("filtered_umap_plot"),
        fluidRow(
          column(width = 3,
                 br(),
                 uiOutput("dataset_filter_ui")

                           ),
          column(width = 3,
                 br(),
                 numericInput("filtered_umap_plot_dpi", "Resolution", value = 300, min = 72, step = 72)
          ),
          column(width = 3,
                 br(),
                 br(),
                 checkboxInput("show_labels_merge", "Show labels", value = TRUE)          ),
          column(width = 3,
                 br(),
                 br(),
                 downloadButton("download_filtered_umap_plot", "Download UMAP plot")
          )),
        box(  title = "Compares one cluster with all others", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
              column(width = 4,
                     selectInput("selected_cluster", label = "Select a cluster for comparison", choices = NULL),
                     actionButton("calculate_DE", "Start analysis")
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

                     actionButton("compare_clusters_button", "Start analysis between those Clusters")
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
                   actionButton("compare_datasets_button", "start analysis  between those datasets")


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
        box(title = "Compares a cluster between two datasets ", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
            actionButton("generate_cluster_table", "Generate Cluster Table"),

            DTOutput("cluster_table")
        ),
      ),

      ############################## Multiple/Subseting seurat object ##############################
      tabItem(tabName = "subset_merge",
              box(width = 14, plotOutput("global_umap_merge")),

              fluidRow(
                column(width = 6,
                       box(title = "Subset by clusters", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = NULL,
                           selectInput("select_ident_subset_merge", "Select the clusters to include:", choices = NULL, multiple = TRUE),
                           actionButton("apply_subset_merge", "Apply cluster based subset", class = "btn-primary")
                       )
                ),

                column(width = 6,
                       box(title = "Subset by genes expressions", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = NULL,
                           numericInput("expression_threshold_merge", "Expression threshold", value = 1, min = 0),
                           textInput("gene_list_merge", "Enter genes separated by comma (,)", value = ""),
                           numericInput("num_genes_to_express_merge", "Number of genes to be expressed from the list:", value = 1, min = 1),
                           actionButton("apply_gene_subset_merge", "Apply Gene Subsetting", class = "btn-primary")
                       )
                )
              ),

              box(width = 14, plotOutput("subset_umap_merge")),

              box(title = "Save Subset", status = "success", solidHeader = TRUE, width = 14, align = "center",
                  downloadButton("download_subset_merge", "Save subset as .RDS", class = "btn-lg")
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
                box(title = "Analysis Steps", status = "primary", solidHeader = TRUE, width = 12,
                    column(6,
                           fileInput('load_seurat_file_monocle', 'Choose Seurat file (.rds)', accept = ".rds"),
                           actionButton("convertToMonocle", "Convert to Monocle",
                                        class = "btn-primary", style = "width: 100%; margin-bottom: 10px;"),
                           actionButton("constructGraph", "Construct Graph",
                                        class = "btn-success", style = "width: 100%;")
                    ),
                    column(6,
                           numericInput("cluster_k", "Number of clusters (k):", value = 5, min = 1, step = 1),
                           verbatimTextOutput("selected_root_cell"),
                           actionButton("startRootSelection", "Select root cell",
                                        class = "btn-warning", style = "width: 100%;")
                    )
                ),
                                box(title = "Trajectory Plot", status = "primary", solidHeader = TRUE, width = 7,
                    div(class = "text-center", style = "margin-bottom: 10px;",
                        plotlyOutput("trajectoryPlot", height = "500px")
                    )
                ),
                box(title = "Pseudotime Distribution", status = "primary", solidHeader = TRUE, width = 5,
                    div(class = "text-center", style = "margin-bottom: 10px;",
                        plotOutput("pseudotimePlot", height = "500px")
                    )
                )
              )
      )
      ,
      ############################## Trajectory/Differential expressed genes along the trajectory ##############################
      tabItem(
        tabName = "genes_pseudotime",
        # Information Box
        box(title = "About Pseudotime Gene Expression Analysis", status = "info", solidHeader = TRUE, collapsible = TRUE, width = 14,
            p("This analysis identifies genes that change their expression patterns along the trajectory, revealing the molecular progression 
              of cells through biological processes."),
            
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
            
            # Source citation
            div(style = "border-top: 1px solid #ddd; margin-top: 20px; padding-top: 10px; color: #666; font-size: 0.9em;",
                p("For more information, see: ",
                  tags$a(href="http://cole-trapnell-lab.github.io/monocle-release/docs/#differentialgetest-details-and-options", 
                         "Monocle Documentation - Differential Expression Testing",
                         target="_blank")
                )
            )
        ),
        
        fluidRow(
          box(title = "Differential Expression Analysis", status = "primary", solidHeader = TRUE, width = 12,
              column(4, 
                     actionButton("run_diff_gene_pseudotime", "Run Differential Gene Test", class = "btn-primary"),
                     div(style="margin-top: 5px;",
                         "Identifies genes that change significantly over pseudotime")
              ),
              column(4, 
                     downloadButton("download_pseudotime_diff_genes", "Download Results (CSV)", class = "btn-info"),
                     div(style="margin-top: 5px;",
                         "Export complete analysis results")
              ),
              column(4, 
                     numericInput("sig_genes_cutoff", "q-value cutoff:", value = 0.05, min = 0, max = 1, step = 0.01),
                     div(style="margin-top: 5px;",
                         "Adjust significance threshold for differential expression")
              ),
              DTOutput("diffGeneTable")
          ),
          
          box(title = "Gene Expression Visualization", status = "primary", solidHeader = TRUE, width = 12,
              column(6,
                     pickerInput("gene_picker", "Select Genes:", choices = NULL, 
                                 multiple = TRUE, 
                                 options = list(`actions-box` = TRUE, `live-search` = TRUE)),
                     actionButton("visualize_gene_trajectory", "Plot Selected Genes", class = "btn-success")
              ),
              column(6,
                     numericInput("dpi_selection", "Plot Resolution (DPI):", value = 300, min = 72, max = 1200),
                     downloadButton("download_trajectory_plot", "Download Plot", class = "btn-info")
              ),
              plotOutput("geneTrajectoryPlot", height = "500px")
          )
        )
      ),
    ############################## NicheNet UI ##############################
    # UI complète de l'onglet NicheNet avec introduction et explication
    tabItem(
      tabName = "nichenet_load_and_define",

      fluidRow(
        box(title = "What is NicheNet?", status = "info", solidHeader = TRUE, collapsible = TRUE, width = 12,
            p("NicheNet is a computational method that allows you to infer which ligands,
       secreted by interacting cells (senders), have influenced the gene expression
       changes observed in target cells (receivers). This method uses ligand-receptor networks
       to predict which ligands are likely responsible for inducing specific gene expression programs
       in the receiving cell population."),
            p("To run a NicheNet analysis, you will need:"),
            tags$ul(
              tags$li("A Seurat object containing your single-cell RNA-seq data."),
              tags$li("A ligand-receptor network file for ",
                      a("Human", href = "https://zenodo.org/record/10229222", target = "_blank"),
                      " or ",
                      a("Mouse", href = "https://zenodo.org/record/10229222/files/lr_network_mouse_allInfo_30112033.rds", target = "_blank")),
              tags$li(a("A ligand-target matrix file", href = "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds", target = "_blank")),
              tags$li(a("Weighted networks file for ligand-receptor interactions", href = "https://zenodo.org/record/7074291", target = "_blank"))
            ),
            div(style = "border-top: 1px solid #ddd; margin-top: 20px; padding-top: 10px; color: #666; font-size: 0.9em;",
                p("Sources:", 
                  tags$ul(
                    tags$li(
                      tags$a(href="https://github.com/saeyslab/nichenetr", 
                             "NicheNet GitHub Repository", 
                             target="_blank")
                    ),
                    tags$li(
                      "Browaeys et al., NicheNet: modeling intercellular communication by linking ligands to target genes. ",
                      tags$i("Nature Methods"), 
                      " 17, 159–162 (2020). ",
                      tags$a(href="https://doi.org/10.1038/s41592-019-0667-5", 
                             "DOI: 10.1038/s41592-019-0667-5", 
                             target="_blank")
                    )
                  )
                )
            ),
         div(style = "text-align: center; margin-top: 20px;",  
                actionButton("open_modal_nichenet", "Open Data Loading Menu", icon = icon("upload")),
                verbatimTextOutput("load_data_status_nichenet")
            )
        )

      )
    ,
    fluidRow(
      div(class = "m-3",
          box(title = "Step 2: Define Sender and Receiver", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 6,
              p("In this section, you will select the metadata column that contains cell type identities.
            Then, choose which cell population is the 'Sender' (the cells secreting ligands) and which is the 'Receiver'
            (the cells responding to ligands)."),
              selectInput("cell_identity_column", "Select Cell Identity Column", choices = NULL),
              selectizeInput("sender_select", "Select Sender(s)", choices = NULL, multiple = TRUE),
              selectizeInput("receiver_select", "Select Receiver(s)", choices = NULL, multiple = TRUE)
          )
      ),
      div(class = "m-3",
          box(title = "Select Condition Parameters", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 6,
              p("In this section, choose the column that contains experimental conditions (e.g., treatment vs control),
             and then select the condition of interest and its reference (baseline or control condition)."),
              column(6,
              selectInput("condition_column", "Select Condition Column", choices = NULL),
              selectInput("condition_oi_select", "Condition of Interest", choices = NULL),
              selectInput("condition_reference_select", "Condition Reference", choices = NULL),
             ),
              column(6, numericInput("expression_pct", "Expression Percentage Threshold:", value = 0.05, min = 0.01,  max = 1,  step = 0.01),
                     numericInput("top_n_ligands_nichenet", "Number of top ligands to analyze:", value = 50, min = 1,   max = 200, step = 5),
                     actionButton("run_nichenet", "Run NicheNet Analysis")
      )
    )

    )))
    ,

    ############################## NicheNet Results UI ##############################
    tabItem(
      tabName = "nichenet_run_and_view",
      fluidRow(
        box(title = "NicheNet Results", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
            column(4,
                   selectInput("nichenet_output_type", "Select Output Type",
                        choices = c("Ligand Activities" = "ligand_activities",
                                    "Top Ligands" = "top_ligands",
                                    "Top Targets" = "top_targets",
                                    "Ligand Expression Dotplot" = "ligand_expression_dotplot",
                                    "Ligand Differential Expression Heatmap" = "ligand_differential_expression_heatmap",
                                    "Ligand Target Heatmap" = "ligand_target_heatmap",
                                    "Ligand Receptor Heatmap" = "ligand_receptor_heatmap",
                                    "Ligand-Target Matrix" = "ligand_target_matrix",
                                    "Ligand-Target Dataframe" = "ligand_target_df")),
            actionButton("view_nichenet_results", "View Results")),
            column(4,
                   conditionalPanel(
                     condition = "['ligand_expression_dotplot', 'ligand_differential_expression_heatmap', 'ligand_target_heatmap', 'ligand_receptor_heatmap'].includes(input.nichenet_output_type)",
                   )),
                   column(4,
                          numericInput("plot_dpi", "Resolution (DPI):", value = 300, min = 72, max = 1200, step = 50),
                          downloadButton("download_nichenet_result", "Download Result")),
                        uiOutput("nichenet_output_display"),

        ),

      ),

      )


    ,

    ############################## Circos Plot UI ##############################
    tabItem(
      tabName = "circos_plot",

      fluidRow(
        box(title = "Circos Plot Overview", status = "info", solidHeader = TRUE, collapsible = TRUE, width = 12,
            p("The Circos plot visualizes the interactions between ligands and their target genes, showing which ligands from different sender cell types affect specific targets."),
            p("Follow the steps below to configure the Circos plot:"),
            fluidRow(
              box(title = "Step 1: Assign Ligands to Sender Cells", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 6,
                  p("Assign ligands to sender cell types based on their expression levels. Adjust the function for assignment using the standard deviation multiplier."),
                  actionButton("assign_ligands", "Assign Ligands", icon = icon("tasks"))
              ),

              box(title = "Step 2: Define Ligand-Target Links", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 6,
                  p("Define the ligand-target links using the specified target type and a cutoff for interaction strength."),
                  textInput("target_type", "Target Type (e.g., LCMV-DE):", value = "LCMV-DE"),
                  numericInput("link_cutoff", "Link Cutoff (0 - 1):", value = 0.4, min = 0, max = 1, step = 0.05),
                  actionButton("define_links", "Define Links", icon = icon("link"))
              )
        )
      )),


      fluidRow(
        box(title = "Step 3: Customize Colors", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
            column(4,
                   uiOutput("donor_color_inputs"),
                   uiOutput("receiver_color_input")),
            column(4,
                   checkboxInput("transparency_circos", "Enable transparency", value = FALSE),
                   sliderInput("link_size", "Link size:", min = 0.1, max = 2, value = 1, step = 0.1),
                   sliderInput("text_size", "Text size:", min = 0.1, max = 2, value = 0.6, step = 0.1),
                   numericInput("top_n_receptors", "Number of top receptors:", value = 20, min = 1, max = 100)
            ),
            column(4,
                   selectInput("text_position", "Text position:", choices = c("Inside", "Outside", "Both")),
                   selectInput("sort_mode", "Sort mode:", choices = c("Default", "Ascending", "Descending")),
                   numericInput("gap_degree", "Gap between sectors:", value = 1, min = 0, max = 10, step = 0.1),
                   actionButton("draw_circos_plot", "Generate Circos Plot", class = "btn-primary")

        )
      )),

      fluidRow(
        box(title = "Circos Plot", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
            column(6,
                   numericInput("plot_dpi", "Resolution (DPI):", value = 300, min = 72, max = 1200, step = 50),

            ),
            column(6,
                   downloadButton('download_circo_plot', 'Download Plot')

            ),


                   plotOutput("circos_plot_output")

        )
      ))

    ,



    ############################## MultinicheNet/Load and Define Parameters ##############################
    tabItem(
      tabName = "multinichenet_load_and_define",
      fluidRow(
        box(title = "What is MultiNicheNet?", status = "info", solidHeader = TRUE, collapsible = TRUE, width = 12,
            p("MultiNicheNet extends NicheNet's capabilities to analyze ligand-receptor interactions across multiple samples and conditions.
      This allows for the study of changes in intercellular communication patterns between different biological states or experimental conditions."),
            p("To run a MultiNicheNet analysis, you will need:"),
            tags$ul(
              tags$li("A Seurat object or SingleCellExperiment object containing your single-cell RNA-seq data."),
              tags$li(HTML("A ligand-receptor network file for ",
                           "<a href='https://zenodo.org/record/10229222/files/lr_network_human_allInfo_30112033.rds' target='_blank'>Human</a>",
                           " or ",
                           "<a href='https://zenodo.org/record/10229222/files/lr_network_mouse_allInfo_30112033.rds' target='_blank'>Mouse</a>")),
              tags$li(HTML("A ligand-target matrix file for ",
                           "<a href='https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds' target='_blank'>Human</a>",
                           " or ",
                           "<a href='https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds' target='_blank'>Mouse</a>")),
              tags$li(HTML("You can find more resources on the ",
                           "<a href='https://zenodo.org/record/7074291' target='_blank'>Zenodo repository</a>"))
            ),
            div(style = "border-top: 1px solid #ddd; margin-top: 20px; padding-top: 10px; color: #666; font-size: 0.9em;",
                p("Sources:", 
                  tags$ul(
                    tags$li(
                      tags$a(href="https://github.com/saeyslab/multinichenetr", 
                             "MultiNicheNet GitHub Repository", 
                             target="_blank")
                    ),
                    tags$li(
                      "Browaeys et al., MultiNicheNet: a flexible framework for differential cell-cell communication analysis from multi-sample multi-condition single-cell transcriptomics data. ",
                      tags$i("bioRxiv"), 
                      " (2023). ",
                      tags$a(href="https://doi.org/10.1101/2023.06.13.544751", 
                             "DOI: 10.1101/2023.06.13.544751", 
                             target="_blank")
                    )
                  )
                )
            ),
            # Centered button for data loading
            div(style = "text-align: center; margin-top: 20px;",
                actionButton("open_modal", "Open Data Loading Menu", icon = icon("upload")),
                verbatimTextOutput("load_data_status")
            )
        )
      ),
        div(class = "mx-3",
            fluidRow(
              column(6,
                     div(class = "mr-2",
                         box(title = "Define Metadata Parameters", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                             column(6,
                                    selectizeInput("sample_id", "Select Sample ID", choices = NULL, options = list(placeholder = 'Type to search...', searchField = 'text', create = FALSE)),
                                    selectizeInput("group_id", "Select Group ID", choices = NULL, options = list(placeholder = 'Type to search...', searchField = 'text', create = FALSE)),
                                    selectizeInput("celltype_id", "Select Cell Type ID", choices = NULL, options = list(placeholder = 'Type to search...', searchField = 'text', create = FALSE))
                             ),
                             column(6, actionButton("define_params", "Set Metadata Parameters"))
                         )
                     )
              ),
              column(6,
                     div(class = "ml-2",
                         box(title = "Define Contrast", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                             column(6, uiOutput("group_select_ui"), actionButton("save_contrast", "Set Contrast")),
                             column(6, h4("Generated Contrast Formula"), textOutput("generated_contrast_formula"))
                         )
                     )
              )
            )
        )
      )
    ,

    ############################## MultiNicheNet/Run and View Results ##############################
      tabItem(
        tabName = "multinichenet_run_and_view",
        fluidRow(
          box(
            title = "Step 2: Set Analysis Parameters", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
            column(6,
                   numericInput("min_sample_prop", "Minimum sample proportion for gene expression", value = 0.50, step = 0.01, min = 0),
                   numericInput("fraction_cutoff", "Fraction cutoff for gene expression", value = 0.05, step = 0.01, min = 0),
                   numericInput("logFC_threshold", "LogFC threshold", value = 0.50, step = 0.01, min = 0),
                   numericInput("p_val_threshold", "P-value threshold", value = 0.05, step = 0.01, min = 0),
                   actionButton("calculate_abundance_and_frq", "Calculate Abundance & Expression"),
            ),
            column(6,
                   checkboxInput("empirical_pval", "Use empirical p-values?", value = FALSE),
                   checkboxInput("p_val_adj", "Adjust p-values for multiple testing?", value = FALSE),
                   numericInput("top_n_target", "Top N target genes per ligand", value = 250, min = 1),
                   numericInput("n_cores", "Number of cores for parallel processing", value = 8, min = 1),
                   numericInput("min_cells", "Minimul cells", value =10, min = 1),
                   actionButton("run_expression_and_DE", "Run Expression & DE Analysis"),
            ),

        )),
        fluidRow(
          box(title = "Differential Expression P-values Distribution", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
              column(12,
              plotOutput("de_pvals_hist", height = "600px"),
              downloadButton("download_de_pvals_hist", "Download Plot (TIFF)")
          ))
        )




    ),

    ############################## Ligand activities ##############################
    tabItem(
      tabName = "genes_set_ligan_inferance",
      fluidRow(
        box(title = "Geneset and ligand activity analysis Evaluation Settings", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
            column(4,
                   numericInput("logFC_threshold_geneset_multinichenet", "Log FC threshold:", value = 0.50, step = 0.05),
                   numericInput("p_val_threshold_geneset_multinichenet", "P-value threshold:", value = 0.05, step = 0.01),
                   checkboxInput("p_val_adj_geneset_multinichenet", "Use adjusted p-values", value = FALSE)
            ),
            column(4,
                   checkboxInput("verbose_ligand_multinichenet", "Show detailed progress", value = TRUE),
                   checkboxInput("ligand_activity_down_multinichenet", "Consider downregulated ligands", value = FALSE),
                   helpText("If TRUE, consider both up and down regulated genes. If FALSE, focus only on upregulating ligands."),
                            ),
            column(4,
                   numericInput("top_n_target_multinichenet", "Number of top targets:", value = 250, min = 1),
                   numericInput("n_cores_multinichenet", "Number of cores for parallel processing:",
                                value = 1, min = 1, step = 1),
                   actionButton("run_combined_analysis", "Run Geneset & Ligand Analysis"),
                   br(), br(),
                   verbatimTextOutput("analysis_status")
                   ),
        )
      ),

      fluidRow(
        box(title = "Ligand Activities", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
            downloadButton('download_ligand_activities', 'Download Ligand Activities'),
            tableOutput("ligand_activities_table_multinichenet")
        ),
        box(title = "Prioritization Results", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
            tableOutput("prioritization_table_multinichenet"),
            downloadButton('download_prioritization', 'Download Prioritization Table'),
            downloadButton('download_all_results', 'Save Complete Analysis (RDS)'),
            checkboxInput("calculate_correlations_multinichenet", "Calculate ligand-receptor-target correlations", value = FALSE)
        )
      ),

            box(title = "Circos Plot Visualization", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
          fluidRow(
            column(4,
                   numericInput("top_n_interactions_multinichenet", "Number of top interactions:",
                                value = 50, min = 1, step = 5),
                   selectInput("color_palette_multinichenet", "Color Palette:",
                               choices = c("Spectral", "Set1", "Set2", "Set3", "Paired"))
            ),


            column(4,
                   selectInput("selected_group", "Select group to download:",
                               choices = NULL),
                   numericInput("circos_plot_dpi", "Download Resolution (DPI):",
                                value = 300, min = 72, max = 1200, step = 72)
            ),
            column(4,
                   br(),
                   actionButton("generate_circos_multinichenet", "Generate Circos Plots", class = "btn-primary"),
                   br(), br(),
                   downloadButton("download_circos_plot", "Download Plot")
            )
          ),
          fluidRow(
            column(12,
                   uiOutput("circos_plots_ui_multinichenet")
            )
          )
      )
    ),


    ############################## LR Product Activity Plots ##############################

    tabItem(
      tabName = "lr_activity_plots",
      box(title = "LR Activity Plot Settings", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
          column(3,
                 selectInput("selected_group_viz", "Select Group:", choices = NULL),
                 numericInput("top_n_interactions", "Number of top interactions:",
                              value = 50, min = 1, step = 5)
          ),
          column(3,
                 # Optional sender selection
                 selectInput("specific_sender", "Select Sender (Optional):",
                             choices = NULL,
                             multiple = FALSE,
                             selected = NULL),
                 # Optional receiver selection
                 selectInput("specific_receiver", "Select Receiver (Optional):",
                             choices = NULL,
                             multiple = FALSE,
                             selected = NULL)
          ),
          column(3,
                 actionButton("generate_lr_plot", "Generate Plot")
          ),
          column(3,
                 numericInput("lr_plot_dpi", "Image resolution:",
                              value = 300, min = 72, step = 72),
                 br(),
                 downloadButton("download_lr_plot", "Download Plot")
          )
      ),
      box(title = "LR Product Activity Plot", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
          plotOutput("lr_activity_plot", height = "600px")
      )
    ),

############################## Regulatory Network ##############################

tabItem(
  tabName = "regulatory_network",
  box(title = "Network Settings", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
      fluidRow(
        column(4,
               numericInput("top_n_network", "Number of top interactions:",
                            value = 50, min = 1, step = 5),
               checkboxInput("use_correlations", "Filter by correlations", value = FALSE),
               actionButton("generate_network", "Generate Network")

        ),

        column(4,
               selectInput("sender_to_color", "Select Sender to Color:",
                           choices = NULL),
               colourInput("sender_color", "Choose Color:", "pink")
        ),
        column(4,
               br(), br(),
               numericInput("network_plot_dpi", "Image resolution:", value = 300, min = 72, step = 72),

               downloadButton("download_network", "Download Network Plot")
        )
      )
  ),
  # Zone d'affichage du réseau
  box(title = "Regulatory Network Visualization", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
      plotOutput("regulatory_network_plot", height = "800px")
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

                    # Header
                    h1("Acknowledgements & License",
                       class = "text-center",
                       style = "font-size: 36px; margin-bottom: 40px;"
                    ),

                    # Main content
                    div(class = "row",
                        div(class = "col-md-10 col-md-offset-1",
                            # Development credits
                            div(class = "section", style = "margin-bottom: 30px;",
                                h3("Development", style = "color: #7ACFB0;"),
                                p("This application was developed by Gaspard Macaux and is the property of the Neuromuscular Development,
                 Genetics and Physiopathology laboratory directed by Dr. Pascal Maire.")
                            ),

                            # Contributors
                            div(class = "section", style = "margin-bottom: 30px;",
                                h3("Contributors", style = "color: #7ACFB0;"),
                                p("Special thanks to:"),
                                tags$ul(
                                  tags$li("Edgar Jauliac, Léa Delivry and Hugues Escoffier for their expertise in transcriptomic analysis"),
                                  tags$li("Maxime Di Gallo for testing the application")
                                )
                            ),

                            # Technologies
                            div(class = "section", style = "margin-bottom: 30px;",
                                h3("Technologies", style = "color: #7ACFB0;"),
                                tags$ul(
                                  tags$li("Seurat - Comprehensive single-cell analysis toolkit"),
                                  tags$li("Shiny - Web application framework"),
                                  tags$li("R Studio - Development environment"),
                                  tags$li("Monocle - Trajectory analysis"),
                                  tags$li("NicheNet - Cell-cell communication analysis")
                                )
                            ),

                            # License
                            div(class = "section", style = "margin-top: 40px;",
                                h3("License", style = "color: #7ACFB0;"),
                                p("This application is licensed under the GPL3."),
                                tags$a(href = "https://www.gnu.org/licenses/gpl-3.0.html",
                                       "Learn more about GPL3",
                                       style = "color: #7ACFB0;")
                            )
                        )
                    )
                )
            )
    )
    )



    , tags$script(HTML('$(function () { $("[data-toggle=\'popover\']").popover(); });')))



)
