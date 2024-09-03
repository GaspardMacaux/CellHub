
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

      menuItem("Acknowledgement & Licence", tabName = "acknowledgement", icon = icon("heart"))
    )
  ),

  ############################## Single/Introduction ##############################

  dashboardBody(
    tabItems(
      tabItem(
        tabName = "introduction",
        h1("Introduction"),
        p("Welcome to our Single-Cell and Single-Nucleus RNA sequencing (scRNA-seq and snRNA-seq) data analysis application.
                  The data you will be exploring are produced using 10x Genomics technology, a cutting-edge method for genome-wide profiling of single cells and nuclei. This allows us to analyze gene expression at the level of a single cell or nucleus, offering unprecedented resolution for understanding biological mechanisms.
                  This application utilizes the Seurat library, a popular R platform for scRNA-seq and snRNA-seq analysis. Seurat offers a suite of tools for quality, analysis, exploration, and visualization of such data.
                  A key principle of this application is the clustering of cells based on their gene expression profiles. This means we group cells into subsets (clusters) based on the similarity of their gene expression profiles. These clusters can often correspond to different cell types or distinct cellular states.
                  Additionally, we focus on differentially expressed genes, that is, genes that are significantly more or less expressed in one cluster compared to others. These differentially expressed genes provide us with valuable clues about the unique characteristics of cells in each cluster.
                  We hope you will find this application useful and informative. Happy exploring!"),
        div(style = "text-align: center; display: block;",
            tags$img(src = 'pipeline.png', style="max-height: 100vh; max-width: 140vw; height: auto; width: auto;")
        )
      )
      ,

      tabItem(tabName = "intro",
              fluidPage(
                h2("Single Dataset Analysis"),
                p("To load your transcriptomics data, please compress the three 10X files (barcodes.tsv.gz, matrix.mtx.gz and features.tsv.gz) into a single .zip file."),
                box(title = "Options", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                    column(6,
                           radioButtons("species_choice", label = "Select Species:", choices = list("Mouse" = "mouse", "Human" = "human"), selected = "mouse")
                    ),
                    column(6,
                           radioButtons("dataset_type", "Choose Dataset Type", choices = list("snRNA-seq" = "snRNA", "Multiome" = "multiome"))
                    )
                ),
                box(title = "Load Dataset", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
                    column(6,
                           fileInput('file', 'Choose a ZIP file ', accept = c('.zip')),
                    ),
                    column(6,
                           fileInput("load_seurat_file", "Choose Seurat Object", accept = ".rds")
                    )
                )


              )),

      ############################## Single/QC metrics and normalization  ##############################

      tabItem(
        tabName = "qc",

        sidebarLayout(


          sidebarPanel(

            div(style = "display: inline-block; width: 80%;",
                sliderInput("nFeature_range", "Unique genes detected in each cell", min = 0, max = 5000, value = c(200, 2500))
            ),
            div(style = "display: inline-block; width: 18%;",actionButton("exclamation1",label = icon("exclamation-triangle"), `data-toggle` = "popover", `data-html` = "true",`data-content` = "The number of unique genes detected in each cell. <ul><li>Low-quality cells or empty droplets will often have very few genes</li><li>Cell doublets or multiplets may exhibit an aberrantly high gene count.</li><li>Be careful, in most cases it's better to stay between 200 and 2500</li></ul>")
            )
            ,
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
            infoBoxOutput("nuclei_count"),

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

      ############################## Single/Scaling, PCA and elbow plot  ##############################

      tabItem(
        tabName = "perform_linear_dimensional_reduction",
        box(  title = "Select genes to visualize", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
              column(width = 6,
                     div(
                       class = "button-group",
                       actionButton("scale_button", "Scale Data & PCA"),
                       actionButton("exclamation4",label = icon("exclamation-triangle"), `data-toggle` = "popover", `data-html` = "true",
                                    `data-content` = "Next, we apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA. The ScaleData() function:<ul><li>Shifts the expression of each gene, so that the mean expression across cells is 0</li><li>Scales the expression of each gene, so that the variance across cells is 1</li><li>This step gives equal weight in downstream analyses, so that highly-expressed genes do not dominate") ),
              ),
              column(width = 6,
                     div(
                       class = "button-group",
                       actionButton("exclamation5",label = icon("exclamation-triangle"),
                                    `data-toggle` = "popover", `data-html` = "true",
                                    `data-content` = "PCA is a very powerful technique and can have much utility if you can get a grasp of it and what it means. It was initially developed to analyse large volumes of data in order to tease out the differences/relationships between the logical entities being analysed (for example, a data-set consisting of a large number of samples, each with their own data points/variables). It extracts the fundamental structure of the data without the need to build any model to represent it. This ‘summary’ of the data is arrived at through a process of reduction that can transform the large number of variables into a lesser number that are uncorrelated (i.e. the ‘principal components'), whilst at the same time being capable of easy interpretation on the original data.")) ,
              ),
              verbatimTextOutput("pca_results"),
              plotOutput("loading_plot"),
        ),

        box(title = "Select genes to visualize", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
            div(
              class = "button-group",
              actionButton("run_elbow", "Show ElbowPlot"),
              actionButton("exclamation6",label = icon("exclamation-triangle"), `data-toggle` = "popover", `data-html` = "true",`data-content` = "An alternative heuristic method generates an ‘Elbow plot’: a ranking of principle components based on the percentage of variance explained by each one") ),
            plotOutput("elbow_plot")

        ),
      ),

      ############################## Single/Neighbors calculation and clustering ##############################
      tabItem(
        tabName = "cluster_the_cells",
        box(  title = "Cluster cells", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
              column(width = 4,
                     numericInput("dimension_1", "Number of dimensions :", value = 10, min = 1),
                     actionButton("run_neighbors", "Find neighbors"),
                     actionButton("exclamation8",label = icon("exclamation-triangle"), `data-toggle` = "popover", `data-html` = "true",`data-content` = "In single-cell analyses, calculating nearest neighbors involves determining which cells are most similar to each other based on their gene expression profiles. This is vital for several downstream analyses, including clustering and defining cellular trajectories. The choice of dimensions influences how similarities are computed, directly impacting the precision of identified neighbors. Selecting an optimal number of dimensions—often via methods like PCA—ensures that the variability explained is maximized while avoiding overfitting, thus providing a robust basis for identifying true biological relationships among cells."),
                     checkboxInput("remove_axes", "Remove Axes", FALSE)
              ),
              column(width = 4,
                     numericInput("resolution_step1", "Resolution :", min = 0.1, max = 2, step = 0.1, value = 0.5),
                     actionButton("run_clustering", "Find clusters"),
                     actionButton("exclamation9",label = icon("exclamation-triangle"), `data-toggle` = "popover", `data-html` = "true",`data-content` = "Clustering in single-cell analysis groups cells together based on the similarity of their expression profiles, aiming to identify distinct cellular populations or states. The resolution parameter in clustering algorithms, like those used in Seurat,                                                                      controls the granularity of these identified clusters: a higher resolution often yields more, smaller clusters, whereas a lower resolution generates fewer, larger clusters. Selecting an appropriate resolution is pivotal, as it influences the biological interpretations by determining the discernibility of subtle cellular subtypes and the comprehensiveness of identified cell populations."),
                     checkboxInput("remove_legend", "Remove Legend", FALSE)

              ),
              column(width = 4,
                     selectInput("algorithm_select", "Select Algorithm:", choices = list("Original Louvain" = 1, "Louvain with Multilevel Refinement" = 2, "SLM Algorithm" = 3)),
                     numericInput("dpi_umap", "Image resolution:", value = 300, min = 72, step = 72),
                     downloadButton("downloadUMAP", "Download UMAP")
              ) ),
        plotOutput("clustering_plot"))

      ,

      ############################## Single/Calculate differential expressed genes for each cluster ##############################
      tabItem(
        tabName = "gene_cluster",
        box(  title = "Differentially expressed genes", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
              column(width = 4,
                     numericInput("logfc_threshold_all_single", label = "Log2 Fold Change threshold", value = 0.25),
                     numericInput("min_pct_all_single", "Percentage threshold ", value = 0.01, min = 0, max = 1, step = 0.01)

              ),
              column(width = 4,
                     sliderInput("number_genes","Number of genes:",min = 10,max = 1000, value = 10,step = 50),
                     actionButton("run_DE", "Finding differentially expressed genes"),

              ),
              column(width = 4,
                     br(),
                     downloadButton('download_DE', 'Download differentially expressed genes'),
                     br(),
                     br(),
                     br(),
                     downloadButton("save_seurat", "Save Seurat Object"),
              )),
        uiOutput("diff_genes_tables")
      )
      ,



      ############################## Single/Visualization of expressed genes ##############################
      tabItem(
        tabName = "visualizing_marker_expression",
        box(
          title = "Select genes to visualize and customize plots", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
          pickerInput("gene_select", "Select Gene:", choices = NULL, multiple = TRUE, options = list(`actions-box` = TRUE, `live-search` = TRUE)),
          numericInput("text_size", "Text Size:", value = 12, min = 1, max = 200),
          checkboxInput("hide_x_labels", "Hide X axis labels", value = FALSE),
          checkboxInput("hide_y_labels", "Hide Y axis labels", value = FALSE),
          checkboxInput("hide_x_title", "Hide X axis title", value = FALSE),
          checkboxInput("hide_y_title", "Hide Y axis title", value = FALSE),
          checkboxInput("bold_text", "Bold Text", value = FALSE),
          numericInput("dpi_plot", "Images resolution for download:", value = 300, min = 72, step = 72)
        ),
        
        # Visualize with a feature plot
        box(title = "Visualize with a feature plot", status = "primary", collapsed = TRUE, solidHeader = TRUE, collapsible = TRUE, width = 14,
            textInput("gene_list_feature", "Selected genes for FeaturePlot:", value = ""),
            fluidRow(
              column(width = 4,
                     actionButton("show_feature", "Display FeaturePlot"),
                     checkboxInput("add_noaxes_feature", "Remove Axes", FALSE),
                     checkboxInput("show_coexpression", "Show co-expression of genes", value = FALSE)
              ),
              column(width = 4,
                     downloadButton("downloadFeaturePlot", "Download Feature plot"),
                     checkboxInput("add_nolegend_feature", "Remove Legend", FALSE)
              ),
              column(width = 4,
                     selectInput("min_cutoff", "Minimum Cutoff:",
                                 choices = c("None" = NA, "q1" = "q1", "q10" = "q10", "q20" = "q20", "q30" = "q30",
                                             "q40" = "q40", "q50" = "q50", "q60" = "q60", "q70" = "q70", "q80" = "q80",
                                             "q90" = "q90", "q99" = "q99"), selected = NA),
                     selectInput("max_cutoff", "Maximum Cutoff:",
                                 choices = c("None" = NA, "q1" = "q1", "q10" = "q10", "q20" = "q20", "q30" = "q30",
                                             "q40" = "q40", "q50" = "q50", "q60" = "q60", "q70" = "q70", "q80" = "q80",
                                             "q90" = "q90", "q99" = "q99"), selected = NA)
              )
            ),
            plotOutput("feature_plot")
        ),
        
        # Visualize with a violin plot
        box(title = "Visualize with a violin plot", status = "primary", collapsed = TRUE, solidHeader = TRUE, collapsible = TRUE, width = 14,
            textInput("gene_list_vln", "Selected genes for VlnPlot:", value = ""),
            fluidRow(
              column(width = 4,
                     actionButton("show_vln", "Display VlnPlot"),
                     checkboxInput("hide_vln_points", "Hide points in VlnPlot", FALSE)
              ),
              column(width = 4,
                     downloadButton("downloadVlnPlot", "Download Violin Plot")
              )
            ),
            plotOutput("vln_plot")
        ),
        
        # Visualize with a dot plot
        box(title = "Visualize with a dot plot", status = "primary", collapsed = TRUE, solidHeader = TRUE, collapsible = TRUE, width = 14,
            textInput("gene_list_dotplot", "Selected genes for DotPlot:", value = ""),
            fluidRow(
              column(width = 4,
                     actionButton("show_dot", "Display DotPlot"),
                     checkboxInput("add_noaxes_dot", "Remove Axes", FALSE)
              ),
              column(width = 4,
                     downloadButton("downloadDotPlot", "Download DotPlot"),
                     checkboxInput("add_nolegend_dot", "Remove Legend", FALSE)
              )
            ),
            plotOutput("dot_plot")
        ),
        
        # Visualize with a ridge plot
        box(title = "Visualize with a ridge plot", status = "primary", collapsed = TRUE, solidHeader = TRUE, collapsible = TRUE, width = 14,
            textInput("gene_list_ridge_plot", "Selected genes for Ridge Plot:", value = ""),
            fluidRow(
              column(width = 4,
                     actionButton("show_ridge", "Display Ridge Plot")
              ),
              column(width = 4,
                     downloadButton("download_ridge_plot", "Download Ridge Plot")
              )
            ),
            plotOutput("ridge_plot")
        ),
        
        box(
          title = "Visualize cell/nuclei number", status = "primary", collapsed = TRUE, solidHeader = TRUE, collapsible = TRUE, width = 14,
          fluidRow(
            column(
              width = 4,
              pickerInput("gene_select_genes_analysis", "Select Genes:",choices = NULL, multiple = TRUE, options = list(`actions-box` = TRUE, `live-search` = TRUE))
            ),
            column(
              width = 4,
              numericInput("logfc_threshold", "Minimum Log2 Fold Change:",value = 1,  min =0, max = 10),
              actionButton("analyze_btn", "Analyze Expression")
            ),
            column(
              width = 4,
              downloadButton("download_genes_number_expresion", "Download table")
            )
          ),
          dataTableOutput("expression_summary")
        )
      )
      
      
      ,

      ############################## Single/Heatmaps and scatter plots ##############################
      tabItem(
        tabName = "heat_maps",
        fluidRow(
          box(
            title = "Heatmap Settings", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
            fluidRow(
              column(
                width = 4,
                pickerInput("gene_select_heatmap", "Select Genes for Heatmap:", choices = NULL, multiple = TRUE, options = list(`actions-box` = TRUE, `live-search` = TRUE)),
                actionButton("generateCustomHeatmap", "Generate Heatmap")
              ),
              column(
                width = 4,
                textInput("text_clusters", "Specify Clusters (comma-separated):", ""),
                checkboxInput("select_all_clusters", "Select All Clusters", TRUE)
              ),
              column(
                width = 4,
                numericInput("dpi_heatmap_single", "Resolution", value = 300, min = 72, step = 72),
                downloadButton("download_heatmap_single", "Download Heatmap")
              )
            ),
            plotOutput("heatmap_single")
          )
        ),
        fluidRow(
          box(
            title = "Feature Scatter Settings", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
            fluidRow(
              column(
                width = 4,
                pickerInput("feature1_select", "Select the first gene:", choices = NULL, multiple = FALSE, options = list(`actions-box` = TRUE, `live-search` = TRUE)),
                pickerInput("feature2_select", "Select the second gene:", choices = NULL, multiple = FALSE, options = list(`actions-box` = TRUE, `live-search` = TRUE))

              ),
              column(
                width = 4,
                textInput("scatter_text_clusters", "Specify Clusters (comma-separated):", ""),
                checkboxInput("scatter_select_all_clusters", "Select All Clusters", TRUE),
                actionButton("generateFeatureScatter", "Generate Feature Scatter Plot")
              ),
              column(
                width = 4,
                numericInput("dpi_scatter", "Resolution", value = 300, min = 72, step = 72),
                downloadButton("download_scatter_plot", "Download Scatter Plot")
              )
            )   ,
            plotOutput("scatterPlot")

          ),

        ),
      )

      ,

      ############################## Single/Final UMAP ##############################
      tabItem(
        tabName = "assigning_cell_type_identity",
        fluidRow(
          box(  title = "Assign cluster identity", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,

                column(width = 4,
                       selectInput("select_cluster", "Select the cluster to rename:", choices = NULL),
                       textInput("rename_single_cluster", "New name for the cluster:"),
                       actionButton("rename_single_cluster_button", "Rename the selected cluster")
                ),
                column(width = 4,
                       textInput("plot_title", "Plot title:", "UMAP Final"),
                       numericInput("label_font_size", "Label font size", 5, min = 1, max = 20, step = 0.5),
                       numericInput("pt_size", "Point size:", min = 0.1, max = 3, value = 0.3, step = 0.1)
                ),
                column(width = 4,
                       selectInput("cluster_select", "Select Cluster:",  choices = NULL),
                       colourInput("cluster_colour", "Select Colour:", value = "red"),
                       actionButton("update_colour", "Update Colour")
                )
          )),
        plotlyOutput("umap_finale"),
        selectInput("plot_type_select", "Select Plot Type:",
        choices = c("FeaturePlot", "VlnPlot", "DotPlot", "RidgePlot")),
        plotOutput("selected_plot_display")

      ),


      ############################## Single/Find markers for a specific cluster ##############################
      tabItem(
        tabName = "cluster_biomarkers",
        plotOutput("umap_cluster_10"),
        fluidRow(
          column(width = 4,
                 numericInput("dpi_umap_cluster", "Resolution for the image", value = 300, min = 72, step = 72)
          ),
          column(width = 4,
                 br(),
                 downloadButton("downloadUMAPCluster", "Download UMAP plot")
          )
        ),
        box(title = "Compares one cluster with all others", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
            column(width = 4,
                   selectInput("cluster", "Cluster identifier:", choices = NULL),
                   sliderInput("number_genes_10", "Number of differential genes to display:",
                               min = 10, max = 1000, value = 10, step = 50)
            ),
            column(width = 4,
                   numericInput("logfc_threshold_single","Log2 Fold Change threshold", value = 0.1, min = 0, max = 5, step = 0.1),
                   actionButton("find_markers", "Start analysis")
            ),
            column(width = 4,
                   numericInput("min_pct_single", "Minimum percentage threshold", value = 0.01, min = 0, max = 1, step = 0.01),
                   downloadButton("download_markers_single_cluster", "Download table (csv)")
            ),
            uiOutput("gene_tables_10")
        ),
        box(title = "Compares one cluster with one other cluster", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
            column(width = 4,
                   selectInput("cluster1", "Identifier of the first cluster:", choices = NULL),
                   selectInput("cluster2", "Second cluster identifier:", choices = NULL),
                   sliderInput("n_diff_markers", "Number of differential genes to display:",
                               min = 10, max = 1000, value = 10, step = 50)
            ),
            column(width = 4,
                   numericInput("logfc_threshold_comparison", "Log2 Fold Change threshold:", value = 0.1, min = 0, max = 5, step = 0.1),
                   actionButton("compare_markers", "Start analysis")
            ),
            column(width = 4,
                   numericInput("min_pct_comparison", "Minimum percentage threshold", value = 0.01, min = 0, max = 1, step = 0.01),
                   downloadButton('download_markers_multiple_clusters', 'Download table (csv)')
            ),
            uiOutput("gene_tables_new")
        )
      ),
      ############################## Single/Subseting ##############################
      tabItem(
        tabName = "subset",
        plotOutput("global_umap"),
        box(title = "Subset by clusters", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
            selectInput("select_ident_subset", "Select the clusters to include:", choices = NULL, multiple = TRUE),
            actionButton("apply_subset", "Apply cluster based subset")
        ),

        box(title = "Subset by genes expressions", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,

            column(width = 6,
                   numericInput("expression_threshold", "Expression threshold (default 0.1)", value = 0.1),
                   textInput("gene_list_subset", "Enter genes separated by comma (,)", value = ""),
                   
                   actionButton("apply_gene_subset", "Apply Gene Subsetting")


            ),
            column(width = 6,
                   numericInput("num_genes_to_express", "Number of genes to be expressed from the list:", value = 1, min = 1),
                   br(),
                   downloadButton("download_subset_seurat", "Save subset as .RDS")
            )

        ),
        plotOutput("subset_umap"))
      ,
      ############################## Multiple/Loading Data ##############################

      tabItem(
        tabName = "merge_dataset",
        p("To load your transcriptomics data, please compress the three 10X files (barcodes.tsv.gz, matrix.mtx.gz and features.tsv.gz) into a single .zip file. "),


        box(title = "Options", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
            column(6,
                   numericInput('num_datasets', 'Enter number of datasets to upload:', value = 2, min = 2)
            ),
            column(6,
                   radioButtons("species_choice_merge", label = "Select Species:", choices = list("Mouse" = "mouse", "Human" = "human"), selected = "mouse")
            )
        ),


        box(title = "Load Datasets", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
            uiOutput("fileInputs")

        ),


        box(title = "Integration", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
            actionButton('integrate', 'Integrate'),
            actionButton('add_field', 'Add Metadata Field'),
            uiOutput("metadata_inputs"),
            actionButton('add_metadata', 'Add Metadata')
        ),


        fileInput("load_seurat_file_merge", "Load Seurat object", accept = ".rds")

      ),

      ############################## Multiple/Scaling and PCA reduction ##############################

      tabItem(
        tabName = "plot_merge",
        fluidRow(

          box(title = "Scaling & PCA", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
              actionButton("runScalePCA", "Run Scaling, PCA and Elbow Plot"),
              plotOutput("elbow_plot2")
          ),


          box(title = "Cluster cells", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12, height = 650,
              column(width = 4,
                     numericInput("dimension_2", "Number of dimensions:", value = 15, min = 1),
                     actionButton("runFindNeighbors", "Find Neighbors and run UMAP"),
                     checkboxInput("remove_axes_umap_merge", "Remove Axes", FALSE)

              ),
              column(width = 4,
                     numericInput("resolution_step2", "Resolution for clustering:", min = 0.01, step = 0.01, value = 0.5),
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
        box(  title = "Differentially expressed genes", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
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

          box(  title = "Select genes for plots", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
                selectInput("group_by_select", "Group By:", choices = c("dataset", "cluster")),
                pickerInput("geneInput_merge", "Select Genes:", choices = NULL, multiple = TRUE, options = list(`actions-box` = TRUE, `live-search` = TRUE)),
                numericInput("dpi_input_merge", "Images resolution for download:", value = 300, min = 72, step = 72),
                numericInput("text_size", "Text Size:", value = 12, min = 1, max = 200),
                checkboxInput("hide_x_labels", "Hide X axis labels", value = FALSE),
                checkboxInput("hide_y_labels", "Hide Y axis labels", value = FALSE),
                checkboxInput("hide_x_title", "Hide X axis title", value = FALSE),
                checkboxInput("hide_y_title", "Hide Y axis title", value = FALSE),
                checkboxInput("bold_text", "Bold Text", value = FALSE)

          ),

          box(title = "Visualize with a Feature plot", status = "primary", collapsed = TRUE, solidHeader = TRUE, collapsible = TRUE, width = 12,
              textInput("gene_list_feature_merge", "Selected genes for FeaturePlot:", value = ""),
              fluidRow(
                column(4,
                       actionButton("runFeaturePlot", "Run Feature Plot"),
                       checkboxInput("add_nolegend_feature_merge", "Display without legend", value = FALSE)),
                column(4,
                       checkboxInput("add_noaxes_feature_merge", "Display without axes", value = FALSE),
                       checkboxInput("show_coexpression_merge", "Show Co-expression", value = FALSE)),
                column(4, selectInput("min_cutoff_feature_merge", "Minimum Cutoff",   choices = c("None" = NA, "q1" = "q1", "q10" = "q10", "q20" = "q20", "q30" = "q30",
                                                                                                   "q40" = "q40", "q50" = "q50", "q60" = "q60", "q70" = "q70", "q80" = "q80",
                                                                                                   "q90" = "q90", "q99" = "q99"), selected = NA),
                       selectInput("max_cutoff_feature_merge", "Maximum Cutoff",   choices = c("None" = NA, "q1" = "q1", "q10" = "q10", "q20" = "q20", "q30" = "q30",
                                                                                                "q40" = "q40", "q50" = "q50", "q60" = "q60", "q70" = "q70", "q80" = "q80",
                                                                                                "q90" = "q90", "q99" = "q99"), selected = NA),
                       downloadButton("downloadFeaturePlotMerge", "Download FeaturePlot")
                          )

              ),
              plotOutput("FeaturePlot2")

          ),
          box(title = "Visualize with a Violin plot", status = "primary", collapsed = TRUE, solidHeader = TRUE, collapsible = TRUE, width = 12,
              textInput("gene_list_vln_merge", "Selected genes for VlnPlot:", value = ""),
              checkboxInput("hide_vln_points_merge", "Hide points on VlnPlot", value = FALSE),
              actionButton("runVlnPlot", "Run Vln Plot"),
             
              downloadButton("downloadVlnPlotMerge", "Download VlnPlot"),
              plotOutput("VlnPlot2")
          ),
          box(title = "Visualize with a Dot plot", status = "primary", collapsed = TRUE, solidHeader = TRUE, collapsible = TRUE, width = 12,
              textInput("gene_list_dot_merge", "Selected genes for DotPlot:", value = ""),
              actionButton("runDotPlot", "Run Dot Plot"),
              downloadButton("downloadDotPlotMerge", "Download DotPlot"),

              plotOutput("DotPlot2")
          ),
          box(title = "Visualize with a Ridge plot", status = "primary", collapsed = TRUE, solidHeader = TRUE, collapsible = TRUE, width = 12,
              textInput("gene_list_ridge_merge", "Selected genes for Ridge Plot:", value = ""),
              actionButton("runRidgePlot", "Run Ridge Plot"),
              downloadButton("downloadRidgePlotMerge", "Download Ridge Plot"),

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
                     selectInput("dataset_select_heatmap", "Group By:", choices = c("dataset", "cluster")),                              ),
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
                     numericInput("pt_size_merge", "Points size:", value = 0.1, min = 0.01, max = 5, step = 0.1)
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
          column(width = 4,
                 br(),
                 uiOutput("dataset_filter_ui")
                           ),
          column(width = 4,
                 br(),
                 numericInput("filtered_umap_plot_dpi", "Resolution for the image", value = 300, min = 72, step = 72)
          ),
          column(width = 4,
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

      tabItem(
        tabName = "subset_merge",
        plotOutput("global_umap_merge"),
        box(title = "Subset by clusters", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
            selectInput("select_ident_subset_merge", "Select the clusters to include:", choices = NULL, multiple = TRUE),
            actionButton("apply_subset_merge", "Apply Cluster-based subset")
        ),
        box(title = "Subset by genes expressions", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,

            textInput("gene_list_merge", "Enter Gene Names (comma-separated)", value = ""),
            numericInput("expression_threshold_merge", "Expression Threshold", value = 1, min = 0),
            numericInput("num_genes_to_express_merge", "Number of genes to be expressed from the list:", value = 1, min = 1),

            actionButton("apply_gene_subset_merge", "Apply Gene-based Subset")
        ),
        plotOutput("subset_umap_merge"),
        downloadButton("download_subset_merge", "Save subset as .RDS")
      ),


      ############################## Trajectory/Monocle Conversion and trajectory ##############################

      tabItem(
        tabName = "trajectory",
        fluidRow(
          box(title = "Options", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 14,
              column(6,
                     fileInput('load_seurat_file_monocle', 'Choose Seurat file'),
                     actionButton("convertToMonocle", "Convert to Monocle"),
                     actionButton("constructGraph", "Construct Graph"),

              ),

              column(6,
                     numericInput("cluster_k", "Number of clusters (k):", value = 5, min = 1, step = 1),

                     verbatimTextOutput("selected_root_cell"),
                     actionButton("startRootSelection", "Select root cell")
              )),
          plotlyOutput("trajectoryPlot"),
          plotOutput("pseudotimePlot")
        ))
      ,

      ############################## Trajectory/Differential expressed genes along the trajectory ##############################
      tabItem(
        tabName = "genes_pseudotime",
        titlePanel("Analysis of Genes Over Pseudotime"),
        fluidRow(
          box(title = "Options", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,


              column(6,
                     actionButton("run_diff_gene_pseudotime", "Run Differential Gene Test"),
              ),
              column(6,
                     downloadButton('download_pseudotime_diff_genes', 'Download Differential Genes CSV')
              ),
              DTOutput("diffGeneTable"),
          ),

      box(title = "Options", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12, height = 650,

    
    column(6,
                     pickerInput("gene_picker", "Select Genes",  choices = NULL, multiple = TRUE, options = list(`actions-box` = TRUE, `live-search` = TRUE)),
                     actionButton("visualize_gene_trajectory", "Visualize Gene Trajectory"),
             ),
    
               column(6,
                      downloadButton("download_trajectory_plot", "Download Gene Trajectory Plot"),
                     numericInput("dpi_selection", "Select DPI for Download", value = 300, min = 72, max = 1200)
               ),
              plotOutput("geneTrajectoryPlot")
          ),
    #     box(title = "Options", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12,
    #
    #
    #          column(6,
    #                 actionButton("run_analysis_button", "Run Module"),
    #                pickerInput("module_picker", "Select a Module", choices = NULL, multiple = TRUE, options = list(`actions-box` = TRUE, `live-search` = TRUE)),
    #         ),
    #          column(6,
    #                 actionButton("plot_module_button", "Plot Selected Module"),
    #         ),
    #         plotOutput("module_plot")
    #     )
    box(title = "Genes along path", status = "primary", solidHeader = TRUE, collapsible = TRUE, width = 12, height = 650,
        column(6,
               pickerInput("gene_list", "Select Genes for path", choices = NULL, multiple = TRUE, options = list(`actions-box` = TRUE, `live-search` = TRUE)),
               actionButton("plot_gene_path", "Plot Gene Dynamics Along Path")
        ),
        column(6,
               pickerInput("cell_type", "Select Cell Type", choices = NULL, multiple = TRUE, options = list(`actions-box` = TRUE, `live-search` = TRUE)),
               downloadButton("download_gene_path_plot", "Download Plot"),
               numericInput("dpi_selection_path", "Select DPI for Download", value = 300, min = 72, max = 1200)
        ),
  
               plotOutput("genePathPlot")
        )
    )
    )
      ,
      ############################## Acknowlegment and Licence ##############################

      tabItem(
        tabName = "acknowledgement",
        fluidPage( tags$div(style = "position:relative; width:630px; height:800px;",
                            tags$img(src = "muscle.png",
                                     style = "position:absolute; top:0; left:0; width:100%; height:100%;"),
                            tags$div("Acknowledgement & License",
                                     style = "position:absolute; top:5%; left:30%; color:white; font-size:30px;"),
                            tags$div("This application was developed by Gaspard Macaux and is the property of the Neuromuscular Development, Genetics and Physiopathology laboratory directed by Dr. Pascal Maire.",
                                     style = "position:absolute; top:20%; left:4%; color:white; font-size:20px;"),
                            tags$div("Many thanks to Edgar Jauliac, Léa Delivry and Hugues Escoffier for their advices and expertise in transcriptomic analysis",
                                     style = "position:absolute; top:35%; left:4%; color:white; font-size:20px;"),
                            tags$div("Many thanks to Maxime Di Gallo and Benoit Viollet for testing the app",
                                     style = "position:absolute; top:35%; left:4%; color:white; font-size:20px;"),
                            tags$div("Many thanks to Seurat, who has built a very useful and well-documented library.",
                                     style = "position:absolute; top:45%; left:4%; color:white; font-size:20px;"),
                            tags$div("Many thankd to Shiny, which makes it so easy to develop graphical user interfaces.",
                                     style = "position:absolute; top:55%; left:4%; color:white; font-size:20px;"),
                            tags$div("Many thanks to Rstudio.",
                                     style = "position:absolute; top:65%; left:4%; color:white; font-size:20px;"),
                            tags$div("This application is licensed under the GPL3.",
                                     style = "position:absolute; top:75%; left:4%; color:white; font-size:20px;"),
        ),



        ))
    )



    , tags$script(HTML('$(function () { $("[data-toggle=\'popover\']").popover(); });')))



)
