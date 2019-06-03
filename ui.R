#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  
  # Application title
  titlePanel("Seurat Shiny"),
  fluidRow(
    uiOutput('debug_ui'),
    column(3,uiOutput('select_dataset_ui')),
    column(3,uiOutput('reduce_dataset_select_ui')),
    #column(3,tags$h6(htmlOutput('value_names'))),
    column(3,radioButtons('save_rds_rb','Save after every step',c(F,T),F,inline = T)),
    column(3,selectInput('features_of_interest',"Features of interest",uploaded_features,c('IL13','IL17A'),multiple = T)),
    #actionButton('load_data','Load_Data'),
    #selectInput('select_dataset','Select Dataset',datasets,datasets[1]),
    column(12,
           textOutput('data_path_text'),
           textOutput('dataset_text')
    ),
    tabsetPanel(selected = '',
                #selected = 'Visualisation',
                #### UPLOAD ####
                tabPanel('Upload Data',
                         textOutput('test_test'),
                         shinyDirButton('directory', 'Folder select', 'Please select a folder'),
                         column(6,actionButton('reload_data','Create Seurat Object')),
                         tags$h4('### Setup the Seurat Object'),
                         tags$h5('
                                 For this tutorial, we will be analyzing the a dataset of Peripheral Blood Mononuclear Cells (PBMC) freely available from 10X Genomics. There are 2,700 single cells that were sequenced on the Illumina NextSeq 500. The raw data can be found [here](https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz).
                                 
                                 We start by reading in the data. The `Read10X` function reads in the output of the [cellranger](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger) pipeline from 10X, returning a unique molecular identified (UMI) count matrix. The values in this matrix represent the number of molecules for each feature (i.e. gene; row) that are detected in each cell (column).
                                 
                                 We next use the count matrix to create a `Seurat` object. The object serves as a container that contains both data (like the count matrix) and analysis (like PCA, or clustering results) for a single-cell dataset. For a technical discussion of the `Seurat` object structure, check out our [GitHub Wiki](https://github.com/satijalab/seurat/wiki). For example, the count matrix is stored in `pbmc[["RNA"]]@counts`.
                                 '),
                         tags$h4('### Standard pre-processing workflow'),
                         tags$h5('
                                 The steps below encompass the standard pre-processing workflow for scRNA-seq data in Seurat. These represent the selection and filtration of cells based on QC metrics, data normalization and scaling, and the detection of highly variable features.
                                 '),
                         
                         tags$h4('### QC and selecting cells for further analysis'),
                         #tabsetPanel(
                         #  tabPanel('Create Seurat Object',
                         tags$h5('
                                 Seurat allows you to easily explore QC metrics and filter cells based on any user-defined criteria. A few QC metrics [commonly used](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4758103/) by the community include
                                 
                                 * The number of unique genes detected in each cell. 
                                 + Low-quality cells or empty droplets will often have very few genes
                                 + Cell doublets or multiplets may exhibit an aberrantly high gene count
                                 * Similarly, the total number of molecules detected within a cell (correlates strongly with unique genes)
                                 * The percentage of reads that map to the mitochondrial genome
                                 + Low-quality / dying cells often exhibit extensive mitochondrial contamination
                                 + We calculate mitochondrial QC metrics with the `PercentageFeatureSet` function, which calculates the percentage of counts originating from a set of features
                                 + We use the set of all genes starting with `MT-` as a set of mitochondrial genes'),
                         column(6,numericInput('min.cells','Minimum number of cells per gene',step = 1,min = 1,value = 3)),
                         column(6,numericInput('min.genes','Minimum number of genes per cell',step = 1,min = 1,value = 200)),
                         
                         column(6,radioButtons('run_reduce_rb','Reduce Data by removing a marker',c(F,T))),
                         column(12,uiOutput('remove_select_gene_ui'),
                                uiOutput('remove_hist_ui'),
                                uiOutput('remove_threshold_ui'),
                                htmlOutput('text_remaining'),
                                tags$hr())
                         
                         
                         
                         
                         #actionButton('run_reduce','Save Data'),
                         
                         #column(6,actionButton('reload_data','Reload Data'))
                         #column(12,plotOutput('features_plot'))
                         ),
                
                # tabPanel('Upload Data',
                #          column(6,
                #                 
                #          ),
                #          column(6,
                #                 tags$h5('Upload a tar.gz file containing the outs folder generated by CellRanger, the name of the containing folder will be used to name the dataset (i.e immune_hsa)'),
                #                 imageOutput('file_stucture_image',height = 300, width = 300),
                #                 fileInput("file", "Upload Zip file"),
                #                 #shinyDirButton('directory', 'Folder select', 'Please select a folder'),
                #                 actionButton("unzip", "Unzip Files")
                #          )
                #          
                # ),
                #### DATA PROCESSING ####
                tabPanel('Data Processing',
                         tabsetPanel(selected = '',
                                     ### FILTERING AND NORMALISATION ####
                          
                                     #### GENE PLOT ####
                                     tabPanel("Gene Plot",
                                              
                                              
                                              plotOutput('features_plot'),
                                              plotOutput('scatter_plot'),
                                              tags$h5('option to subset by number of features could be added here')
                                              #plotOutput('gene_plot')
                                     ),
                                     #### FILTER ####
                                     tabPanel('Filter, Normalise and Find Variable Genes',
                                              tags$h4('### Normalizing the data'),
                                              
                                              tags$h5('After removing unwanted cells from the dataset, the next step is to normalize the data. By default, we employ a global-scaling normalization method "LogNormalize" that normalizes the feature expression measurements for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result. Normalized values are stored in `pbmc[["RNA"]]@data`.'),
                                              #tags$h5('Seurat calculates highly variable genes and focuses on these for downstream analysis. **`FindVariableGenes`** calculates the average expression and dispersion for each gene, places these genes into bins, and then calculates a z-score for dispersion within each bin. This helps control for the relationship between variability and average expression. This function is unchanged from (Macosko *et al*.), but new methods for variable gene expression identification are coming soon. We suggest that users set these parameters to mark visual outliers on the dispersion plot, but the exact parameter settings may vary based on the data type, heterogeneity in the sample, and normalization strategy. The parameters here identify ~2,000 variable genes, and represent typical parameter settings for UMI data that is normalized to a total of 1e4 molecules.'),
                                              textOutput('FilterCells_text'),
                                              textOutput('NormalizeData_text'),
                                              
                                              tags$h4('### Identification of highly variable features (feature selection)'),
                                              
                                              tags$h5('We next calculate a subset of features that exhibit high cell-to-cell variation in the dataset (i.e, they are highly expressed in some cells, and lowly expressed in others). We and [others](https://www.nature.com/articles/nmeth.2645) have found that focusing on these genes in downstream analysis helps to highlight biological signal in single-cell datasets.
                                                      
                                                      Our procedure in Seurat3 is described in detail [here](https://www.biorxiv.org/content/early/2018/11/02/460147.full.pdf), and improves on previous versions by directly modeling the mean-variance relationship inherent in single-cell data, and is implemented in the `FindVariableFeatures` function. By default, we return 2,000 features per dataset. These will be used in downstream analysis, like PCA.'),
                                              
                                              textOutput('FindVariableGenes_text'),
                                              plotOutput('var_plot')
                                              #))
                                              ),
                                     ### SCALING #####
                                     tabPanel('Scaling',
                                              
                                              tags$h4('### Scaling the data and removing unwanted sources of variation'),
                                              
                                              tags$h5("Your single cell dataset likely contains 'uninteresting' sources of variation. This could include not only technical noise, but batch effects, or even biological sources of variation (cell cycle stage). As suggested in [Buettner *et al*, NBT, 2015](https://www.nature.com/articles/nbt.3102), regressing these signals out of the analysis can improve downstream dimensionality reduction and clustering. To mitigate the effect of these signals, Seurat constructs linear models to predict gene expression based on user-defined variables. The scaled z-scored residuals of these models are stored in the scale.data slot, and  are used for dimensionality reduction and clustering."),
                                              
                                              tags$h5("We can regress out cell-cell variation in gene expression driven by batch (if applicable), cell alignment rate (as provided by Drop-seq tools for Drop-seq data),  the number of detected molecules, and mitochondrial gene expression. For cycling cells, we can also learn a 'cell-cycle' score (see example [here](http://satijalab.org/seurat/cell_cycle_vignette.html)) and regress this out as well. In this simple example here for post-mitotic blood cells, we regress on the number of detected molecules per cell as well as the percentage mitochondrial gene content. "), 
                                              
                                              tags$h5("Seurat v2.0 implements this regression as part of the data scaling process. Therefore, the `RegressOut` function has been deprecated, and replaced with the vars.to.regress argument in `ScaleData`."),
                                              
                                              actionButton('run_scaling','Run Scaling'),
                                              textOutput('ScaleData_text'),
                                              textOutput('scaling_done_text')
                                              
                                     ),
                                     
                                     #tabPanel('Remove Data'
                                     
                                     #         ),
                                     #### PCA ####
                                     tabPanel('PCA',
                                              tags$h4('### Perform linear dimensional reduction'),
                                              
                                              tags$h5("Next we perform PCA on the scaled data. By default, the genes in `object@var.genes` are used as input, but can be defined using pc.genes. We have typically found that running dimensionality reduction on highly variable genes can improve performance. However, with UMI data - particularly after regressing out technical variables, we often see that PCA returns similar (albeit slower) results when run on much larger subsets of genes, including the whole transcriptome."),
                                              #radioButtons('pca_full_dataset','Full Dataset',c(T,F)),
                                              actionButton('run_pca','Run PCA'),
                                              
                                              textOutput('pca_1'),
                                              #verbatimTextOutput('pca_text'),
                                              plotOutput('VizPCA'),
                                              plotOutput('VizPlot')
                                              
                                     ),
                                     tabPanel('PCA Heatmap',
                                              numericInput('pca_heatmap_dims','Dimmensions',2),
                                              plotOutput('DimHeatmap')
                                     ),
                                     # tabPanel('PCAHeatmap',
                                     #          tags$h5("In particular `PCHeatmap` allows for easy exploration of the primary sources of heterogeneity in a dataset, and can be useful when trying to decide which PCs to include for further downstream analyses. Both cells and genes are ordered according to their PCA scores. Setting cells.use to a number plots the 'extreme' cells on both ends of the spectrum, which dramatically speeds plotting for large datasets. Though clearly a supervised analysis, we find this to be a valuable tool for exploring correlated gene sets."),
                                     #          actionButton('ProjectPCA','Run ProjectPCA'),
                                     #          plotOutput('PCHeatmap'),
                                     #          plotOutput('PCHeatmap_2')
                                     #          
                                     # ),
                                     tabPanel('PCA Statistics',
                                              tags$h4("### Determine statistically significant principal components"),
                                              
                                              tags$h5("To overcome the extensive technical noise in any single gene for scRNA-seq data, Seurat clusters cells based on their PCA scores, with each PC essentially representing a 'metagene' that combines information across a correlated gene set. Determining how many PCs to include downstream is therefore an important step."),
                                              
                                              tags$h5("In [Macosko *et al*](http://www.cell.com/abstract/S0092-8674(15)00549-8), we implemented a resampling test inspired by the jackStraw procedure. We randomly permute a subset of the data (1% by default) and rerun PCA, constructing a 'null distribution' of gene scores, and repeat this procedure. We identify 'significant' PCs as those who have a strong enrichment of low p-value genes."),
                                              actionButton('Jacksaw', 'Run Jacksaw'),
                                              
                                              tags$h5("The `JackStrawPlot` function provides a visualization tool for comparing the distribution of p-values for each PC with a uniform distribution (dashed line). 'Significant' PCs will show a strong enrichment of genes with low p-values (solid curve above the dashed line). In this case it appears that PCs 1-10 are significant."),
                                              plotOutput('JackStrawPlot'),
                                              
                                              tags$h5("A more ad hoc method for determining which PCs to use is to look at a plot of the standard deviations of the principle components and draw your cutoff where there is a clear elbow in the graph. This can be done with `PCElbowPlot`. In this example, it looks like the elbow would fall around PC 9."),
                                              plotOutput('PCElbowPlot'),
                                              
                                              tags$h5("PC selection -- identifying the true dimensionality of a dataset -- is an important step for Seurat, but can be challenging/uncertain for the user. We therefore suggest these three approaches to consider. The first is more supervised, exploring PCs to determine relevant sources of heterogeneity, and could be used in conjunction with GSEA for example. The second implements a statistical test based on a random null model, but is time-consuming for large datasets, and may not return a clear PC cutoff. The third is a heuristic that is commonly used, and can be calculated instantly. In this example, all three approaches yielded similar results, but we might have been justified in choosing anything between PC 7-10 as a cutoff. We followed the jackStraw  here, admittedly buoyed by seeing the PCHeatmap returning interpretable signals (including canonical dendritic cell markers) throughout these PCs. Though the results are only subtly affected by small shifts in this cutoff (you can test below), we strongly suggest always explore the PCs they choose to include downstream.")
                                              
                                     ),
                                     ##### CLUSTER #####
                                     tabPanel('Cluster Cells',
                                              tags$h4("### Cluster the cells"),
                                              
                                              tags$h5("Seurat now includes an graph-based clustering approach compared to (Macosko *et al*.). Importantly, the *distance metric* which drives the clustering analysis (based on previously identified PCs) remains the same. However, our approach to partioning the cellular distance matrix into clusters has dramatically improved. Our approach was heavily inspired by recent manuscripts which applied graph-based clustering approaches to scRNA-seq data [[SNN-Cliq, Xu and Su, Bioinformatics, 2015]](http://bioinformatics.oxfordjournals.org/content/early/2015/02/10/bioinformatics.btv088.abstract) and CyTOF data [[PhenoGraph, Levine *et al*., Cell, 2015]](http://www.ncbi.nlm.nih.gov/pubmed/26095251). Briefly, these methods embed cells in a graph structure - for example a K-nearest neighbor (KNN) graph, with edges drawn between cells with similar gene expression patterns, and then attempt to partition this graph into highly interconnected 'quasi-cliques' or 'communities'. As in PhenoGraph, we first construct a KNN graph based on the euclidean distance in PCA space, and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity). To cluster the cells, we apply modularity optimization techniques such as the Louvain algorithm (default) or SLM [[SLM, Blondel *et al*., Journal of Statistical Mechanics]](http://dx.doi.org/10.1088/1742-5468/2008/10/P10008), to iteratively group cells together, with the goal of optimizing the standard modularity function."),
                                              
                                              tags$h5("The `FindClusters` function implements the procedure, and contains a resolution parameter that sets the 'granularity' of the downstream clustering, with increased values leading to a greater number of clusters. We find that setting this parameter between 0.6-1.2 typically returns good results for single cell datasets of around 3K cells. Optimal resolution often increases for larger datasets. The clusters are saved in the `object@ident` slot."),
                                              actionButton('find_clusters','Run FindClusters')
                                     ),
                                     tabPanel('uMap',
                                              actionButton('run_umap','Run uMap'),
                                              plotOutput('umap_plot')
                                     )
                                     # ,
                                     # tabPanel('tSNE',
                                     #          tags$h4("### Run Non-linear dimensional reduction (tSNE)"),
                                     #          
                                     #          tags$h5("Seurat continues to use tSNE as a powerful tool to visualize and explore these datasets. While we no longer advise clustering directly on tSNE components, cells within the graph-based clusters determined above should co-localize on the tSNE plot. This is because the tSNE aims to place cells with similar local neighborhoods in high-dimensional space together in low-dimensional space. As input to the tSNE, we suggest using the same PCs as input to the clustering analysis, although computing the tSNE based on scaled gene expression is also supported using the genes.use argument."),
                                     #          actionButton('run_tsne','Run TSNE'),
                                     #          plotOutput('TSNEPlot')
                                     # )
                                     )
                                              ),
                #### Visualisation #####
                tabPanel('Data Visualisation',
                         tabsetPanel(selected = '',

                           tabPanel('Seurat Plots',
                                    
                                    textOutput('values_text'),
                                    tabsetPanel(selected = 'Assign Cells',
                                                tabPanel('Differential Genes',
                                                         tags$h4("### Finding differentially expressed genes (cluster biomarkers)"),
                                                         
                                                         tags$h5("Seurat can help you find markers that define clusters via differential expression. By default, it identifes positive and negative markers of a single cluster (specified in ident.1), compared to all other cells. FindAllMarkers automates this process for all clusters, but you can also test groups of clusters vs. each other, or against all cells."),
                                                         
                                                         tags$h5("The min.pct argument requires a feature to be detected at a minimum percentage in either of the two groups of cells, and the thresh.test argument requires a feature to be differentially expressed (on average) by some amount between the two groups. You can set both of these to 0, but with a dramatic increase in time - since this will test a large number of features that are unlikely to be highly discriminatory. As another option to speed up these computations, max.cells.per.ident can be set. This will downsample each identity class to have no more cells than whatever this is set to. While there is generally going to be a loss in power, the speed increases can be significiant and the most highly differentially expressed features will likely still rise to the top."),
                                                         tags$h6('Additional options available here, see test.use = "roc"'),
                                                         
                                                         column(12,selectInput('find_single_cluster_test',"Test to Use",FindMarkers_test_list,selected = "Wilcoxon rank sum test (default)")),
                                                         column(4,numericInput('logfc.threshold','logfc.threshold',value = 0.25)),
                                                         column(4,numericInput('min.pct','min.pct',value = 0.1)),
                                                         column(4),
                                                         column(12,
                                                                tabsetPanel(
                                                                  tabPanel('Top markers for All Clusters',
                                                                           
                                                                           actionButton('find_all_clusters','Run All Cluster'),
                                                                           tags$h4('This can take up to half and hour to process'),
                                                                           numericInput('clust_num_display','Number of Markers to Display',value = 2),
                                                                           dataTableOutput('data_markers_table')
                                                                  ),
                                                                  tabPanel('Single Cluster',
                                                                           column(2,numericInput('single_cluster_number','Cluster Number',min = 0, max = 50,step = 1,value = 1)),
                                                                           column(10),
                                                                           column(12,
                                                                                  tabsetPanel(
                                                                                    tabPanel('Single',
                                                                                             column(3,actionButton('find_single_cluster','Run Single Cluster FindMarker')),
                                                                                             
                                                                                             column(12,dataTableOutput('single_markers_table'))
                                                                                             
                                                                                    ),
                                                                                    tabPanel('Gene Unique to Cluster - removed in Seurat 3',
                                                                                             actionButton('find_distinguishing_cluster','Run Distinguish Cluster '),
                                                                                             selectInput('find_distinguishing_cluster_select','Select Clusters',c(1:50),selected = c(2,3),multiple = T),
                                                                                             dataTableOutput('distinguish_markers_table')
                                                                                    )
                                                                                    
                                                                                  )))
                                                                ))
                                                         
                                                         
                                                ),
                                                tabPanel('Visualise Marker Expression',
                                                         tags$h5("We include several tools for visualizing marker expression. VlnPlot (shows expression probability distributions across clusters), and FeaturePlot (visualizes feature expression on a tSNE or PCA plot) are our most commonly used visualizations. We also suggest exploring RidgePlot, CellScatter, and DotPlot as additional methods to view your dataset."),
                                                         
                                                         
                                                         uiOutput('plot_select_genes_ui'),
                                                         tags$h6('VlnPlot(object = pbmc, features = c("MS4A1", "CD79A"))'),
                                                         plotOutput('vlnplot_1'),
                                                         tags$h6('VlnPlot(object = pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)'),
                                                         plotOutput('vlnplot_2'),
                                                         tags$h6('FeaturePlot(object = pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))'),
                                                         plotOutput('feature_plot_1'),
                                                         tags$h6('DoHeatmap(object = data, features = input$plot_select_genes) + NoLegend()'),
                                                         plotOutput('doheatmap_select')
                                                         
                                                         
                                                         
                                                         
                                                         # 
                                                         # ```{r markerplots, fig.height=8, fig.width=15,}
                                                         # VlnPlot(object = pbmc, features.plot = c("MS4A1", "CD79A"))
                                                         # 
                                                         # # you can plot raw UMI counts as well
                                                         # VlnPlot(object = pbmc, features.plot = c("NKG7", "PF4"), use.raw = TRUE, y.log = TRUE)
                                                         # 
                                                         # FeaturePlot(object = pbmc, features.plot = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"), cols.use = c("grey", "blue"), reduction.use = "tsne")
                                                ),
                                                tabPanel('DoHeatmap',
                                                         tags$h5("DoHeatmap generates an expression heatmap for given cells and features. In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster."),
                                                         tags$h6('pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) -> top10
                                                                 DoHeatmap(object = pbmc, features = top10$gene) + NoLegend()'),
                                                         # ```{r clusterHeatmap, fig.height=8, fig.width=15, message=FALSE, warning=FALSE}
                                                         # pbmc.markers %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10
                                                         # # setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name
                                                         # DoHeatmap(object = pbmc, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
                                                         # ```
                                                         plotOutput('doheatmap')
                                                         
                                                         ),
                                                tabPanel('Assign Cells',
                                                         tags$h4("### Assigning cell type identity to clusters"),
                                                         
                                                         #tags$h5("Fortunately in the case of this dataset, we can use canonical markers to easily match the unbiased clustering to known cell types:"),
                                                         column(4,uiOutput('cluster_0_text')),
                                                         column(4,uiOutput('cluster_1_text')),
                                                         column(4,uiOutput('cluster_2_text')),
                                                         column(4,uiOutput('cluster_3_text')),
                                                         column(4,uiOutput('cluster_4_text')),
                                                         column(4,uiOutput('cluster_5_text')),
                                                         column(4,uiOutput('cluster_6_text')),
                                                         column(4,uiOutput('cluster_7_text')),
                                                         column(4,uiOutput('cluster_8_text')),
                                                         column(4,uiOutput('cluster_9_text')),
                                                         column(4,uiOutput('cluster_10_text')),
                                                         column(4,uiOutput('cluster_11_text')),
                                                         column(4,uiOutput('cluster_12_text')),
                                                         column(4,uiOutput('cluster_13_text')),
                                                         
                                                         
                                                         #   Cluster ID | Markers       | Cell Type
                                                         # -----------|---------------|----------
                                                         #   0          | IL7R          | CD4 T cells
                                                         # 1          | CD14, LYZ     | CD14+ Monocytes
                                                         # 2          | MS4A1         | B cells
                                                         # 3          | CD8A          | CD8 T cells
                                                         # 4          | FCGR3A, MS4A7 | FCGR3A+ Monocytes
                                                         # 5          | GNLY, NKG7    | NK cells
                                                         # 6          | FCER1A, CST3  | Dendritic Cells
                                                         # 7          | PPBP          | Megakaryocytes
                                                         # 
                                                         # 
                                                         # ```{r labelplot, fig.height=5, fig.width=9, warning = FALSE}
                                                         # current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)
                                                         # new.cluster.ids <- c("CD4 T cells", "CD14+ Monocytes", "B cells", "CD8 T cells", "FCGR3A+ Monocytes", "NK cells", "Dendritic cells", "Megakaryocytes")
                                                         # pbmc@ident <- plyr::mapvalues(x = pbmc@ident, from = current.cluster.ids, to = new.cluster.ids)
                                                         # TSNEPlot(object = pbmc, do.label = TRUE, pt.size = 0.5)
                                                         # ```
                                                         # 
                                                         column(12,plotOutput('assign_clusters'))
                                                         
                                                         #tags$h5("### Further subdivisions within cell types"),
                                                         # 
                                                         #tags$h4("If you perturb some of our parameter choices above (for example, setting  `resolution=0.8` or changing the number of PCs), you might see the CD4 T cells subdivide into two groups. You can explore this subdivision to find markers separating the two T cell subsets. However, before reclustering (which will overwrite `object@ident`), we can stash our renamed identities to be easily recovered later.")
                                                         # 
                                                         # ```{r investigate_split, fig.width=15}
                                                         # # First lets stash our identities for later
                                                         # pbmc <- StashIdent(object = pbmc, save.name = "ClusterNames_0.6")
                                                         # 
                                                         # # Note that if you set save.snn=T above, you don't need to recalculate the SNN, and can simply put: 
                                                         # # pbmc <- FindClusters(pbmc,resolution = 0.8)
                                                         # pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:10, resolution = 0.8, print.output = FALSE)
                                                         # 
                                                         # # Demonstration of how to plot two tSNE plots side by side, and how to color points based on different criteria
                                                         # plot1 <- TSNEPlot(object = pbmc, do.return = TRUE, no.legend = TRUE, do.label = TRUE)
                                                         # plot2 <- TSNEPlot(object = pbmc, do.return = TRUE, group.by = "ClusterNames_0.6", no.legend = TRUE, do.label = TRUE)
                                                         # plot_grid(plot1, plot2)
                                                         # 
                                                         # # Find discriminating markers
                                                         # tcell.markers <- FindMarkers(object = pbmc, ident.1 = 0, ident.2 = 1)
                                                         # 
                                                         # # Most of the markers tend to be expressed in C1 (i.e. S100A4). However, we can see that CCR7 is upregulated in 
                                                         # # C0, strongly indicating that we can differentiate memory from naive CD4 cells.
                                                         # # cols.use demarcates the color palette from low to high expression
                                                         # FeaturePlot(object = pbmc, features.plot = c("S100A4", "CCR7"), cols.use = c("green", "blue"))
                                                         # ```
                                                         # 
                                                         # The memory/naive split is bit weak, and we would probably benefit from looking at more cells to see if this becomes more convincing. In the meantime, we can restore our old cluster identities for downstream processing.
                                                         # 
                                                         # ```{r restore}
                                                         # pbmc <- SetAllIdent(object = pbmc, id = "ClusterNames_0.6")
                                                         # saveRDS(pbmc, file = paste0(output_dir,output_name'_Seurat_final_output.rds'))
                                                         # ```
                                                         # 
                                                )
                                                
                                    )),
                           #### Features #####
                           tabPanel('Features',
                                    tabsetPanel(selected = '',
                                        
                                      tabPanel('Individual Datasets',
                                               
                                               column(12,
                                                      column(3,uiOutput('original_data_select_features_ui')),
                                                      column(3,numericInput('data_density_y','Density ylim',0)),
                                                      column(3,numericInput('vis_data_value_min','Min Value Cutoff',0)),
                                                      column(3,numericInput('vis_feature_threshold','Threshold',0))
                                               ),
                                               column(12,
                                                      column(6,plotOutput('original_data_boxplot')),
                                                      
                                                      column(6,plotOutput('vis_feature_histogram'))
                                                      
                                               ),
                                               column(12,plotOutput('original_data_density')),
                                               column(12,dataTableOutput('vis_feature_count')),
                                               column(12,dataTableOutput('original_data_table'))
                                      ),
                                      tabPanel('Combined Features',
                                               column(12,
                                                      column(3,uiOutput('select_combined_group_ui')),
                                                      column(3,uiOutput('combined_feature_value_min_ui'))
                                                      ),
                                               column(12,tags$hr()),
                                               column(6,
                                                      column(6,uiOutput('select_combined_feature_1_ui')),
                                                      column(6,uiOutput('combined_feature_1_threshold_ui')),
                                                      column(12,
                                                             plotOutput('combined_feature_1_boxplot'),
                                                             plotOutput('combined_feature_1_hist'),
                                                             dataTableOutput('combined_feature_1_count_df'))
                                                      ),
                                               column(6,
                                                      column(6,uiOutput('select_combined_feature_2_ui')),
                                                      column(6,uiOutput('combined_feature_2_threshold_ui')),
                                                      column(12,
                                                             plotOutput('combined_feature_2_boxplot'),
                                                             plotOutput('combined_feature_2_hist'),
                                                             dataTableOutput('combined_feature_2_count_df')
                                                             )
                                                      )
                                               
                                               
                                               )
                                      
                                      ))
                         ))
                
                                     )
    
                )
  
    ))
