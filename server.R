#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {


  output$debug_ui = renderUI({
    #sysinf <- Sys.info()
    if (!is.null(sysinf)){
      #os <- sysinf['sysname']
      if (os == 'Darwin'){
        actionButton('debug','Debug')
      }
    }
  })
  
  observeEvent(input$debug,{
    browser()
  })
  
  values = reactiveValues()
  #values_save = readRDS('www/values_save.rds')

  read_rds = T
  if(read_rds == T){
    values_save = readRDS('www/values_save.rds')
  }else{
    values_save = reactiveValues()
    saveRDS(values_save, 'www/values_save.rds')
  }
  output$features_plot = renderPlot({
    print('Create Seurat Object')
    if(read_rds == T | input$re_load == T){
      withProgress(message = 'readRDS',{
        Seurat_data = readRDS('www/Seurat_data.rds')
        #values_save = readRDS('www/values_save.rds')
        
      })
    }else{
      print('create')
      #values_save = reactiveValues()
      #saveRDS(values_save, 'www/values_save.rds')
      Seurat_data = CreateSeuratObject(raw.data = data, min.cells = input$min.cells, min.genes = input$min.genes, project = "immune_hsa")
    }
    #Seurat_data = CreateSeuratObject(raw.data = data, min.cells = 3, min.genes = 200, project = "immune_hsa")

    mito.genes <- grep(pattern = "^MT-", x = rownames(x = Seurat_data@data), value = TRUE)
    percent.mito <- Matrix::colSums(Seurat_data@raw.data[mito.genes, ]) / Matrix::colSums(Seurat_data@raw.data)

    # AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats
    Seurat_data <- AddMetaData(object = Seurat_data, metadata = percent.mito, col.name = "percent.mito")
      print('feature_plot')
      values$Seurat_data = Seurat_data

      VlnPlot(object = Seurat_data, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
    })
  #
  #
    output$gene_plot = renderPlot({
      print('gene_plot')
      par(mfrow = c(1, 2))
      
      GenePlot(object = values$Seurat_data, gene1 = "nUMI", gene2 = "percent.mito")
      GenePlot(object = values$Seurat_data, gene1 = "nUMI", gene2 = "nGene")
    })
  #
  #
    output$var_plot = renderPlot({
      Seurat_data = values$Seurat_data
      Seurat_data <- FilterCells(object = Seurat_data, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))


      Seurat_data <- NormalizeData(object = Seurat_data, normalization.method = "LogNormalize", scale.factor = 1e4)

      Seurat_data <- FindVariableGenes(object = Seurat_data, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
      print('var_plot')
      if (os == 'Darwin'){
        saveRDS(Seurat_data,'www/Seurat_data.rds')
      }
    
      values$Seurat_data_norm <- FindVariableGenes(object = Seurat_data, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

    })

  # output$features_plot = renderPlot({
  #   print('feature_plot')
  #
  #   VlnPlot(object = Seurat_data(), features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
  # })
  # 
  # Seurat_data = reactive({
  #   Seurat_data = CreateSeuratObject(raw.data = data, min.cells = input$min.cells, min.genes = input$min.genes, project = "immune_hsa")
  # 
  #   mito.genes <- grep(pattern = "^MT-", x = rownames(x = Seurat_data@data), value = TRUE)
  #   percent.mito <- Matrix::colSums(Seurat_data@raw.data[mito.genes, ]) / Matrix::colSums(Seurat_data@raw.data)
  # 
  #   # AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats
  #   Seurat_data <- AddMetaData(object = Seurat_data, metadata = percent.mito, col.name = "percent.mito")
  # 
  # 
  # 
  #   output$gene_plot = renderPlot({
  #     print('gene_plot')
  #     par(mfrow = c(1, 2))
  #     GenePlot(object = Seurat_data, gene1 = "nUMI", gene2 = "percent.mito")
  #     GenePlot(object = Seurat_data, gene1 = "nUMI", gene2 = "nGene")
  #   })
  # 
  # 
  # 
  #     Seurat_data <- FilterCells(object = Seurat_data, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))
  # 
  # 
  #     Seurat_data <- NormalizeData(object = Seurat_data, normalization.method = "LogNormalize", scale.factor = 1e4)
  # 
  #     Seurat_data <- FindVariableGenes(object = Seurat_data, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
  #   output$var_plot = renderPlot({
  #     print('var_plot')
  # 
  #     Seurat_data <- FindVariableGenes(object = Seurat_data, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
  # 
  #   })
  # 
  #   #svalues$Seurat_data = Seurat_data
  #   Seurat_data
  # })
  # 
  # output$test_test = renderText({
  #   Seurat_data()
  #   print('test')
  # })
  # 
  # 
  # 
  # 
  observeEvent(input$run_scaling,{
    withProgress(message = 'running scaling',{
    Seurat_data <- ScaleData(object = values$Seurat_data_norm, vars.to.regress = c("nUMI", "percent.mito"))
    
    if (os == 'Darwin'){
      saveRDS(Seurat_data,'www/Seurat_data.rds')
    }
    
    values$Seurat_data_scale = Seurat_data
    })
  })
  # 
  # 
  output$scaling_done_text = renderText({
    if(!is.null(values$Seurat_data_scale)){
      print('done')
    }
  })
  #output$pca_1 = renderText({
  #  print('pca_1')
  observeEvent(input$run_pca,{
    print('pca_1')
      if(!is.null(values$Seurat_data_scale)){
        data = values$Seurat_data_scale
      }else{
        withProgress(message = 'readRDS',{
          data = readRDS('www/Seurat_data.rds')
        })
      }
    withProgress(message = 'RunPCA',{
      
      Seurat_data = RunPCA(object = data, pc.genes = data@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
    })
      values$Seurat_data_pca_1 <- Seurat_data
      withProgress(message = 'saveRDS',{
        
        if (os == 'Darwin'){
          saveRDS(Seurat_data,'www/Seurat_data.rds')
        }
      })
      
      print('pca')
      
      #PrintPCA(object = Seurat_data, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
      #VizPCA(object = Seurat_data, pcs.use = 1:2)
      #PCAPlot(object = Seurat_data, dim.1 = 1, dim.2 = 2)
      
      
    #})
  })
  
  output$pca_text = renderText({
    print('pca_text')
    if(!is.null(values$Seurat_data_pca_1)){
      PrintPCA(object = values$Seurat_data_pca_1, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
    }
  })
  
  output$VizPCA = renderPlot({
    print('VizPCA')
    if(!is.null(values$Seurat_data_pca_1)){
      
      VizPCA(object = values$Seurat_data_pca_1, pcs.use = 1:2)
    }
  })
  
  output$VizPlot = renderPlot({
    print('VizPlot')
    if(!is.null(values$Seurat_data_pca_1)){
      
      PCAPlot(object = values$Seurat_data_pca_1, dim.1 = 1, dim.2 = 2)
    }
  })
  
  observeEvent(input$ProjectPCA,{
    if(is.null(values$Seurat_data_pca_1)){
      withProgress(message = 'readRDS',{
        
        data = readRDS('www/Seurat_data.rds')
      })
    }else{
      data = values$Seurat_data_pca_1
    }
    # ProjectPCA scores each gene in the dataset (including genes not included in the PCA) based on their correlation 
    # with the calculated components. Though we don't use this further here, it can be used to identify markers that 
    # are strongly correlated with cellular heterogeneity, but may not have passed through variable gene selection. 
    # The results of the projected PCA can be explored by setting use.full=T in the functions above
    withProgress(message = 'ProjectPCA',{
      
      Seurat_data <- ProjectPCA(object = data, do.print = FALSE)
    })
    values$Seurat_data_ProjectPCA = Seurat_data
    withProgress(message = 'saveRDS',{
      if (os == 'Darwin'){
        saveRDS(Seurat_data,'www/Seurat_data.rds')
      }
    })
      
      
      
  })
  
  output$PCHeatmap = renderPlot({
    if(!is.null(values$Seurat_data_ProjectPCA)){
      PCHeatmap(object = values$Seurat_data_ProjectPCA, pc.use = 1, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
    }
  })
  output$PCHeatmap_2 = renderPlot({
    if(!is.null(values$Seurat_data_ProjectPCA)){
      PCHeatmap(object = values$Seurat_data_ProjectPCA, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
    }
  })
  
  observeEvent(input$Jacksaw,{
    print('JackStraw')
    if(is.null(values$Seurat_data_ProjectPCA)){
      withProgress(message = 'readRDS',{
        
        data = readRDS('www/Seurat_data.rds')
      })
    }else{
      data = values$Seurat_data_pca_1
    }
    withProgress(message = 'JackStraw',{
      Seurat_data <- JackStraw(object = data, num.replicate = 100, display.progress = FALSE)
    })
    values$Seurat_data_jack = Seurat_data
    withProgress(message = 'saveRDS',{
      if (os == 'Darwin'){
        saveRDS(Seurat_data,'www/Seurat_data.rds')
      }
    })
  })
    output$JackStrawPlot = renderPlot({
      if(!is.null(values$Seurat_data_jack)){
        JackStrawPlot(object = values$Seurat_data_jack, PCs = 1:12)
      }
    })
    
    output$PCElbowPlot = renderPlot({
      if(!is.null(values$Seurat_data_jack)){
        
        PCElbowPlot(object = values$Seurat_data_jack)
      }
    })
    
    observeEvent(input$find_clusters,{
      print('FindClusters')
      if(is.null(values$Seurat_data_jack)){
        withProgress(message = 'readRDS',{
          
          data = readRDS('www/Seurat_data.rds')
        })
      }else{
        data = values$Seurat_data_jack
      }
      # save.SNN = T saves the SNN so that the clustering algorithm can be rerun using the same graph
      # but with a different resolution value (see docs for full details)
      withProgress(message = 'FindClusters',{
        
        Seurat_data <- FindClusters(object = data, reduction.type = "pca", dims.use = 1:10, resolution = 0.6, print.output = 0, save.SNN = TRUE)
      })
      values$Seurat_data_clusters = Seurat_data
      if (os == 'Darwin'){
        saveRDS(Seurat_data,'www/Seurat_data.rds')
      }
    })
    
    output$cluster_text = renderText({
      if(!is.null(values$Seurat_data_clusters)){
        PrintFindClustersParams(object = values$Seurat_data_clusters)
      }
    })
    
    observeEvent(input$run_tsne,{
      print('TSNE')
      if(is.null(values$Seurat_data_clusters)){
        withProgress(message = 'readRDS',{
          data = readRDS('www/Seurat_data.rds')
        })
      }else{
        data = values$Seurat_data_clusters
      }
      withProgress(message = 'RunTSNE',{
        Seurat_data <- RunTSNE(object = data, dims.use = 1:10)
      })
      values$Seurat_data_tsne = Seurat_data
      if (os == 'Darwin'){
        saveRDS(Seurat_data,'www/Seurat_data.rds')
      }
    })
    
    output$TSNEPlot = renderPlot({
      if(!is.null(values$Seurat_data_tsne)){
        TSNEPlot(object = values$Seurat_data_tsne)
      }
    })
    ##### Visualisations ########
    
    plot_data = reactive({
      withProgress(message = 'readRDS',{
        data = readRDS('www/Seurat_data.rds')
      })
      
    })
    output$FindMarkers_text = renderPlot({
      data = values$Seurat_data_clusters
      
      cluster1.markers <- FindMarkers(object = data, ident.1 = 1, min.pct = 0.25)
      print(x = head(x = cluster1.markers, n = 5))
      str(cluster1.markers)
      # find all markers distinguishing cluster 5 from clusters 0 and 3
      cluster5.markers <- FindMarkers(object = data, ident.1 = 5, ident.2 = c(0,3), min.pct = 0.25)
      str(cluster5.markers)
      print(x = head(x = cluster5.markers, n = 5))
      # find markers for every cluster compared to all remaining cells, report only the positive ones
      data.markers <- FindAllMarkers(object = data, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
      data.markers = pbmc.markers
      data.markers.top_n = data.markers %>% group_by(cluster) %>% top_n(2, avg_logFC)
      
      str(data.markers.top_n)
      
    })
    observeEvent(input$find_single_cluster,{
      data = plot_data()
      withProgress(message = 'FindMarkers',{
        values_save$cluster.marker <- FindMarkers(object = data, 
                                                  ident.1 = input$single_cluster_number, 
                                                  test.use = input$find_single_cluster_test, 
                                                  min.pct = input$min.pct, 
                                                  logfc.threshold = input$logfc.threshold)
        values_save$cluster.marker$gene = rownames(values_save$cluster.marker)
        values_save$cluster.marker$cluster = input$single_cluster_number
      })
      withProgress(message = 'saveRDS values save',{
        if (os == 'Darwin'){
          saveRDS(values_save,'www/values_save.rds')
        }
        
      })
    })
    
    observeEvent(input$find_single_cluster_roc,{
      withProgress(message = 'readRDS',{
        data = readRDS('www/Seurat_data.rds')
      })
      withProgress(message = 'FindMarkers',{
        values_save$cluster.marker.roc <- FindMarkers(object = data, ident.1 = 0, thresh.use = 0.25, test.use = "roc", only.pos = TRUE)
        values_save$cluster.marker.roc$gene = rownames(values$cluster.marker.roc)
        if (os == 'Darwin'){
          saveRDS(values_save,'www/values_save.rds')
        }
        
      })
    })
    
    observeEvent(input$find_distinguishing_cluster,{
      withProgress(message = 'readRDS',{
        data = readRDS('www/Seurat_data.rds')
      })
      withProgress(message = 'FindMarkers',{
        # find all markers distinguishing cluster 5 from clusters 0 and 3
        values_save$cluster_markers_distinguishing <- FindMarkers(object = data, 
                                                                  ident.1 = input$single_cluster_number, 
                                                                  ident.2 = as.numeric(input$find_distinguishing_cluster_select), 
                                                                  min.pct = input$min.pct, 
                                                                  logfc.threshold = input$logfc.threshold)
        values_save$cluster_markers_distinguishing$gene = rownames(values_save$cluster_markers_distinguishing)
        if (os == 'Darwin'){
          saveRDS(values_save,'www/values_save.rds')
        }
        
      })
    })
      
    
    observeEvent(input$find_all_clusters,{
      withProgress(message = 'readRDS',{
        data = readRDS('www/Seurat_data.rds')
      })
      withProgress(message = 'FindAllMarkers',{
        values_save$data.markers <- FindAllMarkers(object = data, 
                                                   only.pos = TRUE, 
                                                   min.pct = input$min.pct, 
                                                   logfc.threshold = input$logfc.threshold)
        if (os == 'Darwin'){
          saveRDS(values_save,'www/values_save.rds')
        }
      })

    })
  
    output$single_markers_table = renderDataTable({
      if(!is.null(values_save$cluster.marker)){
        values_save$cluster.marker
      }
      
    })
    
    output$single_markers_table_roc = renderDataTable({
      if(!is.null(values_save$cluster.marker)){
        values_save$cluster.marker.roc
      }
      
    })
    
    output$distinguish_markers_table = renderDataTable({
      if(!is.null(values_save$cluster_markers_distinguishing)){
        values_save$cluster_markers_distinguishing
      
        
      }
      
    })
    
    output$data_markers_table = renderDataTable({
      if(!is.null(values_save$data.markers)){
        data.markers.top_n = values_save$data.markers %>% group_by(cluster) %>% top_n(input$clust_num_display, avg_logFC)
      }
      
    })

    
    output$plot_select_genes_ui = renderUI({
      gene_list = values_save$data.markers$gene[1:3]
      gene_list = c("MS4A1", "IGLL1", "JCHAIN","HGB", "HLA-DQA1",
                    "XCL1",'XCL2',
                    "FSCN1",'CCL17','CCl27',
                    'IL1A','CCL7','CXCL5',
                    'LRRN2','CD88','GIMAP4')
      selectInput('plot_select_genes','Select Genes',values_save$data.markers$gene,gene_list,multiple = T)
      
    })
    output$vlnplot_1 = renderPlot({
      data = plot_data()
      VlnPlot(object = data, features.plot = input$plot_select_genes)
    }) 
    output$vlnplot_2 = renderPlot({
      data = plot_data()
      # you can plot raw UMI counts as well
      VlnPlot(object = data, features.plot = input$plot_select_genes, use.raw = TRUE, y.log = TRUE)
    }) 
    output$feature_plot_1 = renderPlot({
      #data = readRDS('www/Seurat_data.rds')
      data = plot_data()
      
      FeaturePlot(object = data, features.plot = input$plot_select_genes, cols.use = c("grey", "blue"), reduction.use = "tsne")
    })
    
    output$doheatmap = renderPlot({
      # ```{r clusterHeatmap, fig.height=8, fig.width=15, message=FALSE, warning=FALSE}
      values_save$data.markers %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10
      # # setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name
      DoHeatmap(object = plot_data(), genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
      # ```
    })
    
    output$assign_clusters = renderPlot({
      data = plot_data()
      #current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)
      #new.cluster.ids <- c("CD4 T cells", "CD14+ Monocytes", "B cells", "CD8 T cells", "FCGR3A+ Monocytes", "NK cells", "Dendritic cells", "Megakaryocytes")
      #data@ident <- plyr::mapvalues(x = data@ident, from = current.cluster.ids, to = new.cluster.ids)
      TSNEPlot(object = data, do.label = TRUE, pt.size = 0.5)
    })

})
    

  
