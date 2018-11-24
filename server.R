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
  
  
  output$select_dataset_ui = renderUI({
    datasets = values$data_list
    selectInput('select_dataset','Select Dataset',datasets,datasets[1])
  })
  observeEvent(input$debug,{
    browser()
  })
  
  observeEvent(input$unzip,{
    withProgress(message = 'untar',{
      untar(input$file$datapath, exdir="./data/")
      values$data_list = list.files('data/')
    })
    #unzip(input$file$datapath, list = TRUE, exdir = 'data/')
  })
  
  values = reactiveValues(data_list = data_list)
  values_save = reactiveValues()
  
  values_path = reactive({
    if(!is.null(input$select_dataset)){
      if(!is.null(input$reduce_select)){
        if(input$reduce_select !=  '_'){
          values_path = paste0('www/values_',input$select_dataset,'.rds',input$reduce_select)
        }else{
          values_path = paste0('www/values_',input$select_dataset,'.rds')
        }
        #values_path = paste0('www/values_',input$select_dataset,'.rds')
        if(input$run_reduce_rb == T){
          
          if(!is.null(input$remove_select_gene)){
            if(input$remove_select_gene != '_'){
              values_path = paste0(values_path,'_remove_',input$remove_select_gene)
            }
          }
        }
        values_path
      }
    }
  })
  #### REDUCE ####
  data_df = reactive({
    if(input$run_reduce_rb == T){
      withProgress(message = 'Read Data',{
      print(data_path())
      data <- Read10X(data.dir = data_path())
      df = as.data.frame(as.matrix(data))
      dim(df)
      df = as.data.frame(as.matrix(data))
      df
      })
    }
  })
  

  output$remove_select_gene_ui = renderUI({
    if(input$run_reduce_rb == T){
      
      selectInput('remove_select_gene','Select Marker to Remove',c('_',rownames(data_df())),'FSCN1')
    }
  })
  
  output$remove_hist = renderPlot({
    if(input$run_reduce_rb == T){
      
      if(!is.null(input$remove_select_gene)){
        withProgress(message = 'histogram',{
          df = data_df()[input$remove_select_gene,]
          #str(df)
          df_l = as.data.frame(t(df))
          ggplot(df_l) +
            geom_histogram(aes_string(input$remove_select_gene))
        })
      }
        
    }
  })
  output$remove_threshold_ui = renderUI({
    if(input$run_reduce_rb == T){
      
      numericInput('remove_threshold','Select Expression Threshold',1)
    }
  })
  
  data_reduce = reactive({
    if(input$run_reduce_rb == T){
      if(!is.null(input$remove_select_gene)){
        withProgress(message = 'data reduce',{
          df = data_df()
          df_r = df[,as.numeric(df[input$remove_select_gene,]) < input$remove_threshold]
          #dim(df_r)
          
          df_r
        })
      }
    }
  })
  
  
  
  output$text_remaining = renderText({
    if(input$run_reduce_rb == T){
      
      print(paste('Total Number of cells =',dim(data_df())[2],'<br>Remaining number of cells = ',dim(data_reduce())[2]))
    }
  })
  
  observeEvent(input$reload_data,{
    withProgress(message = 'Re-CreateSeuratObject',{
      print(data_path())
      if(input$run_reduce_rb == T){
        if(input$remove_select_gene == '_'){
          data <- Read10X(data.dir = data_path())
        }else{
          data = as(as.matrix(data_reduce()),'dgCMatrix')
        }
      }else{
        data <- Read10X(data.dir = data_path())
        
      }
      dim(data)
      Seurat_data = CreateSeuratObject(raw.data = data, min.cells = input$min.cells, min.genes = input$min.genes, project = input$select_dataset)
      values$Seurat_data = Seurat_data
      values_save = reactiveValues()
      
      
    })
    withProgress(message = paste('saveRDS',dataset_path()),{
      
      saveRDS(Seurat_data,dataset_path())
      
      saveRDS(values_save,values_path())
    })
  })
  
  output$reduce_dataset_select_ui = renderUI({
    if(!is.null(input$select_dataset)){
      rds_list = list.files('www')
      reduce_list = grep(input$select_dataset,rds_list,value = T)
      #if(length(reduce_list) > 1){
      reduce_list_sub = gsub(paste0(input$select_dataset,'.rds'),'',reduce_list)
      values_list = grep('values_',reduce_list_sub,value = T)
      reduce_list_sub = reduce_list_sub[!(reduce_list_sub %in% values_list)]
      reduce_list_sub
      selectInput('reduce_select','Select sub datasets',c('_',reduce_list_sub,'_'))
      # }
    }
  })
  
  dataset_path = reactive({
    if(!is.null(input$select_dataset)){
      if(!is.null(input$reduce_select)){
        
        dataset_rds_name = 'Seurat_data.rds'
        if(input$reduce_select !=  '_'){
          dataset_rds_name = paste0(input$select_dataset,'.rds',input$reduce_select)
        }else{
          dataset_rds_name = paste0(input$select_dataset,'.rds')
        }
        dataset_path = paste0('www/',dataset_rds_name)
        dataset_path
        if(input$run_reduce_rb == T){
          
          if(!is.null(input$remove_select_gene)){
            if(input$remove_select_gene != '_'){
              dataset_path = paste0(dataset_path,'_remove_',input$remove_select_gene)
            }
          }
        }
        print(dataset_path)
        dataset_path
      }
    }
  })
  
  observeEvent(input$load_data,{
    
    dataset_path = dataset_path()
    file_list = list.files('www')
    file_list
    dataset_rds_name = 'Seurat_data.rds'
    if(input$reduce_select !=  '_'){
      dataset_rds_name = paste0(input$select_dataset,'.rds',input$reduce_select)
    }else{
      dataset_rds_name = paste0(input$select_dataset,'.rds')
    }
    if(dataset_rds_name %in% file_list){
      withProgress(message = paste('readRDS',dataset_path()),{
        values$Seurat_data = readRDS(dataset_path)
      })
    }else{      
      withProgress(message = 'CreateSeuratObject',{
        
        data_dir = dataset_list[[input$select_dataset]]
        print(data_dir)
        data <- Read10X(data.dir = data_dir)
        Seurat_data = CreateSeuratObject(raw.data = data_path(), min.cells = input$min.cells, min.genes = input$min.genes, project = "immune_hsa")
        values$Seurat_data = Seurat_data
        values_save = reactiveValues()
        
        
      })
      
      withProgress(message = paste('saveRDS',dataset_path()),{
        saveRDS(Seurat_data,dataset_path)
        saveRDS(values_save,values_path())
        
      })
    }
    print(paste0('values_',dataset_rds_name))
    # if(paste0('values_',dataset_rds_name) %in% file_list){
    #   print('readRDS values')
    #   withProgress(message = 'readRDS values',{
    #     values_save = reactiveValues()
    #     names(values_save)
    #     values_load = readRDS(values_path())
    #     names(values_load)
    #     for(name in names(values_load)){
    #       print(name)
    #       values_save[[name]] = values_load[[name]]
    #     }
    #   })
    #   print(names(values_save))
    # }else{      
    #   withProgress(message = 'values_saves',{
    #     values_save = reactiveValues()
    #     saveRDS(values_save,values_path())
    #   })
    # }
    dataset_path
    
  })
  
  data_path = reactive({
    if(!is.null(input$select_dataset)){
      path_list = list.dirs(paste0('data/',input$select_dataset))
      mex_path = grep('mex',path_list,value = TRUE)  
      data_path = paste0(mex_path[nchar(mex_path) == max(nchar(mex_path))],'/')
      data_path
    }
  })
  output$dataset_text = renderText({
    dataset_path()
  })
  
  output$data_path_text = renderText({data_path()})
  
  
  #values_save = readRDS('www/values_save.rds')
  
  #read_rds = T
  #if(read_rds == T){
  #  values_save = readRDS('www/values_save.rds')
  #}else{
  #saveRDS(values_save, 'www/values_save.rds')
  #}
  output$features_plot = renderPlot({
    print('features_plot')
    withProgress(message = paste('readRDS',dataset_path()),{
      Seurat_data = readRDS(dataset_path())
    })
    withProgress(message = 'AddMetaData',{
      mito.genes <- grep(pattern = "^MT-", x = rownames(x = Seurat_data@data), value = TRUE)
      percent.mito <- Matrix::colSums(Seurat_data@raw.data[mito.genes, ]) / Matrix::colSums(Seurat_data@raw.data)
      
      # AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats
      Seurat_data <- AddMetaData(object = Seurat_data, metadata = percent.mito, col.name = "percent.mito")
    })
    if(input$save_rds_rb == T){
      withProgress(message = paste('saveRDS',dataset_path()),{
        saveRDS(Seurat_data,dataset_path())
      })
    }
    
    withProgress(message = 'VlnPlot',{
      print('feature_plot')
      values$Seurat_data = Seurat_data
      VlnPlot(object = Seurat_data, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
    })
  })
  #
  #
  output$gene_plot = renderPlot({
    print('gene_plot')
    par(mfrow = c(1, 2))
    if("percent.mito" %in% names(values$Seurat_data@meta.data)){
      GenePlot(object = values$Seurat_data, gene1 = "nUMI", gene2 = "percent.mito")
    }
    GenePlot(object = values$Seurat_data, gene1 = "nUMI", gene2 = "nGene")
  })
  #
  #
  output$var_plot = renderPlot({
    Seurat_data = values$Seurat_data
    withProgress(message = 'FilterCells',{
      Seurat_data <- FilterCells(object = Seurat_data, subset.names = c("nGene"), low.thresholds = c(200), high.thresholds = c(2500))
    })
    output$FilterCells_text = renderText({
      print("Seurat_data <- FilterCells(object = Seurat_data, subset.names = c('nGene'), low.thresholds = c(200), high.thresholds = c(2500))")
      
    })
    withProgress(message = 'NormalizeData',{
      Seurat_data <- NormalizeData(object = Seurat_data, normalization.method = "LogNormalize", scale.factor = 1e4)
    })
    output$NormalizeData_text = renderText({
      print("Seurat_data <- NormalizeData(object = Seurat_data, normalization.method = 'LogNormalize', scale.factor = 1e4)")
    })
    withProgress(message = 'FindVariableGenes',{
      Seurat_data <- FindVariableGenes(object = Seurat_data, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
    })
    output$FindVariableGenes_text = renderText({
      print("Seurat_data <- FindVariableGenes(object = Seurat_data, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)")
    })
    
    print('var_plot')
    if(input$save_rds_rb == T){
      
      withProgress(message = paste('saveRDS',dataset_path()),{
        saveRDS(Seurat_data,dataset_path())
      })
    }
    
    withProgress(message = 'FindVariableGenes',{
      Seurat_data <- FindVariableGenes(object = Seurat_data, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
    })
    
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
    if(is.null(values$Seurat_data_norm)){
      withProgress(message = paste('readRDS',dataset_path()),{
        print(dataset_path())
        data = readRDS(dataset_path())
      })
    }else{
      data = values$Seurat_data_norm
    }
    withProgress(message = 'running scaling',{
      Seurat_data <- ScaleData(object = data, vars.to.regress = c("nUMI"))
      
      if(input$save_rds_rb == T){
        withProgress(message = paste('saveRDS',dataset_path()),{
          saveRDS(Seurat_data,dataset_path())
        })
      }
      
      values$Seurat_data_scale = Seurat_data
    })
  })
  # 
  # 
  output$ScaleData_text = renderText({
    print("Seurat_data <- ScaleData(object = values$Seurat_data_norm, vars.to.regress = c('nUMI', 'percent.mito'))")
  })
  
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
      withProgress(message = paste('readRDS',dataset_path()),{
        print(dataset_path())
        data = readRDS(dataset_path())
      })
    }
    
    withProgress(message = 'RunPCA',{
      
      Seurat_data = RunPCA(object = data, pc.genes = data@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
    })
    values$Seurat_data_pca_1 <- Seurat_data
    withProgress(message = paste('saveRDS',dataset_path()),{
      
      if(input$save_rds_rb == T){
        withProgress(message = paste('saveRDS',dataset_path()),{
          saveRDS(Seurat_data,dataset_path())
        })
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
      withProgress(message = paste('readRDS',dataset_path()),{
        
        data = readRDS(dataset_path())
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
    if(input$save_rds_rb == T){
      withProgress(message = paste('saveRDS',dataset_path()),{
        saveRDS(Seurat_data,dataset_path())
      
      })
    }
    
    
    
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
      withProgress(message = paste('readRDS',dataset_path()),{
        
        data = readRDS(dataset_path())
      })
    }else{
      data = values$Seurat_data_pca_1
    }
    withProgress(message = 'JackStraw',{
      Seurat_data <- JackStraw(object = data, num.replicate = 100, display.progress = FALSE)
    })
    values$Seurat_data_jack = Seurat_data
    if(input$save_rds_rb == T){
      withProgress(message = paste('saveRDS',dataset_path()),{
        saveRDS(Seurat_data,dataset_path())
      })
    }
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
      withProgress(message = paste('readRDS',dataset_path()),{
        
        data = readRDS(dataset_path())
      })
    }else{
      data = values$Seurat_data_jack
    }
    # save.SNN = T saves the SNN so that the clustering algorithm can be rerun using the same graph
    # but with a different resolution value (see docs for full details)
    withProgress(message = 'FindClusters',{
      
      Seurat_data <- FindClusters(object = data, 
                                  reduction.type = "pca", 
                                  dims.use = 1:10, 
                                  resolution = 0.6, 
                                  print.output = 0, 
                                  save.SNN = TRUE,
                                  force.recalc = TRUE)
    })
    print('done FindClusters')
    values$Seurat_data_clusters = Seurat_data
    if(input$save_rds_rb == T){
      withProgress(message = paste('saveRDS',dataset_path()),{
        saveRDS(Seurat_data,dataset_path())
      })
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
      withProgress(message = paste('readRDS',dataset_path()),{
        data = readRDS(dataset_path())
      })
    }else{
      data = values$Seurat_data_clusters
    }
    withProgress(message = 'RunTSNE',{
      Seurat_data <- RunTSNE(object = data, dims.use = 1:10)
    })
    values$Seurat_data_tsne = Seurat_data
    #if (os == 'Darwin'){
    withProgress(message = paste('saveRDS',dataset_path()),{
      
      saveRDS(Seurat_data,dataset_path())
    })
    #}
  })
  
  output$TSNEPlot = renderPlot({
    if(!is.null(values$Seurat_data_tsne)){
      TSNEPlot(object = values$Seurat_data_tsne)
    }
  })
  ##### Visualisations ########
  
  plot_data = reactive({
    withProgress(message = paste('readRDS',dataset_path()),{
      
      data = readRDS(dataset_path())
      
    })
    
  })
  
  observeEvent(input$find_single_cluster,{
    data = plot_data()
    values_save = readRDS(values_path())
    
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
      #if (os == 'Darwin'){
      saveRDS(values_save,values_path())
      #}
      
    })
  })
  
  
  
  observeEvent(input$find_distinguishing_cluster,{
    data = plot_data()
    values_save = readRDS(values_path())
    
    withProgress(message = 'FindMarkers',{
      # find all markers distinguishing cluster 5 from clusters 0 and 3
      values_save$cluster_markers_distinguishing <- FindMarkers(object = data, 
                                                                ident.1 = input$single_cluster_number, 
                                                                ident.2 = as.numeric(input$find_distinguishing_cluster_select), 
                                                                min.pct = input$min.pct, 
                                                                logfc.threshold = input$logfc.threshold)
      values_save$cluster_markers_distinguishing$gene = rownames(values_save$cluster_markers_distinguishing)
      #if (os == 'Darwin'){
      saveRDS(values_save,values_path())
      #}
      
    })
  })
  
  
  observeEvent(input$find_all_clusters,{
    
    data = plot_data()
    values_save = readRDS(values_path())
    
    withProgress(message = 'FindAllMarkers',{
      values_save$data.markers <- FindAllMarkers(object = data, 
                                                 only.pos = TRUE, 
                                                 min.pct = input$min.pct, 
                                                 logfc.threshold = input$logfc.threshold)
      #if (os == 'Darwin'){
      saveRDS(values_save,values_path())
      #}
    })
    
  })
  
  output$single_markers_table = renderDataTable({
    #values_save = readRDS(values_path())
    
    if(!is.null(values_save$cluster.marker)){
      values_save$cluster.marker
    }
    
  })
  
  
  
  output$distinguish_markers_table = renderDataTable({
    #values_save = readRDS(values_path())
    values_save = readRDS(values_path())
    
    if(!is.null(values_save$cluster_markers_distinguishing)){
      values_save$cluster_markers_distinguishing
      
      
    }
    
  })
  
  observeEvent(input$select_dataset,{
    if(!is.null(values_path())){
      print('load values_save')
      #values_load = readRDS(values_path())
      values_save = reactiveValues()
      
      values_load = readRDS(values_path())
      names(values_load)
      for(name in names(values_load)){
        print(name)
        values_save[[name]] = values_load[[name]]
      }
    }
    
  })
  
  data_markers_table_top_n = reactive({
    if(!is.null(input$select_dataset)){
      values_save = readRDS(values_path())
      names(values_save)
      if(!is.null(values_save$data.markers)){
        data.markers.top_n = values_save$data.markers %>% group_by(cluster) %>% top_n(input$clust_num_display, avg_logFC)
        data.markers.top_n
      }
    }
  })
  
  output$data_markers_table = renderDataTable({
    data_markers_table_top_n()
  })
  
  
  
  output$plot_select_genes_ui = renderUI({
    values_save = readRDS(values_path())
    
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
    values_save = readRDS(values_path())
    
    # ```{r clusterHeatmap, fig.height=8, fig.width=15, message=FALSE, warning=FALSE}
    values_save$data.markers %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10
    # # setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name
    DoHeatmap(object = plot_data(), genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
    # ```
  })
  
  output$cluster_0_text = renderUI({
    if(!is.null(data_markers_table_top_n())){
      data = data_markers_table_top_n()
      dim(data)
      head(data)
      label = paste(unlist(data$gene[data$cluster == 0]),collapse = ', ')
      label
      textInput('cluster_0','Cluster 0',label)
    }
  })
  output$cluster_1_text = renderUI({
    if(!is.null(data_markers_table_top_n())){
      data = data_markers_table_top_n()
      dim(data)
      head(data)
      label = paste(unlist(data$gene[data$cluster == 1]),collapse = ', ')
      label
      textInput('cluster_1','Cluster 1',label)
    }
  })    
  
  output$cluster_2_text = renderUI({
    if(!is.null(data_markers_table_top_n())){
      data = data_markers_table_top_n()
      dim(data)
      head(data)
      label = paste(unlist(data$gene[data$cluster == 2]),collapse = ', ')
      label
      textInput('cluster_2','Cluster 2',label)
    }
  })    
  output$cluster_3_text = renderUI({
    if(!is.null(data_markers_table_top_n())){
      data = data_markers_table_top_n()
      dim(data)
      head(data)
      label = paste(unlist(data$gene[data$cluster == 3]),collapse = ', ')
      label
      textInput('cluster_3','Cluster 3',label)
    }
  })    
  output$cluster_4_text = renderUI({
    if(!is.null(data_markers_table_top_n())){
      data = data_markers_table_top_n()
      dim(data)
      head(data)
      label = paste(unlist(data$gene[data$cluster == 4]),collapse = ', ')
      label
      textInput('cluster_4','Cluster 4',label)
    }
  })    
  output$cluster_5_text = renderUI({
    if(!is.null(data_markers_table_top_n())){
      data = data_markers_table_top_n()
      dim(data)
      head(data)
      label = paste(unlist(data$gene[data$cluster == 5]),collapse = ', ')
      label
      textInput('cluster_5','Cluster 5',label)
    }
  })    
  output$cluster_6_text = renderUI({
    if(!is.null(data_markers_table_top_n())){
      data = data_markers_table_top_n()
      dim(data)
      head(data)
      label = paste(unlist(data$gene[data$cluster == 6]),collapse = ', ')
      label
      textInput('cluster_6','Cluster 6',label)
    }
  })    
  output$cluster_7_text = renderUI({
    if(!is.null(data_markers_table_top_n())){
      data = data_markers_table_top_n()
      dim(data)
      head(data)
      label = paste(unlist(data$gene[data$cluster == 7]),collapse = ', ')
      label
      textInput('cluster_7','Cluster 7',label)
    }
  })
  
  output$cluster_8_text = renderUI({
    if(!is.null(data_markers_table_top_n())){
      data = data_markers_table_top_n()
      dim(data)
      head(data)
      label = paste(unlist(data$gene[data$cluster == 8]),collapse = ', ')
      label
      textInput('cluster_8','Cluster 8',label)
    }
  })
  
  output$cluster_9_text = renderUI({
    if(!is.null(data_markers_table_top_n())){
      data = data_markers_table_top_n()
      dim(data)
      head(data)
      label = paste(unlist(data$gene[data$cluster == 9]),collapse = ', ')
      label
      textInput('cluster_9','Cluster 9',label)
    }
  })
  
  output$cluster_10_text = renderUI({
    if(!is.null(data_markers_table_top_n())){
      data = data_markers_table_top_n()
      dim(data)
      head(data)
      label = paste(unlist(data$gene[data$cluster == 10]),collapse = ', ')
      label
      textInput('cluster_10','Cluster 10',label)
    }
  })
  
  output$cluster_11_text = renderUI({
    if(!is.null(data_markers_table_top_n())){
      data = data_markers_table_top_n()
      dim(data)
      head(data)
      label = paste(unlist(data$gene[data$cluster == 11]),collapse = ', ')
      label
      textInput('cluster_11','Cluster 11',label)
    }
  })
  
  output$cluster_12_text = renderUI({
    if(!is.null(data_markers_table_top_n())){
      data = data_markers_table_top_n()
      dim(data)
      head(data)
      label = paste(unlist(data$gene[data$cluster == 12]),collapse = ', ')
      label
      textInput('cluster_12','Cluster 12',label)
    }
  })
  
  output$cluster_13_text = renderUI({
    if(!is.null(data_markers_table_top_n())){
      data = data_markers_table_top_n()
      dim(data)
      head(data)
      label = paste(unlist(data$gene[data$cluster == 13]),collapse = ', ')
      label
      textInput('cluster_13','Cluster 13',label)
    }
  })
  
  
  
  output$assign_clusters = renderPlot({
    if(!is.null(data_markers_table_top_n())){
      data = plot_data()
      current.cluster.ids = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13)
      length(current.cluster.ids)
      new.cluster.ids = c(input$cluster_0,input$cluster_1,input$cluster_2,input$cluster_3,input$cluster_4,input$cluster_5,
                          input$cluster_6,input$cluster_7,input$cluster_8,input$cluster_9,input$cluster_10,input$cluster_11,
                          input$cluster_12,input$cluster_13)
      length(new.cluster.ids)
      #current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)
      #new.cluster.ids <- c("CD4 T cells", "CD14+ Monocytes", "B cells", "CD8 T cells", "FCGR3A+ Monocytes", "NK cells", "Dendritic cells", "Megakaryocytes")
      data@ident <- plyr::mapvalues(x = data@ident, from = current.cluster.ids, to = new.cluster.ids)
      TSNEPlot(object = data, do.label = TRUE, pt.size = 0.5,label.size = 6)
    }
  })
  
})
