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
shinyServer(function(input, output, session) {
  volumes <- getVolumes()
  shinyDirChoose(input, 'directory', roots=root_path, session=session)
  path_full <- reactive({
    return(print(parseDirPath(root_path, input$directory)))
  })
  project_name = reactive({
    project_name = '_'
    if(!is.null(names(input$directory))){
      project_name = unlist(input$directory$path[2])
      if(!is.null(input$select_dataset)){
        if(input$select_dataset != '_'){
          project_name = gsub('.rds','',input$select_dataset)
        }
      }
    }else{
      if(!is.null(input$select_dataset)){
          project_name = gsub('.rds','',input$select_dataset)
        }
      #project_name = '_'
    }
    project_name  
  })
  
  values = reactiveValues()
  Seurat_values = reactiveValues()
  #values_save = reactiveValues()
  vis_values = reactiveValues()
   
  # observeEvent(input$select_dataset,{
  #   if(!is.null(values_path())){
  #       
  #       #values_load = readRDS(values_path())
  #       #if(!is.null(names(values))){
  #     print('removing values')
  #     #values = reactiveValues()
  #     #Seurate_values = reactiveValues(NULL)
  #         for(entry in names(Seurat_values)){
  #           print(entry)
  #           Seurat_values[[entry]] = NULL
  #         }
  #       #}
  #       if(file.exists(values_path())){
  #         withProgress(message = paste('Loading values from',values_path()),{
  #           print('load values')
  #           values_load = readRDS(values_path())
  #           names(values_load)
  #           for(name in names(values_load)){
  #             print(name)
  #             Seurat_values[[name]] = values_load[[name]]
  #           }
  #         })
  #       }else{
  #         print('values_path() does not exist')
  #         #values[[input$select_dataset]] = list()
  #       }
  #     }
  #   
  # })
  # 
  # observeEvent(input$save_rds_rb,{
  #   if(input$save_rds_rb == T){
  #     withProgress(message = paste('saveRDS(values',values_path()),{
  #       saveRDS(values,values_path())
  #     })
  #   }
  # })
  
  
  output$debug_ui = renderUI({
    #sysinf <- Sys.info()
    if (!is.null(sysinf)){
      #os <- sysinf['sysname']
      if (os == 'Darwin'){
        actionButton('debug','Debug')
      }
    }
  })
  
  output$value_names = renderText({
    #values = readRDS(values_path())
    if(!is.null(values)){
      paste(paste(names(values),collapse = ', '),'<br> <br>',paste(names(Seurat_values),collaspe = ', '))
    }
  })
  
  output$select_dataset_ui = renderUI({
    #datasets = values$data_list
    (file_list = list.files('www'))
    (datasets = grep('rds',file_list,value = T))
    dataset_entry = '_'
    if(!is.null(names(input$directory))){
      project_name = unlist(input$directory$path[2])
      dataset_entry = '_'
    }
    selectInput('select_dataset','Select Dataset',c('_',datasets),dataset_entry)
  })
  observeEvent(input$debug,{
    browser()
  })
  output$file_stucture_image <- renderImage({
 
    outfile <- 'Seurat_file_Structure.png'
    list(src = outfile,
         contentType = 'image/png',
         width = 300,
         height = 300,
         alt = "This is alternate text")
  },deleteFile = F)
  
  observeEvent(input$unzip,{
    withProgress(message = 'untar',{
      untar(input$file$datapath, exdir="./data/")
      #values$data_list = list.files('data/')
    })
    #unzip(input$file$datapath, list = TRUE, exdir = 'data/')
  })
  
  values_path = reactive({
    values_path = paste0('www/values_',project_name(),'.rds')
    values_path
  })
  
  vis_values_path = reactive({
    values_path = paste0('www/vis_values_',project_name(),'.rds')
    values_path
  })
  
  # values_path = reactive({
  #   if(!is.null(input$select_dataset)){
  #     if(!is.null(input$reduce_select)){
  #       if(input$reduce_select !=  '_'){
  #         values_path = paste0('www/values_',input$select_dataset,'.rds',input$reduce_select)
  #       }else{
  #         values_path = paste0('www/values_',input$select_dataset,'.rds')
  #       }
  #       #values_path = paste0('www/values_',input$select_dataset,'.rds')
  #       if(input$run_reduce_rb == T){
  #         
  #         if(!is.null(input$remove_select_gene)){
  #           if(input$remove_select_gene != '_'){
  #             values_path = paste0(values_path,'_remove_',input$remove_select_gene,'_',input$remove_threshold)
  #           }
  #         }
  #       }
  #       values_path
  #     }
  #   }
  # })
  #### REDUCE ####
  data_df = reactive({
    if(input$run_reduce_rb == T){
      withProgress(message = 'Read Data',{
      print(data_path())
      data <- Read10X(data.dir = data_path())
      df = as.data.frame(as.matrix(data))
      dim(df)
      #df = as.data.frame(as.matrix(data))
      df
      })
    }else{
      if(file.exists(dataset_path())){
        withProgress(message = paste('readRDS',dataset_path()),{
          data = readRDS(dataset_path())@assays$RNA@data
        })
          df = as.data.frame(as.matrix(data))
          dim(df)
          df
        }
        
      }
    })
#  })
  

  output$remove_select_gene_ui = renderUI({
    if(input$run_reduce_rb == T){
      
      selectInput('remove_select_gene','Select Marker to Remove',c('_',rownames(data_df())),'FSCN1')
    }
  })
  
  
  output$remove_hist_ui = renderUI({
    if(input$run_reduce_rb == T){
      plotOutput('remove_hist')
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
    }else{
      print(paste('Total Number of cells =',dim(data_df())[2]))
    }
  })
  
  observeEvent(input$reload_data,{
    #withProgress(message = 'Re-CreateSeuratObject',{
      #print(data_path())
    
      data <- Read10X(data.dir = data_path())
      # if(input$run_reduce_rb == T){
      #   if(input$remove_select_gene == '_'){
      #     if(!is.null(path_full())){
      #       data = Read10X(data.dir = path_full())
      #     }else{
      #       data <- Read10X(data.dir = data_path())
      #     }
      #   }else{
      #     data = as(as.matrix(data_reduce()),'dgCMatrix')
      #   }
      # }else{
      #   if(!is.null(path_full())){
      #     data = Read10X(data.dir = path_full())
      #   }else{
      #     data <- Read10X(data.dir = data_path())
      #   }
      #   
      # }
      dim(data)
      Seurat_data <- CreateSeuratObject(counts = data, project = project_name(), min.cells = input$min.cells, min.features = input$min.genes)
      #values = NULL
      #for(entry in names(values)){
      #  print(entry)
      #  values[[entry]] = NULL
      #}
      #values = reactiveValues()
      #Seurat_values = NULL
      values$data_path = data_path()
      values$Seurat_data = Seurat_data
      #values$data_path = path_full()
      
      #Seurat_data = CreateSeuratObject(raw.data = data, min.cells = input$min.cells, min.genes = input$min.genes, project = input$select_dataset)
      values$data = data
      #values$Seurat_data = Seurat_data
      #values_save = reactiveValues()
      
      
    #})
    withProgress(message = paste('saveRDS',dataset_path()),{
      
      saveRDS(Seurat_data,original_dataset_path())
      
      saveRDS(values,values_path())
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
    (dataset_rds_name = paste0(project_name(),'.rds'))
    dataset_path = paste0('www/',dataset_rds_name)
    dataset_path
  })
  original_dataset_path = reactive({
    (dataset_rds_name = paste0(project_name(),'.rds'))
    dataset_path = paste0('www/original_',dataset_rds_name)
    dataset_path
  })
  

  data_path = reactive({
    if(!is.null(path_full())){
      data_path = path_full()
    }else{
      if(!is.null(input$select_dataset)){
        (path_list = list.dirs(paste0('data/',input$select_dataset)))
        (mex_path = grep('mex',path_list,value = TRUE))
        (data_path = paste0(mex_path[nchar(mex_path) == max(nchar(mex_path))],'/'))
        data_path
      }
    }
    #values$data_path = data_path
    data_path
  })
  output$dataset_text = renderText({
    dataset_path()
  })
  
  output$data_path_text = renderText({
    values$data_path
  })
  
  
  #values_save = readRDS('www/values_save.rds')
  
  #read_rds = T
  #if(read_rds == T){
  #  values_save = readRDS('www/values_save.rds')
  #}else{
  #saveRDS(values_save, 'www/values_save.rds')
  #}
  output$data_info = renderPrint({
    data = values$data
    #dataset_path()
    #Seurat_data = readRDS(dataset_path())
    #Seurat_data
    data[c("CD3D","TCL1A","MS4A1"), 1:30]
    dense.size <- object.size(x = as.matrix(x = data))
    dense.size
    sparse.size <- object.size(x = data)
    sparse.size
    dense.size / sparse.size
  })
  
  output$features_plot = renderPlot({
    print('features_plot')
    withProgress(message = paste('readRDS',dataset_path()),{
      Seurat_data = readRDS(dataset_path())
      
    })
    withProgress(message = 'AddMetaData',{
      #mito.genes <- grep(pattern = "^MT-", x = rownames(x = Seurat_data@assays$RNA@data), value = TRUE)
      #percent.mito <- Matrix::colSums(Seurat_data@raw.data[mito.genes, ]) / Matrix::colSums(Seurat_data@raw.data)
      Seurat_data[["percent.mt"]] <- PercentageFeatureSet(object = Seurat_data, pattern = "^MT-")
      
      
      # AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats
      #Seurat_data <- AddMetaData(object = Seurat_data, metadata = percent.mito, col.name = "percent.mito")
    })
    if(input$save_rds_rb == T){
      withProgress(message = paste('saveRDS',dataset_path()),{
        saveRDS(Seurat_data,dataset_path())
        saveRDS(values,values_path())
      })
    }
    
    withProgress(message = 'VlnPlot',{
      print('feature_plot')
      values$Seurat_data = Seurat_data
      
      #VlnPlot(object = Seurat_data, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
      VlnPlot(object = Seurat_data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
      
    })
  })
  
  
  output$scatter_plot = renderPlot({
    Seurat_data = values$Seurat_data
    plot1 <- FeatureScatter(object = Seurat_data, feature1 = "nCount_RNA", feature2 = "percent.mt") 
    plot2 <- FeatureScatter(object = Seurat_data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") 
    CombinePlots(plots = list(plot1,plot2))
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
      Seurat_data <- subset(x = Seurat_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
      
      #Seurat_data <- FilterCells(object = Seurat_data, subset.names = c("nGene"), low.thresholds = c(200), high.thresholds = c(2500))
    })
    output$FilterCells_text = renderText({
      #print("Seurat_data <- FilterCells(object = Seurat_data, subset.names = c('nGene'), low.thresholds = c(200), high.thresholds = c(2500))")
      print('Seurat_data <- subset(x = Seurat_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)')
    })
    withProgress(message = 'NormalizeData',{
      Seurat_data <- NormalizeData(object = Seurat_data, normalization.method = "LogNormalize", scale.factor = 1e4)
      
      #Seurat_data <- NormalizeData(object = Seurat_data, normalization.method = "LogNormalize", scale.factor = 1e4)
    })
    output$NormalizeData_text = renderText({
      print("Seurat_data <- NormalizeData(object = Seurat_data, normalization.method = 'LogNormalize', scale.factor = 1e4)")
    })
    withProgress(message = 'FindVariableFeatures',{
      Seurat_data <- FindVariableFeatures(object = Seurat_data,selection.method = 'vst', nfeatures = 2000)
      
      #Seurat_data <- FindVariableGenes(object = Seurat_data, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
    })
    output$FindVariableGenes_text = renderText({
      #print("Seurat_data <- FindVariableGenes(object = Seurat_data, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)")
      print("Seurat_data <- FindVariableFeatures(object = Seurat_data,selection.method = 'vst', nfeatures = 2000)")
      })
    
    print('var_plot')
    if(input$save_rds_rb == T){
      
      withProgress(message = paste('saveRDS',dataset_path()),{
        saveRDS(Seurat_data,dataset_path())
        saveRDS(values,values_path())
      })
    }
    Seurat_values$Seurat_data_norm = Seurat_data
    
    plot1 <- VariableFeaturePlot(object = Seurat_data)
    top10 <- head(x = VariableFeatures(object = Seurat_data), 20)
    
    plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    CombinePlots(plots = list(plot1, plot2))
    
    #withProgress(message = 'FindVariableGenes',{
    #  Seurat_data <- FindVariableGenes(object = Seurat_data, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
    #})
    #values$Seurat_data_norm = Seurat_data
    
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
  
  #obsereEvent(input$subset{
  #  Seurat_data = readRDS()
  #  Seurat_data <- subset(x = Seurat_data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
  #  
  #})
  
  observeEvent(input$run_scaling,{
    if(is.null(Seurat_values$Seurat_data_norm)){
      withProgress(message = paste('readRDS',dataset_path()),{
        print(dataset_path())
        Seurat_data = readRDS(dataset_path())
      })
    }else{
      Seurat_data = Seurat_values$Seurat_data_norm
    }
    withProgress(message = 'running scaling',{
      all.genes <- rownames(x = Seurat_data)
      length(all.genes)
      Seurat_data <- ScaleData(object = Seurat_data)
      #Seurat_data <- ScaleData(object = Seurat_data, vars.to.regress = c("nUMI"))
      
      if(input$save_rds_rb == T){
        withProgress(message = paste('saveRDS',dataset_path()),{
          saveRDS(Seurat_data,dataset_path())
          saveRDS(values,values_path())
        })
      }
      
      Seurat_values$Seurat_data_scale = Seurat_data
    })
  })
  # 
  # 
  output$ScaleData_text = renderText({
    print("Seurat_data <- ScaleData(object = Seurat_values$Seurat_data_norm, vars.to.regress = c('nUMI'))")
  })
  
  output$scaling_done_text = renderText({
    if(!is.null(Seurat_values$Seurat_data_scale)){
      print('done')
    }
  })
  #output$pca_1 = renderText({
  #  print('pca_1')
  

  
  
  
  observeEvent(input$run_pca,{ 
    print('pca_1')
    
    if(!is.null(Seurat_values$Seurat_data_scale)){
      Seurat_data = Seurat_values$Seurat_data_scale
    }else{
      withProgress(message = paste('readRDS',dataset_path()),{
        print(dataset_path())
        Seurat_data = readRDS(dataset_path())
      })
    }
    
    withProgress(message = 'RunPCA',{
      Seurat_data <- RunPCA(object = Seurat_data, features = VariableFeatures(object = Seurat_data))
      #Seurat_data <- RunPCA(object = Seurat_data)
      #Seurat_data = RunPCA(object = Seurat_data, pc.genes = Seurat_data@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
      #Seurat_data = RunPCA(object = data, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
      
    })
    Seurat_values$Seurat_data_pca_1 <- Seurat_data
    withProgress(message = paste('saveRDS',dataset_path()),{
      
      if(input$save_rds_rb == T){
        withProgress(message = paste('saveRDS',dataset_path()),{
          saveRDS(Seurat_data,dataset_path())
          saveRDS(values,values_path())
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
    if(!is.null(Seurat_values$Seurat_data_pca_1)){
      PrintPCA(object = Seurat_values$Seurat_data_pca_1, pcs.print = 1:5, genes.print = 5, use.full = FALSE)
    }
  })
  
  output$VizPCA = renderPlot({
    print('VizPCA')
    if(!is.null(Seurat_values$Seurat_data_pca_1)){
      
      #VizPCA(object = Seurat_values$Seurat_data_pca_1, pcs.use = 1:2)
      Seurat_data = Seurat_values$Seurat_data_pca_1
      print(x = Seurat_data[['pca']], dims = 1:5, nfeatures = 5)
      VizDimLoadings(object = Seurat_data, dims = 1:2, reduction = 'pca')
      #DimPlot(object = Seurat_data, reduction = 'pca')
    }
  })
  
  output$VizPlot = renderPlot({
    print('VizPlot')
    if(!is.null(Seurat_values$Seurat_data_pca_1)){
      DimPlot(object = Seurat_values$Seurat_data_pca_1, reduction = 'pca')
      
      #PCAPlot(object = Seurat_values$Seurat_data_pca_1, dim.1 = 1, dim.2 = 2)
    }
  })
  
  observeEvent(input$ProjectPCA,{
    if(is.null(Seurat_values$Seurat_data_pca_1)){
      withProgress(message = paste('readRDS',dataset_path()),{
        
        data = readRDS(dataset_path())
      })
    }else{
      data = Seurat_values$Seurat_data_pca_1
    }
    # ProjectPCA scores each gene in the dataset (including genes not included in the PCA) based on their correlation 
    # with the calculated components. Though we don't use this further here, it can be used to identify markers that 
    # are strongly correlated with cellular heterogeneity, but may not have passed through variable gene selection. 
    # The results of the projected PCA can be explored by setting use.full=T in the functions above
    withProgress(message = 'ProjectPCA',{
      
      Seurat_data <- ProjectPCA(object = data, do.print = FALSE)
    })
    Seurat_values$Seurat_data_ProjectPCA = Seurat_data
    if(input$save_rds_rb == T){
      withProgress(message = paste('saveRDS',dataset_path()),{
        saveRDS(Seurat_data,dataset_path())
        saveRDS(values,values_path())
      })
    }
    
    
    
  })
  
  output$PCHeatmap = renderPlot({
    if(!is.null(Seurat_values$Seurat_data_ProjectPCA)){
      PCHeatmap(object = Seurat_values$Seurat_data_ProjectPCA, pc.use = 1:input$pca_heatmap_dims, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
    }
  })
  output$PCHeatmap_2 = renderPlot({
    if(!is.null(Seurat_values$Seurat_data_ProjectPCA)){
      PCHeatmap(object = Seurat_values$Seurat_data_ProjectPCA, pc.use = 1:12, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
    }
  })
  
  
  output$DimHeatmap = renderPlot({ 
    Seurat_data = Seurat_values$Seurat_data_pca_1
    Seurat_data
    DimHeatmap(object = Seurat_data, dims = 1:input$pca_heatmap_dims, cells = 500, balanced = TRUE)
  })
  
  observeEvent(input$Jacksaw,{
    print('JackStraw')
    if(is.null(Seurat_values$Seurat_data_pca_1)){
      withProgress(message = paste('readRDS',dataset_path()),{
        
        Seurat_data = readRDS(dataset_path())
      })
    }else{
      Seurat_data = Seurat_values$Seurat_data_pca_1
    }
    withProgress(message = 'JackStraw',{
      Seurat_data <- JackStraw(object = Seurat_data, num.replicate = 100)
      Seurat_data <- ScoreJackStraw(object = Seurat_data, dims = 1:20)
      #Seurat_data <- JackStraw(object = data, num.replicate = 100, display.progress = FALSE)
    })
    Seurat_values$Seurat_data_jack = Seurat_data
    if(input$save_rds_rb == T){
      withProgress(message = paste('saveRDS',dataset_path()),{
        saveRDS(Seurat_data,dataset_path())
        saveRDS(values,values_path())
      })
    }
  })
  output$JackStrawPlot = renderPlot({
    if(!is.null(Seurat_values$Seurat_data_jack)){
      JackStrawPlot(object = Seurat_values$Seurat_data_jack, dims = 1:15)
      #JackStrawPlot(object = Seurat_values$Seurat_data_jack, PCs = 1:12)
    }
  })
  
  output$PCElbowPlot = renderPlot({
    if(!is.null(Seurat_values$Seurat_data_jack)){
      ElbowPlot(object = Seurat_values$Seurat_data_jack)
      #PCElbowPlot(object = Seurat_values$Seurat_data_jack)
    }
  })
  
  observeEvent(input$find_clusters,{
    print('FindClusters')
    if(is.null(Seurat_values$Seurat_data_jack)){
      withProgress(message = paste('readRDS',dataset_path()),{
        
        Seurat_data = readRDS(dataset_path())
      })
    }else{
      Seurat_data = Seurat_values$Seurat_data_jack
    } 
    
    
    # save.SNN = T saves the SNN so that the clustering algorithm can be rerun using the same graph
    # but with a different resolution value (see docs for full details)
    withProgress(message = 'FindClusters',{
      Seurat_data <- FindNeighbors(object = Seurat_data, dims = 1:10)
      Seurat_data <- FindClusters(object = Seurat_data, resolution = 0.5)
      
      # Seurat_data <- FindClusters(object = data, 
      #                             reduction.type = "pca", 
      #                             dims.use = 1:10, 
      #                             resolution = 0.6, 
      #                             print.output = 0, 
      #                             save.SNN = TRUE,
      #                             force.recalc = TRUE)
    })
    print('done FindClusters')
    Seurat_values$Seurat_data_clusters = Seurat_data
    if(input$save_rds_rb == T){
      withProgress(message = paste('saveRDS',dataset_path()),{
        saveRDS(Seurat_data,dataset_path())
        saveRDS(values,values_path())
    })
    }
  })
  
  output$cluster_text = renderText({
    if(!is.null(Seurat_values$Seurat_data_clusters)){
      head(x = Idents(object = Seurat_values$Seurat_data_clusters), 5)
      DimPlot(object = pbmc, reduction = 'umap')
      #PrintFindClustersParams(object = Seurat_values$Seurat_data_clusters)
    }
  })
  
  observeEvent(input$run_umap,{
    print('FindClusters')
    if(is.null(Seurat_values$Seurat_data_jack)){
      withProgress(message = paste('readRDS',dataset_path()),{
        
        Seurat_data = readRDS(dataset_path())
      })
    }else{
      Seurat_data = Seurat_values$Seurat_data_clusters
    } 

    withProgress(message = 'uMAP',{
      Seurat_data <- RunUMAP(object = Seurat_data, dims = 1:10)
      
    })
    print('done RunUMAP')
    Seurat_values$Seurat_umap = Seurat_data
    #if(input$save_rds_rb == T){
    withProgress(message = paste('saveRDS',dataset_path()),{
        saveRDS(Seurat_data,dataset_path())
        saveRDS(values,values_path())
    })
    #}
  })
  
  output$umap_plot = renderPlot({
    if(!is.null(Seurat_values$Seurat_umap)){
      DimPlot(object = Seurat_values$Seurat_umap, reduction = 'umap')
    }
  }) 
  
  observeEvent(input$run_tsne,{
    print('TSNE')
    if(is.null(Seurat_values$Seurat_data_clusters)){
      withProgress(message = paste('readRDS',dataset_path()),{
        data = readRDS(dataset_path())
      })
    }else{
      data = Seurat_values$Seurat_data_clusters
    }
    withProgress(message = 'RunTSNE',{
      Seurat_data <- RunTSNE(object = data, dims.use = 1:10)
    })
    Seurat_values$Seurat_data_tsne = Seurat_data
    #if (os == 'Darwin'){
    withProgress(message = paste('saveRDS',dataset_path()),{
      
      saveRDS(Seurat_data,dataset_path())
    })
    #}
  })
  
  output$TSNEPlot = renderPlot({
    if(!is.null(Seurat_values$Seurat_data_tsne)){
      TSNEPlot(object = Seurat_values$Seurat_data_tsne)
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
    vis_values = readRDS(values_vis_path())
     
    withProgress(message = 'FindMarkers',{
      vis_values$cluster.marker <- FindMarkers(object = data, 
                                      ident.1 = input$single_cluster_number, 
                                      min.pct = input$min.pct)
      #cluster5.markers <- FindMarkers(object = data, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
      
      #values_save$cluster.marker <- FindMarkers(object = data, 
                                                #ident.1 = input$single_cluster_number, 
                                                #test.use = input$find_single_cluster_test, 
                                                #min.pct = input$min.pct, 
                                                #logfc.threshold = input$logfc.threshold)
      vis_values$cluster.marker$gene = rownames(vis_values$cluster.marker)
      vis_values$cluster.marker$cluster = input$single_cluster_number
      names(vis_values)
    })
    if(input$save_rds_rb == T){
      withProgress(message = 'saveRDS vis_values save',{
        #if (os == 'Darwin'){
        saveRDS(vis_values,vis_values_path())
      
      })
    }
  })
  
  
  
  observeEvent(input$find_distinguishing_cluster,{
    data = plot_data()
    vis_values = readRDS(vis_values_path())
    
    withProgress(message = 'FindMarkers',{
      # find all markers distinguishing cluster 5 from clusters 0 and 3
      vis_values$cluster_markers_distinguishing <- FindMarkers(object = data, 
                                                                ident.1 = input$single_cluster_number, 
                                                                ident.2 = as.numeric(input$find_distinguishing_cluster_select), 
                                                                min.pct = input$min.pct, 
                                                                logfc.threshold = input$logfc.threshold)
      vis_values$cluster_markers_distinguishing$gene = rownames(vis_values$cluster_markers_distinguishing)
      #if (os == 'Darwin'){
      #if(input$sav)
      saveRDS(vis_values,vis_values_path())
      #}
      
    })
  })
  
  
  observeEvent(input$find_all_clusters,{
    
    data = plot_data()
    vis_values = readRDS(vis_values_path())
    
    withProgress(message = 'FindAllMarkers',{
      vis_values$data.markers <- FindAllMarkers(object = data, only.pos = TRUE, min.pct = input$min.pct, logfc.threshold = input$logfc.threshold)
      #values$data.markers.group %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
      #values_save$data.markers <- FindAllMarkers(object = data, 
                                                 # only.pos = TRUE, 
                                                 # min.pct = input$min.pct, 
                                                 # logfc.threshold = input$logfc.threshold)
      #if (os == 'Darwin'){
      print('Done FineAllMarkers')
      if(input$save_rds_rb == T){
        withProgress(message = paste('saveRDS',vis_values_path()),{
          saveRDS(vis_values,vis_values_path())
        })
      }
      #}
    })
    
  })
  
  output$single_markers_table = renderDataTable({
    vis_values = readRDS(vis_values_path())
    
    if(!is.null(vis_values$cluster.marker)){
      vis_values$cluster.marker
    }
    
  })
  
  
  
  output$distinguish_markers_table = renderDataTable({
    #values_save = readRDS(values_path())
    #values_save = readRDS(values_path())
    
    if(!is.null(vis_values$cluster_markers_distinguishing)){
      vis_values$cluster_markers_distinguishing
      
      
    }
    
  })
  

  
  data_markers_table_top_n = reactive({
    #if(!is.null(input$select_dataset)){
    #  values_save = readRDS(values_path())
    #  names(values_save)
      
      if(!is.null(vis_values$data.markers)){
        data.markers.top_n = vis_values$data.markers %>% group_by(cluster) %>% top_n(input$clust_num_display, avg_logFC)
        data.markers.top_n
      }
    #}
  })
  
  output$data_markers_table = renderDataTable({
    data_markers_table_top_n()
  })
  
  
  
  output$plot_select_genes_ui = renderUI({
    vis_values = readRDS(vis_values_path())
    full_gene_list = vis_values$data.markers$gene
    length(full_gene_list)
    gene_list = vis_values$data.markers$gene[1:3]
    length(gene_list)
    gene_search = c('IL13','IL17')
    gene_list = c()
    for(gene in gene_search){
      gene_list = c(gene_list,grep(gene,full_gene_list,value = T))
    }
    gene_list = input$features_of_interest
    selectInput('plot_select_genes','Select Genes',full_gene_list,gene_list,multiple = T)
    
  })
  output$vlnplot_1 = renderPlot({
    data = plot_data()
    VlnPlot(object = data, features = input$plot_select_genes)
    
    #VlnPlot(object = data, features.plot = input$plot_select_genes)
  }) 
  output$vlnplot_2 = renderPlot({
    data = plot_data()
    VlnPlot(object = data, features = input$plot_select_genes, slot = 'counts', log = TRUE)
    
    # you can plot raw UMI counts as well
    #VlnPlot(object = data, features.plot = input$plot_select_genes, use.raw = TRUE, y.log = TRUE)
  }) 
  output$feature_plot_1 = renderPlot({
    #data = readRDS('www/Seurat_data.rds')
    data = plot_data()
    FeaturePlot(object = data, features = input$plot_select_genes)
    
    #FeaturePlot(object = data, features.plot = input$plot_select_genes, cols.use = c("grey", "blue"), reduction.use = "tsne")
  })
  
  output$doheatmap_select = renderPlot({
    data = plot_data()
    DoHeatmap(object = data, features = input$plot_select_genes) + NoLegend()
  })
  
  output$doheatmap = renderPlot({
    
    data = plot_data()
    vis_values = readRDS(vis_values_path())
    
    # ```{r clusterHeatmap, fig.height=8, fig.width=15, message=FALSE, warning=FALSE}
    #values_save$data.markers %>% group_by(cluster) %>% top_n(10, avg_logFC) -> top10
    # # setting slim.col.label to TRUE will print just the cluster IDS instead of every cell name
    #DoHeatmap(object = plot_data(), genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
    # ```
    top10 = vis_values$data.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
    top10$gene[1:10]
    DoHeatmap(object = data, features = top10$gene[1:10]) + NoLegend()
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
      num  =  length(levels(x = data))
      #new.cluster.ids <- c("Memory CD4 T", "CD14+ Mono", "Naive CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Mk")[1:num]
      new.cluster.ids = c(input$cluster_0,input$cluster_1,input$cluster_2,input$cluster_3,input$cluster_4,input$cluster_5,
                          input$cluster_6,input$cluster_7,input$cluster_8,input$cluster_9,input$cluster_10,input$cluster_11,
                          input$cluster_12,input$cluster_13)[1:num]
      new.cluster.ids
      names(x = new.cluster.ids) <- levels(x = data)
      new.cluster.ids
      data <- RenameIdents(object = data, new.cluster.ids)
      DimPlot(object = data, reduction = 'umap', label = TRUE, pt.size = 0.5) + NoLegend()
      
      #current.cluster.ids = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13)
      #length(current.cluster.ids)
      #new.cluster.ids = c(input$cluster_0,input$cluster_1,input$cluster_2,input$cluster_3,input$cluster_4,input$cluster_5,
      #                    input$cluster_6,input$cluster_7,input$cluster_8,input$cluster_9,input$cluster_10,input$cluster_11,
      #                    input$cluster_12,input$cluster_13)
      #length(new.cluster.ids)
      #current.cluster.ids <- c(0, 1, 2, 3, 4, 5, 6, 7)
      #new.cluster.ids <- c("CD4 T cells", "CD14+ Monocytes", "B cells", "CD8 T cells", "FCGR3A+ Monocytes", "NK cells", "Dendritic cells", "Megakaryocytes")
      #data@ident <- plyr::mapvalues(x = data@ident, from = current.cluster.ids, to = new.cluster.ids)
      #TSNEPlot(object = data, do.label = TRUE, pt.size = 0.5,label.size = 6)
    }
  })
  
  
  ##### FEATURES ######
  
  #### Individual Datasets ####
  original_data_df = reactive({
    Seurat_data =  plot_data()
    data = as.data.frame(Seurat_data@assays$RNA@data)
    dim(data)
    data$features = rownames(data)
    data = data %>% select(features,everything())
    data
    saveRDS(data,paste0('www/data/',project_name(),'.rds'))
    write.csv2(data,paste0('www/data_csv/',project_name(),'.csv'))
    data_l = data %>%  
      tidyr::gather(variables,value,c(2:dim(data)[2]))
    dim(data_l)
    as.tbl(data_l)
    data_l = data_l %>% filter(value > 0)
    dim(data_l)
    saveRDS(data_l,paste0('www/data_long/',project_name(),'.rds'))
    data
  })
  
  # observeEvent(input$save_data,{
  #   saveRDS(original_data_df(),paste0('www/data/',project_name(),'.rds'))
  #   data_l = original_data_df() %>%  
  #     tidyr::gather(variables,value,c(2:dim(original_data_df())[2]))
  #   dim(data_l)
  #   as.tbl(data_l)
  #   data_l = data_l %>% filter(value > 0)
  #   dim(data_l)
  #   saveRDS(data_l,paste0('www/data_long/',project_name(),'.rds'))
  #   
  # })
  
  output$original_data_select_features_ui = renderUI({
    features = unique(original_data_df()$features)
    length(features)
    saveRDS(unique(c(uploaded_features,features)),'www/saved_features.rds')
    selected_features = c('IL13','IL16')
    gene_search = c('IL13','IL17')
    gene_list = c()
    for(gene in gene_search){
      gene_list = c(gene_list,grep(gene,features,value = T))
    }
    gene_list = input$features_of_interest
    selectInput('original_data_features_select','Select Features',features,gene_list,multiple = T)
  })
  
  
  
  original_data_plot = reactive({
    if(!is.null(input$original_data_features_select)){
      if(length(input$original_data_features_select) > 0){
        
        selected_features = input$original_data_features_select
        data = original_data_df()
        plot_data = data %>% filter(features %in% selected_features)
        dim(plot_data)
        #plot_data
        plot_data_l = plot_data %>%  
          tidyr::gather(variables,value,c(2:dim(plot_data)[2]))
        dim(plot_data_l)
        as.tbl(plot_data_l)
        #plot_name = paste0('www/feature_data/',project_name(),'__',paste(selected_features,collapse = '_'),'.rds')
        plot_name = paste0('www/feature_data/',project_name(),'.rds')
        
        plot_name
        withProgress(message = paste('save rds : ',plot_name),{ 
          saveRDS(plot_data_l,plot_name)
        })
        #plot_data
        as.tbl(plot_data_l)
        plot_data_l_min = plot_data_l %>% filter(value >= input$vis_data_value_min)
        as.tbl(plot_data_l)
        boxplot = ggplot(plot_data_l_min, aes(y = value, x = features)) +
          geom_boxplot(aes(col = features)) + 
          theme(axis.text.x = element_text(angle = 90)) + 
          ggtitle(project_name()) 
        
        ggsave(paste0('www/plots/boxplot_',project_name(),'.png'),boxplot,device = 'png')
        
        density = ggplot(plot_data_l_min,aes(x = value, col = features)) +
          geom_density() +
          ggtitle(project_name())
        
        ggsave(paste0('www/plots/density_',project_name(),'.png'),density,device = 'png')
        
        reduced_data = plot_data_l %>% filter(value >= input$vis_feature_threshold)
        
        (data_threshold = reduced_data %>% 
            group_by(features) %>% 
            summarise(count = n()) %>% 
            ungroup())
        
        
        #threshold = 0.5
        g_hist = ggplot(plot_data_l_min,aes(x = value,fill = features)) + 
          geom_histogram(position = 'dodge',col = 'black',alpha =0.5) +
          ggtitle(project_name())
        if(input$vis_feature_threshold > 0){
          g_hist = g_hist + geom_vline(xintercept = input$vis_feature_threshold, size = 2,col = 'red')
          g_hist = g_hist + annotate('text',x = Inf,y = Inf, label = paste('value >=',input$vis_feature_threshold), vjust = 1, hjust = 1)
        }
        for(i in c(1:dim(data_threshold)[1])){
          print(i)
          text = paste(data_threshold[i,],collapse = ' = ')
          text
          print(text)
          g_hist = g_hist + annotate('text',x = Inf,y = Inf, label = text, vjust = i*1.25+1.5, hjust = 1)
        }
        
        
        g_hist
        
        ggsave(paste0('www/plots/histogram_',project_name(),'.png'),g_hist,device = 'png')
        
        
        
        
        
        
        
        list(plot_data = plot_data,boxplot = boxplot, density = density,g_hist = g_hist,data_threshold = data_threshold, reduced_data = reduced_data)
      }
    }
  })
  
  output$vis_feature_histogram = renderPlot({
    if(!is.null(original_data_plot())){
      original_data_plot()$g_hist
    }
  })
  
  output$vis_feature_count = renderDataTable({
    if(!is.null(original_data_plot())){
      original_data_plot()$data_threshold
    }
  })
  
  output$original_data_boxplot = renderPlot({
    if(!is.null(original_data_plot())){
      original_data_plot()$boxplot
    }
  })
  output$original_data_density = renderPlot({
    if(!is.null(original_data_plot())){
      
      p = original_data_plot()$density
      
      if(input$data_density_y != 0){
        p = p + ylim(0,input$data_density_y)
      }
      p
    }
  })
  
  
  output$original_data_table = renderDataTable({
    if(!is.null(original_data_plot())){
      
      plot_data = original_data_plot()$reduced_data
      plot_data
    }
    
  })
  
  
  #### Combined ###
  
  combine_df_l_function = function(data_path){
    #data_path = 'www/data_long/'
    file_list = list.files(data_path)
    file_list
    entry = file_list[1]
    file_path = paste0(data_path,entry)
    print(file_path)
    df = readRDS(file_path)
    (project_name = gsub('.rds','',entry))
    
    #(project_name = unlist(strsplit(entry,'__'))[1])
    df$project = project_name
    dim(df)
    as.tbl(df)
    
    for(entry in file_list[2:length(file_list)]){
      file_path = paste0(data_path,entry)
      print(file_path)
      df_n = readRDS(file_path)
      print(dim(df_n))
      if(dim(df_n)[1] > 0){
        (project_name = gsub('.rds','',entry))
        #(project_name = unlist(strsplit(entry,'__'))[1])
        df_n$project = project_name
        df = rbind(df,df_n)
      }
    }
    
    dim(df)
    as.tbl(df)
    return(df)
  }
  
  combined_features = reactive({
    data_path = 'www/feature_data/'
    df = combine_df_l_function(data_path)
    
    dim(df)
    as.tbl(df)
    
    project_list = unique(df$project)
    print(project_list)
    
    psoriasis_grouping = c("19006-04-02","19006-07-02")
    non_psoriasis_grouping = c("19006-01","19006-04-01","abd154","19006-06-02","19006-07-01","19006-07-03","19006-05-1")
    
    spongiotic_grouping = c("19006-07-03","19006-07-01")
    non_spongiotic_grouping = c("19006-01","19006-04-01","abd154","19006-06-02","19006-04-02","19006-07-02","19006-05-1")
    
    
    df$psoriasis = NA
    df$spongiotic = NA
    
    df$psoriasis[df$project %in% psoriasis_grouping] = T
    df$psoriasis[df$project %in% non_psoriasis_grouping] = F
    unique(df$project[is.na(df$psoriasis)])
    
    df$spongiotic[df$project %in% spongiotic_grouping] = T
    df$spongiotic[df$project %in% non_spongiotic_grouping] = F
    unique(df$project[is.na(df$spongiotic)])
    
    (features = unique(df$features))
    colnames(df)
    (groups = colnames(df)[5:length(colnames(df))])
    
    list(df = df,features = features,groups = groups)
  })
  
  output$select_combined_feature_1_ui = renderUI({
    selectInput('combined_feature_1_selected_feature','Select Feature 1',combined_features()$features,combined_features()$features[1])
  })
  output$select_combined_feature_2_ui = renderUI({
    selectInput('combined_feature_2_selected_feature','Select Feature 2',combined_features()$features,combined_features()$features[2])
  })
  
  output$select_combined_group_ui = renderUI({
    selectInput('combined_feature_group','Select Group',combined_features()$groups,combined_features()$groups[1])
  })
  
  output$combined_feature_value_min_ui = renderUI({
    numericInput('combined_feature_value_min', 'Value Minimim Cutoff',0.0001)
  })
  
  output$combined_feature_1_threshold_ui = renderUI({
    numericInput('combined_feature_1_threshold', 'Feature 1 expression threshold',0.5)
  })
  output$combined_feature_2_threshold_ui = renderUI({
    numericInput('combined_feature_2_threshold', 'Feature 2 expression threshold',0.5)
  })
  
  
  combined_histogram_function = function(df,value_min = 0.01,
                                         selected_feature = 'IL13',
                                         group = 'psoriasis',
                                         threshold = 1){
    #value_min = 0.01
    #selected_feature = 'IL13'
    #group = 'psoriasis'
    var <- rlang::parse_quosures(group)[[1]] 
    var
   
    
    #threshold = 1
    dim(df)
    feature_data = df %>% filter(features == selected_feature, 
                                 value >= value_min,
                                 !is.na(!!var))
    dim(feature_data)
    as.tbl(feature_data)
    df_count = feature_data %>% 
      filter(value >= threshold) %>% 
      group_by(project,!!var) %>% 
        summarise(count = n()) %>% 
      ungroup()
    df_count
    gg_hist = ggplot(feature_data, aes(x = value,fill = project), aes_string(group = group)) +
      geom_histogram(position = 'dodge',binwidth = 0.2) + 
      facet_grid(. ~ features) +
      geom_vline(xintercept = threshold,col = 'red',size = 2) +
      ggtitle(group)
    if(threshold > 0){
      gg_hist = gg_hist + geom_vline(xintercept = threshold,col = 'red',size = 2)
      gg_hist = gg_hist + annotate('text',x = Inf,y = Inf, label = paste('value >=',input$vis_feature_threshold), vjust = 1, hjust = 1)
    }
    gg_hist
    text_df = df_count %>% filter(!!var)
    text_df
    for(i in c(1:dim(text_df)[1])){
      print(i)
      text = paste(text_df[i,c(1,3)],collapse = ' = ')
      text
      print(text)
      gg_hist = gg_hist + annotate('text',x = Inf,y = Inf, label = 'TRUE', vjust = 2.5, hjust = 8)
      gg_hist = gg_hist + annotate('text',x = Inf,y = Inf, label = text, vjust = i*1.25+3, hjust = 3)
    }
    text_df = df_count %>% filter(!(!!var))
    text_df
    for(i in c(1:dim(text_df)[1])){
      print(i)
      text = paste(text_df[i,c(1,3)],collapse = ' = ')
      text
      print(text)
      gg_hist = gg_hist + annotate('text',x = Inf,y = Inf, label = 'FALSE', vjust = 2.5, hjust = 1.5)
      gg_hist = gg_hist + annotate('text',x = Inf,y = Inf, label = text, vjust = i*1.25+3, hjust = 1)
    }
    gg_hist
    
    gg_boxplot  = ggplot(feature_data,aes(x = !!var,y = value,fill = project)) + 
      geom_violin() +
      geom_dotplot(aes(y = value, group = project, x = !!var, fill = project),
                   binaxis = "y", stackdir = "center", position = 'dodge',size = 0.1) +
      #geom_point(aes(col = project,group = project)) +
      ggtitle(group)
    gg_boxplot
    if(threshold > 0){
      gg_boxplot = gg_boxplot + 
        facet_grid(. ~ features) +
        geom_hline(yintercept = threshold, col = 'red', size = 2)
    }
    gg_boxplot
    ggsave(paste0('www/plots/',selected_feature,'_',group,'_combined_boxplot.png'),gg_boxplot,device = 'png')
    ggsave(paste0('www/plots/',selected_feature,'_',group,'_combined_hist.png'),gg_hist,device = 'png')
    df_count
    (min_true = min((df_count %>% filter((!!var)) %>% pull(count)),na.rm = T))
    (max_false = max((df_count %>% filter(!(!!var)) %>% pull(count)),na.rm = T))
    if(min_true > 15 & max_false <= 5){
      print(paste('Saving ',selected_feature,'positive cell data'))
      save_positive_cells_function(feature_data,df_count,group,selected_feature)
    }else{
      print(paste(selected_feature,'DOES NOT MEET THE CRITERIA'))
    }
  
    
    return(list(gg_hist = gg_hist,gg_boxplot = gg_boxplot,df_count = df_count, feature_data = feature_data))
  }
  
  
  
  combined_feature_1 = reactive({
    value_min = 0.01
    selected_feature = 'IL17A'
    group = 'psoriasis'
    threshold = 1
    
    df = combined_features()$df
    group = input$combined_feature_group
    value_min = input$combined_feature_value_min
    
    selected_feature = input$combined_feature_1_selected_feature
    threshold = input$combined_feature_1_threshold
    
    plot_list = combined_histogram_function(df,value_min,selected_feature,group,threshold)
    plot_list
  })
  output$combined_feature_1_boxplot = renderPlot({
    combined_feature_1()$gg_boxplot
  })
  output$combined_feature_1_hist = renderPlot({
    combined_feature_1()$gg_hist
  })
  output$combined_feature_1_count_df = renderDataTable({
    combined_feature_1()$count_df
  })
  
  combined_feature_2 = reactive({
    value_min = 0.01
    selected_feature = 'IL13'
    group = 'psoriasis'
    threshold = 1
    
    df = combined_features()$df
    group = input$combined_feature_group
    value_min = input$combined_feature_value_min
    
    selected_feature = input$combined_feature_2_selected_feature
    threshold = input$combined_feature_2_threshold
    
    plot_list = combined_histogram_function(df,value_min,selected_feature,group,threshold)
    plot_list
  })
  output$combined_feature_2_boxplot = renderPlot({
    combined_feature_2()$gg_boxplot
  })
  output$combined_feature_2_hist = renderPlot({
    combined_feature_2()$gg_hist
  })
  


    
      
  save_positive_cells_function = function(feature_data,df_count,group,selected_feature){
      var <- rlang::parse_quosures(group)[[1]] 
      var
      #saveRDS(feature_data,positive_list_name)
      (positive_samples = unique(feature_data %>% filter(!!var) %>% pull(project)))
      sample = positive_samples[1]
      for(sample in positive_samples){
        df = readRDS(paste0('www/data/',sample,'.rds'))
        (positive_cells = feature_data %>% filter(project == sample) %>% pull(variables)) 
        df_reduced = df[,c('features',positive_cells)]
        dim(df_reduced)
        (group_path = paste0('www/positive_cells/',group,'/'))
        system(paste('mkdir',group_path))
        (feature_path = paste0(group_path,selected_feature,'/'))
        (cmd = paste0('mkdir ',feature_path))
        system(cmd)
        saveRDS(df_reduced,paste0(feature_path,sample,'.rds'))
        write.csv2(df_reduced,paste0(feature_path,sample,'.csv'))
        df_reduced[df_reduced == 0] = NA
        as.tbl(df_reduced)
        df_reduced_c = df_reduced[complete.cases(df_reduced),]
        dim(df_reduced_c)
        df_reduced_c$features
        write.csv2(df_reduced_c,paste0(feature_path,sample,'_expressed_in_all.csv'))
      
      }
    }

  
  
  output$combined_feature_2_count_df = renderDataTable({
    combined_feature_2()$count_df
  })

  #### Extract Genes ######
  
  full_long_data_combined = reactive({
    data_path = 'www/data_long/'
    df = combine_df_l_function(data_path)
    
  })
  
  extract_genes = reactive({
    
    positive_cells = feature_data %>% filter(!!var) %>% 
    pull(variables)
    positive_cells
    (min_true = min((df_count %>% filter((!!var)) %>% pull(count)),na.rm = T))
    (max_false = max((df_count %>% filter(!(!!var)) %>% pull(count)),na.rm = T))
    (positive_list_name = paste0('www/positive_cells/',selected_feature,'_',threshold,'_',length(positive_cells),'.rds'))
    if(min_true > 15 & max_false <= 5){
      saveRDS(positive_cells,positive_list_name)
    }
    data_path = 'www/data/'
    file_list = list.files(data_path)
    file_list
    file_name = file_list[1]
    for(file_name in file_list){
      print(file_name)
      df = readRDS(paste0(data_path,file_name))
      dim(df)
      colnames(df)
      df$rownames = df$features
      df$rownames
      
      df$features = NULL
      colnames(df)
      df_reduce = df[,colnames(df)[as.numeric(df[1,]) > 0]]
      df_reduce > 0
      dim(df_reduce)
      df_reduce
    }
  })
  
  
})
