library(Seurat)
library(dplyr)
library(shiny)
options(shiny.maxRequestSize=3000*1024^2) 


# Load the PBMC dataset
sysinf <- Sys.info()
os <- sysinf['sysname']
os = 'Darwin'
#setwd('/srv/shiny-server/sample-apps/Seural_shiny/')


data_list = list.files('data/')
datasets = c('immune_hsa','immune_mmu')
dataset_list = list('immune_hsa' = 'data/immune_hsa/outs/filtered_gene_bc_matrices_mex/GRCh38/',
                    'immune_mmu' = 'data/immune_mmu/outs/filtered_gene_bc_matrices_mex/mm10/')

#options(rsconnect.http = "curl")
output_dir = '../output/'
output_name = 'immune_hsa'
#output_name = 'immune_mmu'
#pbmc.data <- Read10X(data.dir = '../Seurat/filtered_gene_bc_matrices/hg19/')
#data <- Read10X(data.dir = 'data/GRCh38/')
#data <- Read10X(data.dir = 'data/aggr_for_Ray/immune_hsa/filtered_gene_bc_matrices_mex/GRCh38/')
#data <- Read10X(data.dir = 'data/aggr_for_Ray/immune_mmu/filtered_gene_bc_matrices_mex/mm10/)

#pbmc.data <- Read10X(data.dir = '../data/aggr_for_Ray/immune_mmu/outs/filtered_gene_bc_matrices_mex/mm10/')



# Examine the memory savings between regular and sparse matrices
#dense.size <- object.size(x = as.matrix(x = data))
#dense.size
#sparse.size <- object.size(x = data)
#sparse.size
#dense.size / sparse.size

#values = readRDS('www/values.rds')


FindMarkers_test_list = list(
  "Wilcoxon rank sum test (default)" =  "wilcox",
  "Likelihood-ratio test for single cell gene expression, (McDavid et al., Bioinformatics, 2013)" =  "wilcox",
  "Standard AUC classifier" = "roc",
  "Student's t-test" = "t",
  "Tobit-test for differential gene expression (Trapnell et al., Nature Biotech, 2014)" = "tobit",
  "Likelihood ratio test assuming an underlying poisson distribution. Use only for UMI-based datasets" = "poisson",
  "Likelihood ratio test assuming an underlying negative binomial distribution. Use only for UMI-based datasets" = "negbinom",
  "GLM-framework that treates cellular detection rate as a covariate (Finak et al, Genome Biology, 2015)" =  "MAST",
  "DE based on a model using the negative binomial distribution (Love et al, Genome Biology, 2014)" = "DESeq2"
)

process_data = F
if(process_data == T){
  print("Seurat data")
  print('read')
  Seurat_data <- CreateSeuratObject(raw.data = data, min.cells = 3, min.genes = 200, project = "immune_hsa")
  
  
  # The number of genes and UMIs (nGene and nUMI) are automatically calculated for every object by Seurat.
  # For non-UMI data, nUMI represents the sum of the non-normalized values within a cell
  # We calculate the percentage of mitochondrial genes here and store it in percent.mito using AddMetaData.
  # We use object@raw.data since this represents non-transformed and non-log-normalized counts
  # The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = Seurat_data@data), value = TRUE)
  percent.mito <- Matrix::colSums(Seurat_data@raw.data[mito.genes, ]) / Matrix::colSums(Seurat_data@raw.data)
  
  # AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats
  Seurat_data <- AddMetaData(object = Seurat_data, metadata = percent.mito, col.name = "percent.mito")
  
  VlnPlot(object = Seurat_data, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
  
  
  # GenePlot is typically used to visualize gene-gene relationships, but can be used for anything 
  # calculated by the object, i.e. columns in object@meta.data, PC scores etc.
  # Since there is a rare subset of cells with an outlier level of high mitochondrial percentage
  # and also low UMI content, we filter these as well
  par(mfrow = c(1, 2))
  GenePlot(object = Seurat_data, gene1 = "nUMI", gene2 = "percent.mito")
  GenePlot(object = Seurat_data, gene1 = "nUMI", gene2 = "nGene")
  
  
  # We filter out cells that have unique gene counts over 2,500 or less than 200
  # Note that low.thresholds and high.thresholds are used to define a 'gate'.
  # -Inf and Inf should be used if you don't want a lower or upper threshold.
  par(mfrow = c(1, 1))
  
  Seurat_data <- FilterCells(object = Seurat_data, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))
  
  
  Seurat_data <- NormalizeData(object = Seurat_data, normalization.method = "LogNormalize", scale.factor = 1e4)
  
  Seurat_data <- FindVariableGenes(object = Seurat_data, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
  
  
  
  #Seurat_data
  
  
  
  print(length(x = Seurat_data()@var.genes))
  
  Seurat_data <- ScaleData(object = Seurat_data, vars.to.regress = c("nUMI", "percent.mito"))
  
  Seurat_data <- RunPCA(object = Seurat_data, pc.genes = Seurat_data@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
  
}
