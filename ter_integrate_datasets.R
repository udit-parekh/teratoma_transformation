library(Seurat)
library(ggplot2)
library(cowplot)

ter.list = list()
ter.names <- c("M1", "M2", "M3", "M4")

for (i in seq_along(ter.names)) {
  #Read in doublets
  df.doublets <- read.csv(paste0(ter.names[i],"/",ter.names[i],"_doublets.csv"))
  rownames(df.doublets) <- df.doublets$X
  df.doublets$X <- NULL
  
  #Read in cellranger metadata
  metadata <- read.csv(paste0(ter.names[i],"/outs/analysis/gem_classification.csv"), row.names = 1, header = TRUE)
  
  #Add sample name column
  metadata$sample <- rep(paste0(ter.names[i],"-"), nrow(metadata))
  
  #Remove hg19 and mm10 columns
  metadata$hg19 <- NULL
  metadata$mm10 <- NULL
  
  #Remove "-1" from rownames, retaining only cell barcode sequence
  rownames(metadata) = gsub(pattern = "-1", replacement = "", x = rownames(metadata))
  
  #Read in 10X data
  cancer_orfs.data <- Read10X(data.dir = paste0(ter.names[i],"/outs/filtered_feature_bc_matrix/"))
  colnames(cancer_orfs.data) = gsub(pattern = "-1", replacement = "", x = colnames(cancer_orfs.data))
  
  #Subset human genes only
  human.genes <- grep(pattern = "hg19", x = rownames(x = cancer_orfs.data), value = TRUE)
  human.cells <- subset(metadata, call == "hg19")
  cancer_orfs.data.human <- cancer_orfs.data[human.genes,]
  cancer_orfs.data.human <- cancer_orfs.data.human[,colnames(cancer_orfs.data.human) %in% rownames(human.cells)]
  print(ncol(cancer_orfs.data.human))
  
  #Create seurat object, add sample name to cell names and retain only human cells and singlets 
  print(round(ncol(cancer_orfs.data.human)*0.001))
  cancer_orfs_human<-CreateSeuratObject(counts = cancer_orfs.data.human,
                                        min.cells = round(ncol(cancer_orfs.data.human)*0.001),
                                        min.features = 200)
  cancer_orfs_human <- RenameCells(cancer_orfs_human, add.cell.id = ter.names[i])
  #cancer_orfs_human <- subset(x=cancer_orfs.human, subset = call == "hg19")
  cancer_orfs_human <- AddMetaData(cancer_orfs_human, df.doublets, col.name = "df.doublets")
  
  ter.list[[i]] <- subset(x = cancer_orfs_human, subset = df.doublets =="Singlet")
  
  #Filter out cells with high mitochondrial gene fraction and low feature counts
  mito.features <- grep(pattern = "^hg19-MT-", x = rownames(x = ter.list[[i]]), value = TRUE)
  percent.mito <- Matrix::colSums(x = GetAssayData(object = ter.list[[i]], slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = ter.list[[i]], slot = 'counts'))
  ter.list[[i]][['percent.mito']] <- percent.mito
  ter.list[[i]] <- subset(x = ter.list[[i]], subset = nFeature_RNA > 200 & percent.mito < 0.2)
  
  #Normalize data using log normalization and total feature counts
  ter.list[[i]] <- NormalizeData(object=ter.list[[i]], normalization.method = "LogNormalize", scale.factor = 1e4)
  
  #Find variable features
  ter.list[[i]] <- FindVariableFeatures(object=ter.list[[i]], selection.method = 'vst', mean.cutoff = c(0, Inf), dispersion.cutoff = c(0.5, Inf), nfeatures = 4000)
  print(length(x = VariableFeatures(object = ter.list[[i]])))
}

#saveRDS(ter.list, file = "teratoma_list.rds")

#Find anchors, batch correct and teratoma integrate datasets
teratoma.anchors <- FindIntegrationAnchors(object.list = ter.list, anchor.features = 4000, dims = 1:50)
teratoma_orf.integrated <- IntegrateData(anchorset = teratoma.anchors, dims = 1:50)
saveRDS(teratoma_orf.integrated, file = "teratoma_integrated.rds")

DefaultAssay(teratoma_orf.integrated) <- "integrated"

#Scale data and regress out library depth
teratoma_orf.integrated <- ScaleData(object = teratoma_orf.integrated, features = rownames(x = teratoma_orf.integrated), vars.to.regress = c("nCount_RNA", "percent.mito"))

#Run PCA
teratoma_orf.integrated <- RunPCA(object = teratoma_orf.integrated, features = VariableFeatures(object = teratoma_orf.integrated),  npcs = 50, verbose = FALSE)
ElbowPlot(object = teratoma_orf.integrated, ndims = 50)

#Generate k-nearest neighbours graph and find clusters
teratoma_orf.integrated <- FindNeighbors(object = teratoma_orf.integrated, dims = 1:40, k.param = 10)

clusters.resolution <- 0.8
teratoma_orf.integrated <- FindClusters(object = teratoma_orf.integrated, resolution = clusters.resolution)

#Plot clusters with UMAP
teratoma_orf.integrated <- RunUMAP(object = teratoma_orf.integrated, reduction = "pca", dims = 1:40)
pdf("teratoma_orf_clusters.pdf")
DimPlot(object = teratoma_orf.integrated, reduction = "umap", label = TRUE)
dev.off()

saveRDS(teratoma_orf.integrated, file = "~/Teratoma_ORF/teratoma_integrated.rds")
