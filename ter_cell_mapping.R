library(Seurat)
library(swne)
library(perturbLM)
library(glmnet)
library(stringr)
library(plyr)
library(dplyr)

ReadGenotypes <- function (file.name) {
  geneset.to.use <- as.list(readLines(file.name))
  geneset.to.use <- lapply(geneset.to.use, function(v) strsplit(v, ",")[[1]])
  geneset.names <- unlist(lapply(geneset.to.use, function(x) x[[1]]))
  geneset.to.use <- lapply(geneset.to.use, function(v) v[2:length(v)])
  names(geneset.to.use) <- geneset.names
  geneset.to.use <- lapply(geneset.to.use, function(v) gsub('"',"", v))
  return(geneset.to.use)
}

## Map teratoma screen cells to reference clusters
input.seurat.file <- "teratoma_integrated.rds"
ref.seurat.file <- "ref_ter.Robj"

## Load data
ref <- readRDS(ref.seurat.file)
ter <- readRDS(input.seurat.file)

DefaultAssay(ref) <- "integrated"
DefaultAssay(ter) <- "integrated"

## Remove genome identifier from gene names
new.rownames <- sapply(rownames(ter), function(x) gsub("hg19-","",x))
names(new.rownames) <- NULL
rownames(ter@assays$integrated@data) <- new.rownames

head(rownames(ter))
length(intersect(rownames(ref),rownames(ter)))

## Load metadata
meta.data <- read.table("~/cellmapper/ref_ter_metadata.tsv", header = T, sep = "\t")
ref.ident <- meta.data$cluster; names(ref.ident) <- meta.data$cell;
table(ref.ident)

## Classify cells with kNN classifier
ter.anchors <- FindTransferAnchors(reference = ref, query = ter, 
                                        dims = 1:20)
predictions <- TransferData(anchorset = ter.anchors, refdata = ref.ident, 
                            dims = 1:20)

#Prediction score assignment
orig.ident <- ter$integrated_snn_res.0.8
pred.ident.scores <- as.matrix(predictions[,grepl("prediction.score", colnames(predictions))])
pred.ident.scores <- pred.ident.scores[,1:(ncol(pred.ident.scores) - 1)]
pred.ident.scores <- pred.ident.scores[names(orig.ident),]

cluster.pred.ident.scores <- apply(pred.ident.scores, 2, function(x) {
  tapply(x, orig.ident, mean)
})


#Cell type assignment
cluster.celltype <- colnames(cluster.pred.ident.scores)[apply(cluster.pred.ident.scores,1,which.max)]
cluster.celltype <- gsub("prediction.score.","",cluster.celltype)
current.cluster.ids <- levels(ter$integrated_snn_res.0.8)
ter$cluster_celltype <- ter$integrated_snn_res.0.8
ter$cluster_celltype <- plyr::mapvalues(x = ter$cluster_celltype, from = current.cluster.ids, to = cluster.celltype)

## Assign predicted clusters
pred.ident.clusters <- predictions$predicted.id;
names(pred.ident.clusters) <- rownames(predictions)
ter.ident <- ter$cluster_celltype
ter <- AddMetaData(ter, metadata = ter.ident, col.name = "cluster_celltype")
Idents(ter) <- "cluster_celltype"

pdf("~teratoma_orf_celltypes.pdf", width = 8, height = 6)
DimPlot(object = ter, reduction = "umap", label = TRUE, repel = TRUE) & NoLegend() & NoAxes()
dev.off()

saveRDS(ter, file = "teratoma_integrated.rds")
