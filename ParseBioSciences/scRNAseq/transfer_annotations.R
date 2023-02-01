library(scran)
library(Seurat)
library(ggplot2)

path2data    <- '/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/all-well/DGE_unfiltered/'
path2data10x <- '/data1/ivanir/Ilaria2021/UCSC_data/'

sce_pb  <- readRDS(paste0(path2data, 'sce.rds'))
colnames(sce_pb) <- paste0(colnames(sce_pb),'_pb')
rownames(sce_pb) <- rowData(sce_pb)$gene_name
sce_10x <- readRDS(paste0(path2data10x, 'sce2023.rds'))
colnames(sce_10x) <- paste0(colnames(sce_10x),'_10x')

shared_genes <- intersect(rownames(sce_pb), rownames(sce_10x))

sce <- SingleCellExperiment(list(counts=cbind(counts(sce_pb[shared_genes,]), counts(sce_10x[shared_genes,]))))
rownames(sce) <- shared_genes

batch <- c(rep("query", ncol(sce_pb)), rep("reference", ncol(sce_10x))) 

metadata <- data.frame(cbind(cell=c(colnames(sce_pb), colnames(sce_10x)), batch))
rownames(metadata) <- metadata$cell

### Seurat integration
seurat_integ <- CreateSeuratObject(cbind(counts(sce_pb[shared_genes,]), counts(sce_10x[shared_genes,])), meta.data = metadata)
seurat_list  <- SplitObject(seurat_integ, split.by = "batch")

for (i in 1:length(x = seurat_list)) {
   seurat_list[[i]] <- NormalizeData(object = seurat_list[[i]], verbose = TRUE)
   seurat_list[[i]] <- FindVariableFeatures(object = seurat_list[[i]], selection.method = "vst", 
        nfeatures = 2000, verbose = TRUE)
}

reference_dataset <- which(names(seurat_list) == "reference")

seurat_query     <- seurat_list[["query"]]
seurat_reference <- seurat_list[["reference"]]
seurat_anchors <- FindTransferAnchors(reference = seurat_reference, query = seurat_query, 
    dims = 1:30)
predictions    <- TransferData(anchorset = seurat_anchors, refdata = colData(sce_10x)$celltype.manuscript, 
    dims = 1:30)
seurat_query   <- AddMetaData(object = seurat_query, metadata = predictions)

meta <- cbind(colData(sce_pb), seurat_prediction=seurat_query$predicted.id, seurat_max.score=seurat_query$prediction.score.max)

saveRDS(meta, "transferred_annot_meta.rds")


