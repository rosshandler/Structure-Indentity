library(scran)
library(irlba)
library(Seurat)
library(ggplot2)

#library(umap)
#library(reticulate)
#use_condaenv(condaenv="scanpy-p3.9")

#umap = import('umap')

condition_colours <- c(
  "Kit"="#e12f2f",
  "Fast"="#6889ff",
  "Diss"="#9a9100",
  "Unembed"="#00733a"
)

sample_order  <- c(
 "B6Kit", "B7Kit", "B8Kit",
 "B6Fast", "B7Fast", "B8Fast",
 "B6diss", "B7diss", "B8diss",
 "B6Unembed", "B8Unembed_A", "B8Unembed_B")

setwd("/data1/ivanir/Ilaria2021/data/")
sce <- readRDS("sce.rds")

sce <- sce[calculateAverage(sce)>0.01, colData(sce)$scDblFinder.class == "singlet"]

sce <- logNormCounts(sce)

write.table(colData(sce), file="metadata_qc_doublet_removal.tab", sep="\t", row.names = FALSE)
write.table(as.matrix(t(logcounts(sce))), file="normalised_counts_qc_doublet_removal.tab", row.names=FALSE, col.names=FALSE, sep="\t")
writeLines(colnames(sce), "barcodes_qc_doublet_removal.txt")
writeLines(colData(sce)$cell, "barcodes_qc_doublet_removal_loom.txt")
writeLines(rownames(sce), "gene_ids_loom_filtered.txt")

#decomp  <- modelGeneVar(sce)
#hvgs    <- rownames(decomp)[decomp$FDR < 0.5]
#pca     <- prcomp_irlba(t(logcounts(sce[hvgs,])), n = 30)
#layout  <- umap(pca$x, method="umap-learn", umap_learn_args=c("n_neighbors", "n_epochs"), n_neighbors=15)

meta <- data.frame(colData(sce))

meta$sample <- factor(meta$sample, levels=sample_order)

#########################
### Seurat integration
setwd('/data1/ivanir/Ilaria2021/data')

seurat_integ <- CreateSeuratObject(as.matrix(counts(sce)), meta.data = meta)
seurat_list  <- SplitObject(seurat_integ, split.by = "sample")

# normalize and identify variable features for each dataset independently
seurat_list <- lapply(X = seurat_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration run PCA on each
# dataset using these features
features <- SelectIntegrationFeatures(object.list = seurat_list)
seurat_list <- lapply(X = seurat_list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors  <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features, reduction = "rpca")

combined <- IntegrateData(anchorset = anchors)

DefaultAssay(combined) <- "integrated"

combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)

pca_corrected <- combined@reductions$pca@cell.embeddings
colnames(pca_corrected) <- paste0("PC",1:30)
write.table(pca_corrected, file="pca_corrected_seurat.csv", row.names=FALSE, col.names=TRUE, sep=",")

umap <- combined@reductions$umap@cell.embeddings
colnames(umap) <- c("UMAP1", "UMAP2")
write.table(umap, file="umap_seurat.csv", row.names=FALSE, col.names=TRUE, sep=",")

p1 <- DimPlot(combined, reduction = "umap", group.by = "sequencing.round")
p2 <- DimPlot(combined, reduction = "umap", group.by = "sample")
p3 <- DimPlot(combined, reduction = "umap", group.by = "condition", cols = condition_colours)

p1 + p2 + p3

FeaturePlot(combined, features = c("TTR","EOMES","CLIC6"))

