library(scran)
library(irlba)
library(scater)
library(Matrix)
library(ggplot2)

library(umap)
library(leiden)
library(reticulate)
use_condaenv(condaenv="scanpy-p3.9")

sce  <- readRDS("/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/all-well/DGE_unfiltered/sce_transferred_annot.rds")
sce  <- logNormCounts(sce)

sce_hvgs <- sce[calculateAverage(sce)>0.01,]
decomp   <- modelGeneVar(sce_hvgs)
hvgs     <- rownames(decomp)[decomp$FDR < 0.0005]
writeLines(hvgs,"/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/all-well/DGE_unfiltered/scranHVGsFDR.01.txt")
markers <- readLines("/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/all-well/DGE_unfiltered/geneMarkers.txt")
markers <- intersect(markers, rownames(sce))
hvgs    <- union(hvgs, markers)
pca     <- prcomp_irlba(t(logcounts(sce[hvgs,])), n = 30)
rownames(pca$x) <- colnames(sce)
layout  <- umap(pca$x, method="umap-learn", umap_learn_args=c("n_neighbors", "n_epochs", "min_dist"), n_neighbors=30, min_dist=.25)

graph <- buildSNNGraph(pca$x, d = NA, transposed = TRUE)
set.seed(42)
clusters <- leiden(graph, resolution_parameter = 2)
names(clusters) <-  colData(sce)$cell

df_plot <- data.frame(
 colData(sce),
 leiden = clusters,
 UMAP1 = layout$layout[,1],
 UMAP2 = layout$layout[,2] 
)

writeMM(t(logcounts(sce[hvgs,])), "/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/sctour/norm_counts_hvgs.mtx", sep="\t")
writeLines(colnames(sce), "/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/sctour/cells_norm_hvgs.txt")
writeLines(hvgs, "/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/sctour/genes_norm_hvgs.txt")
