library(miloR)
library(scran)
library(scater)
library(Matrix)
library(ggplot2)

library(patchwork)

setwd("/data1/ivanir/Ilaria2021/data")

layout_umap <- read.csv("umap_layout_batch_corrected_contamination_cleaned.csv", header = TRUE)

meta_clus   <- read.csv("meta_annotated_updated.csv", header = TRUE)

sce <- readRDS("sce.rds")
sce <- sce[calculateAverage(sce)>0.05, colData(sce)$scDblFinder.class == "singlet"]
sce <- logNormCounts(sce)

pca <- read.csv("pca_corrected_seurat.csv")
rownames(pca) <- colnames(sce)

sce  <- sce[,meta_clus$celltype_annotation != "High metabolism/protein translation"]
pca  <- pca[meta_clus$celltype_annotation != "High metabolism/protein translation",]
meta <- meta_clus[meta_clus$celltype_annotation != "High metabolism/protein translation",]

pca_morph <- readRDS("morphology_pca.rds")
rownames(pca_morph$x) <- c(
  "B6Kit", "B7Kit", "B8Kit",
  "B6Fast", "B7Fast", "B8Fast",
  "B6Unembed", "B8Unembed_A", "B8Unembed_B",
  "B6diss", "B7diss", "B8diss")

set.seed(14)
clustering <- kmeans(pca_morph$x[,1:2], centers = 2, nstart = 25) # same result if three PCs are used

pdf("kmeans2_morphology.pdf")
fviz_cluster(clustering, data = pca$x[,1:2],
             palette = c("#2E9FDF", "#E7B800"), 
             geom = "point",
             ellipse.type = "convex", 
             ggtheme = theme_bw()
)
dev.off()

meta$morphology_clus <- clustering$cluster[as.character(meta$sample)]

meta <- data.frame(
 meta,
 UMAP1 = layout_umap[,1],
 UMAP2 = layout_umap[,2] 
)

setwd("/data1/ivanir/Ilaria2021/data/milo_plots")

reducedDim(sce, "PCA") <- as.matrix(pca)
reducedDims(sce)

milo.obj <- Milo(sce)

milo.obj <- buildGraph(milo.obj, k=20, d=30)
milo.obj <- makeNhoods(milo.obj, k=20, d=30, refined=TRUE, prop=0.2)

pdf("nhood_size_hist.pdf")
plotNhoodSizeHist(milo.obj)
dev.off()

milo.obj <- calcNhoodDistance(milo.obj, d=30)
milo.obj <- countCells(milo.obj, samples="sample", meta.data=meta)

milo.design <- as.data.frame(xtabs(~ condition + sample, data=meta))
milo.design <- milo.design[milo.design$Freq > 0, ]
rownames(milo.design) <- milo.design$sample

milo.res1 <- testNhoods(milo.obj, design=~condition+0, design.df=milo.design,
  model.contrasts="conditionFast - conditionKit")
head(milo.res1[order(milo.res1$PValue, decreasing=FALSE),])

milo.res2 <- testNhoods(milo.obj, design=~condition+0, design.df=milo.design,
  model.contrasts="conditionUnembed - conditionKit")
head(milo.res2[order(milo.res2$PValue, decreasing=FALSE),])

milo.res3 <- testNhoods(milo.obj, design=~condition+0, design.df=milo.design,
  model.contrasts="conditionDiss - conditionKit")
head(milo.res3[order(milo.res3$PValue, decreasing=FALSE),])

milo.obj <- buildNhoodGraph(milo.obj)

layout <- meta[, c("UMAP1","UMAP2")]
rownames(layout) <- colnames(sce)

reducedDim(milo.obj, "UMAP") <- as.matrix(layout) 

pdf("milo_noMG.pdf",width=12,height=6)
#plotUMAP(milo.obj) + plotNhoodGraphDA(milo.obj, milo.res2, alpha=0.05) +
plotNhoodGraphDA(milo.obj, milo.res2, alpha=0.5) +
  plot_layout(guides="collect")
dev.off()
