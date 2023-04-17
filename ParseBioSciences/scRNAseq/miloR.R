library(miloR)
library(scran)
library(irlba)
library(scater)
library(Matrix)
library(ggplot2)

library(patchwork)

path2data   <- '/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/all-well/DGE_unfiltered/'
sce <- readRDS(paste0(path2data,"sce.rds"))
sce <- sce[,colData(sce)$doublet_class == "singlet"]

setwd("/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/sctour")

umap <- read.csv("umap_layout.csv", header = FALSE)
colnames(umap) <- c("UMAP1","UMAP2")
rownames(umap) <- colnames(sce)
meta <- read.csv("metadata_scanpy.csv", header = TRUE)[,-1]

sce_filt <- sce[calculateAverage(sce)>0.01,]
sce_filt <- logNormCounts(sce_filt)

decomp  <- modelGeneVar(sce_filt)
hvgs    <- rownames(decomp)[decomp$FDR < 0.1]
length(hvgs)
#[1] 783

pca <- prcomp_irlba(t(logcounts(sce_filt[hvgs,])), n = 30)
rownames(pca$x) <- colnames(sce_filt)
reducedDim(sce_filt, "PCA") <- as.matrix(pca$x)

setwd("/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/miloR")

milo.obj <- Milo(sce_filt)

milo.obj <- buildGraph(milo.obj, k=20, d=30)
milo.obj <- makeNhoods(milo.obj, k=20, d=30, refined=TRUE, prop=0.2)
plotNhoodSizeHist(milo.obj)

milo.obj <- calcNhoodDistance(milo.obj, d=30)
milo.obj <- countCells(milo.obj, samples="sample_name", meta.data=meta)

reducedDim(milo.obj, "UMAP") <- as.matrix(umap) 
saveRDS(milo.obj,"MiloRscRNASeqObject.rds")

milo.design <- as.data.frame(xtabs(~ condition + sample_name, data=meta))
milo.design <- milo.design[milo.design$Freq > 0, ]
rownames(milo.design) <- milo.design$sample_name

table(meta$condition)
#CTRL_45 CTRL_55 CTRL_70 DISS_48 DISS_55 DISS_70   
#   7274    7291    2042    1460    8824    7312    4875

milo.res0 <- testNhoods(milo.obj, design=~condition+0, design.df=milo.design,
  model.contrasts="conditionCTRL_55 - conditionCTRL_45")
head(milo.res0[order(milo.res0$PValue, decreasing=FALSE),])

milo.res1 <- testNhoods(milo.obj, design=~condition+0, design.df=milo.design,
  model.contrasts="conditionDISS_48 - conditionCTRL_45")

milo.res2 <- testNhoods(milo.obj, design=~condition+0, design.df=milo.design,
  model.contrasts="conditionDISS_48 - conditionCTRL_55")

milo.res3 <- testNhoods(milo.obj, design=~condition+0, design.df=milo.design,
  model.contrasts="conditionDISS_55 - conditionCTRL_45")

milo.res4 <- testNhoods(milo.obj, design=~condition+0, design.df=milo.design,
  model.contrasts="conditionDISS_55 - conditionCTRL_55")

milo.res5 <- testNhoods(milo.obj, design=~condition+0, design.df=milo.design,
  model.contrasts="conditionEMB_55 - conditionCTRL_55")

milo.res6 <- testNhoods(milo.obj, design=~condition+0, design.df=milo.design,
  model.contrasts="conditionDISS_48 - conditionDISS_55")

milo.res7 <- testNhoods(milo.obj, design=~condition+0, design.df=milo.design,
  model.contrasts="conditionCTRL_70 - conditionCTRL_55")

milo.res8 <- testNhoods(milo.obj, design=~condition+0, design.df=milo.design,
  model.contrasts="conditionDISS_70 - conditionDISS_48")

milo.res9 <- testNhoods(milo.obj, design=~condition+0, design.df=milo.design,
  model.contrasts="conditionDISS_70 - conditionDISS_55")

milo.res10 <- testNhoods(milo.obj, design=~condition+0, design.df=milo.design,
  model.contrasts="conditionDISS_70 - conditionCTRL_70")

milo.obj <- buildNhoodGraph(milo.obj)


p0 <- plotNhoodGraphDA(milo.obj, milo.res0, alpha=0.5) +
  plot_layout(guides="collect") + ggtitle("CTRL_55 - CTRL_45")

p1 <- plotNhoodGraphDA(milo.obj, milo.res1, alpha=0.5) +
  plot_layout(guides="collect") + ggtitle("DISS_48 - CTRL_45")

p2 <- plotNhoodGraphDA(milo.obj, milo.res2, alpha=0.5) +
  plot_layout(guides="collect") + ggtitle("DISS_48 - CTRL_55")

p3 <- plotNhoodGraphDA(milo.obj, milo.res3, alpha=0.5) +
  plot_layout(guides="collect") + ggtitle("DISS_55 - CTRL_45")

p4 <- plotNhoodGraphDA(milo.obj, milo.res4, alpha=0.5) +
  plot_layout(guides="collect") + ggtitle("DISS_55 - CTRL_55")

p5 <- plotNhoodGraphDA(milo.obj, milo.res5, alpha=0.5) +
  plot_layout(guides="collect") + ggtitle("EMB_55 - CTRL_55")

p6 <- plotNhoodGraphDA(milo.obj, milo.res6, alpha=0.5) +
  plot_layout(guides="collect") + ggtitle("DISS_48 - DISS_55")

p7 <- plotNhoodGraphDA(milo.obj, milo.res7, alpha=0.5) +
  plot_layout(guides="collect") + ggtitle("CTRL_70 - CTRL_55")

p8 <- plotNhoodGraphDA(milo.obj, milo.res8, alpha=0.5) +
  plot_layout(guides="collect") + ggtitle("DISS_70 - DISS_48")

p9 <- plotNhoodGraphDA(milo.obj, milo.res9, alpha=0.5) +
  plot_layout(guides="collect") + ggtitle("DISS_70 - DISS_55")

p10 <- plotNhoodGraphDA(milo.obj, milo.res10, alpha=0.5) +
  plot_layout(guides="collect") + ggtitle("DISS_70 - CTRL_70")

pdf("miloR.pdf",width=24,height=16)
egg::ggarrange(p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, nrow = 3, ncol = 4)
dev.off()