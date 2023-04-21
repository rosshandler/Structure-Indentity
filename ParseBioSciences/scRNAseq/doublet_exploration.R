library(scran)
library(scran)
library(irlba)
library(Rtsne)
library(Matrix)
library(ggplot2)
library(biomaRt)
library(viridis)
library(scDblFinder)

library(umap)
library(reticulate)
use_condaenv(condaenv="scanpy-p3.9")

umap = import('umap')

############################################
######## Doublet exploration ###############
path2data <- "/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/all-well/DGE_unfiltered/"

sce <- readRDS(paste0(path2data,"sce.rds"))
sce <- readRDS(paste0(path2data,"sce_stringent.rds"))

hist(colData(sce)$doublet_score)

hist(colData(sce[,colData(sce)$doublet_class == "singlet"])$doublet_score)

dim(sce)
#[1] 62704 41576
dim(sce[,colData(sce)$doublet_class == "singlet"])
#[1] 62704 39078
sum(colData(sce[,colData(sce)$doublet_class == "singlet"])$doublet_score < .8)
#[1] 38780
sum(colData(sce[,colData(sce)$doublet_class == "singlet"])$doublet_score < .6)
#[1] 37200
sum(colData(sce[,colData(sce)$doublet_class == "singlet"])$doublet_score < .5)
#[1] 36228

sce <- sce[,colData(sce)$doublet_class == "singlet"]

sce_filt <- sce[calculateAverage(sce)>0.01,colData(sce[,
  colData(sce)$doublet_class == "singlet"])$doublet_score < .1]

sce_filt <- logNormCounts(sce_filt)
gexp <- as.matrix(logcounts(sce_filt))
rownames(gexp) <- rowData(sce_filt)$gene_name

decomp  <- modelGeneVar(sce_filt)
hvgs    <- rownames(decomp)[decomp$FDR < 0.1]
length(hvgs)
#[1] 499
pca     <- prcomp_irlba(t(logcounts(sce_filt[hvgs,])), n = 30)
rownames(pca$x) <- colnames(sce_filt)
layout  <- umap(pca$x, method="umap-learn", umap_learn_args=c("n_neighbors", "n_epochs", "min_dist"), n_neighbors=30, min_dist=.25)

df_plot <- data.frame(
 colData(sce_filt),
 UMAP1 = layout$layout[,1],
 UMAP2 = layout$layout[,2] 
)

plotLayoutExpression <- function(gene="TTR"){
  require(Matrix)
  require(ggplot2)
    logcounts <- as.vector(as.matrix(gexp[gene,]))
    if (sum(logcounts)>0){
        df_tmp    <- data.frame(cbind(df_plot, logcounts))
        plot.index  <- order(df_tmp$logcounts)
        ggplot(df_tmp[plot.index,], aes(x = UMAP1, y = UMAP2, colour = logcounts)) + 
          geom_point(size = 1) +
          scale_color_gradient(low="gray", high="darkgreen") +
          labs(color = paste0(gene,"\nlog(counts)")) +
          theme_minimal() + 
          theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
          theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
          xlab("Dimension 1") + ylab("Dimension 2")
    }else{
    message(gene," was not detected in the expression matrix")
    }
}

#Chp
plotLayoutExpression("TTR")
plotLayoutExpression("CLIC6")
#oRG
plotLayoutExpression("TNC")
plotLayoutExpression("MOXD1")
plotLayoutExpression("PTPRZ1")
plotLayoutExpression("IL6ST")
plotLayoutExpression("FAM107A")
plotLayoutExpression("SLC1A3")
#IPCs
plotLayoutExpression("EOMES")
plotLayoutExpression("NEUROD4")
#INs
plotLayoutExpression("DLX6-AS1")
plotLayoutExpression("GAD1")
plotLayoutExpression("GAD2")
plotLayoutExpression("SST")
plotLayoutExpression("PVALB")
plotLayoutExpression("DLX2")
#UL
plotLayoutExpression("CUX1")
plotLayoutExpression("CUX2")
plotLayoutExpression("GLRA3")
plotLayoutExpression("NRF1")
plotLayoutExpression("RORB")



