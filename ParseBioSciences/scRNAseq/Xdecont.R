library(scran)
library(celda)
library(irlba)
library(scater)
library(Matrix)
library(ggplot2)
library(biomaRt)
library(viridis)

library(umap)
library(leiden)
library(reticulate)
use_condaenv(condaenv="scanpy-p3.9")

umap = import('umap')

############################################
######## Xdecont exploration ###############
path2data <- "/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/all-well/DGE_unfiltered/"
markers   <- readLines("/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/all-well/DGE_unfiltered/geneMarkers.txt")

sce <- readRDS(paste0(path2data,"sce.rds"))
sce <- logNormCounts(sce)

sce_hvgs <- sce[calculateAverage(sce)>0.01,]
decomp   <- modelGeneVar(sce_hvgs)
hvgs     <- rownames(decomp)[decomp$FDR < 0.01]
markers <- intersect(markers, rownames(sce))
hvgs    <- union(hvgs, markers)
pca     <- prcomp_irlba(t(logcounts(sce[hvgs,])), n = 30)
rownames(pca$x) <- colnames(sce)

graph <- buildSNNGraph(pca$x, d = NA, transposed = TRUE)
set.seed(42)
clusters <- leiden(graph, resolution_parameter = 2)
names(clusters) <-  colData(sce)$cell

sce <- decontX(sce, z=clusters)

saveRDS(sce, paste0(path2data,"sce_decontX.rds"))

setwd('/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/QC')

pdf("hist_fraction_contamination.pdf")
hist(colData(sce)$decontX_contamination, main="", xlab="Fraction of contamination")
abline(v=median(colData(sce)$decontX_contamination), col="red")
dev.off()

pdf("hist_doublet_score_before_qc.pdf")
hist(colData(sce)$doublet_score, main="", xlab="Doublet score")
dev.off()

pdf("hist_doublet_score_after_qc.pdf")
hist(colData(sce[,colData(sce)$doublet_class == "singlet"])$doublet_score, main="", xlab="Doublet score")
dev.off()

sce <- logNormCounts(sce, exprs_values = "decontXcounts", name = "decontXlogcounts")
    
sce <- sce[,colData(sce)$doublet_class == "singlet"]

gexp <- as.matrix(decontXlogcounts(sce))
rownames(gexp) <- rowData(sce)$gene_name

calculateAverage(sce,assay.type = "decontXcounts")

sce_filt <- sce[rowMeans(decontXcounts(sce))>0.01, colData(sce)$doublet_class == "singlet"]

decomp  <- modelGeneVar(sce_filt, assay.type = "decontXlogcounts")
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



