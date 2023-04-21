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
######## Decontamination with Xdecont ######
setwd('/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/QC')

path2data <- "/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/all-well/DGE_unfiltered/"
markers   <- readLines("/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/all-well/DGE_unfiltered/geneMarkers.txt")

sce <- readRDS(paste0(path2data,"sce.rds"))
sce <- logNormCounts(sce)

pdf("hist_doublet_score_before_qc.pdf")
hist(colData(sce)$doublet_score, main="", xlab="Doublet score", ylab="Cells")
dev.off()

pdf("hist_doublet_score_after_qc.pdf")
hist(colData(sce[,colData(sce)$doublet_class == "singlet"])$doublet_score, main="", xlab="Doublet score", ylab="Cells")
dev.off()

sce_hvgs <- sce[calculateAverage(sce)>0.01,]
rownames(sce_hvgs) <- rowData(sce_hvgs)$gene_name
decomp   <- modelGeneVar(sce_hvgs)
hvgs     <- rownames(decomp)[decomp$FDR < 0.01]
markers  <- intersect(markers, rownames(sce_hvgs))
hvgs     <- union(hvgs, markers)
pca      <- prcomp_irlba(t(logcounts(sce_hvgs[hvgs,])), n = 30)
rownames(pca$x) <- colnames(sce_hvgs)

graph <- buildSNNGraph(pca$x, d = NA, transposed = TRUE)
set.seed(42)
clusters <- leiden(graph, resolution_parameter = 2)
names(clusters) <-  colData(sce_hvgs)$cell

sce <- decontX(sce, z=clusters)

saveRDS(sce, paste0(path2data,"sce_decontX.rds"))
metadata(sce)$decontX$estimates$all_cells$delta
#[1] 1.2352165 0.9965378

sce <- decontX(sce, z=clusters, delta = c(1, 10), estimateDelta = FALSE)
saveRDS(sce, paste0(path2data,"sce_decontX_stringent.rds"))

pdf("hist_fraction_contamination.pdf")
hist(colData(sce)$decontX_contamination, main="decontX", xlab="Fraction of contamination", ylab="Cells")
abline(v=median(colData(sce)$decontX_contamination), col="red")
dev.off()

lib.sizes <- colSums(decontXcounts(sce))
ngenes <- colSums(decontXcounts(sce) > 0)
sum(lib.sizes > 700)
#[1] 14198
sum(ngenes > 500)
#[1] 34318


sce_tmp <- sce[,colData(sce)$doublet_class == "singlet" & lib.sizes > 700 & ngenes > 500]
sce_tmp <- logNormCounts(sce_tmp, exprs_values = "decontXcounts", name = "decontXlogcounts")

saveRDS(sce_tmp, paste0(path2data,"sce_hardQC.rds"))

###########################################################################
### Inspect markers to check whether cell type separation has improved ####

gexp <- as.matrix(logcounts(sce_tmp,assay.type = "decontXcounts"))
rownames(gexp) <- rowData(sce_tmp)$gene_name

sce_filt <- sce_tmp[calculateAverage(sce_tmp, assay.type = "decontXcounts")>0.01,]
rownames(sce_filt) <- rowData(sce_filt)$gene_name

decomp  <- modelGeneVar(sce_filt, assay.type = "decontXlogcounts")
hvgs    <- rownames(decomp)[decomp$FDR < 0.01]
markers <- intersect(markers, rownames(sce_filt))
hvgs    <- union(hvgs, markers)
pca     <- prcomp_irlba(t(logcounts(sce_filt,assay.type = "decontXcounts")[hvgs,]), n = 30)
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


###################################################
######## Generate output for scanpy/sctour ########
#sce_tmp <- readRDS(paste0(path2data,"sce_hardQC.rds"))

setwd("/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/sctour")

writeMM(t(decontXcounts(sce_tmp)), "raw_counts.mtx")
writeLines(colnames(sce_tmp), "cells.txt")
writeLines(rownames(sce_tmp), "genes.txt")

meta <- cbind(cell=rownames(colData(sce_tmp)), colData(sce_tmp))
write.table(data.frame(meta), file="metadata.tab", sep="\t", row.names=FALSE, quote=FALSE)

setwd("/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/scanpy")

writeMM(t(logcounts(sce_filt,assay.type = "decontXcounts")), "norm_counts.mtx", sep="\t")
writeLines(colnames(sce_filt), "cells_norm.txt")
writeLines(rowData(sce_filt)$gene_name, "genes_norm.txt")

meta <- cbind(cell=rownames(colData(sce_filt)), colData(sce_filt))
write.table(data.frame(meta), file="metadata.tab", sep="\t", row.names=FALSE, quote=FALSE)

write.table(df_plot[,"UMAP1","UMAP2"], file="UMAP.tab", sep="\t", row.names=FALSE, quote=FALSE)

