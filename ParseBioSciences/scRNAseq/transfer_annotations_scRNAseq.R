library(scran)
library(Seurat)
library(ggplot2)

library(irlba)
library(Rtsne)

library(umap)
library(reticulate)
use_condaenv(condaenv="scanpy-p3.9")

path2data    <- '/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/all-well/DGE_unfiltered/'
path2data10x <- '/data1/ivanir/Ilaria2021/UCSC_data/'

sce_pb <- readRDS(paste0(path2data, 'sce.rds'))
sce_pb <- sce_pb[,colData(sce_pb)$doublet_class == "singlet"]
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
        nfeatures = 500, verbose = TRUE)
}

reference_dataset <- which(names(seurat_list) == "reference")

seurat_query     <- seurat_list[["query"]]
seurat_reference <- seurat_list[["reference"]]
seurat_anchors <- FindTransferAnchors(reference = seurat_reference, query = seurat_query, 
    dims = 1:30)
predictions    <- TransferData(anchorset = seurat_anchors, refdata = colData(sce_10x)$celltype.manuscript, 
    dims = 1:30)
seurat_query   <- AddMetaData(object = seurat_query, metadata = predictions)

setwd("/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/sctour")

meta_scanpy <- read.csv("metadata_scanpy.csv", header = TRUE)[,-1]
colData(sce_pb) <- DataFrame(meta_scanpy)

meta <- cbind(colData(sce_pb), seurat_prediction=seurat_query$predicted.id, seurat_max.score=seurat_query$prediction.score.max)
colData(sce_pb) <- DataFrame(meta)

sce_filt <- sce_pb[calculateAverage(sce_pb)>0.01,]
sce_filt <- logNormCounts(sce_filt)

decomp  <- modelGeneVar(sce_filt)
hvgs    <- rownames(decomp)[decomp$FDR < 0.1]
length(hvgs)
#[1] 783

pca <- prcomp_irlba(t(logcounts(sce_filt[hvgs,])), n = 30)
rownames(pca$x) <- colnames(sce_filt)

meta_scanpy_dmap <- read.csv("metadata_scanpy_dmap_res2.5.csv", header = TRUE)[,-1]
colData(sce_pb)$leiden_dmap <- meta_scanpy_dmap$leiden_dmap

umap <- read.csv("umap_layout_dmap_res2.5.csv", header = FALSE)
colnames(umap) <- c("UMAP1","UMAP2")
rownames(umap) <- colnames(sce_filt)

reducedDim(sce_pb, "PCA")  <- as.matrix(pca$x)
reducedDim(sce_pb, "UMAP") <- umap

colnames(sce_pb) <- colData(sce_pb)$cell
saveRDS(sce_pb, paste0(path2data,'sce_transferred_annot.rds'))

##
## Plotting
##

df_plot <- data.frame(
 colData(sce_pb),umap
)

####################
setwd('/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/plots')

population_colours <- c(
"Mitotic RG" = "#005dd8",
"Cycling RG" = "#4aadd6",
"Differentiating RG" = "#6083d1",
"RG" = "#02d7ec",
"Glycolytic RG" = "#0086cb",
"High trancription RG" = "#4aac8d",
"Chp" = "#935de6",
"Cortical hem" = "#9e4d70",
"Cycling IPCs" = "#d74897",
"IPCs"   = "#ff8ba6",
"IN"="#228B22",
"Mature IN" = "#90b600",
"CR cells" = "#00ff00",
"Migrating excitatory neurons" = "#a6023e",
"UL neurons"   = "#fa2274",
"DL neurons"   = "#9e5d56",
"Mature excitatory neurons"="#d3b000")


pdf("transfer_label.pdf", width=12, height=8)
ggplot(df_plot, aes(x = UMAP1, y = UMAP2, col = factor(seurat_prediction))) +
  geom_point(size = 1) +
  #geom_point(size = 3, aes(alpha = seurat_max.score)) + scale_alpha("Mapping score") +
  scale_color_manual(values=population_colours, name = "Cell population mapped", labels = names(population_colours)) +
  theme_minimal() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=7)))
dev.off()

ggplot(df_plot, aes(x = UMAP1, y = UMAP2, col = factor(seurat_prediction))) +
  geom_point(size = 1) +
  #geom_point(size = 3, aes(alpha = seurat_max.score)) + scale_alpha("Mapping score") +
  scale_color_manual(values=population_colours, name = "Cell population mapped", labels = names(population_colours)) +
  theme_minimal() + theme(legend.position = "none") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=7))) +
  facet_wrap(~seurat_prediction)
ggsave("umap_split_transfer_label.pdf")

df_plot$seurat_prediction_cutoff <- rep("-",nrow(df_plot))
df_plot$seurat_prediction_cutoff[df_plot$seurat_max.score > .5] <- as.character(df_plot$seurat_prediction[df_plot$seurat_max.score > .5])
plot.index <- order(df_plot$seurat_prediction_cutoff)

pdf("transfer_label_mapscore_thr.pdf", width=12, height=8)
ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = factor(seurat_prediction_cutoff))) +
geom_point(size = 1) +
scale_color_manual(values=population_colours, name = "Cell population mapped", labels = names(population_colours)) +
theme_minimal() +
theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
guides(colour = guide_legend(override.aes = list(size=7)))
dev.off()

plot.index <- sample(1:nrow(df_plot))
ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = day)) +
geom_point(size = 1) +
#scale_color_manual(values=population_colours, name = "Cell population mapped", labels = names(population_colours)) +
theme_minimal() +
theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
guides(colour = guide_legend(override.aes = list(size=7)))
ggsave("umap_day.pdf")

plot.index <- sample(1:nrow(df_plot))
ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = day)) +
geom_point(size = 1) +
#scale_color_manual(values=population_colours, name = "Cell population mapped", labels = names(population_colours)) +
theme_minimal() + theme(legend.position = "none") +
theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
guides(colour = guide_legend(override.aes = list(size=7))) +
facet_wrap(~day)
ggsave("umap_split_day.pdf")
                      
plot.index <- sample(1:nrow(df_plot))
ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = condition)) +
geom_point(size = 1) +
#scale_color_manual(values=population_colours, name = "Cell population mapped", labels = names(population_colours)) +
theme_minimal() +
theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
guides(colour = guide_legend(override.aes = list(size=7)))
ggsave("umap_condition.pdf")

plot.index <- sample(1:nrow(df_plot))
ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = condition)) +
geom_point(size = 1) +
#scale_color_manual(values=population_colours, name = "Cell population mapped", labels = names(population_colours)) +
theme_minimal() + theme(legend.position = "none") +
theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
guides(colour = guide_legend(override.aes = list(size=7))) +
facet_wrap(~condition)
ggsave("umap_split_condition.pdf")

plot.index <- order(df_plot$mt.fraction)
ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = mt.fraction)) +
geom_point(size = 1) +
scale_color_gradient(low="gray", high="red") +
theme_minimal() +
theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
ggsave("umap_mtfraction.pdf")

plot.index <- order(df_plot$doublet_score)
ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = doublet_score)) +
geom_point(size = 1) +
scale_color_gradient(low="gray", high="black") +
theme_minimal() +
theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank())
ggsave("umap_doublet_score.pdf")






#### Deprecated-to be deleted
setwd('/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/scanpy/')

sce_pb_qc <- sce_pb[,colData(sce_pb)$doublet_class == "singlet"]
meta_qc   <- meta[colData(sce_pb)$doublet_class    == "singlet",]
write.table(as.matrix(logcounts(sce_pb_qc)),"normalised_counts_qc.tab", sep="\t", row.names=FALSE, quote=FALSE)
writeLines(colnames(sce_pb_qc),"cells_qc.txt")
writeLines(rownames(sce_pb_qc),"genes_qc.txt")
write.table(meta_qc,"cell_metadata_qc.tab", sep="\t", row.names=FALSE, quote=FALSE)

sce_ctrl  <- sce_pb[,grep("CTRL",meta$condition)]
meta_ctrl <- meta[grep("CTRL",meta$condition),]
meta_ctrl <- meta_ctrl[colData(sce_ctrl)$doublet_class == "singlet",]
sce_ctrl  <- sce_ctrl[,colData(sce_ctrl)$doublet_class == "singlet"]
write.table(as.matrix(logcounts(sce_ctrl)),"normalised_counts_ctrl.tab", sep="\t", row.names=FALSE, quote=FALSE)
writeLines(colnames(sce_ctrl),"cells_ctrl.txt")
writeLines(rownames(sce_ctrl),"genes_ctrl.txt")
write.table(meta_ctrl,"cell_metadata_ctrl.tab", sep="\t", row.names=FALSE, quote=FALSE)

sce_diss  <- sce_pb[,grep("DISS",meta$condition)]
sce_diss  <- sce_diss[,colData(sce_diss)$doublet_class == "singlet"]
meta_diss <- meta[grep("DISS",meta$condition),]
meta_diss <- meta_diss[colData(sce_diss)$doublet_class == "singlet",]
write.table(as.matrix(logcounts(sce_diss)),"normalised_counts_diss.tab", sep="\t", row.names=FALSE, quote=FALSE)
writeLines(colnames(sce_diss),"cells_diss.txt")
writeLines(rownames(sce_diss),"genes_diss.txt")
write.table(meta_diss,"cell_metadata_diss.tab", sep="\t", row.names=FALSE, quote=FALSE)

sce_emb  <- sce_pb[,grep("EMB",meta$condition)]
meta_emb <- meta[grep("EMB",meta$condition),]
meta_emb <- meta_emb[colData(sce_emb)$doublet_class == "singlet",]
sce_emb  <- sce_emb[,colData(sce_emb)$doublet_class == "singlet"]
write.table(as.matrix(logcounts(sce_emb)),"normalised_counts_emb.tab", sep="\t", row.names=FALSE, quote=FALSE)
writeLines(colnames(sce_emb),"cells_emb.txt")
writeLines(rownames(sce_emb),"genes_emb.txt")
write.table(meta_emb,"cell_metadata_emb.tab", sep="\t", row.names=FALSE, quote=FALSE)

######
sce_ctrl45vsdiss48  <- sce_pb[,
  colData(sce_pb)$doublet_class == "singlet" & (colData(sce_pb)$condition == "CTRL_45" | colData(sce_pb)$condition == "DISS_48")]
meta_ctrl45vsdiss48 <- meta[
  colData(sce_pb)$doublet_class == "singlet" & (colData(sce_pb)$condition == "CTRL_45" | colData(sce_pb)$condition == "DISS_48"),]

