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

sce_pb  <- readRDS(paste0(path2data, 'sce_stringent.rds'))
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

meta <- cbind(colData(sce_pb), seurat_prediction=seurat_query$predicted.id, seurat_max.score=seurat_query$prediction.score.max)

saveRDS(meta, paste0(path2data,'transferred_annot_meta_stringent.rds'))

sce_pb <- logNormCounts(sce_pb)

sce_pb <- sce_pb[calculateAverage(sce_pb)>0.01,]

decomp  <- modelGeneVar(sce_pb)
hvgs    <- rownames(decomp)[decomp$FDR < 0.5]
pca     <- prcomp_irlba(t(logcounts(sce_pb[hvgs,])), n = 30)
rownames(pca$x) <- colnames(sce_pb)
tsne    <- Rtsne(pca$x, pca = FALSE, check_duplicates = FALSE, num_threads=30)
layout  <- umap(pca$x, method="umap-learn", umap_learn_args=c("n_neighbors", "n_epochs", "min_dist"), n_neighbors=30, min_dist=.25)

df_plot <- data.frame(
 meta,
 tSNE1    = tsne$Y[, 1],
 tSNE2    = tsne$Y[, 2], 
 UMAP1 = layout$layout[,1],
 UMAP2 = layout$layout[,2] 
)

####################
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

setwd('/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/plots')

pdf("transfer_label_stringent.pdf", width=12, height=8)
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
ggsave("umap_split_transfer_label_stringent.pdf")

df_plot$seurat_prediction_cutoff <- rep("-",nrow(df_plot))
df_plot$seurat_prediction_cutoff[df_plot$seurat_max.score > .5] <- as.character(df_plot$seurat_prediction[df_plot$seurat_max.score > .5])
plot.index <- order(df_plot$seurat_prediction_cutoff)

pdf("transfer_label_mapscore_thr_stringent.pdf", width=12, height=8)
ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = factor(seurat_prediction_cutoff))) +
geom_point(size = 1) +
scale_color_manual(values=population_colours, name = "Cell population mapped", labels = names(population_colours)) +
theme_minimal() +
theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
guides(colour = guide_legend(override.aes = list(size=7)))
dev.off()


setwd('/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/scanpy/')


sce_ctrl  <- sce_pb[,grep("CTRL",meta$condition)]
meta_ctrl <- meta[grep("CTRL",meta$condition),]
meta_ctrl <- meta_ctrl[colData(sce_ctrl)$doublet_class == "singlet",]
sce_ctrl  <- sce_ctrl[,colData(sce_ctrl)$doublet_class == "singlet"]
write.table(as.matrix(logcounts(sce_ctrl)),"normalised_counts_ctrl_stringent.tab", sep="\t", row.names=FALSE, quote=FALSE)
writeLines(colnames(sce_ctrl),"cells_ctrl_stringent.txt")
writeLines(rownames(sce_ctrl),"genes_ctrl_stringent.txt")
write.table(meta_ctrl,"cell_metadata_ctrl_stringent.tab", sep="\t", row.names=FALSE, quote=FALSE)
