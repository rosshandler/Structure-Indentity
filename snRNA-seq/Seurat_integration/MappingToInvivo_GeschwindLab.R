library(irlba)
library(scran)
library(mclust)
library(Seurat)
library(Matrix)
library(ggplot2)

population_colours <- c(
  "Mic"   = "#935de6",
  "InMGE" = "#90b600",
  "PgG2M" = "#005dd8",
  "ExDp1" = "#d3b000",
  "PgS"   = "#5484ff",
  "InCGE" = "#00c571",
  "ExM-U" = "#fa2274",
  "oRG"   = "#02d7ec",
  "ExM"   = "#ff6500",
  "OPC"   = "#02acf4",
  "End"   = "#bd6200",
  "vRG"   = "#ff7bd8",
  "ExDp2" = "#9e5d56",
  "Per"   = "#ada7de",
  "ExN"   = "#a6023e",
  "IP"    = "#ff8ba6",
  "-"     = "gray"
)

population_labels <- c(
  "Mic"   = "Microglia",
  "InMGE" = "Interneuron MGE",
  "PgG2M" = "Cyclin prog. G2/M phase",
  "ExDp1" = "Excitatory deep layer 1",
  "PgS"   = "Cyclin prog. S phase",
  "InCGE" = "Interneuron CGE",
  "ExM-U" = "Mature excitatory upper enriched",
  "oRG"   = "oRG",
  "ExM"   = "Mature excitatory",
  "OPC"   = "OPC",
  "End"   = "Endothelial",
  "vRG"   = "VRG",
  "ExDp2" = "Excitatory deep layer 2",
  "Per"   = "Pericyte",
  "ExN"   = "Migrating excitatory",
  "IP"    = "IP",
  "-"     = "Bellow threshold"
)
 
setwd("/data1/ivanir/HumaFetalNeocortex2019/data")
sce_invivo  <- readRDS("sce.rds")
meta_invivo <- readRDS("metadata_annotation_refinement.rds")

filter_index <- which(
  meta_invivo$celltype.mapped == "End" |
  meta_invivo$celltype.mapped == "Per" |
  meta_invivo$celltype.mapped == "Mic")

meta_invivo <- meta_invivo[-filter_index,]
sce_invivo  <- sce_invivo[, -filter_index]

sce_invivo  <- sce_invivo[scater::calculateAverage(sce_invivo)>0.05,]

setwd("/data1/ivanir/Ilaria2021/data")
sce_ilaria  <- readRDS("sce.rds")
sce_ilaria  <- sce_ilaria[, colData(sce_ilaria)$scDblFinder.class == "singlet"]
sce_ilaria  <- sce_ilaria[scater::calculateAverage(sce_ilaria)>0.05,]
meta_ilaria <- read.csv("structure-identity_clustering_batch_corrected.csv")

shared_genes <- intersect(rownames(sce_invivo), rownames(sce_ilaria))

sce <- SingleCellExperiment(list(counts=cbind(counts(sce_ilaria[shared_genes,]), counts(sce_invivo[shared_genes,]))))
rownames(sce) <- shared_genes

batch <- c(rep("query", ncol(sce_ilaria)), rep("atlas", ncol(sce_invivo))) 

metadata <- data.frame(cbind(cell=c(meta_ilaria$cell, meta_invivo$Cell), batch))
rownames(metadata)   <- c(colnames(sce_ilaria),colnames(sce_invivo))

### Seurat integration
seurat_integ <- CreateSeuratObject(cbind(counts(sce_ilaria[shared_genes,]), counts(sce_invivo[shared_genes,])), meta.data = metadata)
seurat_list  <- SplitObject(seurat_integ, split.by = "batch")

### SCT transformation
for (i in names(seurat_list)) {
    seurat_list[[i]] <- SCTransform(seurat_list[[i]], verbose = FALSE)
}

seurat_features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3249)
seurat_list     <- PrepSCTIntegration(object.list = seurat_list, anchor.features = seurat_features)

reference_dataset <- which(names(seurat_list) == "atlas")

seurat_anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", 
    anchor.features = seurat_features, reference = reference_dataset)
    
seurat_integrated <- IntegrateData(anchorset = seurat_anchors, normalization.method = "SCT")

## Label transfer
for (i in 1:length(x = seurat_list)) {
   seurat_list[[i]] <- NormalizeData(object = seurat_list[[i]], verbose = TRUE)
   seurat_list[[i]] <- FindVariableFeatures(object = seurat_list[[i]], selection.method = "vst", 
        nfeatures = 2000, verbose = TRUE)
}

reference_dataset <- which(names(seurat_list) == "atlas")

#seurat_anchors    <- FindIntegrationAnchors(object.list = seurat_list, reference = reference_dataset, dims = 1:30)
#seurat_integrated <- IntegrateData(anchorset = seurat_anchors, dims = 1:30)

seurat_query <- seurat_list[["query"]]
seurat_atlas <- seurat_list[["atlas"]]
seurat_anchors <- FindTransferAnchors(reference = seurat_atlas, query = seurat_query, 
    dims = 1:30)
predictions    <- TransferData(anchorset = seurat_anchors, refdata = meta_invivo$celltype.mapped, 
    dims = 1:30)
seurat_query   <- AddMetaData(object = seurat_query, metadata = predictions)

meta_ilaria <- cbind(meta_ilaria, seurat_prediction=seurat_query$predicted.id, seurat_max.score=seurat_query$prediction.score.max)

saveRDS(meta_ilaria, "meta_seurat_mapping.rds")

meta_ilaria_tmp <- readRDS("meta_seurat_mapping.rds")

layout1  <- read.csv("umap_layout_x.csv", header = FALSE)
layout2  <- read.csv("fa2_layout_x.csv", header = FALSE)

layout1  <- read.csv("umap_layout_batch_corrected.csv", header = FALSE)
layout2  <- read.csv("fa2_layout_batch_corrected.csv", header = FALSE)

colnames(layout1) <- c("UMAP1_scanpy", "UMAP2_scanpy")
colnames(layout2)  <- c("Fa1", "Fa2")

df_plot <- cbind(meta_ilaria, layout1, layout2)

population_colours<-population_colours[c("ExDp1", "ExDp2", "ExM-U", "ExN", "InCGE", "InMGE", "IP", "OPC", "oRG", "vRG", "PgG2M", "PgS")]
population_labels<-population_labels[c("ExDp1", "ExDp2", "ExM-U", "ExN", "InCGE", "InMGE", "IP", "OPC", "oRG", "vRG", "PgG2M", "PgS")]

pdf("transfer_label.pdf", width=12, height=8)
ggplot(df_plot, aes(x = UMAP1_scanpy, y = UMAP2_scanpy, col = factor(seurat_prediction))) +
  geom_point(size = 1) +
  #geom_point(size = 3, aes(alpha = seurat_max.score)) + scale_alpha("Mapping score") +
  scale_color_manual(values=population_colours, name = "Cell population mapped", labels = population_labels) +
  theme_minimal() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=7)))
dev.off()

pdf("transfer_label_threshold.pdf", width=12, height=8)
df_plot$seurat_prediction_cutoff <- rep("-",nrow(df_plot))
df_plot$seurat_prediction_cutoff[meta_ilaria$seurat_max.score > .5] <- as.character(df_plot$seurat_prediction[meta_ilaria$seurat_max.score > .5])
plot.index <- order(df_plot$seurat_prediction_cutoff)
ggplot(df_plot[plot.index,], aes(x = UMAP1_scanpy, y = UMAP2_scanpy, col = factor(seurat_prediction_cutoff))) +
geom_point(size = 1) +
scale_color_manual(values=population_colours, name = "Cell population mapped", labels = population_labels) +
theme_minimal() +
theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
guides(colour = guide_legend(override.aes = list(size=7)))
dev.off()
