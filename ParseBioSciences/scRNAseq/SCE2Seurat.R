library(scran)
library(Matrix)
library(Seurat)

sce  <- readRDS("/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/all-well/DGE_unfiltered/sce_resubmission_v1.rds")
rownames(sce) <- rowData(sce)$gene_id

seurat <- as.Seurat(sce, counts = "decontXcounts", data = "decontXlogcounts")

saveRDS(seurat,"/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/all-well/DGE_unfiltered/seurat_resubmission_v1.rds")

