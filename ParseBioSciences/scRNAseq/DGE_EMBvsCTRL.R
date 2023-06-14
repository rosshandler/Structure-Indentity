library(scran)
library(scater)

sce  <- readRDS("/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/all-well/DGE_unfiltered/sce_resubmission_v1.rds")

#####
ctrl_index <- grep("CTRL",colData(sce)$condition)
diss_index <- grep("DISS",colData(sce)$condition)
emb_index  <- grep("EMB", colData(sce)$condition)

condition_full <- rep(NA,nrow(colData(sce)))

condition_full[ctrl_index] <- "CTRL"
condition_full[diss_index] <- "DISS"
condition_full[emb_index]  <- "EMB"

sce_subset <- sce[,c(emb_index,ctrl_index)]
condition  <- condition_full[c(emb_index,ctrl_index)]

markers <- findMarkers(sce_subset, groups = condition, direction = "up")
markers_clus <- markers[[2]]
markers_clus <- cbind(symbol=rownames(markers_clus), markers_clus)
setwd('/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/all-well/DGE_unfiltered/')
write.table(markers_clus, file = paste0("DiffExpEMB.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
