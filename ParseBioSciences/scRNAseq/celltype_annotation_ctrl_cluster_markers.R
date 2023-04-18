library(scran)

setwd('/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/celltype_annotation/')
load('plots_Ilaria_PB_17Apr2023.RData')

markers <- findMarkers(sce, groups = df_plot$leiden_dmap, direction = "up")

for (clus in as.character(1:17)){
  markers_clus <- markers[[clus]]
  markers_clus <- cbind(symbol=rownames(markers_clus), markers_clus)
  write.table(markers_clus, file = paste0("DiffExp_Cluster",clus,".txt"), row.names = FALSE, sep = "\t", quote = FALSE)
  print(clus)
}
