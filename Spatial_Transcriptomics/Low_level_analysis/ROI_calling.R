library(RImageJROI)

setwd('/data1/ivanir/Ilaria2021/SpatialData/MolecularCartography/20230_segmentation')

slide1_A1 <- read.table("20230_slide1_A1_DAPI.tiff_measurements.csv", sep = "\t", header = TRUE)
slide1_A2 <- read.table("20230_slide1_A2_DAPI.tiff_measurements.csv", sep = "\t", header = TRUE)
slide6_A2 <- read.table("20230_slide6_A2_DAPI.tiff_measurements.csv", sep = "\t", header = TRUE)
slide6_C1 <- read.table("20230_slide6_C1_DAPI.tiff_measurements.csv", sep = "\t", header = TRUE)
slide6_D1 <- read.table("20230_slide6_D1_DAPI.tiff_measurements.csv", sep = "\t", header = TRUE)

sample_info <- read.table("sample_info.txt", header = TRUE)

setwd('/data1/ivanir/Ilaria2021/SpatialData/MolecularCartography/20230_segmentation/roi_files')
samples <- list.files()

i=1
  setwd(paste0("/data1/ivanir/Ilaria2021/SpatialData/MolecularCartography/20230_segmentation/roi_files/", samples)[i])
  roi <- lapply(list.files()[-1], function(x){read.ijroi(x, verbose = FALSE)})
  names(roi) <- gsub(".roi","",list.files()[-1])
  centroids_tmp <- lapply(roi, function(x){c(median(x$coords[,'x']), median(x$coords[,'y']))})
  centroids <- rlist::list.rbind(centroids_tmp)
  colnames(centroids) <- c("x","y")
  cells <- intersect(colnames(slide1_A1), rownames(centroids))
  centroids <- centroids[cells, ]
  counts    <- slide1_A1[, cells]
  rownames(counts) <- slide1_A1$Gene.ID

meta <- data.frame(
 cell   =    colnames(counts),
 sample =    rep("slide1_A1", ncol(counts)),
 condition = rep("Control", ncol(counts)),
 organoid  = rep(0, ncol(counts))
)
setwd('/data1/ivanir/Ilaria2021/SpatialData/MolecularCartography/20230_segmentation')
saveRDS(meta, "slide1_A1_meta.rds")
saveRDS(counts, "slide1_A1_counts.rds")
saveRDS(centroids, "slide1_A1_centroids.rds")

i=2
  setwd(paste0("/data1/ivanir/Ilaria2021/SpatialData/MolecularCartography/20230_segmentation/roi_files/", samples)[i])
  roi <- lapply(list.files()[-1], function(x){read.ijroi(x, verbose = FALSE)})
  names(roi) <- gsub(".roi","",list.files()[-1])
  centroids_tmp <- lapply(roi, function(x){c(median(x$coords[,'x']), median(x$coords[,'y']))})
  centroids <- rlist::list.rbind(centroids_tmp)
  colnames(centroids) <- c("x","y")
  cells <- intersect(colnames(slide1_A2), rownames(centroids))
  centroids <- centroids[cells, ]
  counts    <- slide1_A2[, cells]
  rownames(counts) <- slide1_A2$Gene.ID

meta <- data.frame(
 cell   =    colnames(counts),
 sample =    rep("slide1_A2", ncol(counts)),
 condition = rep("Control", ncol(counts)),
 organoid  = rep(0, ncol(counts))
)
setwd('/data1/ivanir/Ilaria2021/SpatialData/MolecularCartography/20230_segmentation')
saveRDS(meta, "slide1_A2_meta.rds")
saveRDS(counts, "slide1_A2_counts.rds")
saveRDS(centroids, "slide1_A2_centroids.rds")

i=3
  setwd(paste0("/data1/ivanir/Ilaria2021/SpatialData/MolecularCartography/20230_segmentation/roi_files/", samples)[i])
  roi <- lapply(list.files()[-1], function(x){read.ijroi(x, verbose = FALSE)})
  names(roi) <- gsub(".roi","",list.files()[-1])
  centroids_tmp <- lapply(roi, function(x){c(median(x$coords[,'x']), median(x$coords[,'y']))})
  centroids <- rlist::list.rbind(centroids_tmp)
  colnames(centroids) <- c("x","y")
  cells <- intersect(colnames(slide6_A2), rownames(centroids))
  centroids <- centroids[cells, ]
  counts    <- slide6_A2[, cells]
  rownames(counts) <- slide6_A2$Gene.ID

meta <- data.frame(
 cell   =    colnames(counts),
 sample =    rep("slide6_A2", ncol(counts)),
 condition = rep("No_MG", ncol(counts)),
 organoid  = rep(1, ncol(counts))
)
setwd('/data1/ivanir/Ilaria2021/SpatialData/MolecularCartography/20230_segmentation')
saveRDS(meta, "slide6_A2_meta.rds")
saveRDS(counts, "slide6_A2_counts.rds")
saveRDS(centroids, "slide6_A2_centroids.rds")

i=4
  setwd(paste0("/data1/ivanir/Ilaria2021/SpatialData/MolecularCartography/20230_segmentation/roi_files/", samples)[i])
  roi <- lapply(list.files()[-1], function(x){read.ijroi(x, verbose = FALSE)})
  names(roi) <- gsub(".roi","",list.files()[-1])
  centroids_tmp <- lapply(roi, function(x){c(median(x$coords[,'x']), median(x$coords[,'y']))})
  centroids <- rlist::list.rbind(centroids_tmp)
  colnames(centroids) <- c("x","y")
  cells <- intersect(colnames(slide6_C1), rownames(centroids))
  centroids <- centroids[cells, ]
  counts    <- slide6_C1[, cells]
  rownames(counts) <- slide6_C1$Gene.ID

meta <- data.frame(
 cell   =    colnames(counts),
 sample =    rep("slide6_C1", ncol(counts)),
 condition = rep("No_MG", ncol(counts)),
 organoid  = rep(1, ncol(counts))
)
setwd('/data1/ivanir/Ilaria2021/SpatialData/MolecularCartography/20230_segmentation')
saveRDS(meta, "slide6_C1_meta.rds")
saveRDS(counts, "slide6_C1_counts.rds")
saveRDS(centroids, "slide6_C1_centroids.rds")

i=5
  setwd(paste0("/data1/ivanir/Ilaria2021/SpatialData/MolecularCartography/20230_segmentation/roi_files/", samples)[i])
  roi <- lapply(list.files()[-1], function(x){read.ijroi(x, verbose = FALSE)})
  names(roi) <- gsub(".roi","",list.files()[-1])
  centroids_tmp <- lapply(roi, function(x){c(median(x$coords[,'x']), median(x$coords[,'y']))})
  centroids <- rlist::list.rbind(centroids_tmp)
  colnames(centroids) <- c("x","y")
  cells <- intersect(colnames(slide6_D1), rownames(centroids))
  centroids <- centroids[cells, ]
  counts    <- slide6_D1[, cells]
  rownames(counts) <- slide6_D1$Gene.ID

meta <- data.frame(
 cell   =    colnames(counts),
 sample =    rep("slide6_D1", ncol(counts)),
 condition = rep("No_MG", ncol(counts)),
 organoid  = rep(2, ncol(counts))
)
setwd('/data1/ivanir/Ilaria2021/SpatialData/MolecularCartography/20230_segmentation')
saveRDS(meta, "slide6_D1_meta.rds")
saveRDS(counts, "slide6_D1_counts.rds")
saveRDS(centroids, "slide6_D1_centroids.rds")

