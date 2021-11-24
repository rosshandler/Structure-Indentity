library(scran)
library(scater)
library(igraph)
library(Matrix)

path2data  <- "/data1/ivanir/Ilaria2021/data/"

counts <- readMM(paste0(path2data, "raw_counts_qc.mtx"))
meta   <- read.table(paste0(path2data, "meta_qc.tab"), header = TRUE)
genes  <- gsub(" ","",readLines(paste0(path2data, "gene_ids_loom.csv")))

sce <- SingleCellExperiment(list(counts=counts))
rownames(sce) <- genes

lib.sizes <- colSums(counts(sce))
sce       <- sce[calculateAverage(sce)>0.1,]

clusts <- as.numeric(quickCluster(sce, method = "igraph", min.size = 100))

min.clust <- min(table(clusts))/2
new_sizes <- c(floor(min.clust/3), floor(min.clust/2), floor(min.clust))
sce <- computeSumFactors(sce, clusters = clusts, sizes = new_sizes, max.cluster.size = 3000)

ggplot(data = data.frame(X = lib.sizes, Y = sizeFactors(sce)), mapping = aes(x = X, y = Y)) +
  geom_point() +
  scale_x_log10(breaks = c(500, 2000, 5000, 10000, 30000), labels = c("5,00", "2,000", "5,000", "10,000", "30,000") ) +
  scale_y_log10(breaks = c(0.2, 1, 5)) +
  theme_minimal() +
  theme(text = element_text(size=20))  +
  labs(x = "Number of UMIs", y = "Size Factor")
ggsave("/data1/ivanir/Ilaria2021/QC/sizefactors.pdf")

write.table(sizeFactors(sce), quote = F, col.names = F, row.names = F,
  file = paste0(path2data, "sizefactors.tab"))
  
