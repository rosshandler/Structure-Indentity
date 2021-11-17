library(irlba)
library(Rtsne)
library(scran)
library(Matrix)
library(igraph)
library(scater)
library(ggplot2)
library(scDblFinder)
library(matrixStats)
library(BiocParallel)
library(SingleCellExperiment)

path2data  <- "/data1/ivanir/Ilaria2021/data/"

counts <- readMM(paste0(path2data, "raw_counts_qc.mtx"))
meta   <- read.table(paste0(path2data, "meta_qc.tab"), header = TRUE)
genes  <- gsub(" ","",readLines(paste0(path2data, "gene_ids_loom.csv")))

condition <- rep(NA,nrow(meta))
condition[grep("Kit",meta$sample)]  <- "Kit"
condition[grep("Fast",meta$sample)] <- "Fast"
condition[grep("diss",meta$sample)] <- "Diss"
condition[grep("Unembed",meta$sample)]  <- "Unembed"

meta$condition <- condition

sce <- SingleCellExperiment(list(counts=counts),colData=DataFrame(meta))
rownames(sce) <- genes

bp <- MulticoreParam(12, RNGseed=1234)
bpstart(bp)
sce <- scDblFinder(sce, samples="sample", BPPARAM=bp)
bpstop(bp)
table(sce$scDblFinder.class)
#singlet doublet 
#30366    1576
  
sf <- readLines(paste0(path2data, "sizefactors.tab"))
sizeFactors(sce) <- as.numeric(sf)

saveRDS(sce, paste0(path2data, "sce.rds"))

lib.size <- colSums(counts(sce))

df_plot <- data.frame(
 lib.size = log10(lib.size),
 doublet  = colData(sce)$scDblFinder.class
)

ggplot(df_plot, aes(x=doublet, y=lib.size)) + 
  geom_boxplot() +
  labs(x = "", y = "log10(Library size)") 
ggsave("/data1/ivanir/Ilaria2021/QC/library_size_doublets.pdf")

sce <- sce[calculateAverage(sce)>0.1,]

sce <- logNormCounts(sce)

decomp  <- modelGeneVar(sce)
hvgs    <- rownames(decomp)[decomp$FDR < 0.1]
pca     <- prcomp_irlba(t(logcounts(sce[hvgs,])), n = 30)
rownames(pca$x) <- colnames(sce)
tsne <- Rtsne(pca$x, pca = FALSE, check_duplicates = FALSE)

df_plot <- data.frame(
 Dim1    = tsne$Y[, 1],
 Dim2    = tsne$Y[, 2], 
 doublet = colData(sce)$scDblFinder.class
)

plot.index <- order(df_plot$doublet)
ggplot(df_plot[plot.index,], aes(x = Dim1, y = Dim2, col = factor(doublet))) +
  geom_point(size = 0.4) +
  scale_color_manual(values=c("gray","#0169c1"), name = "") +
  labs(x = "Dim 1", y = "Dim 2") +
  theme_minimal() + #theme(legend.position = "none") +
  guides(colour = guide_legend(override.aes = list(size=7))) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=7)))
ggsave("/data1/ivanir/Ilaria2021/QC/tsne_doublets.pdf")
