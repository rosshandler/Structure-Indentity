library(Matrix)
library(biomaRt)
library(ggplot2)

path2data  <- "/data1/ivanir/Ilaria2021/data/"
spliced    <- readMM(paste0(path2data, "spliced_matrix.mtx"))
unspliced  <- readMM(paste0(path2data, "unspliced_matrix.mtx"))
genes      <- readLines(paste0(path2data, "gene_ids_loom.csv"))
barcodes   <- readLines(paste0(path2data, "barcodes_loom.csv"))
counts <- spliced + unspliced

## Sample and sequencing round metadata

sample_names <- c(
  "B6Fast"      = 1,
  "B6Kit"       = 2,
  "B7diss"      = 3,
  "B7Fast"      = 4,
  "B7Kit"       = 5,
  "B8diss"      = 6,
  "B8Kit"       = 7,
  "B8Unembed_A" = 8,
  "B6diss"      = 9,
  "B6Unembed"   = 10,
  "B8Fast"      = 11,
  "B8Unembed_B" = 12
)

samples   <- unlist(lapply(strsplit(barcodes, split = ":"), function(x) {x[1]}))
samples[21436:23082] <- paste0(samples[21436:23082],"_A")
samples[39618:43478] <- paste0(samples[39618:43478],"_B")

sequencing.round <- rep("SLX-19865", length(samples))
sequencing.round[samples %in% names(sample_names[9:12])] <- "SLX-20646"

## Library size

lib.sizes <- colSums(counts)

## Library complexity

ngenes    <- colSums(counts > 0)

## Mitochondrial gene expression

ensembl <- useEnsembl(biomart = "ensembl",  dataset = "hsapiens_gene_ensembl", mirror = "useast")

gene_map  <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "chromosome_name"),
  filters = "hgnc_symbol", values = genes, mart = ensembl)

mt.index    <- gene_map$chromosome_name == "MT"
mt.counts   <- counts[which(genes %in% gene_map$hgnc_symbol[mt.index]), ]
mt.fraction <- colSums(mt.counts)/lib.sizes

plot_df <- data.frame(
  sample   = samples,
  batch    = sequencing.round,
  lib.size = lib.sizes,
  n.genes  = ngenes,
  mt.fraction
)

# UMIs by sequencing round
pdf("/data1/ivanir/Ilaria2021/QC/UMIsByBatch.pdf")
ggplot(plot_df, aes (x = factor(batch), y = lib.size)) +
  geom_boxplot() +
  theme_bw() + 
  labs(x = "Batch", y = "Number of UMIs") +
  scale_y_log10(breaks = c(1000, 5000, 10000, 50000, 100000),
    labels = c("1,000", "5,000", "10,000", "50,000", "100,000"))
dev.off()

# Cell complexity thresholding
pdf("/data1/ivanir/Ilaria2021/QC/cell_complexity.pdf")
qplot(lib.sizes, ngenes, col = ifelse(ngenes < 500, "drop", "keep")) +
  scale_x_log10() +
  scale_y_log10() +
  theme_minimal() + 
  theme(text = element_text(size=20),legend.position = "none")  +
  labs(x = "UMI count", y = "Number of expressed genes") +
  scale_color_manual(values = c("drop" = "grey50", "keep" = "black"), name = "")
dev.off()

#####
##### Filtering step 1

## Library complexity

filter.cell.1 <- ngenes > 500

counts   <- counts[, filter.cell.1]
barcodes <- barcodes[filter.cell.1]
samples  <- samples[filter.cell.1]
sequencing.round <- sequencing.round[filter.cell.1]

lib.sizes <- colSums(counts)
ngenes    <- colSums(counts > 0)

#####
##### Filtering step 2

## Mitochondrial gene expression

mt.index    <- gene_map$chromosome_name == "MT"
mt.counts   <- counts[which(genes %in% gene_map$hgnc_symbol[mt.index]), ]
mt.fraction <- colSums(mt.counts)/lib.sizes

mt.p   <- pnorm(mt.fraction, mean = median(mt.fraction), sd = mad(mt.fraction), lower.tail = FALSE)
mt.lim <- min(mt.fraction[which(p.adjust(mt.p, method = "fdr") < 0.05)])

# Threhdold
mt.lim
[1] 0.004136505
            
qplot(lib.sizes, mt.fraction, col = ifelse(mt.fraction>mt.lim, "drop", "keep")) +
  scale_x_log10() +
  labs(x = "UMI count", y = "MT read fraction") +
  theme_minimal() + 
  theme(text = element_text(size=20),legend.position = "none")  +
  scale_color_manual(values = c("drop" = "grey50", "keep" = "black"), name = "")
ggsave("/data1/ivanir/Ilaria2021/QC/mtreadfraction.pdf")

filter.cell.2 <- mt.fraction > mt.lim

plot_df_filtered <- data.frame(
  sample = samples,
  lib.size = lib.sizes,
  mt.fraction,
  filter.cell = filter.cell.2
)

correct_sample_order  <- c(
 "B6Kit", "B7Kit", "B8Kit",
 "B6Fast", "B7Fast", "B8Fast",
 "B6diss", "B7diss", "B8diss",
 "B6Unembed", "B8Unembed_A", "B8Unembed_B")
 
plot_df_filtered$sample <- factor(plot_df_filtered$sample, levels = correct_sample_order)

ggplot(plot_df_filtered, aes(fill=filter.cell, y=1, x=sample)) +
    geom_bar(position="fill", stat="identity") +
    scale_fill_viridis(discrete = T, option = "E", name = "Filter cell") + 
    ylab("Fraction of cells removed by stress signals") + xlab("Sample")
ggsave("/data1/ivanir/Ilaria2021/QC/mt_fraction_reomval.pdf")

counts      <- counts[, mt.fraction < mt.lim]
barcodes    <- barcodes[mt.fraction < mt.lim]
samples     <- samples[mt.fraction  < mt.lim]
sequencing.round <- sequencing.round[mt.fraction < mt.lim]

mt.fraction <- mt.fraction[mt.fraction < mt.lim]
lib.sizes   <- colSums(counts)
n.genes     <- colSums(counts>0)

plot_df_filtered <- data.frame(sample = samples, lib.size = lib.sizes, batch = sequencing.round)
 
plot_df_filtered$sample <- factor(plot_df_filtered$sample, levels = correct_sample_order)

ggplot(plot_df_filtered, aes(x = factor(sample), y = lib.size)) +
  geom_boxplot() +
  theme_bw() +
  scale_y_log10() +
  labs(x = "Sample", y = "Number of UMIs")
ggsave("/data1/ivanir/Ilaria2021/QC/UMIsBySample.pdf")

#####
##### Export post-qc data
meta <- data.frame(cell = barcodes, sequencing.round = sequencing.round, sample = samples, mt.fraction = mt.fraction,  lib.size = lib.sizes, n.genes
)

writeMM(counts, file = paste0(path2data, "raw_counts_qc.mtx"))

write.table(meta, file = paste0(path2data, "meta_qc.tab"),
  row.names = F, col.names = T, quote = F, sep = "\t")
 
#####
##### Export pre-qc data

writeMM(spliced + unspliced, file = paste0(path2data, "raw_counts_preqc.mtx"))

cells   <- readLines(paste0(path2data, "barcodes_loom.csv"))
qc.filter.pass <- cells %in% meta$cell
plot_df <- cbind(cell = cells, plot_df, qc.filter.pass)
write.table(plot_df, file = paste0(path2data, "meta_preqc.tab"),
  row.names = F, col.names = T, quote = F, sep = "\t")

