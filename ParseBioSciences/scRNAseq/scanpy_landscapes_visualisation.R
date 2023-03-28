library(ggplot2)
library(Seurat)

setwd("/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/scanpy")

meta <- read.csv("structure-identity_pb.csv", header = TRUE)

umap <- read.csv("umap_layout.csv", header = FALSE)
colnames(umap) <- c("UMAP1", "UMAP2")

colnames(meta)[1] <- "cell"
df_plot <- cbind(meta, umap)

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

pdf("umap_scanpy_transfer_celltype.pdf", width=12, height=8)
ggplot(df_plot, aes(x = UMAP1, y = UMAP2, col = factor(seurat_prediction))) +
  geom_point(size = 1) +
  scale_color_manual(values=population_colours, name = "Cell population mapped", labels = names(population_colours)) +
  theme_minimal() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=7)))
dev.off()

plot.index <- sample(1:nrow(df_plot))
ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = as.factor(day))) +
geom_point(size = 1) +
#scale_color_manual(values=population_colours, name = "Cell population mapped", labels = names(population_colours)) +
theme_minimal() +
theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
guides(colour = guide_legend(override.aes = list(size=7)))
ggsave("umap_scanpy_day.pdf")

plot.index <- sample(1:nrow(df_plot))
ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = condition)) +
geom_point(size = 1) +
#scale_color_manual(values=population_colours, name = "Cell population mapped", labels = names(population_colours)) +
theme_minimal() +
theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
guides(colour = guide_legend(override.aes = list(size=7)))
ggsave("umap_scanpy_condition.pdf")

ggplot(df_plot, aes(x = UMAP1, y = UMAP2, col = factor(seurat_prediction))) +
  geom_point(size = 1) +
  #geom_point(size = 3, aes(alpha = seurat_max.score)) + scale_alpha("Mapping score") +
  scale_color_manual(values=population_colours, name = "Cell population mapped", labels = names(population_colours)) +
  theme_minimal() + theme(legend.position = "none") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=7))) +
  facet_wrap(~seurat_prediction)
ggsave("umap_scanpy_split_transfer_label.pdf")

plot.index <- sample(1:nrow(df_plot))
ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = as.factor(day))) +
geom_point(size = 1) +
#scale_color_manual(values=population_colours, name = "Cell population mapped", labels = names(population_colours)) +
theme_minimal() + theme(legend.position = "none") +
theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
guides(colour = guide_legend(override.aes = list(size=7))) +
facet_wrap(~day)
ggsave("umap_scanpy_split_day.pdf")


plot.index <- sample(1:nrow(df_plot))
ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = condition)) +
geom_point(size = 1) +
#scale_color_manual(values=population_colours, name = "Cell population mapped", labels = names(population_colours)) +
theme_minimal() + theme(legend.position = "none") +
theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
guides(colour = guide_legend(override.aes = list(size=7))) +
facet_wrap(~condition)
ggsave("umap_scanpy_split_condition.pdf")

ggplot(df_plot, aes(x = UMAP1, y = UMAP2, col = factor(leiden_pca))) +
  geom_point(size = 1) +
  #scale_color_manual(values=population_colours, name = "Cell population mapped", labels = names(population_colours)) +
  theme_minimal() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=7)))
  ggsave("umap_scanpy_leiden.pdf")

ggplot(df_plot, aes(x = UMAP1, y = UMAP2, col = phase)) +
  geom_point(size = 1) +
  #scale_color_manual(values=population_colours, name = "Cell population mapped", labels = names(population_colours)) +
  theme_minimal() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=7)))
  ggsave("umap_scanpy_phase.pdf")


##################
##################

library(ggplot2)

setwd("/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/scanpy")

meta <- read.csv("structure-identity_pb_ctrl.csv", header = TRUE)

umap <- read.csv("umap_layout_ctrl.csv", header = FALSE)
colnames(umap) <- c("UMAP1", "UMAP2")

colnames(meta)[1] <- "cell"
df_plot <- cbind(meta, umap)

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

pdf("umap_scanpy_transfer_celltype_ctrl.pdf", width=12, height=8)
ggplot(df_plot, aes(x = UMAP1, y = UMAP2, col = factor(seurat_prediction))) +
  geom_point(size = 1) +
  scale_color_manual(values=population_colours, name = "Cell population mapped", labels = names(population_colours)) +
  theme_minimal() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=7)))
dev.off()

plot.index <- sample(1:nrow(df_plot))
ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = as.factor(day))) +
geom_point(size = 1) +
#scale_color_manual(values=population_colours, name = "Cell population mapped", labels = names(population_colours)) +
theme_minimal() +
theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
guides(colour = guide_legend(override.aes = list(size=7)))
ggsave("umap_scanpy_day_ctrl.pdf")

plot.index <- sample(1:nrow(df_plot))
ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = condition)) +
geom_point(size = 1) +
#scale_color_manual(values=population_colours, name = "Cell population mapped", labels = names(population_colours)) +
theme_minimal() +
theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
guides(colour = guide_legend(override.aes = list(size=7)))
ggsave("umap_scanpy_condition_ctrl.pdf")

ggplot(df_plot, aes(x = UMAP1, y = UMAP2, col = factor(seurat_prediction))) +
  geom_point(size = 1) +
  #geom_point(size = 3, aes(alpha = seurat_max.score)) + scale_alpha("Mapping score") +
  scale_color_manual(values=population_colours, name = "Cell population mapped", labels = names(population_colours)) +
  theme_minimal() + theme(legend.position = "none") +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=7))) +
  facet_wrap(~seurat_prediction)
ggsave("umap_scanpy_split_transfer_label_ctrl.pdf")

plot.index <- sample(1:nrow(df_plot))
ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = as.factor(day))) +
geom_point(size = 1) +
#scale_color_manual(values=population_colours, name = "Cell population mapped", labels = names(population_colours)) +
theme_minimal() + theme(legend.position = "none") +
theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
guides(colour = guide_legend(override.aes = list(size=7))) +
facet_wrap(~day)
ggsave("umap_scanpy_split_day_ctrl.pdf")

plot.index <- sample(1:nrow(df_plot))
ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = condition)) +
geom_point(size = 1) +
#scale_color_manual(values=population_colours, name = "Cell population mapped", labels = names(population_colours)) +
theme_minimal() + theme(legend.position = "none") +
theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
guides(colour = guide_legend(override.aes = list(size=7))) +
facet_wrap(~condition)
ggsave("umap_scanpy_split_condition_ctrl.pdf")

ggplot(df_plot, aes(x = UMAP1, y = UMAP2, col = factor(leiden_pca))) +
  geom_point(size = 1) +
  #scale_color_manual(values=population_colours, name = "Cell population mapped", labels = names(population_colours)) +
  theme_minimal() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=7)))
  ggsave("umap_scanpy_leiden_ctrl.pdf")


cluster_colours <- c("#8ea900",
"#db64e9",
"#5bde5d",
"#bf008a",
"#2a9600",
"#74339d",
"#99d767",
"#931a80",
"#8c8800",
"#4199ff",
"#d07000",
"#02bbe9",
"#db1834",
"#02b399",
"#ff4863",
"#006328",
"#ff73bb",
"#475a02",
"#b8a0ff",
"#972e04",
"#00588b",
"#b50020",
"#78bfa0",
"#962a42",
"#acd289",
"#fdaafa",
"#feb86c",
"#ab6880")
names(cluster_colours)<-0:27

ggplot(df_plot, aes(x = UMAP1, y = UMAP2, col = factor(leiden_dmap))) +
  geom_point(size = 1) +
  scale_color_manual(values=cluster_colours, name = "Leiden", labels = names(cluster_colours)) +
  theme_minimal() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=7)))
  ggsave("umap_scanpy_leiden_dmap_ctrl.pdf")


 ggplot(df_plot, aes(x = UMAP1, y = UMAP2, col = phase)) +
  geom_point(size = 1, shape = 16) +
  #scale_color_manual(values=population_colours, name = "Cell population mapped", labels = names(population_colours)) +
  theme_minimal() + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  guides(colour = guide_legend(override.aes = list(size=7)))
  ggsave("umap_scanpy_phase_ctrl.pdf")
