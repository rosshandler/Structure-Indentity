library(scran)
library(irlba)
library(scater)
library(Matrix)
library(ggplot2)

library(umap)
library(leiden)
library(reticulate)
use_condaenv(condaenv="scanpy-p3.9")

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

leiden_colours <- c(
 "#985088",
"#67b649",
"#8d5dcd",
"#537b35",
"#c850ba",
"#aaa947",
"#0054b9",
"#ce5d33",
"#50b48f",
"#d03c49",
"#697dc8",
"#dc9a38",
"#59a9dd",
"#a5743c",
"#d38fd3",
"#ce7276",
"#d24683")
names(cluster_colours)<-1:17

day_colours <- c(
"#0054b9",
"#ff506d",
"#905700")
names(day_colours) <-c("48","55","70")

setwd("/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/scanpy")
meta <- read.csv("structure-identity_pb_ctrl.csv", header = TRUE)
colnames(meta)[1] <- "cell"

#umap <- read.csv("umap_layout_ctrl.csv", header = FALSE)
#colnames(umap) <- c("UMAP1", "UMAP2")
#df_plot <- cbind(meta, umap)

#df_plot$sample <- factor(df_plot$sample, levels=names(sample_colours))
#df_plot$condition <- factor(df_plot$condition, levels=names(condition_colours))
#df_plot$celltype_annotation <- factor(df_plot$celltype_annotation, levels=names(population_colours))
#df_plot$celltype_manuscript <- factor(df_plot$celltype_manuscript, levels=names(population_colours_manuscript))

sce_ctrl <- readRDS(
"/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/all-well/DGE_unfiltered/sce.rds")[,gsub("_pb","",as.character(meta$cell))]
sce_ctrl <- logNormCounts(sce_ctrl)
rownames(sce_ctrl) <- rowData(sce_ctrl)$gene_name
gexp <- logcounts(sce_ctrl)

sce_ctrl <- sce_ctrl[calculateAverage(sce_ctrl)>0.01,]
decomp  <- modelGeneVar(sce_ctrl)
hvgs    <- rownames(decomp)[decomp$FDR < 0.01]
length(hvgs)
#[1] 930
pca     <- prcomp_irlba(t(logcounts(sce_ctrl[hvgs,])), n = 30)
rownames(pca$x) <- colnames(sce_ctrl)
#tsne    <- Rtsne(pca$x, pca = FALSE, check_duplicates = FALSE, num_threads=30)
layout  <- umap(pca$x, method="umap-learn", umap_learn_args=c("n_neighbors", "n_epochs", "min_dist"), n_neighbors=30, min_dist=.25)

graph <- buildSNNGraph(pca$x, d = NA, transposed = TRUE)
set.seed(42)
clusters <- leiden(graph, resolution_parameter = 2)
names(clusters) <-  meta$cell

df_plot <- data.frame(
 meta,
 leiden = clusters,
 #tSNE1    = tsne$Y[, 1],
 #tSNE2    = tsne$Y[, 2], 
 UMAP1 = layout$layout[,1],
 UMAP2 = layout$layout[,2] 
)
#saveRDS(df_plot, ""/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/all-well/DGE_unfiltered/umap_sce.rds"")

plotLayoutExpression <- function(gene="TTR"){
  require(Matrix)
  require(ggplot2)
    logcounts <- as.vector(as.matrix(gexp[gene,]))
    if (sum(logcounts)>0){
        df_tmp    <- data.frame(cbind(df_plot, logcounts))
        plot.index  <- order(df_tmp$logcounts)
        ggplot(df_tmp[plot.index,], aes(x = UMAP1, y = UMAP2, colour = logcounts)) + 
          geom_point(size = 1) +
          scale_color_gradient(low="gray", high="darkgreen") +
          labs(color = paste0(gene,"\nlog(counts)")) +
          theme_minimal() + 
          theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
          theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
          xlab("Dimension 1") + ylab("Dimension 2")
    }else{
    message(gene," was not detected in the expression matrix")
    }
}

plotLayoutCelltypeMapped <- function(layout="UMAP"){
  require(ggplot2)
    ggplot(df_plot, aes(x = UMAP1, y = UMAP2, col = factor(seurat_prediction))) +
      geom_point(size = 1) +
      scale_color_manual(values=population_colours, name = "Cell population mapped", labels = names(population_colours)) +
      theme_minimal() + 
      labs(col="Leiden") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7)))  
}

plotLayoutLeiden <- function(layout="UMAP"){
  require(ggplot2)
    ggplot(df_plot, aes(x = UMAP1, y = UMAP2, col = factor(leiden))) +
      geom_point(size = 1) +
      scale_color_manual(values=leiden_colours, name = "Leiden") +
      theme_minimal() + 
      labs(col="Leiden") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7)))  
}

plotLayoutCondition <- function(layout="UMAP"){
  require(ggplot2)
  plot.index  <- sample(nrow(df_plot))
    ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = factor(condition))) +
      geom_point(size = 1) +        
      #scale_color_manual(values=condition_colours, name = "Condition", labels=condition_labels) +
      theme_minimal() + 
      labs(col="Condition/Day") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7)))
}

plotLayoutSample <- function(layout="UMAP"){
  require(ggplot2)
  plot.index  <- sample(nrow(df_plot))
    ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = factor(sample_name))) +
      geom_point(size = 1) +
      theme_minimal() +
      labs(col="Sample") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7)))
}

plotLayoutPhase <- function(layout="UMAP"){
  require(ggplot2)
  plot.index  <- sample(nrow(df_plot))
     ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = factor(phase))) +
      geom_point(size = 1) +
      theme_minimal() +
      labs(col="Phase") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7)))
}

plotViolinExpressionLeiden <- function(gene="TTR"){
    require(Matrix)
    require(ggplot2)
    logcounts <- as.vector(as.matrix(gexp[gene,]))
    if (sum(logcounts)>0){
        leiden <- df_plot$leiden
        ggplot(mapping =  aes(x = leiden, 
                              y = logcounts, 
                              fill = factor(leiden))) +
        geom_violin(scale = "width") +
        geom_boxplot(width=0.1) +
        scale_fill_manual(values=leiden_colours, name = "Leiden") +
        labs(y = "Log2 normalised count", x = "Leiden clusters") + 
        ggtitle(paste0(gene," expression across leiden clusters")) +
        theme_minimal() + theme(legend.position = "none") +
        theme(axis.title = element_text(face = "bold", size = 12),
              axis.text.y = element_text(size = 12, face = "bold"),
              axis.text.x = element_blank(),
              axis.title.x = element_blank()) 
    }else{
      message(gene," was not detected in the expression matrix")
    }
}

save.image(file='/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/all-well/DGE_unfiltered/plots_Ilaria_PB_28March2023.RData')
