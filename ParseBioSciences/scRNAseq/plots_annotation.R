library(scran)
library(scater)
library(Matrix)
library(ggplot2)

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

day_colours <- c("#0054b9",
"#ff506d",
"#905700")
names(day_colours) <-c("48","55","70")

setwd("/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/scanpy")

meta <- read.csv("structure-identity_pb_ctrl.csv", header = TRUE)

umap <- read.csv("umap_layout_ctrl.csv", header = FALSE)
colnames(umap) <- c("UMAP1", "UMAP2")

colnames(meta)[1] <- "cell"
df_plot <- cbind(meta, umap)

#df_plot$sample <- factor(df_plot$sample, levels=names(sample_colours))
#df_plot$condition <- factor(df_plot$condition, levels=names(condition_colours))
#df_plot$celltype_annotation <- factor(df_plot$celltype_annotation, levels=names(population_colours))
#df_plot$celltype_manuscript <- factor(df_plot$celltype_manuscript, levels=names(population_colours_manuscript))

sce_ctrl <- readRDS(
"/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/all-well/DGE_unfiltered/sce.rds")[,gsub("_pb","",as.character(meta$cell))]
sce_ctrl <- logNormCounts(sce_ctrl)
rownames(sce_ctrl) <- rowData(sce_ctrl)$gene_name
gexp <- logcounts(sce_ctrl)

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

sce_ctrl <- sce_ctrl[calculateAverage(sce_ctrl)>0.01,]
decomp  <- modelGeneVar(sce_ctrl)
hvgs    <- rownames(decomp)[decomp$FDR < 0.01]
length(hvgs)
#
library(irlba)
pca     <- prcomp_irlba(t(logcounts(sce_ctrl[hvgs,])), n = 30)
rownames(pca$x) <- colnames(sce_ctrl)
#tsne    <- Rtsne(pca$x, pca = FALSE, check_duplicates = FALSE, num_threads=30)
library(umap)
library(reticulate)
use_condaenv(condaenv="scanpy-p3.9")
layout  <- umap(pca$x, method="umap-learn", umap_learn_args=c("n_neighbors", "n_epochs", "min_dist"), n_neighbors=30, min_dist=.25)

df_plot <- data.frame(
 meta,
 #tSNE1    = tsne$Y[, 1],
 #tSNE2    = tsne$Y[, 2], 
 UMAP1 = layout$layout[,1],
 UMAP2 = layout$layout[,2] 
)

plotLayoutExpression()
saveRDS(df_plot, ""/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/all-well/DGE_unfiltered/umap_sce.rds"")







##to be edited
plotLayoutLeiden <- function(layout="UMAP"){
  require(ggplot2)
  if (layout=="FA"){ 
    ggplot(df_plot, aes(x = FA1, y = FA2, col = factor(leiden))) +
      geom_point(size = 1) +        
      scale_color_manual(values=leiden_colours, name = "Leiden") +
      theme_minimal() + 
      labs(col="Leiden") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7)))
  }else if(layout =="UMAP"){
    ggplot(df_plot, aes(x = UMAP1, y = UMAP2, col = factor(leiden))) +
      geom_point(size = 1) +
      scale_color_manual(values=leiden_colours, name = "Leiden") +
      theme_minimal() + 
      labs(col="Leiden") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7)))  
  }else{
    message("Layout name not found")
  }
}

plotLayoutCondition <- function(layout="UMAP"){
  require(ggplot2)
  plot.index  <- sample(nrow(df_plot))
  if (layout=="FA"){ 
    ggplot(df_plot[plot.index,], aes(x = FA1, y = FA2, col = factor(condition))) +
      geom_point(size = 1) +        
      scale_color_manual(values=condition_colours, name = "Condition", labels=condition_labels) +
      theme_minimal() + 
      labs(col="Condition") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7)))
  }else if(layout =="UMAP"){
    ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = factor(condition))) +
      geom_point(size = 1) +
      scale_color_manual(values=condition_colours, name = "Condition", labels=condition_labels) +
      theme_minimal() + 
      labs(col="Condition") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7)))  
  }else{
    message("Layout name not found")
  }
}

plotLayoutMorphology <- function(layout="UMAP"){
  require(ggplot2)
  plot.index  <- sample(nrow(df_plot))
  if (layout=="FA"){
    ggplot(df_plot[plot.index,], aes(x = FA1, y = FA2, col = factor(morphology_clus))) +
      geom_point(size = 1) +
      scale_color_manual(values=morphology_colours, name = "Morphology", labels=morphology_labels) +
      theme_minimal() +
      labs(col="Morphology") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7)))
  }else if(layout =="UMAP"){
    ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = factor(morphology_clus))) +
      geom_point(size = 1) +
      scale_color_manual(values=morphology_colours, name = "Morphology", labels=morphology_labels) +
      theme_minimal() +
      labs(col="Morphology") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7)))
  }else{
    message("Layout name not found")
  }
}

plotLayoutPopulation <- function(layout="UMAP"){
  require(ggplot2)
  plot.index  <- sample(nrow(df_plot))
  if (layout=="FA"){ 
    ggplot(df_plot[plot.index,], aes(x = FA1, y = FA2, col = factor(celltype_annotation))) +
      geom_point(size = 1) +        
      scale_color_manual(values=population_colours, name = "Cell type", labels = population_labels) +
      theme_minimal() + 
      labs(col="Cell type") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7)))
  }else if(layout =="UMAP"){
    ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = factor(celltype_annotation))) +
      geom_point(size = 1) +
      scale_color_manual(values=population_colours, name = "Cell type", labels = population_labels) +
      theme_minimal() + 
      labs(col="Cell type") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7)))  
  }else{
    message("Layout name not found")
  }
}

plotLayoutPopulationManuscript <- function(layout="UMAP"){
  require(ggplot2)
  plot.index  <- sample(nrow(df_plot))
  if (layout=="FA"){ 
    ggplot(df_plot[plot.index,], aes(x = FA1, y = FA2, col = factor(celltype_annotation))) +
      geom_point(size = 1) +        
      scale_color_manual(values=population_colours_manuscript, name = "Cell type", labels = population_labels_manuscript) +
      theme_minimal() + 
      labs(col="Cell type") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7)))
  }else if(layout =="UMAP"){
    ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = factor(celltype_annotation))) +
      geom_point(size = 1) +
      scale_color_manual(values=population_colours, name = "Cell type", labels = population_labels) +
      theme_minimal() + 
      labs(col="Cell type") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7)))  
  }else{
    message("Layout name not found")
  }
}

plotLayoutPopulationManuscript <- function(layout="UMAP"){
  require(ggplot2)
  plot.index  <- sample(nrow(df_plot))
  if (layout=="FA"){ 
    ggplot(df_plot[plot.index,], aes(x = FA1, y = FA2, col = factor(celltype_manuscript))) +
      geom_point(size = 1) +        
      scale_color_manual(values=population_colours_manuscript, name = "Cell type") +
      theme_minimal() + 
      labs(col="Cell type") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7)))
  }else if(layout =="UMAP"){
    ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = factor(celltype_manuscript))) +
      geom_point(size = 1) +
      scale_color_manual(values=population_colours_manuscript, name = "Cell type") +
      theme_minimal() + 
      labs(col="Cell type") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7)))  
  }else{
    message("Layout name not found")
  }
}

plotLayoutSample <- function(layout="UMAP"){
  require(ggplot2)
  plot.index  <- sample(nrow(df_plot))
  if (layout=="FA"){
    ggplot(df_plot[plot.index,], aes(x = FA1, y = FA2, col = factor(sample))) +
      geom_point(size = 1) +
      theme_minimal() +
      labs(col="Sample") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7)))
  }else if(layout =="UMAP"){
    ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = factor(sample))) +
      geom_point(size = 1) +
      theme_minimal() +
      labs(col="Sample") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7)))
  }else{
    message("Layout name not found")
  }
}

plotLayoutPhase <- function(layout="UMAP"){
  require(ggplot2)
  plot.index  <- sample(nrow(df_plot))
  if (layout=="FA"){
    ggplot(df_plot[plot.index,], aes(x = FA1, y = FA2, col = factor(phase))) +
      geom_point(size = 1) +
      theme_minimal() +
      labs(col="Sample") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7)))
  }else if(layout =="UMAP"){
    ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = factor(phase))) +
      geom_point(size = 1) +
      theme_minimal() +
      labs(col="Sample") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7)))
  }else{
    message("Layout name not found")
  }
}

plotLayoutSeqRound <- function(layout="UMAP"){
  require(ggplot2)
  plot.index  <- sample(nrow(df_plot))
  if (layout=="FA"){
    ggplot(df_plot[plot.index,], aes(x = FA1, y = FA2, col = factor(sequencing.round))) +
      scale_color_manual(values=c("#bf0024","#0169c1"), name = "Sequencing round") +
      geom_point(size = 1) +
      theme_minimal() +
      labs(col="Sequencing round") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7)))
  }else if(layout =="UMAP"){
    ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = factor(sequencing.round))) +
      scale_color_manual(values=c("#bf0024","#0169c1"), name = "Sequencing round") +
      geom_point(size = 1) +
      theme_minimal() +
      labs(col="Sequencing round") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7)))
  }else{
    message("Layout name not found")
  }
}


plotViolinExpressionSample <- function(gene="MALAT1"){
    require(Matrix)
    require(ggplot2)
    logcounts <- as.vector(as.matrix(gexp[gene,]))
    if (sum(logcounts)>0){
        sample <- df_plot$sample
        ggplot(mapping =  aes(x = sample, 
                              y = logcounts, 
                              fill = factor(sample))) +
        geom_violin(scale = "width") +
        geom_boxplot(width=0.1) +
        scale_fill_manual(values=sample_colours, name = "Sample") +
        labs(y = "Log2 normalised count") + 
        ggtitle(paste0(gene," expression across samples")) +
        theme_minimal() +
        theme(axis.title = element_text(face = "bold", size = 12),
              axis.text.y = element_text(size = 12, face = "bold"),
              axis.text.x = element_text(size = 12, face = "bold", angle = 90, hjust = 1, vjust = 0.5),
              legend.position = "none",
              axis.title.x = element_blank())
    }else{
      message(gene," was not detected in the expression matrix")
    }
}

plotViolinExpressionLeiden <- function(gene="MALAT1"){
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

plotViolinExpressionCondition <- function(gene="MALAT1"){
    require(Matrix)
    require(ggplot2)
    logcounts <- as.vector(as.matrix(gexp[gene,]))
    if (sum(logcounts)>0){
        condition <- df_plot$condition
        ggplot(mapping =  aes(x = condition, 
                              y = logcounts, 
                              fill = factor(condition))) +
        geom_violin(scale = "width") +
        geom_boxplot(width=0.1) +
        scale_fill_manual(values=condition_colours, name = "Condition") +
        labs(y = "Log2 normalised count") + 
        ggtitle(paste0(gene," expression across conditions")) +
        theme_minimal() + theme(legend.position = "none") +
        theme(axis.title = element_text(face = "bold", size = 12),
              axis.text.y = element_text(size = 12, face = "bold"),
              axis.text.x = element_blank(),
              legend.position = "none",
              axis.title.x = element_blank())
    }else{
      message(gene," was not detected in the expression matrix")
    }
}

plotViolinExpressionPopulation <- function(gene="MALAT1"){
    require(Matrix)
    require(ggplot2)
    logcounts <- as.vector(as.matrix(gexp[gene,]))
    if (sum(logcounts)>0){
        population <- df_plot$celltype_annotation
        ggplot(mapping =  aes(x = population, 
                              y = logcounts, 
                              fill = factor(population))) +
        geom_violin(scale = "width") +
        geom_boxplot(width=0.1) +
        scale_fill_manual(values=population_colours, name = "Cell type") +
        labs(y = "Log2 normalised count") + 
        ggtitle(paste0(gene," expression across populations")) +
        theme_minimal() + theme(legend.position = "none") +
        theme(axis.title = element_text(face = "bold", size = 12),
              axis.text.y = element_text(size = 12, face = "bold"),
              axis.text.x = element_blank(),
              legend.position = "none",
              axis.title.x = element_blank())
    }else{
      message(gene," was not detected in the expression matrix")
    }
}

plotViolinExpressionPopulationManuscript <- function(gene="MALAT1"){
    require(Matrix)
    require(ggplot2)
    logcounts <- as.vector(as.matrix(gexp[gene,]))
    if (sum(logcounts)>0){
        population <- df_plot$celltype_manuscript
        ggplot(mapping =  aes(x = population, 
                              y = logcounts, 
                              fill = factor(population))) +
        geom_violin(scale = "width") +
        geom_boxplot(width=0.1) +
        scale_fill_manual(values=population_colours_manuscript, name = "Cell type") +
        labs(y = "Log2 normalised count") + 
        ggtitle(paste0(gene," expression across populations")) +
        theme_minimal() + theme(legend.position = "none") +
        theme(axis.title = element_text(face = "bold", size = 12),
              axis.text.y = element_text(size = 12, face = "bold"),
              axis.text.x = element_blank(),
              legend.position = "none",
              axis.title.x = element_blank())
    }else{
      message(gene," was not detected in the expression matrix")
    }
}

save.image(file='plots_Ilaria_20Nov2021.RData')
