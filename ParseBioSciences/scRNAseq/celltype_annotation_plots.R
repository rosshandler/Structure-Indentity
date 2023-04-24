library(scran)
library(irlba)
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

leiden_colours <- c(
"#d54d92",
"#5cc556",
"#b254bf",
"#9cb735",
"#6c69ca",
"#cda937",
"#5d8fcb",
"#dc5b31",
"#42bdc0",
"#dd4663",
"#56993f",
"#c78acc",
"#61bf8c",
"#b63f37",
"#407a48",
"#9e4a6b",
"#b5ae68",
"#e18880",
"#757327",
"#d48a3c",
"#9c5e32")
names(leiden_colours)<-1:21

day_colours <- c(
"#0054b9",
"#ff506d",
"#905700")
names(day_colours) <-c("48","55","70")

sce  <- readRDS("/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/all-well/DGE_unfiltered/sce_transferred_annot.rds")
gexp <- as.matrix(logcounts(sce,assay.type = "decontXcounts"))
rownames(gexp) <- rowData(sce)$gene_name

df_plot <- data.frame(colData(sce), reducedDim(sce, "UMAP"))

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
    ggplot(df_plot, aes(x = UMAP1, y = UMAP2, col = factor(leiden_dmap))) +
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

save.image(file='/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/celltype_annotation/plots_Ilaria_PB_24Apr2023.RData')
