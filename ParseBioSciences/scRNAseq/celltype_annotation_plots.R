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

leiden_annot_colours <- c(
"Differentiating RG"="#d54d92",
"Mixed Identity RG/Neu 1"="#5cc556",
"Mixed Identity RG/Neu 2"="#b254bf",
"Chp 1"="#9cb735",
"Glicolytic Neuronal"="#6c69ca",
"RG 1"="#cda937",
"IPCs"="#5d8fcb",
"Migrating Excitatory Neurons"="#dc5b31",
"Mitotic RG 1"="#42bdc0",
"Mixed Identity RG/Neu 3"="#dd4663",
"UL Neurons"="#56993f",
"DL Neurons 1"="#c78acc",
"Mixed Identity RG/Neu 4"="#61bf8c",
"RG 2"="#b63f37",
"Chp 2"="#407a48",
"Inhibitory Neurons"="#9e4a6b",
"DL Neurons 2"="#b5ae68",
"CR Cells"="#e18880",
"Mitotic RG 2"="#757327",
"Mixed Identity Chp/Neu"="#d48a3c",
"Chemochine Signaling"="#9c5e32")

leiden_colours <- leiden_annot_colours
names(leiden_colours)<- 0:20

day_colours <- c(
"#0054b9",
"#ff506d",
"#905700")
names(day_colours) <-c("48","55","70")

sce  <- readRDS("/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/all-well/DGE_unfiltered/sce_resubmission_v1.rds")
gexp <- as.matrix(logcounts(sce,assay.type = "decontXcounts"))
rownames(gexp) <- rowData(sce)$gene_name

df_plot <- data.frame(colData(sce), reducedDim(sce, "UMAP"))

leiden_annot_factor_order <- c(
  "RG 1","Differentiating RG","RG 2",
  "Mitotic RG 1","Mitotic RG 2",
  "IPCs",
  "Glicolytic Neuronal","Migrating Excitatory Neurons","DL Neurons 1","DL Neurons 2","UL Neurons","Inhibitory Neurons","CR Cells",
  "Mixed Identity RG/Neu 1","Mixed Identity RG/Neu 2","Mixed Identity RG/Neu 3","Mixed Identity RG/Neu 4",
  "Chp 1","Chp 2","Mixed Identity Chp/Neu",
  "Chemochine Signaling"
)
df_plot$leiden_annot <- factor(df_plot$leiden_annot, levels=leiden_annot_factor_order)

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

plotLayoutPseudotime <- function(layout="UMAP"){
  require(ggplot2)
  plot.index  <- order(df_plot$pt_monocle3)
  ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, colour = pt_monocle3)) + 
    geom_point(size = 1) +
    scale_color_gradient(low="gold", high="darkgreen") +
    theme_minimal() + 
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
    xlab("UMAP 1") + ylab("UMAP 2")
}

plotLayoutCelltypeMapped <- function(layout="UMAP"){
  require(ggplot2)
    ggplot(df_plot, aes(x = UMAP1, y = UMAP2, col = factor(seurat_prediction))) +
      geom_point(size = 1) +
      scale_color_manual(values=population_colours, name = "Cell population mapped", labels = names(population_colours)) +
      theme_minimal() + 
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
      guides(colour = guide_legend(override.aes = list(size=7),nrow=21))  
}

plotLayoutLeidenAnnot <- function(layout="UMAP"){
  require(ggplot2)
    ggplot(df_plot, aes(x = UMAP1, y = UMAP2, col = factor(leiden_annot))) +
      geom_point(size = 1) +
      scale_color_manual(values=leiden_annot_colours[leiden_annot_factor_order], name = "Leiden Annotated") +
      theme_minimal() + 
      labs(col="Leiden") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7),nrow=21))  
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

save.image(file='/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/celltype_annotation/plots_Ilaria_PB_16May2023.RData')

