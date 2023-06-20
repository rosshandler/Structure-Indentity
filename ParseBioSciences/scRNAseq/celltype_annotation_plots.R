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

setwd('/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/scanpy/')
sce_ctrl  <- sce[,grep("CTRL", colData(sce)$condition)]
colData(sce_ctrl)$day <- paste0(colData(sce_ctrl)$day," day")

write.table(as.matrix(logcounts(sce_ctrl)),"normalised_counts_ctrl.tab", sep="\t", row.names=FALSE, quote=FALSE)
writeLines(colnames(sce_ctrl),"cells_ctrl.txt")
writeLines(rownames(sce_ctrl),"genes_ctrl.txt")
write.table(data.frame(colData(sce_ctrl)),"cell_metadata_ctrl.tab", sep="\t", row.names=FALSE, quote=FALSE)

sce_diss  <- sce[,grep("DISS", colData(sce)$condition)]

colData(sce_diss)$day <- paste0(colData(sce_diss)$day," day")

write.table(as.matrix(logcounts(sce_diss)),"normalised_counts_diss.tab", sep="\t", row.names=FALSE, quote=FALSE)
writeLines(colnames(sce_diss),"cells_diss.txt")
writeLines(rownames(sce_diss),"genes_diss.txt")
write.table(data.frame(colData(sce_diss)),"cell_metadata_diss.tab", sep="\t", row.names=FALSE, quote=FALSE)

##############
sample_subset <- c("6/07_d45_B", "6/07_d48_B", "6/07_d55_DISSB","6/07_d70_DISSB")

df_subset <- df_plot[df_plot$sample_name %in% sample_subset,]

plotLayoutSampleSub <- function(layout="UMAP"){
  require(ggplot2)
  plot.index  <- sample(nrow(df_subset))
    ggplot(df_subset[plot.index,], aes(x = UMAP1, y = UMAP2, col = factor(sample_name))) +
      geom_point(size = 2) +
      theme_minimal() +
      labs(col="Sample") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7)))
}

p <- plotLayoutSampleSub()
p + facet_wrap(~sample_name)

#####
ctrl_index <- grep("CTRL",df_plot$condition)
diss_index <- grep("DISS",df_plot$condition)
emb_index  <- grep("EMB", df_plot$condition)

condition_full <- rep(NA,nrow(df_plot))

condition_full[ctrl_index] <- "CTRL"
condition_full[diss_index] <- "DISS"
condition_full[emb_index]  <- "EMB"

df_plot$condition_full <- condition_full

plotLayoutConditionExtra <- function(layout="UMAP"){
  require(ggplot2)
  plot.index  <- sample(nrow(df_plot))
    ggplot(df_plot[plot.index,], aes(x = UMAP1, y = UMAP2, col = factor(condition_full))) +
      geom_point(size = 1) +        
      theme_minimal() + 
      labs(col="Condition") +
      theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
      theme(axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
      guides(colour = guide_legend(override.aes = list(size=7)))
}

plotLayoutConditionExtra() + facet_wrap(~condition_full)

###### Heatmap mixed identity
df_plot$leiden_annot <- gsub("DL Neurons 2","Mixed Identity RG/Neu 5", df_plot$leiden_annot)

RG_markers <- c("DACH1","GLI3","MEIS2","CREB5","SHROOM3","PAX6","NPAS3","ZFHX4","SLC1A3")

N_markers  <- c("CTNNA2","ZFPM2","GRIA2","NRXN1","SLC24A2","DCC","BCL11B","PTPRD","KCNQ3")

markers <- c(RG_markers, N_markers)

dat1 <- apply(gexp[markers, df_plot$leiden_annot == "Mixed Identity RG/Neu 1"],1,mean)
dat2 <- apply(gexp[markers, df_plot$leiden_annot == "Mixed Identity RG/Neu 2"],1,mean)
dat3 <- apply(gexp[markers, df_plot$leiden_annot == "Mixed Identity RG/Neu 3"],1,mean)
dat4 <- apply(gexp[markers, df_plot$leiden_annot == "Mixed Identity RG/Neu 4"],1,mean)
dat5 <- apply(gexp[markers, df_plot$leiden_annot == "Mixed Identity RG/Neu 5"],1,mean)
dat6 <- apply(gexp[markers, df_plot$leiden_annot == "RG 1"],1,mean)
dat7 <- apply(gexp[markers, df_plot$leiden_annot == "Migrating Excitatory Neurons"],1,mean)

dat <- cbind(dat6,dat7,dat1,dat2,dat3,dat4,dat5)

pheatmap::pheatmap(dat,scale="row",cluster_cols=FALSE,cluster_rows=FALSE)

colnames(dat) <- c("RG 1","Migrating Excitatory Neurons","Mixed Identity RG/Neu 1","Mixed Identity RG/Neu 2","Mixed Identity RG/Neu 3","Mixed Identity RG/Neu 4","Mixed Identity RG/Neu 5")


############
###### Heatmap mixed identity
library(viridis)

df_plot$leiden_annot <- gsub("DL Neurons 2","Mixed Identity RG/Neu 5", df_plot$leiden_annot)

RG_markers <- c("DACH1","GLI3","MEIS2","CREB5","SHROOM3","PAX6","NPAS3","ZFHX4","SLC1A3")

#N_markers  <- c("CTNNA2","ZFPM2","GRIA2","NRXN1","SLC24A2","DCC","BCL11B","PTPRD","KCNQ3")
N_markers  <- c("KCND2","KCNQ3","GRIA2","GRIA1","NRXN1","SLC24A2","GRM1","GRIN1")
N_markers  <- c("KCND2","KCNQ3","GRIA2","GRIA1","NRXN1","SLC24A2")

markers <- c(RG_markers, N_markers)

dat1 <- gexp[markers, df_plot$leiden_annot == "Mixed Identity RG/Neu 1"]
dat2 <- gexp[markers, df_plot$leiden_annot == "Mixed Identity RG/Neu 2"]
dat3 <- gexp[markers, df_plot$leiden_annot == "Mixed Identity RG/Neu 3"]
dat4 <- gexp[markers, df_plot$leiden_annot == "Mixed Identity RG/Neu 4"]
dat5 <- gexp[markers, df_plot$leiden_annot == "Mixed Identity RG/Neu 5"]
dat6 <- gexp[markers, df_plot$leiden_annot == "RG 1"]
dat7 <- gexp[markers, df_plot$leiden_annot == "Migrating Excitatory Neurons"]

dat <- cbind(dat6,dat7,dat1,dat2,dat3,dat4,dat5)

dat_annot <- data.frame(c(
  rep("RG",ncol(dat6)),
  rep("Migrating Excitatory Neurons",ncol(dat7)),
  rep("Mixed Identity RG/Neu 1",ncol(dat1)),
  rep("Mixed Identity RG/Neu 2",ncol(dat2)),
  rep("Mixed Identity RG/Neu 3",ncol(dat3)),
  rep("Mixed Identity RG/Neu 4",ncol(dat4)),
  rep("Mixed Identity RG/Neu 5",ncol(dat5))
))
rownames(dat_annot) <- colnames(dat)
colnames(dat_annot) <- "Clusters"
dat_annot$Clusters <- factor(dat_annot$Clusters, levels=c("RG","Migrating Excitatory Neurons","Mixed Identity RG/Neu 1","Mixed Identity RG/Neu 2","Mixed Identity RG/Neu 3","Mixed Identity RG/Neu 4","Mixed Identity RG/Neu 5"))

dat_annotation_colors <- list(c(
"RG"="#cda937",
"Migrating Excitatory Neurons"="#dc5b31",
"Mixed Identity RG/Neu 1"="#5cc556",
"Mixed Identity RG/Neu 2"="#b254bf",
"Mixed Identity RG/Neu 3"="#dd4663",
"Mixed Identity RG/Neu 4"="#61bf8c",
"Mixed Identity RG/Neu 5"="#b5ae68"
))
names(dat_annotation_colors) <- "Clusters"

setwd('/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/plots')
pheatmap::pheatmap(dat, scale="row", cluster_cols=FALSE, cluster_rows=FALSE, annotation_col=dat_annot, show_colnames = FALSE, color=viridis_pal()(20), file="hmap1.pdf")
pheatmap::pheatmap(dat, scale="row", cluster_cols=FALSE, cluster_rows=FALSE, annotation_col=dat_annot, show_colnames = FALSE, annotation_colors=dat_annotation_colors, file="hmap2.pdf")









