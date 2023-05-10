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

sce  <- readRDS("/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/all-well/DGE_unfiltered/sce_transferred_annot.rds")
gexp <- as.matrix(logcounts(sce,assay.type = "decontXcounts"))
rownames(gexp) <- rowData(sce)$gene_name

df_plot <- data.frame(colData(sce), reducedDim(sce, "UMAP"))

leiden_annot <- as.character(df_plot$leiden)
leiden_annot <- gsub("^0$","Differentiating RG",leiden_annot)
leiden_annot <- gsub("^1$","Mixed Identity RG/Neu 1",leiden_annot)
leiden_annot <- gsub("^2$","Mixed Identity RG/Neu 2",leiden_annot)
leiden_annot <- gsub("^3$","Chp 1",leiden_annot)
leiden_annot <- gsub("^4$","Glicolytic Neuronal",leiden_annot)
leiden_annot <- gsub("^5$","RG 1",leiden_annot)
leiden_annot <- gsub("^6$","IPCs",leiden_annot)
leiden_annot <- gsub("^7$","Migrating Excitatory Neurons",leiden_annot)
leiden_annot <- gsub("^8$","Mitotic RG 1",leiden_annot)
leiden_annot <- gsub("^9$","Mixed Identity RG/Neu 3",leiden_annot)
leiden_annot <- gsub("^10$","UL Neurons",leiden_annot)
leiden_annot <- gsub("^11$","DL Neurons 1",leiden_annot)
leiden_annot <- gsub("^12$","Mixed Identity RG/Neu 4",leiden_annot)
leiden_annot <- gsub("^13$","RG 2",leiden_annot)
leiden_annot <- gsub("^14$","Chp 2",leiden_annot)
leiden_annot <- gsub("^15$","Inhibitory Neurons",leiden_annot)
leiden_annot <- gsub("^16$","DL Neurons 2",leiden_annot)
leiden_annot <- gsub("^17$","CR Cells",leiden_annot)
leiden_annot <- gsub("^18$","Mitotic RG 2",leiden_annot)
leiden_annot <- gsub("^19$","Mixed Identity Chp/Neu",leiden_annot)
leiden_annot <- gsub("^20$","Chemochine Signaling",leiden_annot)

leiden_annot_factor_order <- c(
  "RG 1","Differentiating RG","RG 2",
  "Mitotic RG 1","Mitotic RG 2",
  "IPCs",
  "Glicolytic Neuronal","Migrating Excitatory Neurons","DL Neurons 1","DL Neurons 2","UL Neurons","Inhibitory Neurons","CR Cells",
  "Mixed Identity RG/Neu 1","Mixed Identity RG/Neu 2","Mixed Identity RG/Neu 3","Mixed Identity RG/Neu 4",
  "Chp 1","Chp 2","Mixed Identity Chp/Neu",
  "Chemochine Signaling"
)
df_plot$leiden_annot <- factor(leiden_annot, levels=leiden_annot_factor_order)
