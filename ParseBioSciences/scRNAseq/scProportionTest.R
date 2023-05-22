library(miloR)
library(scran)
library(ggplot2)

library(scProportionTest)

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

setwd("/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/miloR")

path2data   <- '/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/all-well/DGE_unfiltered/'

sce  <- readRDS(paste0(path2data,"sce_resubmission_v1.rds"))
