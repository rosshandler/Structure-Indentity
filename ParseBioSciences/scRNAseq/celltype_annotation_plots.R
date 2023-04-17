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

sce <- readRDS("/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/all-well/DGE_unfiltered/sce_transferred_annot.rds")
sce <- logNormCounts(sce)
