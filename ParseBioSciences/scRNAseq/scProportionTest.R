library(miloR)
library(scran)
library(ggplot2)

library(scProportionTest)

setwd("/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/miloR")

path2data   <- '/data1/ivanir/Ilaria2023/ParseBS/newvolume/analysis/sCell/combined/all-well/DGE_unfiltered/'

sce  <- readRDS(paste0(path2data,"sce_resubmission_v1.rds"))
