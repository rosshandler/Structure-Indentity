library(dtwclust)

library("reshape2")
library("ggplot2")

library("ape")

setwd('/data1/ivanir/Ilaria2021/data')
hsm <- read.delim("HSMneu_brain_atlas.txt",header=TRUE,sep="\t")
lsm <- read.delim("LSMneu_brain_atlas.txt",header=TRUE,sep="\t")
fcd <- read.delim("FCD_brain_atlas.txt",header=TRUE,sep="\t")
sri <- read.delim("Srivastava_brain_atlas.txt",header=TRUE,sep="\t")

labels <- factor(c(rep("hsm", nrow(hsm)), rep("lsm", nrow(lsm)), rep("fcd", nrow(fcd)), rep("sri", nrow(sri))))
colors <- c(rep("red", nrow(hsm)), rep("blue", nrow(lsm)), rep("pink", nrow(fcd)), rep("green", nrow(sri)))

series <- rbind(hsm,lsm,fcd,sri)

pc.l2 <- tsclust(series[,-1], k = 4L,
                 distance = "L2", centroid = "pam",
                 seed = 3247, trace = TRUE,
                 control = partitional_control(nrep = 10L))

hc.l2 <- tsclust(series[,-1], type = "hierarchical",
                 k = 4L, trace = TRUE,
                 control = hierarchical_control(method = "all",
                                                distmat = pc.l2[[1L]]@distmat))

# Plot the best dendrogram according to variation of information
plot(hc.l2[[which.min(sapply(hc.l2, cvi, b = labels, type = "VI"))]])

plot(as.phylo(hc.l2[[which.min(sapply(hc.l2, cvi, b = labels, type = "VI"))]]), type = "fan", tip.color = colors,cex = 0.7, pch = 21)

clusters <- cutree(hc.l2[[which.min(sapply(hc.l2, cvi, b = labels, type = "VI"))]], k = 100)

clus1 <- which(clusters == 1)
clus2 <- which(clusters == 2)
clus3 <- which(clusters == 6)

series1 <- series[clus1,]
series2 <- series[clus2,]
series3 <- series[clus3,]

series1_tmp <- t(series1[,-1])
colnames(series1_tmp) <- series1$Gene
PCW <- as.numeric(gsub("PCW","",gsub("X","",colnames(series1[,-1]))))
test_data1 <- data.frame(cbind(series1_tmp, PCW))

test_data_long1 <- melt(test_data1, id="PCW")  # convert to long format

ggplot(data=test_data_long1,
       aes(x=PCW, y=value, colour=variable)) +
       geom_line()
       
labels1 <- factor(labels[clus1])
test_data2 <- data.frame(cbind(t(series1_tmp),gene_set = as.character(labels1)))

test_data_long2 <- melt(test_data2, id="gene_set")  # convert to long format

test_data_long <- data.frame(cbind(test_data_long1,gene_set=test_data_long2$gene_set))

test_data_long$gene_set <- factor(test_data_long$gene_set, levels=c("hsm","lsm","fcd","sri"))
col_key <- c("#F8766D","#00BFC4","#7CAE00","#C77CFF")
names(col_key) <- levels(test_data_long$gene_set)

ggplot(data=test_data_long,
       aes(x=PCW, y=value, group=variable, colour=gene_set)) +
       geom_line(size=1.5) +
       scale_color_manual(values=col_key)
       
#library("ggdendro")
#ggdendrogram(hc.l2[[which.min(sapply(hc.l2, cvi, b = labels, type = "VI"))]], rotate = TRUE, theme_dendro = FALSE)
