
setwd('/data1/ivanir/Ilaria2021/data')

hsm=read.csv('hsm_ltime.csv', header=TRUE)
hsm=hsm[-c(1:2),]

lsm=read.csv('lsm_ltime.csv', header=TRUE)

morphology <- c(rep("hsm",nrow(hsm)),rep("lsm",nrow(lsm)))

df_plot <- data.frame(t(rbind(hsm,lsm)))
genes   <- df_plot[1,]
df_plot <- df_plot[-1,]
colnames(df_plot) <- genes
df_plot1 <- data.frame(cbind(latent_time=gsub("X","",rownames(df_plot)),df_plot))
df_plot2 <- data.frame(cbind(morphology,t(df_plot)))

library("reshape2")
library("ggplot2")

g = ggplot(df_plot1,aes(x = as.numeric(latent_time), y = as.numeric(MALAT1)))
g + geom_line()

test_data_long_time <- melt(df_plot1, id="latent_time")  # convert to long format

test_data_long_morph <- melt(df_plot2, id="morphology")  # convert to long format

test_data_long_time_ext <- cbind(test_data_long_time,morphology=test_data_long_morph$morphology)

test_data_long <- melt(df_plot1[,c("latent_time","RBFOX1", "MALAT1", "TIMP3", "ITSN2")], id="latent_time")  # convert to long format

ggplot(data=test_data_long,
       aes(x=as.numeric(latent_time), y=as.numeric(value), colour=variable)) +
       geom_line()
       
pdf("test.pdf")
ggplot(data=test_data_long_time_ext,
       aes(x=as.numeric(latent_time), y=as.numeric(value), colour=morphology)) +
       geom_line()
 dev.off()      
