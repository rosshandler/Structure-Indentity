library("ggplot2")
library("reshape2")

setwd('/data1/ivanir/Ilaria2021/data')

hsm <- read.csv('hsm_ltime.csv', header=TRUE)
hsm <- hsm[-c(1:2),]

lsm=read.csv('lsm_ltime.csv', header=TRUE)

hsm_mean <- apply(hsm[,-1],2,mean)
lsm_mean <- apply(lsm[,-1],2,mean)

hsm_sd <- apply(hsm[,-1],2,sd)
lsm_sd <- apply(lsm[,-1],2,sd)

latent_time=as.numeric(gsub("X","",colnames(hsm[,-1])))

df <- data.frame(cbind(latent_time,hsm_mean,lsm_mean))
df_long1 <- melt(df, id="latent_time")
colnames(df_long1) <- c("latent_time","morphology","mean")

df <- data.frame(cbind(latent_time,hsm_sd,lsm_sd))
df_long2 <- melt(df, id="latent_time")
colnames(df_long2) <- c("latent_time","morphology","sd")

df_long <- data.frame(df_long1,sd=df_long2$sd)

pdf("HSMvsLSM_Trend.pdf")
ggplot(data=df_long,
       aes(x=as.numeric(latent_time), y=as.numeric(mean), colour=morphology)) +
       geom_line() + theme_minimal() + theme(legend.position = "none") +
       geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = morphology), alpha = .2, colour = NA) +
       xlab("Cells in latent time") + ylab("Mean pooled gene expression")
dev.off()

hsm <- read.csv('hsm_ltime_prog.csv', header=TRUE)

lsm <- read.csv('lsm_ltime_prog.csv', header=TRUE)

hsm_mean <- apply(hsm[,-1],2,mean)
lsm_mean <- apply(lsm[,-1],2,mean)

hsm_sd <- apply(hsm[,-1],2,sd)
lsm_sd <- apply(lsm[,-1],2,sd)

latent_time=as.numeric(gsub("X","",colnames(hsm[,-1])))

df <- data.frame(cbind(latent_time,hsm_mean,lsm_mean))
df_long1 <- melt(df, id="latent_time")
colnames(df_long1) <- c("latent_time","morphology","mean")

df <- data.frame(cbind(latent_time,hsm_sd,lsm_sd))
df_long2 <- melt(df, id="latent_time")
colnames(df_long2) <- c("latent_time","morphology","sd")

df_long <- data.frame(df_long1,sd=df_long2$sd)

pdf("HSMvsLSM_Trend_prog.pdf")
ggplot(data=df_long,
       aes(x=as.numeric(latent_time), y=as.numeric(mean), colour=morphology)) +
       geom_line() + theme_minimal() + theme(legend.position = "none") +
       geom_ribbon(aes(y = mean, ymin = mean - sd, ymax = mean + sd, fill = morphology), alpha = .2, colour = NA) +
       xlab("Cells in latent time") + ylab("Mean pooled gene expression")
dev.off()


