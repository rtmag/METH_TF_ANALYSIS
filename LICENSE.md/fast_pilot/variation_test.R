library(vegan)
library(beeswarm)
library(graphics)
library(gplots)
library(ggplot2)

meth_mat <- readRDS("HyperTFBS_methMatrix.rds")
unmeth_mat <- readRDS("hypoTFBS_methMatrix.rds")

meth <- meth_mat[,4:dim(meth_mat)[2]]
unmeth <- unmeth_mat[,4:dim(unmeth_mat)[2]]

######################################################################################################
meth.sd<-apply(meth,1,sd)
unmeth.sd<-apply(unmeth,1,sd)

meth.shannon<-apply(meth,1,function(x){ diversity(x[!is.na(x)], index = "shannon") } )
unmeth.shannon<-apply(unmeth,1,function(x){ diversity(x[!is.na(x)], index = "shannon") } )

sd_obj<-rbind(data.frame(cluster="Methylated",sd=meth.sd),data.frame(cluster="Unmethylated",sd=unmeth.sd))
colnames(sd_obj) <- c("cluster","sd")
sd_obj<-sd_obj[!is.na(sd_obj$sd),]

shannon_obj<-rbind(data.frame(cluster="Methylated",shannon=meth.shannon),
                   data.frame(cluster="Unmethylated",shannon=unmeth.shannon))
colnames(shannon_obj) <- c("cluster","shannon")
shannon_obj<-shannon_obj[!is.na(shannon_obj$shannon),]

######################################################################################################
#meth vs non-meth
boxplot(meth.sd,unmeth.sd)
boxplot(meth.shannon,unmeth.shannon)

par(mfrow=c(2,1))
stripchart(sd ~ cluster, vertical = TRUE, data = sd_obj,jitter = 0.4,
           ylab = expression('SD of Beta-scores in TF clusters across 16 cell lines'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.003),cex = 2,las=1)
boxplot(sd ~ cluster,add = TRUE, vertical = TRUE, data = sd_obj,boxlwd = 3,las=1,outline=FALSE)

stripchart(shannon ~ cluster, vertical = TRUE, data = shannon_obj,jitter = 0.4,
           ylab = expression('Shannon index of Beta-scores in TF clusters across 16 cell lines'),
    method = "jitter", pch = 20, col = alpha(colour='red',alpha=.003),cex = 2,las=1)
boxplot(shannon ~ cluster,add = TRUE, vertical = TRUE, data = shannon_obj,boxlwd = 3,las=1,outline=FALSE)

pdf("variation_boxplots_nojitter.pdf",width=4,height=8)
par(mfrow=c(2,1))
boxplot(sd ~ cluster,add = F, vertical = TRUE, data = sd_obj,boxlwd = 3,las=1,outline=FALSE,
        ylab = 'SD of B-val in TF clusters across 16 cell lines')
boxplot(shannon ~ cluster,add = F, vertical = TRUE, data = shannon_obj,boxlwd = 3,las=1,outline=FALSE,
        ylab = 'Shannon index of B-val in TF clusters across 16 cell lines')
dev.off()
