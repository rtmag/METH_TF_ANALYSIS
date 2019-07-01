library(ggplot2)
library(dplyr)
library(reshape2)
library(vegan)

file.list <- list.files("./",pattern="*TFclusters_beta.bed",recursive=TRUE,full.names = FALSE)
#pdf(paste0("HYPER_TFcluster_BetaValues.pdf"),height=2)
for(i in 1:16){

beta = read.table(file.list[i],sep="\t",header=TRUE)
name <- gsub("_TFclusters_beta.bed","",file.list[i])
name <- gsub("-",".",name)

beta_m = beta[,7:dim(beta)[2]]
colnames(beta_m) <- gsub("_beta","",colnames(beta_m))

j = which(colnames(beta_m)==name)

hyper = beta_m[beta_m[,j]>80,]
hypo = beta_m[beta_m[,j]<20,]
table(beta[beta_m[,j]>80,4])

hyper_long <- melt(hyper)
hyper_long$condition <- "1"
hyper_long[hyper_long[,1]==name, 'condition'] <- '2'
colnames(hyper_long) <- c("CellLine","Beta","condition")
hyper_hypo <- melt(hypo)

x=ggplot(data=hyper_long, aes(x=CellLine, y=Beta)) +
geom_violin(aes(fill=condition, color=condition)) + ggtitle(name) +
geom_boxplot(width=.1, outlier.shape=NA) +
theme_minimal()+theme(legend.position="none",axis.text.x=element_blank())
#print(x)
}
#dev.off()


file.list <- list.files("./",pattern="*TFclusters_beta.bed",recursive=TRUE,full.names = FALSE)
for(i in 1:16){

beta = read.table(file.list[i],sep="\t",header=TRUE)
name <- gsub("_TFclusters_beta.bed","",file.list[i])
name <- gsub("-",".",name)

beta_m = beta[,7:dim(beta)[2]]
colnames(beta_m) <- gsub("_beta","",colnames(beta_m))

j = which(colnames(beta_m)==name)

x=as.character(beta[beta_m[,j]>80,4])
xx=sort(table(unlist(strsplit( x, ',' ) )))
write.table(xx, paste0(name,'TF_on_Hypermethylated_regions_table.txt'), sep='\t', col.names=F, row.names=F, quote=F)
}
#####################################################################################################################
#####################################################################################################################
#####################################################################################################################
# 1)Meth Variance in  methylated regions bound by any TF compared to regions bound by TF but not menthylated for each cell line….?
file.list <- list.files("./",pattern="*TFclusters_beta.bed",recursive=TRUE,full.names = FALSE)
pdf("Shannon_index_across_16Cells.pdf")
par(mfrow=c(4,4))
for(i in 1:16){
beta = read.table(file.list[i],sep="\t",header=TRUE)
name <- gsub("_TFclusters_beta.bed","",file.list[i])
name <- gsub("-",".",name)
cellname <- gsub("_.+","",file.list[i],perl=TRUE)

beta_m = beta[,7:dim(beta)[2]]
colnames(beta_m) <- gsub("_beta","",colnames(beta_m))

j = which(colnames(beta_m)==name)
thr_cov <- (beta[,6]/beta[,5])>10

hyper = beta_m[beta_m[,j]>80 & thr_cov,]
hypo = beta_m[beta_m[,j]<20 & thr_cov,]

hyper.shannon<-apply(hyper,1,function(x){ diversity(x[!is.na(x)], index = "shannon") } )
hypo.shannon<-apply(hypo,1,function(x){ diversity(x[!is.na(x)], index = "shannon") } )

shannon_obj<-rbind(data.frame(cluster="Hyper",shannon=hyper.shannon),
                   data.frame(cluster="Hypo",shannon=hypo.shannon))
colnames(shannon_obj) <- c("cluster","shannon")
shannon_obj<-shannon_obj[!is.na(shannon_obj$shannon),]
boxplot(shannon ~ cluster,add = F, vertical = TRUE, data = shannon_obj,boxlwd = 3,las=1,outline=FALSE,
        ylab = 'Beta Shannon index',main=cellname)
}
dev.off()
