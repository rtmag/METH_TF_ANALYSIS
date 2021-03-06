####################################################################################################################
library(data.table)
library(rtracklayer)
options(scipen=999)

file.list <- c("/home/rtm/methmotif_cov/WGBS_MethMotif/A549/A549.bismark.destranded.pooled.cov",
               "/home/rtm/methmotif_cov/WGBS_MethMotif/GM12878/GM12878.bismark.destranded.pooled.cov",
               "/home/rtm/methmotif_cov/WGBS_MethMotif/H1-hESC/H1-hESC.bismark.destranded.pooled.cov",
               "/home/rtm/methmotif_cov/WGBS_MethMotif/H9/H9.bismark.destranded.pooled.cov",
               "/home/rtm/methmotif_cov/WGBS_MethMotif/HCT116/combine_hct116_WGBS_1_bismark_bt2_pe.bismark.destranded.cov",
               "/home/rtm/methmotif_cov/WGBS_MethMotif/HEK293/GSM1254259_r1_bismark_bt2_pe.bismark.destranded.cov",
               "/home/rtm/methmotif_cov/WGBS_MethMotif/HEK293T/HEK239T_WGBS_r1_bismark_bt2_pe.bismark.destranded.cov",
               "/home/rtm/methmotif_cov/WGBS_MethMotif/HeLa-S3/SRR3574510_1_val_1_bismark_bt2_pe.bismark.destranded.cov",
               "/home/rtm/methmotif_cov/WGBS_MethMotif/HepG2/HepG2.bismark.destranded.pooled.cov",
               "/home/rtm/methmotif_cov/WGBS_MethMotif/HUES64/HUES64.bismark.destranded.pooled.cov",
               "/home/rtm/methmotif_cov/WGBS_MethMotif/IMR-90/ENCFF937OSM.bismark.destranded.cov",
               "/home/rtm/methmotif_cov/WGBS_MethMotif/K562/K562.bismark.destranded.pooled.cov",
               "/home/rtm/methmotif_cov/WGBS_MethMotif/LNCaP/SRR5219980_1_bismark_bt2_pe.bismark.destranded.cov",
               "/home/rtm/methmotif_cov/WGBS_MethMotif/MCF-7/SRR1169729_1_bismark_bt2_pe.bismark.destranded.cov",
               "/home/rtm/methmotif_cov/WGBS_MethMotif/SK-N-SH/SK-N-SH.bismark.destranded.pooled.cov",
               "/home/rtm/methmotif_cov/WGBS_MethMotif/SNU398/SNU398_WGBS_r1_bismark_bt2_pe.bismark.destranded.cov")

wgbs = list()
wgbs.gr = list()

for(i in 1:length(file.list) ){
  cell.name <- gsub("/home/rtm/methmotif_cov/WGBS_MethMotif/","",file.list[i])
  cell.name <- gsub("\\/.+cov","",cell.name,perl=TRUE)
  print(cell.name)
  covfile <- fread( file.list[i],
                  sep="\t",
                  data.table=TRUE, header=FALSE)
  covfile[,2] = covfile[,2]-1
  wgbs[[cell.name]] <- covfile
  
  wgbs_df <- as.data.frame(covfile[,1:3])
  colnames(wgbs_df) <- c("chr","start","end")
  wgbs_i.gr <- makeGRangesFromDataFrame(wgbs_df)
  wgbs.gr[[cell.name]] <- wgbs_i.gr
  }

cells <- unlist(file.list)
cells <- gsub("/home/rtm/methmotif_cov/WGBS_MethMotif/","",cells)
cells <- gsub("/.+","",cells,perl=TRUE)
##################################################################################################################################
hypo_merged <- read.table("hypometh_regions_merged.bed",stringsAsFactors=FALSE,sep="\t")
colnames(hypo_merged) <- c("chr","start","end")
hypo_merged.gr <- makeGRangesFromDataFrame(hypo_merged) 
#######redo

  cpg_cell_merged <- hypo_merged

for( i in 1:length(cells)){
      print(i)

    prev_names = colnames(cpg_cell_merged)
    cpg_cell_merged <- cbind(cpg_cell_merged,newbeta=NA)

    hits <- findOverlaps(hypo_merged.gr, wgbs.gr[[cells[i]]])
    hits.df <- as.data.frame(hits)
    if(is.unsorted(hits.df[,1])){ print("hits1 is unsorted") }
    if(is.unsorted(hits.df[,2])){ print("hits2 is unsorted") }
#
    cpgs_in_bed=wgbs[[cells[i]]][hits.df[,2]]
#
  
    cpgs_in_bed.dt = cpgs_in_bed[, .(CpGnum = .N,ReadNum=(sum(V5)+sum(V6)),beta=round(sum(V5)*100/(sum(V5)+sum(V6))) ), by = hits.df[,1] ]
    cpgs_in_bed.dt[(cpgs_in_bed.dt$ReadNum/cpgs_in_bed.dt$CpGnum)<5]$beta = NA #Filter by coverage
    cpg_cell_merged$newbeta[cpgs_in_bed.dt$hits.df] <- cpgs_in_bed.dt$beta
    colnames(cpg_cell_merged) <- c(prev_names,paste(cells[i],"_beta",sep="") )
}


options(scipen=999)
library(gplots)
library(factoextra)
library(RColorBrewer)
options(bitmapType="cairo")

colors <- rev(colorRampPalette( (brewer.pal(9, "RdBu")) )(9))

beta_matrix <- cpg_cell_merged[,4:dim(cpg_cell_merged)[2]]

#all.meth.norm = beta_matrix[complete.cases(beta_matrix),]
all.meth.norm = beta_matrix
colnames(all.meth.norm) <- gsub("_beta","",colnames(all.meth.norm))

all.meth.norm <- all.meth.norm[rowSums(is.na(all.meth.norm)) != ncol(all.meth.norm),]
all.meth.norm = all.meth.norm[apply(all.meth.norm,1,sd)!=0,]

png("heatmap_hypometh_cells.png",width= 3.25,
  height= 3.25,units="in",
  res=1200,pointsize=4)
heatmap.2(as.matrix(all.meth.norm),col=colors,scale="none", trace="none",srtCol=90,
labRow = FALSE,xlab="", ylab="CpGs",key.title="Methylation lvl",na.color="grey",na.rm=TRUE)
dev.off()
###
#

dist_no_na <- function(mat) {
    edist <- dist(mat)
    edist[which(is.na(edist))] <- max(edist, na.rm=TRUE) * 1.1 
    return(edist)
}

png("heatmap_hypometh_cells.png",width= 3.25,
  height= 3.25,units="in",
  res=1200,pointsize=4)
heatmap.2(as.matrix(all.meth.norm),col=colors,scale="none", trace="none",srtCol=90,distfun=dist_no_na,
labRow = FALSE,xlab="", ylab="CpGs",key.title="Methylation lvl",na.color="grey",na.rm=TRUE)
dev.off()
###

