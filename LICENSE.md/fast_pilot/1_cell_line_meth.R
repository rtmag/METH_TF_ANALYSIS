
library(methylKit) #t
options(scipen=999)

setwd("/home/rtm/methmotif_cov/tfregulomeR/mm_tf_matrix/")
############
# WGBS
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


file.list <- as.list(file.list)


cells <- gsub(".+WGBS_MethMotif\\/","",file.list,perl=TRUE)
cells <- gsub("\\/.+","",cells,perl=TRUE)

cell_sample_number <- as.numeric(factor(cells))


myobj=methRead(file.list,
           sample.id = as.list(cells),
           assembly = "hg38",
           treatment = cell_sample_number,
           context = "CpG",
           pipeline = "bismarkCoverage",
           header = FALSE,
           mincov=5)
######################################

file.list <- list.files("/home/rtm/methmotif_cov/tfregulomeR/mm_tf_matrix/",pattern="_all.txt*",recursive=TRUE,full.names = TRUE)
file.id <- data.frame( do.call( rbind, strsplit( file.list, '//' ) ) )
file.id <- data.frame( do.call( rbind, strsplit( as.character(file.id[,2]), '_' ) ) )

cells <- as.character(unique(file.id[,1]))

for( i in 1:length(cells)){
file.cell <- list.files("/home/rtm/methmotif_cov/tfregulomeR/mm_tf_matrix/",pattern=paste("*",cells[i],"*",sep=""))
command = paste("cat",paste(file.cell,collapse=" "),"| sort -k1,1 -k2,2n|mergeBed -i -")
cell_merged <- read.table(pipe(command),stringsAsFactors=FALSE,sep="\t")

#wgbs.dir <- paste("/home/rtm/methmotif_cov/WGBS_MethMotif/",cells[i],"/",sep="")
#wgbs.path <- list.files(wgbs.dir,pattern="*.pooled.cov",recursive=TRUE,full.names = TRUE)
#if( identical(wgbs.path, character(0)) ){ wgbs.path <- list.files(wgbs.dir,pattern="*destranded.cov",recursive=TRUE,full.names = TRUE) }

  
}


bedt
