library(data.table)
library(rtracklayer)
options(scipen=999)

setwd("/home/rtm/methmotif_cov/tfregulomeR/mm_tf_matrix/")
############
# WGBS bed cov file
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

######################################
# TF PEAKS
file.list <- list.files("/home/rtm/methmotif_cov/tfregulomeR/mm_tf_matrix/",pattern="_all.txt*",recursive=TRUE,full.names = TRUE)
file.id <- data.frame( do.call( rbind, strsplit( file.list, '//' ) ) )
file.id <- data.frame( do.call( rbind, strsplit( as.character(file.id[,2]), '_' ) ) )

cells <- as.character(unique(file.id[,1]))

for( i in 1:length(cells)){
file.cell <- list.files("/home/rtm/methmotif_cov/tfregulomeR/mm_tf_matrix/",pattern=paste("*",cells[i],"*",sep=""))
command = paste("cat",paste(file.cell,collapse=" "),"| sort -k1,1 -k2,2n|mergeBed -i -")
cell_merged <- read.table(pipe(command),stringsAsFactors=FALSE,sep="\t")
colnames(cell_merged) <- c("chr","start","end")
cell_merged.gr <- makeGRangesFromDataFrame(cell_merged) 

hits <- findOverlaps(cell_merged.gr, wgbs.gr[[cells[i]]])
hits.df <- as.data.frame(hits)
if(is.unsorted(hits.df[,1])){ print("hits1 is unsorted") }
if(is.unsorted(hits.df[,2])){ print("hits2 is unsorted") }

cpgs_in_bed=wgbs[[cells[i]]][hits.df[,2]]
x = cpgs_in_bed[, .(cpgNum = .N,beta=round(sum(V5)*100/(sum(V5)+sum(V6))) ), by = hits.df[,1] ]
  
x= head(cpgs_in_bed,500)
x1 = head(hits.df,500)
t1=x[, .(cpgNum = .N,beta=round(sum(V5)*100/(sum(V5)+sum(V6))) ), by = x1[,1] ]
#data table aggregate fast
dat[, .(count = .N, var = sum(VAR)), by = MNTH]

# ENDGAME
# cell_merged: chr,start,end,cpgNum,betaScore1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,b15,b16
}

