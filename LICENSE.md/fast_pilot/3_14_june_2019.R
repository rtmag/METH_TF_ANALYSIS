 # 1)Meth Variance in  methylated regions bound by any TF compared to regions but by TF but not menthylated
library(data.table)
library(rtracklayer)
options(scipen=999)

system("mkdir /home/rtm/methmotif_cov/cell_cluster_methylation")

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


file.list <- list.files("/home/rtm/methmotif_cov/tfregulomeR/mm_tf_matrix/",pattern="_all.txt*",recursive=TRUE,full.names = TRUE)
file.id <- data.frame( do.call( rbind, strsplit( file.list, '//' ) ) )
file.id <- data.frame( do.call( rbind, strsplit( as.character(file.id[,2]), '_' ) ) )

cells <- as.character(unique(file.id[,1]))

for( i in 1:length(cells)){
  file.cell <- list.files("/home/rtm/methmotif_cov/tfregulomeR/mm_tf_matrix/",pattern=paste("*",cells[i],"*",sep=""))
  command = paste("cat",paste(file.cell,collapse=" "),"| sort -k1,1 -k2,2n|mergeBed -i - -c 4 -o distinct")
  cell_merged <- read.table(pipe(command),stringsAsFactors=FALSE,sep="\t")
  colnames(cell_merged) <- c("chr","start","end","TF")
  cell_merged.gr <- makeGRangesFromDataFrame(cell_merged[,1:3]) 
  
  hits <- findOverlaps(cell_merged.gr, wgbs.gr[[cells[i]]])
  hits.df <- as.data.frame(hits)
  
  #save bed clusters without CpGs in its own WGBS
  no_cpg_cell_merged <- cell_merged[which(!1:dim(cell_merged)[1] %in% hits.df[,1]),]
  write.table(no_cpg_cell_merged,
              paste("/home/rtm/methmotif_cov/cell_cluster_methylation_2nd/",cells[i],"_TFclusters_withoutCpG.bed",sep=""),
             sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
  #update merge bed clusters with cpgs in its WGBS
  cpg_cell_merged_ori <- cell_merged[unique(hits.df[,1]),]
  cell_merged.gr <- makeGRangesFromDataFrame(cpg_cell_merged_ori)
  
  # generate pre append object
  cpg_cell_merged <- cbind(cpg_cell_merged_ori,CpGnum=0,ReadNum=0)
  
    print(i)
    
    prev_names = colnames(cpg_cell_merged)
    cpg_cell_merged <- cbind(cpg_cell_merged,newbeta=NA)

    hits <- findOverlaps(cell_merged.gr, wgbs.gr[[cells[i]]])
    hits.df <- as.data.frame(hits)
    if(is.unsorted(hits.df[,1])){ print("hits1 is unsorted") }
    if(is.unsorted(hits.df[,2])){ print("hits2 is unsorted") }
#
    cpgs_in_bed=wgbs[[cells[i]]][hits.df[,2]]
#
   cpgs_in_bed.dt = cpgs_in_bed[, .(CpGnum = .N,ReadNum=(sum(V5)+sum(V6)),beta=round(sum(V5)*100/(sum(V5)+sum(V6))) ), by = hits.df[,1]]
#
 cpg_cell_merged$CpGnum <- cpgs_in_bed.dt$CpGnum ; 
    cpg_cell_merged$ReadNum <- cpgs_in_bed.dt$ReadNum 
    cpg_cell_merged$newbeta[cpgs_in_bed.dt$hits.df] <- cpgs_in_bed.dt$beta
    colnames(cpg_cell_merged) <- c(prev_names,paste(cells[i],"_beta",sep="") )
  
  
  write.table(cpg_cell_merged,
              paste("/home/rtm/methmotif_cov/cell_cluster_methylation_2nd/",cells[i],"_TFclusters_beta.bed",sep=""),
             sep="\t", row.names = FALSE, col.names = TRUE, quote=FALSE)
}


file.list <- list.files("/home/rtm/methmotif_cov/cell_cluster_methylation_2nd/",pattern="*TFclusters_beta.bed",
                        recursive=TRUE,full.names = TRUE)

file.id <- data.frame( do.call( rbind, strsplit( file.list, '//' ) ) )
file.id <- data.frame( do.call( rbind, strsplit( as.character(file.id[,2]), '_' ) ) )

cells <- as.character(unique(file.id[,1]))

# merge all TF clusters
file.cell <- list.files("/home/rtm/methmotif_cov/cell_cluster_methylation_2nd/",pattern=paste(cells[i],"_TFclusters_beta.bed",sep=""))
 
command = "cat /home/rtm/methmotif_cov/cell_cluster_methylation_2nd/*TFclusters_beta.bed|awk -F'\t' '{if($7>79){print $0}}'|grep -v 'CpGnum'| sort -k1,1 -k2,2n|mergeBed -i - "
hyper_merged <- read.table(pipe(command),stringsAsFactors=FALSE,sep="\t")
colnames(hyper_merged) <- c("chr","start","end")
hyper_merged.gr <- makeGRangesFromDataFrame(hyper_merged[,1:3]) 
hyper_cpg_merged <- cbind(hyper_merged)

# get the WGBS for the merged final bed one
  for( j in 1:length(cells) ){
    print(j)
    
    prev_names = colnames(hyper_cpg_merged)
    hyper_cpg_merged <- cbind(hyper_cpg_merged,newbeta=NA)

    hits <- findOverlaps(hyper_merged.gr, wgbs.gr[[cells[j]]])
    hits.df <- as.data.frame(hits)
    if(is.unsorted(hits.df[,1])){ print("hits1 is unsorted") }
    if(is.unsorted(hits.df[,2])){ print("hits2 is unsorted") }

    cpgs_in_bed=wgbs[[cells[j]]][hits.df[,2]]
    cpgs_in_bed.dt = cpgs_in_bed[, .(beta=round(sum(V5)*100/(sum(V5)+sum(V6))) ), by = hits.df[,1] ]
    hyper_cpg_merged$newbeta[cpgs_in_bed.dt$hits.df] <- cpgs_in_bed.dt$beta
    colnames(hyper_cpg_merged) <- c(prev_names,paste(cells[j],"_beta",sep="") )
  }

saveRDS(hyper_cpg_merged,"/home/rtm/methmotif_cov/cell_cluster_methylation_2nd/HyperTFBS_methMatrix.rds")



