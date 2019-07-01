library(methylKit) #t
options(scipen=999)

setwd("/home/rtm/methmotif_cov/methylkit")
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


           
meth=unite(myobj, destrand=FALSE,mc.cores=20,allow.cartesian=TRUE)
mydiff=calculateDiffMeth(meth,num.cores=22)
saveRDS(mydiff,"mydiff.rds")
################################################
tiles=tileMethylCounts(myobj,win.size=1000,step.size=1000)
meth_tile=unite(tiles, destrand=FALSE,mc.cores=20,allow.cartesian=TRUE)
mydiff_tile=calculateDiffMeth(meth_tile,num.cores=22)
saveRDS(mydiff_tile,"mydiff_tile.rds")
################################################################################################
################################################################################################
################################################################################################


i_matrix = combn(1:16,2)

for(ix in 1:dim(i_matrix)[2]){
  meth_tmp = reorganize(meth_tile,
                        sample.ids=meth_tile@sample.ids[ c(i_matrix[1,ix],i_matrix[2,ix]) ],
                        treatment=c(i_matrix[1,ix],i_matrix[2,ix]) )
  
  mydiff_tmp = calculateDiffMeth(meth_tmp,num.cores=22)
  
  filename = paste(meth_tile@sample.ids[ i_matrix[1,ix] ],
                   "_VS_",
                   meth_tile@sample.ids[ i_matrix[2,ix] ],
                   ".tiles.mydiff.txt",sep="")
  write.table(mydiff_tmp,filename,sep="\t", row.names = FALSE, col.names = FALSE)
}

####

testdata = tiles[[1]][ as.character(getData(tiles[[1]])[,1]) %in% c(paste("chr",1:22,sep=""),"chrX","chrY"), ]

table(as.character(getData(tiles[[1]])[,1]))

pdf("test_plots_methseg.pdf")
 test_methSeg=methSeg(testdata ,diagnostic.plot=TRUE,maxInt=100,minSeg=10)
dev.off()

methseg = list()
for(i in 1:16){
  tiles_mainchr <- tiles[[i]][ as.character(getData(tiles[[i]])[,1]) %in% c(paste("chr",1:22,sep=""),"chrX","chrY"), ]
  filename <- paste(tiles[[i]]@sample.id,"_Segmentation_plots.pdf",sep="")
  pdf(filename)
  tiles_methSeg=methSeg(tiles_mainchr ,diagnostic.plot=TRUE,maxInt=100,minSeg=10,G=1:2)
  dev.off()
  methseg[[tiles[[i]]@sample.id]] <- tiles_methSeg
  }
  
methseg_4gp = list()
for(i in 1:16){
  tiles_mainchr <- tiles[[i]][ as.character(getData(tiles[[i]])[,1]) %in% c(paste("chr",1:22,sep=""),"chrX","chrY"), ]
  filename <- paste(tiles[[i]]@sample.id,"_Segmentation_plots_4gp.pdf",sep="")
  pdf(filename)
  tiles_methSeg=methSeg(tiles_mainchr ,diagnostic.plot=TRUE,maxInt=100,minSeg=10,G=1:4)
  dev.off()
  methseg[[tiles[[i]]@sample.id]] <- tiles_methSeg
  }



methseg_4gp_cpg = list()
for(i in 1:16){
  tiles_mainchr <- myobj[[i]][ as.character(getData(myobj[[i]])[,1]) %in% c(paste("chr",1:22,sep=""),"chrX","chrY"), ]
  filename <- paste("No_gp_constrain_",tiles[[i]]@sample.id,"_Segmentation_plots_cpg.pdf",sep="")
  pdf(filename)
  tiles_methSeg=methSeg(tiles_mainchr ,diagnostic.plot=TRUE,maxInt=100,minSeg=10,G=1:4)
  dev.off()
  methseg_4gp_cpg[[tiles[[i]]@sample.id]] <- tiles_methSeg
  }



methseg_gp_cpg = list()
for(i in 1:16){
  tiles_mainchr <- myobj[[i]][ as.character(getData(myobj[[i]])[,1]) %in% c(paste("chr",1:22,sep=""),"chrX","chrY"), ]
  filename <- paste("No_gp_constrain_",tiles[[i]]@sample.id,"_Segmentation_plots_cpg.pdf",sep="")
  pdf(filename)
  tiles_methSeg=methSeg(tiles_mainchr ,diagnostic.plot=TRUE,maxInt=100,minSeg=10)
  dev.off()
  methseg_gp_cpg[[tiles[[i]]@sample.id]] <- tiles_methSeg
  }
####################################################################################################################
####################################################################################################################
methseg_5gp_cpg = list()
for(i in 1:16){
  tiles_mainchr <- myobj[[i]][ as.character(getData(myobj[[i]])[,1]) %in% c(paste("chr",1:22,sep=""),"chrX","chrY"), ]
  filename <- paste("cpg_gp_1to5_",tiles[[i]]@sample.id,"_Segmentation_plots_cpg.pdf",sep="")
  pdf(filename)
  tiles_methSeg=methSeg(tiles_mainchr ,diagnostic.plot=TRUE,maxInt=100,minSeg=10,G=1:5)
  dev.off()
  methseg_5gp_cpg[[tiles[[i]]@sample.id]] <- tiles_methSeg
  }
####################################################################################################################
####################################################################################################################
test <- makeGRangesFromDataFrame(methseg_5gp_cpg[[1]])
#
hypo_list = list()
for(i in 1:16){
      hypo_list[[tiles[[i]]@sample.id]] <- methseg_5gp_cpg[[i]][methseg_5gp_cpg[[i]]$seg.mean<10,]
}

for(i in 1:16){print(length(hypo_list[[i]])) }

hypo_uniq_list = hypo_list
for(i in 1:16){
  for(j in 1:16){
    if(i!=j){
        hits <- findOverlaps(hypo_list[[i]], hypo_list[[j]])
        hits.df <- as.data.frame(hits)
        matched.pos <- unique(hits.df[,1])
        hypo_uniq_list[[i]] <- hypo_uniq_list[[i]][-matched.pos,]
      }
  }
}

for(i in 1:16){print(length(hypo_uniq_list[[i]])) }

##
x=0
for(i in 1:16){
  for(j in (i+1):16){
    print(c(tiles[[i]]@sample.id,tiles[[j]]@sample.id))
    x=x+1
    }
  }
####################################################################################################################
####################################################################################################################
for(i in 1:16){ 
  gr <- hypo_list[[i]]
  df <- data.frame(seqnames=seqnames(gr),
  starts=start(gr)-1,
  ends=end(gr),
  names=c(rep(".", length(gr))),
  scores=c(rep(".", length(gr))),
  strands=strand(gr))

  file <- paste0("/home/rtm/methmotif_cov/methylkit/hypometh_regions/",tiles[[i]]@sample.id,"_hypomethylated_regions.bed")
  write.table(df, file, quote=F, sep="\t", row.names=F, col.names=F)
}




