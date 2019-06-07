library(methylKit) #t
library(data.table)

setwd("/home/rtm/methmotif_cov/methylkit")

file.list <- list.files("/home/rtm/methmotif_cov/WGBS_MethMotif/",pattern=".destranded.cov",recursive=TRUE,full.names = TRUE)

dir_file <- gsub("/home/rtm/methmotif_cov/WGBS_MethMotif//","",file.list)
dir_file <- gsub("/.+","",dir_file,perl=TRUE)
dir_file_dup <- duplicated(dir_file)
dir_file_dup[which(dir_file_dup)-1] = TRUE
dir_file_dup_ix <- which(dir_file_dup)

for( i in 1:(length(dir_file_dup_ix)/2)){
 
f1 <- file.list[dir_file_dup_ix[(i*2)-1]]
f2 <- file.list[dir_file_dup_ix[(i*2)]]
print(f1)
print(f2)
print("############")

covf1 <- fread( f1,sep="\t", colClasses=c("character","integer","integer","integer","integer","integer"),
               data.table=FALSE, header=FALSE)
covf2 <- fread( f2,sep="\t", colClasses=c("character","integer","integer","integer","integer","integer"),
               data.table=FALSE, header=FALSE)

bed1 <- paste(covf1[,1],covf1[,2],covf1[,3],sep="_")
bed2 <- paste(covf2[,1],covf2[,2],covf2[,3],sep="_")

ix <- match(bed1,bed2)

# Add cov2 intersection counts to file cov1
covf1[!is.na(ix),5] = covf1[!is.na(ix),5] + covf2[ix[!is.na(ix)],5] 
covf1[!is.na(ix),6] = covf1[!is.na(ix),6] + covf2[ix[!is.na(ix)],6] 

# recalculate beta scores
covf1[ !is.na(ix) , 4 ] <- round((covf1[ !is.na(ix) , 5 ] / 
                                  (covf1[ !is.na(ix) , 5 ] + covf1[ !is.na(ix) , 6 ]))*100)
  
# Add cov2 outsection counts to file cov1
jx <- !bed2 %in% bed1 
pooled <- rbind(covf1,covf2[jx,])

# write tmp file
newfile <- gsub(".+//","",f1,perl=TRUE)
newfile <- gsub("/.+","",newfile,perl=TRUE)
newfile <- paste("/home/rtm/methmotif_cov/WGBS_MethMotif/",newfile,"/",newfile,".bismark.destranded.pooled.cov.TMPFILE",sep="")
write.table(pooled, newfile, sep="\t", row.names = FALSE, col.names = FALSE)
# sort cov1 file
newfile_sorted <- gsub(".TMPFILE","",newfile)
command <- paste("sort -k1,1 -k2,2n",newfile,">",newfile_sorted)
system(command)
# delete tmp file
command <- paste("rm",newfile)
system(command)
}
###############################################################################################################################
###############################################################################################################################
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


cells <- gsub("\\/home\\/rtm\\/methmotif_cov\\/WGBS_MethMotif\\/\\/","",file.list)
cells <- gsub("\\/.+","",cells,perl=TRUE)

cell_sample_number <- as.numeric(factor(cells))

cells <- make.unique(cells)
cells <- gsub("\\.1","\\.2",cells)
cells[grep("\\.2",cells)-1] <- paste(cells[grep("\\.2",cells)-1],".1",sep="")

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

###############################################################################################################################
###############################################################################################################################
x=list(myobj[[1]],myobj[[2]])
myMethRawList=new("methylRawList",x,treatment=c(1,1))



