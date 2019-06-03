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

covf1 <- fread( f1,sep="\t", colClasses=c("character","integer","integer","integer","integer","integer"),
               data.table=FALSE, header=FALSE)
covf2 <- fread( f2,sep="\t", colClasses=c("character","integer","integer","integer","integer","integer"),
               data.table=FALSE, header=FALSE)

bed1 <- paste(covf1[,1],covf1[,2],covf1[,3],sep="_")
bed2 <- paste(covf2[,1],covf2[,2],covf2[,3],sep="_")

ix = match(bed1,bed2)
 
covf1[!is.na(ix),] 
covf2[ix[!is.na(ix)],] 

 
 
 
!bed2 %in% bed1

t1 = c(1,3,4,5,6,7,2,8,9,10)
t2 = c(5,6,8,9,10,20,7)
tix= match(t1,t2)
# for first file
t1[!is.na(tix)]
# for second file
t2[tix[!is.na(tix)]]

}
###############################################################################################################################
###############################################################################################################################

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
           


