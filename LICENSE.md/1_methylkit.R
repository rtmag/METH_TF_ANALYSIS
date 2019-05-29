library(methylKit)
 
setwd("/home/rtm/methmotif_cov/methylkit")

file.list <- list.files("/home/rtm/methmotif_cov/WGBS_MethMotif/",pattern=".destranded.cov",recursive=TRUE,full.names = TRUE)
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
    
           
