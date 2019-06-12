
mm <- read.table("mm_tf_matrix.bed",header=TRUE,stringsAsFactors=FALSE)
coord <- mm[,1:3]

# system("cut -f1,2,3 mm_tf_matrix.bed|bedtools merge -i - > mm_tf_matrix_merged.bed ")
merged <- read.table("mm_tf_matrix_merged.bed",sep="\t",header=FALSE,stringsAsFactors=FALSE)
colnames(merged) <- c("chr","start","end")

library(rtracklayer)
coord.gr <- makeGRangesFromDataFrame(coord) 
merged.gr <- makeGRangesFromDataFrame(merged) 
hits <- findOverlaps(coord.gr, merged.gr)
hits.df <- as.data.frame(hits)

library(data.table)
mm.dt <- as.data.table(mm)

merged_data <- data.frame(mm[1,],stringsAsFactors=FALSE)
merged_data[,1] <- as.character(merged_data[,1])
merged_data[,5] <- as.character(merged_data[,5])

colnames(merged) <- c("chrom","start","end")

for( i in 1:max(tail(hits.df[,2])) ){
print(i)
mm_i = mm.dt[hits.df[,2]==i,]
mm_i <- as.data.frame(mm_i)
TF <- unique(sort(unlist(strsplit(as.character(mm_i[,5]), ","))))
tfmat <- as.numeric(colSums(mm_i[,6:dim(mm_i)[2]])>0)
tfmat <- rbind(tfmat)
colnames(tfmat) = colnames(merged_data)[6:dim(merged_data)[2]]
rownames(tfmat) = ""
d1.tmp <- data.frame( merged[i,],num=length(TF),list=paste(TF,collapse=",") ,stringsAsFactors=FALSE)
d2.tmp <- cbind(d1.tmp,tfmat)
merged_data <- rbind(merged_data,d2.tmp,stringsAsFactors=FALSE)
}

