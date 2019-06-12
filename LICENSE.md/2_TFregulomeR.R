library(TFregulomeR)
all_record <- TFBSBrowser()

# download all peaks with motif
for (i in all_record$ID){
peak_i = loadPeaks(id=i, includeMotifOnly=T)
write.table(peak_i, paste0(i,'_peaks_with_motif.txt'), sep='\t', col.names=F, row.names=F, quote=F)}

# Downlload all peaks
for (i in all_record$ID){
peak_i = loadPeaks(id=i, includeMotifOnly=F)
write.table(peak_i, paste0(i,'_peaks_all.txt'), sep='\t', col.names=F, row.names=F, quote=F)}

####################################################################################################

setwd("/home/rtm/methmotif_cov/tfregulomeR/mm_tf_matrix/")

file.list <- list.files("/home/rtm/methmotif_cov/tfregulomeR/all_peaks/",pattern="MM1_HSA_*",recursive=TRUE,full.names = TRUE)
file.id <- data.frame( do.call( rbind, strsplit( file.list, '_' ) ) )

# Create sorted and +/- 100 BP files
new_file <- paste(file.id[,5],file.id[,6],file.id[,8],sep="_")

command = paste("awk -F\"\t\" \'{if(($2-100)<1){print $1\"\t\"1\"\t\"$3+100\"\t\"$4}else{print $1\"\t\"$2-100\"\t\"$3+100\"\t\"$4}}\' ",
    file.list,
    "|perl -pe \'s/\\_all\\_peaks.+//g\'|perl -pe \'s/MM1.+\\_//g\' |sort -k1,1 -k2,2n|bedtools merge -i - -c 4 -o distinct > ",
     "/home/rtm/methmotif_cov/tfregulomeR/mm_tf_matrix/",
     new_file,
     sep="")

for( i in 1:length(command) ){ system(command[i]) }

##############################################################
file.list <- list.files("/home/rtm/methmotif_cov/tfregulomeR/mm_tf_matrix/",recursive=TRUE,full.names = TRUE)
file.id <- data.frame( do.call( rbind, strsplit( file.list, '_' ) ) )
names_id <- paste(file.id[,4],file.id[,5],sep="_")
names_id <- gsub("matrix//","",names_id)
names_id <- paste(names_id,collapse=" ")

command <- paste("multiIntersectBed -i", 
                 paste(file.list,collapse=" "),
                " -names ",
                 names_id,
                 " -header > /home/rtm/methmotif_cov/tfregulomeR/mm_tf_matrix.bed" )
                 
system(command)

command <- paste("multiIntersectBed -i", 
                 paste(file.list,collapse=" "),
                " -names ",
                 names_id,
                 "-cluster -header > /home/rtm/methmotif_cov/tfregulomeR/mm_tf_matrix_cluster.bed" )
                 
system(command)
#########################################################################################################
#########################################################################################################
#########################################################################################################
#mm <- read.table(pipe("cut -f1-5 mm_tf_matrix_cluster_methTF_quy.bed"),header=TRUE)
mm_header <- read.table(pipe("head -n1 mm_tf_matrix_cluster_methTF_quy.bed"),header=FALSE)
mm_header <- mm_header[1,]
mm_header <- as.character(unlist(mm_header))
quy <- read.table("methylated_TFBS_quy_inMotif.txt")
quy <- as.character(sort(quy[,1]))
meth_tf_columns <- which(mm_header %in% quy)
meth_tf_columns <- paste(meth_tf_columns,collapse=",")

command = paste("cut -f1-5,",meth_tf_columns ," mm_tf_matrix_cluster_methTF_quy.bed",sep="")
mm <- read.table(pipe(command),header=TRUE)

mm_ov1 <- mm[mm[,4]>1,]

methtf_counts <- rowSums(mm_ov1[,6:dim(mm_ov1)[2]])
#########################################################################################################
# all TB question
(mm[1,5])
