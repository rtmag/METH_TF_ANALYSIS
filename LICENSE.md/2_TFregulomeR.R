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

command = paste("awk -F\"\t\" \'{if(($2-100)<1){print $1\"\t\"1\"\t\"$3+100\"\t\"$4}else{print $1\"\t\"$2\"\t\"$3+100\"\t\"$4}}\' ",
    file.list,
    " |sort -k1,1 -k2,2n|bedtools merge -i - > ",
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