
mm <- read.table("mm_tf_matrix.bed",header=TRUE)
coord <- mm[,1:3]

# system("cut -f1,2,3 mm_tf_matrix.bed|bedtools merge -i - > mm_tf_matrix_merged.bed ")
merged <- read.table("mm_tf_matrix_merged.bed",sep="\t",header=FALSE)
