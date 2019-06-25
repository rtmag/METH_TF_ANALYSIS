# 1)Meth Variance in  methylated regions bound by any TF compared to regions but by TF but not menthylated
library(data.table)
library(rtracklayer)
options(scipen=999)

setwd("/home/rtm/methmotif_cov/methdiff_tiles_hypo_hyper")

file.list <- list.files("/home/rtm/methmotif_cov/methylkit/",pattern="*tiles.mydiff.txt",recursive=TRUE,full.names = TRUE)

diff_file <- read.table(file.list[1])

for( i in 1:file


