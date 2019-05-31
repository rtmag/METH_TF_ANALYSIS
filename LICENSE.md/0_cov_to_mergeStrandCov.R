library(data.table)
library(bedr)
library(stringr)

file.list <- list.files("/home/rtm/methmotif_cov/WGBS_MethMotif/",pattern="*cov.gz",recursive=TRUE,full.names = TRUE)

for(i in 1:length(file.list)){

print(i)
print(file.list[i])
command <- paste("zcat", file.list[i], "| sort -k1,1 -k2,2n")
covfile <- fread( command,
                  sep="\t", colClasses=c("character","integer","integer","integer","integer","integer"),
                  data.table=FALSE, header=FALSE)
                        
bed <- covfile[,1:3]
bed[,2] <- bed[,2]-1

print("FASTA")
fasta <- get.fasta(bed,fasta = "/home/references/hg38/hg38.fasta",bed12 = FALSE,strand = FALSE,output.fasta = FALSE,
use.name.field = FALSE,check.zero.based = FALSE,check.chr = FALSE,check.valid = FALSE,
check.sort = FALSE,check.merge = FALSE,verbose = TRUE)

fasta_seq <- paste(fasta$sequence,collapse="")
fasta_seq <-tolower(fasta_seq)
# identify the ones that need to be pooled "CG"
fasta_cg_ix <- str_locate_all(pattern ='cg', fasta_seq)
fasta_cg_ix <- fasta_cg_ix[[1]]

covfile[ fasta_cg_ix[,1] , 5 ] <- covfile[ fasta_cg_ix[,1] , 5 ] + covfile[ fasta_cg_ix[,2] , 5 ]
covfile[ fasta_cg_ix[,1] , 6 ] <- covfile[ fasta_cg_ix[,1] , 6 ] + covfile[ fasta_cg_ix[,2] , 6 ]

covfile[ fasta_cg_ix[,1] , 4 ] <- round((covfile[ fasta_cg_ix[,1] , 5 ] / 
                                        (covfile[ fasta_cg_ix[,1] , 5 ] + covfile[ fasta_cg_ix[,1] , 6 ]))*100)
# identify the ones that need to be coordinate substracted "G"
fasta_g_ix <- str_locate_all(pattern ='g', head(fasta_seq))
fasta_g_ix <- fasta_g_ix[[1]]
fasta_g_ix <- fasta_g_ix[ !fasta_g_ix[,2] %in% fasta_cg_ix[,2], ]

covfile[ fasta_g_ix[,2] , 2 ] <- covfile[ fasta_g_ix[,2] , 2 ]-1
covfile[ fasta_g_ix[,2] , 3 ] <- covfile[ fasta_g_ix[,2] , 3 ]-1

# delete the Gs that have been pooled into "CG"
covfile = covfile[ -fasta_cg_ix[,2] ,]

print("writting")
newfile <- file.list[i]
newfile <- gsub("\\.cov\\.gz","\\.destranded\\.cov",newfile)
fwrite(covfile, newfile, sep="\t", row.names = FALSE, col.names = FALSE,buffMB=1000,nThread=12)

# Remove the 
rm(list = c("covfile","bed","fasta","fasta_seq","fasta_cg_ix","fasta_g_ix","newfile"))
}

