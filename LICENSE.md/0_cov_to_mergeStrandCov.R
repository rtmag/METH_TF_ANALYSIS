library(data.table)
library(bedr)

covfile <- fread( "zcat /home/rtm/methmotif_cov/WGBS_MethMotif/A549/ENCFF003JVR.bismark.cov.gz| sort -k1,1 -k2,2n",
                        sep="\t", colClasses=c("character","integer","integer","integer","integer","integer"),
                        data.table=FALSE, header=FALSE)
                        
bed <- covfile[,1:3]
bed[,2] <- bed[,2]-1
                        
fasta <- get.fasta(bed,fasta = "/home/references/hg38/hg38.fasta",bed12 = FALSE,strand = FALSE,output.fasta = FALSE,
use.name.field = FALSE,check.zero.based = FALSE,check.chr = FALSE,check.valid = FALSE,
check.sort = FALSE,check.merge = FALSE,verbose = TRUE)

fasta_seq <- paste(fasta$sequence,collapse="")

fasta$sequence == "g"

# identify the ones that need to be substracted

# identify the ones that need to be pooled

# delete the Gs that have been pooled

