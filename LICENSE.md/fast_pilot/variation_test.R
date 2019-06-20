meth_mat <- readRDS("HyperTFBS_methMatrix.rds")
unmeth_mat <- readRDS("hypoTFBS_methMatrix.rds")

meth <- meth_mat[,4:dim(meth_mat)[2]]
unmeth <- unmeth_mat[,4:dim(unmeth_mat)[2]]

meth.sd<-apply(meth,1,sd)
meth.sd<-apply(meth,1,sd)
meth.sd<-apply(meth,1,sd)

meth vs non-meth
