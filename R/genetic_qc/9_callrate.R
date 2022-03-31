args<- commandArgs(T)
wkdir <- as.character(args[1])
filehead <- as.character(args[2])
outputfile <- as.character(args[3])
setwd(wkdir)

sample_miss<-read.table(paste0(filehead,".imiss"),head=T)
var_miss<-read.table(paste0(filehead,".lmiss"),head=T)
frq<-read.table(paste0(filehead,".frq"),head=T)

var<-merge(frq[,c(2:6)],var_miss,by="SNP")
var.1<-var[which(var$MAF<0.01 & var$F_MISS>0.02),]
dim(var.1)

write.table(var.1[,"SNP",drop=F],outputfile,col.names=F,row.names=F,quote=F,sep="\t")


