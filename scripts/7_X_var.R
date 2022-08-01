args<- commandArgs(T)
hwe_p <- as.double(args[1])
filehead <- as.character(args[2])
outputfile<- as.character(args[3])

bim<-read.table(paste0(filehead,".bim"),sep="\t",head=F)
table(bim$V1)

hwe<-read.table(paste0(filehead,".hwe"),head=T)

#this is commented since we only have cases
#hwe<-hwe[which(hwe$TEST=="UNAFF"),]

frq<-read.table(paste0(filehead,".frq"),sep="",head=T)
var_miss<-read.table(paste0(filehead,".lmiss"),sep="",head=T)

all<-merge(hwe[,c("SNP","P")],frq[,c("SNP","MAF")],by="SNP",sort=F)
all<-merge(all,var_miss[,c("SNP","F_MISS")],by="SNP",sort=F)

all<-all[which(all$MAF>0.05 & all$F_MISS<0.01 & all$P>hwe_p),]
dim(all)

write.table(all,outputfile,col.names=F,row.names=F,quote=F,sep="\t")

