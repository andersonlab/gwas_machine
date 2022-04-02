args<- commandArgs(T)
wkdir <- as.character(args[1])
filehead <- as.character(args[2])
dosagefile <- as.character(args[3])
outputfile<- as.character(args[4])

setwd(wkdir)

frq<-read.table(paste0(filehead,".frq"),head=T)
var_miss<-read.table(paste0(filehead,".lmiss",head=T)

all<-merge(frq[,c("SNP","MAF","NCHROBS")],var_miss[,c("SNP","F_MISS")],by="SNP",sort=F)

## females, find variants with less % of calls:
ped<-read.table(dosagefile,head=T,check.names=F)

dat<-matrix(nrow=nrow(all),ncol=2)
dat<-as.data.frame(dat)
colnames(dat)<-c("variant","percentage_NA")

for (i in 1:nrow(dat)) {
  tmp<-ped[,6+i,drop=F]
  dat$variant[i]<-gsub("_[A-Z]{1}$","",colnames(tmp))
  dat$percentage_NA[i]<-nrow(tmp[which(is.na(tmp)),,drop=F])/nrow(tmp)
}

all<-all[which(all$SNP %in% dat$variant[which(dat$percentage_NA>0.95)]),]


# keep variants with large number of calls in males
all<-all[which(all$NCHROBS>=max(all$NCHROBS,na.rm=T)-(max(all$NCHROBS,na.rm=T)*0.005)),]
dim(all)

write.table(all,outputfile,col.names=F,row.names=F,quote=F,sep="\t")

