args<- commandArgs(T)
cohorts <- as.character(args[1])
outputfile1 <- as.character(args[2])
outputfile2 <- as.character(args[3])

cohorts <- unlist(strsplit(cohorts,spli=" "))

for (i in 1:length(cohorts)){
  if(i==1){
    tmp<-read.table(paste("./data/tmp_",cohorts[i],".lmiss",sep=""),head=T)
    a1<-as.data.frame(table(cut(tmp$F_MISS,breaks=c(-1,0.05,0.1,0.2,0.3,0.4),labels=c("0-0.05","0.05-0.1","0.1-0.2","0.2-0.3","0.3-0.4"))) )
    colnames(a1)<-c("numbers",cohorts[i])
  }else{
    tmp<-read.table(paste("./data/tmp_",cohorts[i],".lmiss",sep=""),head=T)
    a2<-as.data.frame(table(cut(tmp$F_MISS,breaks=c(-1,0.05,0.1,0.2,0.3,0.4),labels=c("0-0.05","0.05-0.1","0.1-0.2","0.2-0.3","0.3-0.4"))) )
    colnames(a2)<-c("numbers",cohorts[i])
    a1<-merge(a1,a2,by="numbers",sort=F)
  }
}

write.table(a1,outputfile1,col.names=T,row.names=F,quote=F,sep=",")


for (i in 1:length(cohorts)){
  if (i==1) {
    tmp<-read.table(paste("./data/tmp_",cohorts[i],".lmiss",sep=""),head=T)
    tmp$study<-cohorts[i]
    tmp<-tmp[,c("SNP","F_MISS","study")]
    all<-tmp
  } else {
    tmp<-read.table(paste("./data/tmp_",cohorts[i],".lmiss",sep=""),head=T)
    tmp$study<-cohorts[i]
    tmp<-tmp[,c("SNP","F_MISS","study")]
    all<-rbind(all,tmp)
  }
}

table(all$study)

variants<-all[which(all$F_MISS>0.1),"SNP",drop=F]
variants<-variants[!duplicated(variants$SNP),,drop=F]
dim(variants)

write.table(variants,outputfile2,col.names=F,row.names=F,quote=F)
