args<-commandArgs(T)
filehead <- as.character(args[1])

bim<-read.table(paste0(filehead,".bim"),sep="",head=F)

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]

otherchr <- bim[ !(bim$V1 %in% c(1:22,"X","Y")),]

all_remove<-rbind(indels,otherchr)
all_remove<-all_remove[!duplicated(all_remove),]

write.table(all_remove[,"V2"],"./data/snps_to_keep", col.names=F,row.names=F,quote=F,sep="\t")

