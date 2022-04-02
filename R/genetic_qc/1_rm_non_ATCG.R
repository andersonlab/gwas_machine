args<-commandArgs(T)
wkdir <- as.character(args[1])
batch <- as.character(args[2])

setwd(wkdir)

bim<-read.table(paste0("data/",batch,"/",batch,".bim"),sep="",head=F)

indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]

#mt<-bim[which(bim$V1 %in% c(0,26)),]
#chry<-bim[which(bim$V1 %in% c(24)),]

otherchr <- bim[ !(bim$V1 %in% c(1:22,"X","Y")),]

all_remove<-rbind(indels,otherchr)
all_remove<-all_remove[!duplicated(all_remove),]

write.table(all_remove[,"V2"],paste0("data/",batch,"/list_indel_var_exclude_",batch), col.names=F,row.names=F,quote=F,sep="\t")

