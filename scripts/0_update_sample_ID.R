args<-commandArgs(T)
filehead <- as.character(args[1])
outputfile<- as.character(args[2])

fam<-read.table(paste0(filehead,".fam"),sep="",head=F)

fidnew <- gsub("([1-9]|1[0-1]):", "", fam$V1)
iidnew <- fam$V2
iidnew[nchar(fam$V2)==17] <- 2

fid_replace <- fidnew[nchar(fam$V2)==17]

iidnew[fidnew%in%fid_replace & nchar(iidnew)==14] <- 1

idmap<- cbind(fam[,c("V1","V2")],fidnew,iidnew)

write.table(idmap,paste0("data/",outputfile), col.names=F,row.names=F,quote=F,sep="\t")


