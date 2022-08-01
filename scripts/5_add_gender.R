library(data.table)

args<-commandArgs(T)
filehead <- as.character(args[1])
phenofile <- as.character(args[2])

fread(phenofile,stringsAsFactor=F,sep="~",header=T)->data
data<- as.data.frame(data)

read.table(paste0(filehead,".fam"),stringsAsFactor=F)->fam
rownames(data)<- data$SPid_1

fam$V5 <- data[fam$V1,"gender"]

## 1: male, 2: female
fam[is.na(fam$V5),"V5"]<- 0

table(fam$V5)

write.table(fam,file=paste0(filehead,"_gender"),row.names=F,quote=F,col.names=F)
