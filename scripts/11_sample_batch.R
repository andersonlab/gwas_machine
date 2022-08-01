args<- commandArgs(T)
inputfile1 <- as.character(args[1])
inputfile2 <- as.character(args[2])

oldbatch <- read.table(inputfile1,stringsAsFactor=F)
newbatch <- read.table(inputfile2,stringsAsFactor=F)

colnames(oldbatch)[2]<- "V3"
colnames(newbatch)[2]<- "V3"

oldbatch$V2 <- oldbatch$V1
newbatch$V2 <- newbatch$V1

oldbatch$V1 <- substr(oldbatch$V1,0,14)
newbatch$V1 <- substr(newbatch$V1,0,14)

ids <- oldbatch[nchar(oldbatch$V1)!=nchar(oldbatch$V2),"V2"]
ids <- unlist(strsplit(ids,split="_"))
oldbatch[nchar(oldbatch$V1)!=nchar(oldbatch$V2),"V2"] <- ids[seq(2,length(ids),2)]

ids <- newbatch[nchar(newbatch$V1)!=nchar(newbatch$V2),"V2"]
ids <- unlist(strsplit(ids,split="_"))
newbatch[nchar(newbatch$V1)!=nchar(newbatch$V2),"V2"] <- ids[seq(2,length(ids),2)]

oldbatch[ oldbatch$V1 %in% newbatch$V1,"V2"]<- 1
newbatch[ newbatch$V1 %in% oldbatch$V1,"V2"]<- 2

batch <- rbind(oldbatch,newbatch)

for(i in c("b04","b06","b15","b19","b20")){
	write.table(batch[batch$V3==i,c("V1","V2")],file=paste0("./data/", i,".sample"),col.names=F,quote=F,row.names=F)
}
