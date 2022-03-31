args<- commandArgs(T)
wkdir <- as.character(args[1])
inputfile<- as.character(args[2])
outputfile<- as.character(args[3])

setwd(wkdir)
data <- read.table(paste0("data/",inputfile),stringsAsFactor=F,header=F)
data$idx <- 1:nrow(data)
rmIDX <- c()

IDcount <- table(data$V1)
dupID <- names(IDcount[IDcount>1])

for(id in dupID){
	datatmp <- data[ data$V1==id,]
	datatmpRM <- datatmp[ datatmp$V3>min(datatmp$V3),]
	rmIDX <- c(rmIDX,datatmpRM$idx)
}

data <- data[!data$idx%in%rmIDX,]
sum(table(data[,1])>1)
#[1] 0

write.table(data[,1:2],file=paste0("data/",outputfile),row.names=F,quote=F,col.names=F)

