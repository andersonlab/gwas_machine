args<-commandArgs(T)
wkdir <- as.character(args[1])
filehead <- as.character(args[2])

setwd(wkdir)

read.table("/lustre/scratch123/hgi/projects/ibdgwas_bioresource/raw/b04_b06_b15_June2021/data_transfer_20210629/IBD_BioRes_phenotypes_updated_20210715.txt",stringsAsFactor=F,sep="~",header=T)->data
rownames(data) <- data$SPid_1

read.table(paste0("data/",filehead,".fam"),stringsAsFactor=F)->fam

rownames(fam)<- fam$V1
intername<- intersect(rownames(fam),rownames(data))
data <- data[intername,]
fam <- fam[intername,]

fam$V5 <- data$gender
write.table(fam,file=paste0("data/",filehead,"_pheno_gender"),row.names=F,quote=F,col.names=F)
