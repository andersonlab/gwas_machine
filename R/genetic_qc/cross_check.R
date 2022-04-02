setwd("/lustre/scratch123/hgi/projects/ibdgwas_bioresource/qz2/07_03_2022/data")
batches <- c("b04","b06","b08","b09","b10","b12","b15","b17","b18","b19")

sample<-c()
for(i in batches){
	tmp <- read.table(paste0(i,"/",i,".fam"),stringsAsFactor=F)
	sample <- rbind(sample,tmp)
}

read.table("/lustre/scratch123/hgi/projects/ibdgwas_bioresource/raw/b04_b06_b15_June2021/data_transfer_20210629/IBD_BioRes_phenotypes_updated_20210715.txt",stringsAsFactor=F,sep="~",header=T)->data
rownames(data) <- data$SPid_1

sum(unique(sample$V1) %in% rownames(data))
#[1] 12413

read.csv("/lustre/scratch123/hgi/projects/ibdgwas_bioresource/raw/data_transfer_20220228/V3_vcfs_b04-b19__ids_by_batch_20220228.csv",stringsAsFactor=F)->vcfID

setdiff(vcfID$SPid,sample$V1)


