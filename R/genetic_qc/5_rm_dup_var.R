args<-commandArgs(T)
wkdir <- as.character(args[1])
batch<- as.character(args[2])

setwd(paste0(wkdir,"/data/",batch))

bim<-read.table(paste0(batch,"_ATCG_aligned.bim"),sep="\t",head=F)
ids<-read.table(paste0(batch,"_ATCG_aligned"),header=F,stringsAsFactor=F) #skip line number might be differ

varmiss<-read.table(paste0(batch,"_ATCG_aligned.lmiss"),sep="",head=T)

colnames(ids)[2]<-"ids"

table(bim$V2==ids$V1)

#the major allele is set to A2 by default by Plink, keep ids with real ref/alt as in vcf using the ids file
bim.1<-cbind(bim,ids[,"ids",drop=F])

bim.1<-merge(bim.1,varmiss[,c("SNP","F_MISS")],by.x="V2",by.y="SNP",sort=F)
identical(bim.1$V2,bim$V2)
#[1] TRUE

bim.1$ids<-as.character(bim.1$ids)

# identify duplicated variants (same chr position ref and alt)
dups<-bim.1[which(duplicated(bim.1$ids)),"ids"]
dups <- unique(dups)
#for (i in 1:length(dups)){
#  tmp<-bim.1[which(bim.1$ids %in% dups[i]),]
#  keep<-tmp[which(tmp$F_MISS==min(tmp$F_MISS)),]
#  if (nrow(keep)>1){
#    keep<-keep[1,]
#  }
#  exclude<-tmp[which(!tmp$V2 %in% keep$V2),]
#  bim.1$ids[which(bim.1$V2 %in% exclude$V2)]<-paste(bim.1$ids[which(bim.1$V2 %in% exclude$V2)],"_rm",sep="")
#  if(i==1000){
#  	print(i)
#  }
#}


#Qian's version, the above one is time-comsuming
bim.1$V2 <- as.character(bim.1$V2)
bim.1$idx <- 1:nrow(bim.1)

rm_dup <- function(bim,dupvar){
	tmp<-bim[bim$ids == dupvar,]
  	keep<-tmp[tmp$F_MISS==min(tmp$F_MISS),]
  	#if there are many duplicated variants with same missing rates, select the first one
  	exclude<-tmp[tmp$V2!= keep[1,"V2"],"idx"]
  	return(exclude)
}

bim.1.dup <- bim.1[bim.1$ids%in%dups,]
excludes <- lapply(dups,function(x) rm_dup(bim.1.dup,x))
excludes <- unlist(excludes)

bim.1$ids[which(bim.1$idx %in% excludes)]<-paste(bim.1$ids[which(bim.1$idx %in% excludes)],"_rm",sep="")

nrow(bim.1[which(duplicated(bim.1$ids)),])

duplicated_variants<-bim.1[grep("_rm",bim.1$ids),"ids"]
length(duplicated_variants)

write.table(duplicated_variants,paste0("list_duplicated_var_exclude_",batch),col.names=F,row.names=F,quote=F,sep="\t")
write.table(bim.1[,c(2,7,3:6)],paste0(batch,"_ATCG_aligned_edited.bim"),col.names=F,row.names=F,quote=F,sep="\t")

