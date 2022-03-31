args<- commandArgs(T)
wkdir <- as.character(args[1])
kingfile <- as.character(args[2])
filehead <- as.character(args[3])
outputfile <- as.character(args[4])

setwd(wkdir)

library(gtools)

kin<-read.table(kingfile,head=T)
table(kin$InfType)
# Duplicated or Monozygotic twin  Dup/MZ
# Parent–offspring                PO
# Full sib                        FS
# 2nd Degree                      2nd

table(cut(kin$Kinship,breaks=c(1,0.354,0.177)))

#kinship coefficient range
# >0.354 duplicate/MZ twin 
# [0.177, 0.354] 1st-degree 
# [0.0884, 0.177] 2nd-degree

dup<-kin[which(kin$Kinship>0.354),c("FID1","FID2","Kinship")]
summary(dup$Kinship)

dup_ids<-c(as.character(dup$FID1),as.character(dup$FID2))
length(dup_ids)
table(length(dup_ids)==(nrow(dup)*2))

dup_ids<-dup_ids[!duplicated(dup_ids)]
length(dup_ids)


#REMOVE ONE IN DUPLICATES WITH LOWER CALL RATES

fam<-read.table(paste0(filehead,".fam"),head=F)
colnames(fam)[5:6]<-c("sex","pheno")
sample_miss<-read.table(paste0(filehead,".imiss"),head=T)

all<-merge(fam[,c(1,5,6)],sample_miss[,c("FID","F_MISS")],by.x="V1",by.y="FID",all.x=T,sort=F)

tmp<-all[which(all$V1 %in% dup_ids),]
dup_ids<-dup_ids[match(tmp$V1,dup_ids)]
rm(tmp)

for (i in 1:length(dup_ids)) {
  tmp1<-dup[which(dup$FID1==dup_ids[i]),]
  tmp2<-dup[which(dup$FID2==dup_ids[i]),]
  colnames(tmp2)[1:2]<-colnames(tmp2)[2:1]
  tmp<-rbind(tmp1,tmp2)
  ids_tmp<-c(as.character(tmp$FID1),as.character(tmp$FID2))
  ids_tmp<-ids_tmp[!duplicated(ids_tmp)]
  #keep same order as in dup_ids and in all
  ids_tmp<-ids_tmp[match(dup_ids[which(dup_ids %in% ids_tmp)],ids_tmp)]
  # get number of possible combinations
  n_possible_combinations<-nrow(permutations(length(ids_tmp), 2))/2
  if ( nrow(dup[which((dup$FID1 %in% ids_tmp) | (dup$FID2 %in% ids_tmp)),])==n_possible_combinations ) {
    # for duplicated samples that have same phenotype and sex, remove the ones with smaller call rate
    data<-as.data.frame(matrix(ncol=1,nrow=length(ids_tmp)))
    data$V1<-ids_tmp
    data<-merge(data,all,by="V1",all.x=T,sort=F)
    if ( (dim(table(data$sex))==1 & dim(table(data$pheno))==1) ) {
      if(!exists("data_remove")) {
        keep_sample<-data$V1[which(data$F_MISS==min(data$F_MISS))][1]
        data_remove<-data[which(!data$V1 %in% keep_sample),]
      } else {
        keep_sample<-data$V1[which(data$F_MISS==min(data$F_MISS))][1]
        data_remove<-rbind(data_remove,data[which(!data$V1 %in% keep_sample),])
      }
    } else {
      data<-as.data.frame(matrix(ncol=1,nrow=length(ids_tmp)))
      data$V1<-ids_tmp
      data<-merge(data,all,by="V1",all.x=T,sort=F)
      # for duplicated samples that do not have same pheno and sex, remove all
      # print("Samples with different sex/pheno:")
      # print(data)
      if(!exists("data_remove")) {
        data_remove<-data
      } else {
        data_remove<-rbind(data_remove,data)
      }
      if(!exists("data_inconsist")){
        
        jj<-1
        data_inconsist<-data
        data_inconsist$group<-jj
        
      } else {
        jj<-jj+1
        data$group<-jj
        data_inconsist<-rbind(data_inconsist,data)
        
      }
    }
    
  } else {
    
    # not all combinations of duplicated pairs found, show list of IDs, to manually inspect issues
    
    print(paste("Number of expected combinations: ",n_possible_combinations,sep=""))
    print(paste("Number of observed combinations: ",nrow(dup[which((dup$FID1 %in% ids_tmp) | (dup$FID2 %in% ids_tmp)),]),sep=""))
    print(ids_tmp)
    
  }
  
}

  
data_remove<-data_remove[!duplicated(data_remove$V1),]
dim(data_remove)


data_remove<-data_remove[,c(1,1)]
colnames(data_remove)<-c("FID","IID")

write.table(data_remove,outputfile,col.names=T,row.names=F,quote=F,sep="\t")



