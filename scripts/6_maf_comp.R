args<- commandArgs(T)
frq1kg <- as.character(args[1])
studyfrq <- as.character(args[2])

gp<-read.table(frq1kg,head=T)
g1<-read.table(studyfrq,head=T)

dim(g1)
dim(gp)

dim(g1[which(!g1$SNP %in% gp$SNP),])

colnames(gp)[3:6]<-paste(colnames(gp)[3:6],"_gp",sep="")
colnames(g1)[3:6]<-paste(colnames(g1)[3:6],"_g1",sep="")

all<-merge(g1,gp[,2:6],by="SNP",all=T)

# variants with different minor allele were checked
check<-all[which(all$A1_g1!=all$A1_gp),]
dim(check)

# Keep only A/T C/G
check<-check[which( (check$A1_g1=="G" & check$A2_g1=="C") | (check$A1_g1=="C" & check$A2_g1=="G") | (check$A1_g1=="A" & check$A2_g1=="T") | (check$A1_g1=="T" & check$A2_g1=="A")),]
dim(check)

# LIST OF VARIANTS TO REMOVE, WE CANNOT REALLY BE SURE WHETHER THERE IS STRAND ISSUE OR NOT
remove<-check[which(check$MAF_g1>=0.45),]
dim(remove)


# LIST OF VARIANTS TO FLIP:
flip<-check[which(check$MAF_g1<0.45),]
dim(flip)

flip<-flip[order(flip$MAF_g1,decreasing=T),]

remove_2<-flip[which(flip$MAF_g1>0.2 & flip$MAF_gp<0.1),]
remove_3<-flip[which(flip$MAF_g1<0.1 & flip$MAF_gp>0.2),]

remove<-rbind(remove,remove_2,remove_3)
dim(remove)

write.table(remove[,"SNP"],"./data/list_variants_to_remove_AT_CG",col.names=F,row.names=F,quote=F,sep="\t")

flip<-flip[which(!flip$SNP %in% remove$SNP),]
dim(flip)

write.table(flip[,"SNP"],"./data/list_variants_to_flip_AT_CG",col.names=F,row.names=F,quote=F,sep="\t")


