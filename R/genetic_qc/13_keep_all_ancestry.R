args<- commandArgs(T)
wkdir <- as.character(args[1])
pop1kg<- as.character(args[2])
filehead <- as.character(args[3])
pcahead <- as.character(args[4])

setwd(wkdir)

read.table(pop1kg,header=T,stringsAsFactors=F)->samples1000
read.table(paste0(filehead,".fam"),header=F,stringsAsFactor=F)->fam

colnames(fam)[1]<-c("sample")
fam$super_pop<- NA
fam$pop<- NA
fam<- fam[,c(1,7,8)]

samples1000<-samples1000[,c("sample","super_pop","pop")]
fam<-rbind(fam,samples1000)

pca<-read.table(paste0(pcahead,".eigenvec"),head=F)
colnames(pca)<-c("FID","IID",paste0("PC",c(1:20)))

pca<-merge(pca,fam,by.x="IID",by.y="sample",all.x=T)

# TOTAL VARIANCE EXPLAINED BY EACH PC:

eigenval<-read.table(paste0(pcahead,".eigenval"),head=F)

eigenval$var_exp<-NA
for (i in 1:nrow(eigenval)){
  eigenval$var_exp[i]<-(eigenval$V1[i] / sum(eigenval$V1))*100
}

library(ggplot2)
library(cowplot)
library(gridExtra)

eigenval$PC<- 1:20
eigenval$cumval<- cumsum(eigenval$var_exp)
p1<-ggplot(data=eigenval,aes(x=PC,y=var_exp))+geom_bar(stat="identity")+theme_bw()+theme(axis.title=element_text(size=14,face="bold"),axis.text=element_text(size=14,face="bold"))+xlab("Principal Component (PC)")+ylab("Proportion of variance explained by each PC (%)")

p2<-ggplot(data=eigenval,aes(x=PC,y=cumval))+geom_point(shape=16)+geom_line()+theme_bw()+theme(axis.title=element_text(size=14,face="bold"),axis.text=element_text(size=14,face="bold"))+xlab("Principal Component (PC)")+ylab("Cumulative proportion of variance explained by top PCs (%)")

p3<-arrangeGrob(p1, p2,ncol=2, nrow = 1, widths = c(2.7, 3), heights = 2.3)

ggsave("../figures/pca_var.png",dpi=300,width=13,height=6,p3)


# INFER POPULATIONS:

pop<-c("EUR","AMR","AFR","EAS","SAS")

pca$inferred_population<-NA

for (i in 1:length(pop)) {
  
  pca$inferred_population[which(pca$PC1>=min(pca$PC1[which(pca$super_pop==pop[i])]) & pca$PC1<=max(pca$PC1[which(pca$super_pop==pop[i])]) &
                                  pca$PC2>=min(pca$PC2[which(pca$super_pop==pop[i])]) & pca$PC2<=max(pca$PC2[which(pca$super_pop==pop[i])]) &
                                  pca$PC3>=min(pca$PC3[which(pca$super_pop==pop[i])]) & pca$PC3<=max(pca$PC3[which(pca$super_pop==pop[i])]) &
                                  pca$PC4>=min(pca$PC4[which(pca$super_pop==pop[i])]) & pca$PC4<=max(pca$PC4[which(pca$super_pop==pop[i])]) &
                                  pca$PC5>=min(pca$PC5[which(pca$super_pop==pop[i])]) & pca$PC5<=max(pca$PC5[which(pca$super_pop==pop[i])]) &
                                  pca$PC6>=min(pca$PC6[which(pca$super_pop==pop[i])]) & pca$PC6<=max(pca$PC6[which(pca$super_pop==pop[i])]) &
                                  pca$PC7>=min(pca$PC7[which(pca$super_pop==pop[i])]) & pca$PC7<=max(pca$PC7[which(pca$super_pop==pop[i])]) &
                                  pca$PC8>=min(pca$PC8[which(pca$super_pop==pop[i])]) & pca$PC8<=max(pca$PC8[which(pca$super_pop==pop[i])]) &
                                  pca$PC9>=min(pca$PC9[which(pca$super_pop==pop[i])]) & pca$PC9<=max(pca$PC9[which(pca$super_pop==pop[i])]) &
                                  pca$PC10>=min(pca$PC10[which(pca$super_pop==pop[i])]) & pca$PC10<=max(pca$PC10[which(pca$super_pop==pop[i])]) )]<-pop[i]
  
}


table(pca$super_pop,pca$inferred_population,useNA="ifany")

# EXPAND LIMITS FOR EUR SUBSET OF SAMPLES

pc1max<-max(pca$PC1[which(pca$super_pop=="EUR")])+((max(pca$PC1[which(pca$super_pop=="EUR")])-min(pca$PC1[which(pca$super_pop=="EUR")]))*0.3)
pc1min<-min(pca$PC1[which(pca$super_pop=="EUR")])-((max(pca$PC1[which(pca$super_pop=="EUR")])-min(pca$PC1[which(pca$super_pop=="EUR")]))*0.3)

pc2max<-max(pca$PC2[which(pca$super_pop=="EUR")])+((max(pca$PC2[which(pca$super_pop=="EUR")])-min(pca$PC2[which(pca$super_pop=="EUR")]))*0.3)
pc2min<-min(pca$PC2[which(pca$super_pop=="EUR")])-((max(pca$PC2[which(pca$super_pop=="EUR")])-min(pca$PC2[which(pca$super_pop=="EUR")]))*0.3)

pc3max<-max(pca$PC3[which(pca$super_pop=="EUR")])+((max(pca$PC3[which(pca$super_pop=="EUR")])-min(pca$PC3[which(pca$super_pop=="EUR")]))*0.3)
pc3min<-min(pca$PC3[which(pca$super_pop=="EUR")])-((max(pca$PC3[which(pca$super_pop=="EUR")])-min(pca$PC3[which(pca$super_pop=="EUR")]))*0.3)

pcaplot<- pca
colnames(pcaplot)[23]<-"Population"
pcaplot[is.na(pcaplot$Population),"Population"]<- "This study"

p1<-ggplot(pcaplot,aes(x=PC1,y=PC2,col=Population))+geom_point()+theme_bw()+theme(axis.title=element_text(size=12,face="bold"),axis.text=element_text(size=12,face="bold"),legend.position="none")+geom_rect(xmin=pc1min,xmax=pc1max,ymin=pc2min,ymax=pc2max,col="red",fill = "transparent")
p2<-ggplot(pcaplot,aes(x=PC2,y=PC3,col=Population))+geom_point()+theme_bw()+theme(axis.title=element_text(size=12,face="bold"),axis.text=element_text(size=12,face="bold"),legend.position="right",legend.title=element_text(size=12,face="bold"),legend.text=element_text(size=12,face="bold"))+geom_rect(xmin=pc2min,xmax=pc2max,ymin=pc3min,ymax=pc3max,col="red",fill = "transparent")

grid.arrange(p1, p2,ncol=2, nrow = 1, widths = c(2.7, 3), heights = 2.7)

p3<-arrangeGrob(p1, p2,ncol=2, nrow = 1, widths = c(2.7, 3), heights = 2.7)

ggsave("../figures/PC123.png",dpi=300,width=13,height=6,p3)


#### reclasify

pca$inferred_population[which(is.na(pca$inferred_population) & pca$PC1<=pc1max & pca$PC1>=pc1min & pca$PC2<=pc2max & pca$PC2>=pc2min& pca$PC3<=pc3max & pca$PC3>=pc3min)]<-"EUR"


# EXPAND LIMITS FOR SAS SUBSET OF SAMPLES

pc1max_sas<-max(pca$PC1[which(pca$super_pop=="SAS")])+((max(pca$PC1[which(pca$super_pop=="SAS")])-min(pca$PC1[which(pca$super_pop=="SAS")]))*0.3)
pc1min_sas<-min(pca$PC1[which(pca$super_pop=="SAS")])-((max(pca$PC1[which(pca$super_pop=="SAS")])-min(pca$PC1[which(pca$super_pop=="SAS")]))*0.3)

pc2max_sas<-max(pca$PC2[which(pca$super_pop=="SAS")])+((max(pca$PC2[which(pca$super_pop=="SAS")])-min(pca$PC2[which(pca$super_pop=="SAS")]))*0.3)
pc2min_sas<-min(pca$PC2[which(pca$super_pop=="SAS")])-((max(pca$PC2[which(pca$super_pop=="SAS")])-min(pca$PC2[which(pca$super_pop=="SAS")]))*0.3)

pc3max_sas<-max(pca$PC3[which(pca$super_pop=="SAS")])+((max(pca$PC3[which(pca$super_pop=="SAS")])-min(pca$PC3[which(pca$super_pop=="SAS")]))*0.3)
pc3min_sas<-min(pca$PC3[which(pca$super_pop=="SAS")])-((max(pca$PC3[which(pca$super_pop=="SAS")])-min(pca$PC3[which(pca$super_pop=="SAS")]))*0.3)


pca$inferred_population[which(is.na(pca$inferred_population) & pca$PC1<=pc1max_sas & pca$PC1>=pc1min_sas & pca$PC2<=pc2max_sas & pca$PC2>=pc2min_sas
                              & pca$PC3<=pc3max_sas & pca$PC3>=pc3min_sas)]<-"SAS"

####################

# EXPAND LIMITS FOR EAS SUBSET OF SAMPLES

pc1max_eas<-max(pca$PC1[which(pca$super_pop=="EAS")])+((max(pca$PC1[which(pca$super_pop=="EAS")])-min(pca$PC1[which(pca$super_pop=="EAS")]))*0.3)
pc1min_eas<-min(pca$PC1[which(pca$super_pop=="EAS")])-((max(pca$PC1[which(pca$super_pop=="EAS")])-min(pca$PC1[which(pca$super_pop=="EAS")]))*0.3)

pc2max_eas<-max(pca$PC2[which(pca$super_pop=="EAS")])+((max(pca$PC2[which(pca$super_pop=="EAS")])-min(pca$PC2[which(pca$super_pop=="EAS")]))*0.3)
pc2min_eas<-min(pca$PC2[which(pca$super_pop=="EAS")])-((max(pca$PC2[which(pca$super_pop=="EAS")])-min(pca$PC2[which(pca$super_pop=="EAS")]))*0.3)

pc3max_eas<-max(pca$PC3[which(pca$super_pop=="EAS")])+((max(pca$PC3[which(pca$super_pop=="EAS")])-min(pca$PC3[which(pca$super_pop=="EAS")]))*0.3)
pc3min_eas<-min(pca$PC3[which(pca$super_pop=="EAS")])-((max(pca$PC3[which(pca$super_pop=="EAS")])-min(pca$PC3[which(pca$super_pop=="EAS")]))*0.3)



pca$inferred_population[which(is.na(pca$inferred_population) & pca$PC1<=pc1max_eas & pca$PC1>=pc1min_eas & pca$PC2<=pc2max_eas & pca$PC2>=pc2min_eas
                              & pca$PC3<=pc3max_eas & pca$PC3>=pc3min_eas)]<-"EAS"


####################

# EXPAND LIMITS FOR AFR SUBSET OF SAMPLES

pc1max_afr<-max(pca$PC1[which(pca$super_pop=="AFR")])+((max(pca$PC1[which(pca$super_pop=="AFR")])-min(pca$PC1[which(pca$super_pop=="AFR")]))*0.3)
pc1min_afr<-min(pca$PC1[which(pca$super_pop=="AFR")])-((max(pca$PC1[which(pca$super_pop=="AFR")])-min(pca$PC1[which(pca$super_pop=="AFR")]))*0.3)

pc2max_afr<-max(pca$PC2[which(pca$super_pop=="AFR")])+((max(pca$PC2[which(pca$super_pop=="AFR")])-min(pca$PC2[which(pca$super_pop=="AFR")]))*0.3)
pc2min_afr<-min(pca$PC2[which(pca$super_pop=="AFR")])-((max(pca$PC2[which(pca$super_pop=="AFR")])-min(pca$PC2[which(pca$super_pop=="AFR")]))*0.3)

pc3max_afr<-max(pca$PC3[which(pca$super_pop=="AFR")])+((max(pca$PC3[which(pca$super_pop=="AFR")])-min(pca$PC3[which(pca$super_pop=="AFR")]))*0.3)
pc3min_afr<-min(pca$PC3[which(pca$super_pop=="AFR")])-((max(pca$PC3[which(pca$super_pop=="AFR")])-min(pca$PC3[which(pca$super_pop=="AFR")]))*0.3)



pca$inferred_population[which(is.na(pca$inferred_population) & pca$PC1<=pc1max_afr & pca$PC1>=pc1min_afr & pca$PC2<=pc2max_afr & pca$PC2>=pc2min_afr
                              & pca$PC3<=pc3max_afr & pca$PC3>=pc3min_afr)]<-"AFR"

####################

# EXPAND LIMITS FOR AMR SUBSET OF SAMPLES

pc1max_amr<-max(pca$PC1[which(pca$super_pop=="AMR")])+((max(pca$PC1[which(pca$super_pop=="AMR")])-min(pca$PC1[which(pca$super_pop=="AMR")]))*0.3)
pc1min_amr<-min(pca$PC1[which(pca$super_pop=="AMR")])-((max(pca$PC1[which(pca$super_pop=="AMR")])-min(pca$PC1[which(pca$super_pop=="AMR")]))*0.3)

pc2max_amr<-max(pca$PC2[which(pca$super_pop=="AMR")])+((max(pca$PC2[which(pca$super_pop=="AMR")])-min(pca$PC2[which(pca$super_pop=="AMR")]))*0.3)
pc2min_amr<-min(pca$PC2[which(pca$super_pop=="AMR")])-((max(pca$PC2[which(pca$super_pop=="AMR")])-min(pca$PC2[which(pca$super_pop=="AMR")]))*0.3)

pc3max_amr<-max(pca$PC3[which(pca$super_pop=="AMR")])+((max(pca$PC3[which(pca$super_pop=="AMR")])-min(pca$PC3[which(pca$super_pop=="AMR")]))*0.3)
pc3min_amr<-min(pca$PC3[which(pca$super_pop=="AMR")])-((max(pca$PC3[which(pca$super_pop=="AMR")])-min(pca$PC3[which(pca$super_pop=="AMR")]))*0.3)

pca$inferred_population[which(is.na(pca$inferred_population) & pca$PC1<=pc1max_amr & pca$PC1>=pc1min_amr & pca$PC2<=pc2max_amr & pca$PC2>=pc2min_amr
                              & pca$PC3<=pc3max_amr & pca$PC3>=pc3min_amr)]<-"AMR"


pca$inferred_population<-as.factor(pca$inferred_population)
pca$inferred_population<-factor(pca$inferred_population, levels=c("EUR","AFR","EAS","SAS","AMR"))


#####################
## some inferred AFR in AMR, reclassify:
#
#table(pca[which(pca$PC3< -0.05),"super_pop"])
#table(pca[which(pca$PC3< -0.05),"inferred_population"])
#
#pca$inferred_population[which((pca$inferred_population=="AFR") & (pca$PC3< -0.05) )]<-"AMR"
#table(pca$super_pop,pca$inferred_population,useNA="ifany")
#
#
#table(pca[which(pca$PC1< -0.05),"super_pop"])
#table(pca[which(pca$PC1> -0.05),"super_pop"])
#
#table(pca[which(pca$PC1< -0.05),"inferred_population"])
#
#pca$inferred_population[which((pca$inferred_population=="AMR") & (pca$PC1< -0.05) )]<-"AFR"
#table(pca$super_pop,pca$inferred_population,useNA="ifany")
#
#
#pca$inferred_population[which((pca$inferred_population=="AFR") & (pca$PC1> -0.05) )]<-"AMR"
#table(pca$super_pop,pca$inferred_population,useNA="ifany")
#

####################
fam<-read.table(paste0(filehead,".fam"), head=F)

pca_nodup<-pca[which(pca$IID %in% fam$V2),]

write.table(pca_nodup[which(pca_nodup$inferred_population=="EUR"),c("FID","IID")],"list_eur_ancestry_samples_noDuplicates",col.names=F,row.names=F,quote=F,sep="\t")

write.table(pca_nodup[which(pca_nodup$inferred_population=="AFR"),c("FID","IID")],"list_afr_ancestry_samples_noDuplicates",col.names=F,row.names=F,quote=F,sep="\t")

write.table(pca_nodup[which(pca_nodup$inferred_population=="EAS"),c("FID","IID")],"list_eas_ancestry_samples_noDuplicates",col.names=F,row.names=F,quote=F,sep="\t")

write.table(pca_nodup[which(pca_nodup$inferred_population=="SAS"),c("FID","IID")],"list_sas_ancestry_samples_noDuplicates",col.names=F,row.names=F,quote=F,sep="\t")

write.table(pca_nodup[which(pca_nodup$inferred_population=="AMR"),c("FID","IID")],"list_amr_ancestry_samples_noDuplicates",col.names=F,row.names=F,quote=F,sep="\t")


