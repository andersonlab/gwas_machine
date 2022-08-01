args<- commandArgs(T)
goodxy <- as.character(args[1])
goodx <- as.character(args[2])
filehead <- as.character(args[3])

library(ggplot2)
library(ggpubr)

sex<-read.table(goodxy,sep="",head=T)
#change FID to FID_IID since sometimes samples could have same FID
sex$FID <- paste0(sex$FID,"_",sex$IID)

sex$PEDSEX<-as.factor(sex$PEDSEX)


# X thresholds:

fmin_male<-mean(sex[which(sex$PEDSEX==1),"F"])-(4*sd(sex[which(sex$PEDSEX==1),"F"]))
fmin_male
fmax_female<-mean(sex[which(sex$PEDSEX==2),"F"])+(4*sd(sex[which(sex$PEDSEX==2),"F"]))
fmax_female

# # Y thresholds:
ymin_male<-mean(sex[which(sex$PEDSEX==1),"YCOUNT"])-(4*sd(sex[which(sex$PEDSEX==1),"YCOUNT"]))
ymin_male
ymax_female<-mean(sex[which(sex$PEDSEX==2),"YCOUNT"])+(4*sd(sex[which(sex$PEDSEX==2),"YCOUNT"]))
ymax_female

p1n<-ggplot(sex[which(sex$PEDSEX==1),],aes(x=F,y=YCOUNT)) + geom_point(colour = "#00AFBB") +
	theme_bw() +
	xlim(min(sex$F),max(sex$F)) +
	ylim(0,max(sex$YCOUNT)) +
	geom_hline(yintercept=ymin_male, linetype="dashed", color = "#00AFBB")+
	geom_vline(xintercept=fmin_male, linetype="dashed", color = "#00AFBB") +
	theme(axis.text= element_text(size=14,face="bold"),axis.title=element_text(size=14,face="bold"))

p2n<-ggplot(sex[which(sex$PEDSEX==2),],aes(x=F,y=YCOUNT)) + geom_point(colour = "#E7B800") +
	theme_bw() +
	xlim(min(sex$F),max(sex$F)) +
	ylim(0,max(sex$YCOUNT)) +
	geom_hline(yintercept=ymax_female, linetype="dashed", color = "#E7B800")+
	geom_vline(xintercept=fmax_female, linetype="dashed", color = "#E7B800") +
	theme(axis.text= element_text(size=14,face="bold"), axis.title=element_text(size=14,face="bold"))

# + geom_vline(xintercept=fmin_female, linetype="dashed", color = "#E7B800")

#p3n<-ggplot(sex[which(sex$PEDSEX==0),],aes(x=F,y=YCOUNT)) + geom_point(colour = "#FC4E07") +
#  xlim(min(sex$F),max(sex$F)) +
#  ylim(0,max(sex$YCOUNT)) +
#  geom_hline(yintercept=ymin_male, linetype="dashed", color = "#00AFBB")+
#  geom_vline(xintercept=fmin_male, linetype="dashed", color = "#00AFBB")+
#  geom_hline(yintercept=ymax_female, linetype="dashed", color = "#E7B800")+
#  geom_vline(xintercept=fmax_female, linetype="dashed", color = "#E7B800")
## + geom_vline(xintercept=fmin_female, linetype="dashed", color = "#E7B800")

#pn<-ggarrange(p1n,p2n,p3n,ncol=3,labels=c("Male","Female","NA"))
#dev.off()
#
#pdf("study_histogram_homozigosity_chrXY_with_newlimits.pdf",width =21)
#pn
#dev.off()

pn<-ggarrange(p1n,p2n,ncol=2,labels=c("Male","Female"))

#ggsave(file="./figures/F_ycount.png",width=12,height=6,dpi=300)

## PLOT DENSITY AS EXAMPLE
#sex$PEDSEX<-as.factor(sex$PEDSEX)
#
#p1<-ggscatterhist(sex[which(sex$PEDSEX==1),], x = "F", y = "YCOUNT",
#                  color = "PEDSEX", # comment out this and last line to remove the split by species
#                  # margin.plot = "histogram", # I'd suggest removing this line to get density plots
#                  # xlim = c(-2, 2),
#                  palette = c("#00AFBB"),
#                  margin.params = list(fill = "PEDSEX", color = "black", size = 0.2, xlim = c(-2, 2))
#)
#
#p2<-ggscatterhist(sex[which(sex$PEDSEX==2),], x = "F", y = "YCOUNT",
#                  color = "PEDSEX", # comment out this and last line to remove the split by species
#                  # margin.plot = "histogram", # I'd suggest removing this line to get density plots
#                  # xlim = c(-2, 2),
#                  palette = c("#E7B800"),
#                  margin.params = list(fill = "PEDSEX", color = "black", size = 0.2, xlim = c(-2, 2))
#)
#
##p3<-ggscatterhist(sex[which(sex$PEDSEX==0),], x = "F", y = "YCOUNT",
##                  color = "PEDSEX", # comment out this and last line to remove the split by species
##                  # margin.plot = "histogram", # I'd suggest removing this line to get density plots
##                  # xlim = c(-2, 2),
##                  palette = c("#FC4E07"),
##                  margin.params = list(fill = "PEDSEX", color = "black", size = 0.2, xlim = c(-2, 2))
##)
##
##p<-ggarrange(p1,p2,p3,ncol=3,labels=c("Male","Female","NA"))
#
#p<-ggarrange(p1,p2,ncol=2,labels=c("Male","Female"))
#
#pdf("study_histogram_homozigosity_chrXY.pdf",width =21)
#ggarrange(p1,p2,p3,ncol=3,labels=c("Male","Female","NA"))
#dev.off()

##############################

sex$SNPSEX2<-0
sex$SNPSEX2[which(sex$F<=fmax_female & sex$YCOUNT<=ymax_female)]<-2
sex$SNPSEX2[which(sex$F>=fmin_male & sex$YCOUNT>=ymin_male)]<-1

sex_chrx<-read.table(goodx,head=T)
#change FID to FID_IID since sometimes samples could have same FID
sex_chrx$FID <- paste0(sex_chrx$FID,"_",sex_chrx$IID)

colnames(sex_chrx)[4]<-"SNPSEX_chrx"
sexm<-merge(sex,sex_chrx[,c("FID","SNPSEX_chrx")],by="FID",sort=F)

table(sexm$PEDSEX,sexm$SNPSEX2)

samples_exclude<-sex[which((sex$PEDSEX==2 & sex$SNPSEX2==1) | (sex$PEDSEX==1 & sex$SNPSEX2==2)),]
table(samples_exclude$PEDSEX,samples_exclude$SNPSEX2)

dim(samples_exclude)

if(nrow(samples_exclude)>0){
	id <- unlist(strsplit(samples_exclude$FID,split="_"))
	samples_exclude$FID <- id[seq(1,length(id),2)]
}

write.table(samples_exclude[,c("FID","IID")],"./data/sex/list_samples_wrong_gender",col.names=T,row.names=F,quote=F,sep="\t")
write.table(samples_exclude[,c("FID","PEDSEX","SNPSEX2","F","YCOUNT")],"./data/sex/study_list_samples_sex_discrepancy",col.names=T,row.names=F,quote=F,sep="\t")

# recode the rest:

fam<-read.table(paste0(filehead,".fam"),sep="",head=F)
fam$FID_IID <- paste0(fam$V1,"_",fam$V2)
fam_ed<-merge(fam,sex[,c("FID","SNPSEX2")],by.x="FID_IID",by.y="FID",sort=F)

table(fam_ed$V2==fam$V2)

fam_ed<-fam_ed[,c("V1","V2","V3","V4","SNPSEX2","V6")]

write.table(fam_ed,paste0(filehead,"_edited.fam"),col.names=F,row.names=F,quote=F,sep="\t")

