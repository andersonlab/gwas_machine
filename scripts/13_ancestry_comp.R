library(ggplot2)
library(ggpubr)
library(gridExtra)

args<- commandArgs(T)
infer_anc<- as.character(args[1])
phenofile <- as.character(args[2])

read.table(infer_anc,header=T,stringsAsFactors=F)->anc_infer
read.table(phenofile,header=T,stringsAsFactors=F,sep="~")->pheno

anc_infer <- anc_infer[substr(anc_infer$FID,0,2)=="SP",]

anc <- merge(anc_infer,pheno[,c("SPid_1","ethnicity_group")],by.x="FID",by.y="SPid_1")
anc[anc$ethnicity_group=="","ethnicity_group"]<- "Not stated"

anc_report <- c("White or White British","Asian or British Asian","Black or Black British","Mixed","Other ethnic","Not stated")
anc_infer <- c("EUR","SAS","EAS","AFR","AMR")

cnt <- c()

for(i in anc_report){
	tmp <- c()
	for(j in anc_infer){
		num <- sum(anc$ethnicity_group==i & anc$inferred_population==j,na.rm=T)
		tmp <- c(tmp,num)
	}
	num <- sum(anc$ethnicity_group==i & is.na(anc$inferred_population))
	tmp <- c(tmp,num)
	cnt <- rbind(cnt,tmp)
}

cnt <- as.data.frame(cnt)
colnames(cnt)<- c(anc_infer,"NA")
rownames(cnt)<- anc_report


p1 <- ggplot(anc[anc$ethnicity_group=="White or White British",],aes(x=PC1,y=PC2,col=inferred_population))+geom_point()+theme_bw()+theme(axis.title=element_text(size=12,face="bold"),strip.text=element_text(size=12,face="bold"), axis.text=element_text(size=12,face="bold"),legend.position="none")

p2 <- ggplot(anc[anc$ethnicity_group=="White or White British",],aes(x=PC2,y=PC3,col=inferred_population))+geom_point()+theme_bw()+theme(axis.title=element_text(size=12,face="bold"),strip.text=element_text(size=12,face="bold"), axis.text=element_text(size=12,face="bold"),legend.position="right",legend.title=element_text(size=12,face="bold"),legend.text=element_text(size=12,face="bold"))

grid.arrange(p1, p2,ncol=2, nrow = 1, widths = c(2.7, 3), heights = 2.7)

p3<-arrangeGrob(p1, p2,ncol=2, nrow = 1, widths = c(2.7, 3), heights = 2.7)


p1 <- ggplot(anc,aes(x=PC1,y=PC2,col=inferred_population))+geom_point()+theme_bw()+theme(axis.title=element_text(size=12,face="bold"),strip.text=element_text(size=12,face="bold"), axis.text=element_text(size=12,face="bold"),legend.position="bottom",legend.title=element_text(size=12,face="bold"),legend.text=element_text(size=12,face="bold"))+facet_wrap(~ethnicity_group)




ggsave("../figures/ancestry_infer_vs_report.png",dpi=300,width=13,height=6,p3)



