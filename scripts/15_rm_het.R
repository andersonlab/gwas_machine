args<-commandArgs(T)
wkdir <- as.character(args[1])
ancestry <- as.character(args[2])
filenodup <- as.character(args[3])

setwd(wkdir)

read.table(paste0("data/ancestry/",ancestry,"/het.het"),header=T,stringsAsFactors=F)->het

het$het<-(het$N.NM.-het$O.HOM)/het$N.NM.

nodup <- read.table(filenodup,header=F,stringsAsFactors=F)
nodup$FID_IID<- paste0(nodup$V1,"_",nodup$V2)
het$FID_IID <- paste0(het$FID,"_",het$IID)

het_nodup <- het[!(het$FID_IID %in% nodup$FID_IID),]

## get the threshold based on non-duplicated samples
limits.lower<-mean(het_nodup$het)-4*sd(het_nodup$het)
limits.upper<-mean(het_nodup$het)+4*sd(het_nodup$het) 
  

library(ggplot2)

ymax <- max(table(cut(het$het,seq(min(het$het),max(het$het),dist(range(het$het))/100))))

ggplot(data=het,aes(het))+geom_histogram(color="darkblue", fill="lightblue",binwidth = dist(range(het$het))/100)+ylab("Number of samples")+xlab("Heterozygosity rate")+theme_bw()+theme(axis.title=element_text(size=14,face="bold"),axis.text=element_text(size=14,face="bold"))+geom_vline(xintercept=limits.lower,lty=2,col="violetred")+geom_vline(xintercept=limits.upper,col="violetred",lty=2)+geom_vline(xintercept=mean(het$het),col="violetred",lty=2)+annotate(geom="text", x=limits.lower, y=ymax*0.95, label="-4SD",color="violetred",size=4,hjust=0)+annotate(geom="text", x=limits.upper, y=ymax*0.95, label="+4SD",color="violetred",size=4,hjust=0)+annotate(geom="text", x=mean(het$het), y=ymax*0.95, label="Average",color="violetred",size=4,hjust=0)

ggsave(paste0("figures/",ancestry,"_het.png"),dpi=300,width=8,height=6)

het.rm<- het[ (het$het > limits.upper)|(het$het < limits.lower),]

write.table(het.rm[,c("FID","IID")],paste0("data/ancestry/",ancestry,"/het_rm.list"),col.names=F,quote=F,row.names=F)
