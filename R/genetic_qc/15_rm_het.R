args<-commandArgs(T)
wkdir <- as.character(args[1])
ancestry <- as.character(args[2])

setwd(paste0(wkdir,"/data/ancestry/",ancestry))

read.table("het.het",header=T,stringsAsFactors=F)->het

het$het<-(het$N.NM.-het$O.HOM)/het$N.NM.

#### CREATE AN HISTOGRAM

# number of bins in histogram
fd=function(x) {
  n=length(x)
  r=IQR(x)
  2*r/n^(1/3)
}

limits.lower<-mean(het$het)-4*sd(het$het)
limits.upper<-mean(het$het)+4*sd(het$het) 
  

library(ggplot2)

ymax <- max(table(cut(het$het,seq(min(het$het),max(het$het),dist(range(het$het))/100))))

ggplot(data=het,aes(het))+geom_histogram(color="darkblue", fill="lightblue",binwidth = dist(range(het$het))/100)+ylab("Number of samples")+xlab("Heterozygosity rate")+theme_bw()+theme(axis.title=element_text(size=14,face="bold"),axis.text=element_text(size=14,face="bold"))+geom_vline(xintercept=limits.lower,lty=2,col="violetred")+geom_vline(xintercept=limits.upper,col="violetred",lty=2)+geom_vline(xintercept=mean(het$het),col="violetred",lty=2)+annotate(geom="text", x=limits.lower, y=ymax*0.95, label="-4SD",color="violetred",size=4,hjust=0)+annotate(geom="text", x=limits.upper, y=ymax*0.95, label="+4SD",color="violetred",size=4,hjust=0)+annotate(geom="text", x=mean(het$het), y=ymax*0.95, label="Average",color="violetred",size=4,hjust=0)

ggsave(paste0(wkdir,"/figures/",ancestry,"_het.png"),dpi=300,width=8,height=6)

het.rm<- het[ (het$het > limits.upper)|(het$het < limits.lower),]

write.table(het.rm[,c("FID","IID")],"het_rm.list",col.names=F,quote=F,row.names=F)
