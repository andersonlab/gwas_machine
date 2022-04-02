args<- commandArgs(T)
wkdir <- as.character(args[1])
filehead <- as.character(args[2])
outfilename1<- as.character(args[3])
outfilename2<- as.character(args[3])
setwd(wkdir)

ca_var_miss<-read.table(paste0(filehead,"_cases.lmiss"),head=T)
ctr_var_miss<-read.table(paste0(filehead,"_ctr.lmiss"), head=T)

ca_var_miss$N_NO_MISS<-ca_var_miss$N_GENO-ca_var_miss$N_MISS
ctr_var_miss$N_NO_MISS<-ctr_var_miss$N_GENO-ctr_var_miss$N_MISS

colnames(ca_var_miss)[c(3,6)]<-paste(colnames(ca_var_miss)[c(3,6)],"_cases",sep="")
colnames(ctr_var_miss)[c(3,6)]<-paste(colnames(ctr_var_miss)[c(3,6)],"_ctr",sep="")

all<-merge(ca_var_miss[,c(2,3,6)],ctr_var_miss[,c(2,3,6)],by="SNP",sort=F)


foo <- function(y){
  # include here as.numeric to be sure that your values are numeric:
  table <-  matrix(as.numeric(c(y[2], y[3], y[4], y[5])), ncol = 2, byrow = TRUE)
  if(any(is.na(table))) p <- "error" else p <- fisher.test(table, alternative="two.sided")$p.value
  p
} 
all$fisher_pvalue <- apply(all, 1, foo)

summary(all$fisher_pvalue)

dim(all[which(all$fisher_pvalue<1E-4),])

pdf("../../figures/hist_missingness_pvalue_ca_ctr.pdf",sep="",height = 5,width = 10)
ggplot(all, aes(x=(-log10(fisher_pvalue)))) + geom_histogram(binwidth=1) + geom_vline(xintercept=(-log10(1E-4)), linetype="dashed",color = "red", size=1)
dev.off()

write.table(all,"comparison_missingness_cases_ctr.lmiss", col.names=T,row.names=F,quote=F,sep="\t")
write.table(all[which(all$fisher_pvalue<1E-4),1,drop=F],outfilename2,col.names=T,row.names=F,quote=F,sep="\t")

