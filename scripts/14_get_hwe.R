args<- commandArgs(T)
wkdir <- as.character(args[1])
filehead <- as.character(args[2])
outputfile <- as.character(args[3])

setwd(wkdir)
read.table(paste0("eur/",filehead,"_eur.hwe"),header=T,stringsAsFactors=F)->hwe
read.table(paste0("eur/",filehead,"_eur.frq"),header=T,stringsAsFactors=F)->frq

read.table(paste0("eur/",filehead,"_eur_CD.hwe"),header=T,stringsAsFactors=F)->hweCD
read.table(paste0("eur/",filehead,"_eur_UC.hwe"),header=T,stringsAsFactors=F)->hweUC


pos<- unlist(strsplit(hwe$SNP,split=":|_"))
pos<- pos[seq(2,length(pos),4)]
hwe$pos<- as.double(pos)
hwe_auto<- hwe[!(hwe$CHR==6 & hwe$pos>=28510120 & hwe$pos<=33480577) & hwe$CHR!=23,]
hwe_auto.sel<- hwe_auto[ hwe_auto$P<1e-12,]
hwe_auto.snp<- hwe_auto.sel$SNP

hweCD.sel <- hweCD[ hweCD$SNP %in% hwe_auto.snp,]
hweUC.sel <- hweUC[ hweUC$SNP %in% hwe_auto.snp,]
frq.sel <- frq[ frq$SNP %in% hwe_auto.snp,]

hweCD.sel$P1 <- hweCD.sel$P
hweCD.sel[ hweCD.sel$P < 1e-300,"P1"]<- 1e-300

hweUC.sel$P1 <- hweUC.sel$P
hweUC.sel[ hweUC.sel$P < 1e-300,"P1"]<- 1e-300


#for SNPs to be removed, check if they also have small HWE P-values (i.e. < 0.05) in CD and UC seperately, if no, that means this SNP could be CD or UC specific, try to not remove them since their small HWE P-value could be caused by MAF difference between UC and CD patients. I chose 0.05 since I want to be conservative (so not too many SNPs would be excluded from the removing list). I also restrict the SNPs for those with MAF > 0.001 since some very rare variants could have high HWE P-value (e.g. 0) in CD/UC, I do not want to exclude these SNPs from removing list.

hwe_auto.snp_torm <- hwe_auto.snp[ !((hweCD.sel$P >= 0.05 | hweUC.sel$P >= 0.05) & frq.sel$MAF > 0.001)]

hwe_sex.snp <- c()
filename <- paste0("eur/",filehead,"_eur_females.hwe")
if(file.exists(filename)){
	read.table(filename,header=T,stringsAsFactors=F)->hwe_sex
	frq <- read.table(paste0("eur/",filehead,"_eur_females.frq"),header=T,stringsAsFactor=F)
	pos<- unlist(strsplit(hwe_sex$SNP,split=":|_"))
	pos<- pos[seq(2,length(pos),4)]
	hwe_sex$pos<- as.double(pos)
	hwe_sex.sel<- hwe_sex[ hwe_sex$P<1e-12,]
	hwe_sex.snp<- hwe_sex.sel$SNP

	read.table(paste0("eur/",filehead,"_eur_females_CD.hwe"),header=T,stringsAsFactors=F)->hwe_sexCD
	read.table(paste0("eur/",filehead,"_eur_females_UC.hwe"),header=T,stringsAsFactors=F)->hwe_sexUC
	
	hwe_sexCD.sel <- hwe_sexCD[ hwe_sexCD$SNP %in% hwe_sex.snp,]
	hwe_sexUC.sel <- hwe_sexUC[ hwe_sexUC$SNP %in% hwe_sex.snp,]
	frq.sel <- frq[ frq$SNP %in% hwe_sex.snp,]
	hwe_sex.snp_torm <- hwe_sex.snp[ !((hwe_sexCD.sel$P >= 0.05 | hwe_sexUC.sel$P >= 0.05) & frq.sel$MAF > 0.001)]

}

hwe.snp <- c(hwe_auto.snp_torm,hwe_sex.snp_torm)
write.table(hwe.snp,file=outputfile,col.names=F,quote=F,row.names=F)


