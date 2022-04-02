args<- commandArgs(T)
wkdir <- as.character(args[1])
filehead <- as.character(args[2])
outputfile <- as.character(args[3])

setwd(wkdir)
read.table(paste0("eur/",filehead,"_eur.hwe"),header=T,stringsAsFactors=F)->hwe

pos<- unlist(strsplit(hwe$SNP,split=":|_"))
pos<- pos[seq(2,length(pos),4)]
hwe$pos<- as.double(pos)
hwe_auto<- hwe[!(hwe$CHR==6 & hwe$pos>=28510120 & hwe$pos<=33480577) & hwe$CHR!=23,]
hwe_auto.sel<- hwe_auto[ hwe_auto$P<1e-12,]
hwe_auto.snp<- hwe_auto.sel$SNP

hwe_sex.snp <- c()
filename <- paste0("eur/",filehead,"_eur_females.hwe")
if(file.exists(filename)){
	read.table(filename,header=T,stringsAsFactors=F)->hwe
	pos<- unlist(strsplit(hwe$SNP,split=":|_"))
	pos<- pos[seq(2,length(pos),4)]
	hwe$pos<- as.double(pos)
	hwe_sex.sel<- hwe_sex[ hwe_sex$P<1e-12,]
	hwe_sex.snp<- hwe_sex.sel$SNP
}

hwe.snp <- c(hwe_auto.snp,hwe_sex.snp)
write.table(hwe.snp,file=outputfile,col.names=F,quote=F,row.names=F)
