args<- commandArgs(T)
wkdir <- as.character(args[1])
filehead <- as.character(args[2])
gwassmyfile <- as.character(args[3])

setwd(wkdir)

library(qqman)
library(ggplot2)
library(data.table)

read.table(paste0("eur/",filehead,"_eur_CD.hwe"),header=T,stringsAsFactors=F)->hweCD
read.table(paste0("eur/",filehead,"_eur_UC.hwe"),header=T,stringsAsFactors=F)->hweUC
read.table(paste0("eur/",filehead,"_eur.hwe"),header=T,stringsAsFactors=F)->hwe

fread(gwassmyfile,sep="\t",stringsAsFactor=F)->smyIBD
smyIBD <- as.data.frame(smyIBD)

smyIBD$id1<- paste0("chr",smyIBD$hm_chrom,":",smyIBD$hm_pos,"_",smyIBD$hm_other_allele,"_",smyIBD$hm_effect_allele)
smyIBD$id2<- paste0("chr",smyIBD$hm_chrom,":",smyIBD$hm_pos,"_",smyIBD$hm_effect_allele,"_",smyIBD$hm_other_allele)

intername1<- intersect(smyIBD$id1, hwe$SNP)
intername2<- intersect(smyIBD$id2, hwe$SNP)

smyIBD1<- smyIBD[smyIBD$id1 %in% intername1,]
smyIBD2<- smyIBD[smyIBD$id2 %in% intername2,]

smyIBD1$id<- smyIBD1$id1
smyIBD2$id<- smyIBD2$id2
smyIBD.sel <- rbind(smyIBD1,smyIBD2)

hwe.sel <- hwe[ hwe$SNP %in%c(intername1,intername2),]

rownames(hwe.sel)<- hwe.sel$SNP
hwe.sel<- hwe.sel[smyIBD.sel$id,]
hwe.sel$hm_pos <- smyIBD.sel$hm_pos
hwe.sel$P1<- hwe.sel$P
hwe.sel[ hwe.sel$P1<1e-300,"P1"]<- 1e-300


png("../../figures/hwe_gwas5e-8.png", units="in", width=16, height=6, res=300)

manhattan(hwe.sel, chr = "CHR", bp = "hm_pos", p = "P1", snp = "SNP",ylim=c(0,300),cex.lab=1.5,font.lab=2,font.axis=2,cex.axis=1.6,las=2,font=2,ylab="-log10(P-value) from HWE test",genomewideline=12,suggestiveline=F,highlight = hwe.sel[smyIBD.sel$p_value<5e-8,"SNP"])

dev.off()

fread("~/proj1/data_06_01_2021/qc_phenotype/hetero/CD_UC_hetero.txt",header=T,stringsAsFactor=F)->hetero
hetero <- as.data.frame(hetero)
hetero$id1<- paste0("chr",hetero$CHR,":",hetero$BP,"_",hetero$A2,"_",hetero$A1)
hetero$id2<- paste0("chr",hetero$CHR,":",hetero$BP,"_",hetero$A1,"_",hetero$A2)

intername1<- intersect(hetero$id1, hwe$SNP)
intername2<- intersect(hetero$id2, hwe$SNP)

hetero1<- hetero[hetero$id1 %in% intername1,]
hetero2<- hetero[hetero$id2 %in% intername2,]

hetero1$id<- hetero1$id1
hetero2$id<- hetero2$id2
hetero.sel <- rbind(hetero1,hetero2)

hwe.sel <- hwe[ hwe$SNP %in%c(intername1,intername2),]

rownames(hwe.sel)<- hwe.sel$SNP
hwe.sel<- hwe.sel[hetero.sel$id,]
hwe.sel$hm_pos <- hetero.sel$BP
hwe.sel$P1<- hwe.sel$P
hwe.sel[ hwe.sel$P1<1e-300,"P1"]<- 1e-300

png("../../figures/hwe_hetero5e-8.png", units="in", width=16, height=6, res=300)

manhattan(hwe.sel, chr = "CHR",bp = "hm_pos", p = "P1", snp = "SNP",ylim=c(0,300),cex.lab=1.5,font.lab=2,font.axis=2,cex.axis=1.6,las=2,font=2,ylab="-log10(P-value) from HWE test",genomewideline=12,suggestiveline=F,highlight = hetero.sel[hetero.sel$P<5e-8  & !is.na(hetero.sel$P),"id"])

dev.off()




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


png("../../figures/hwe_CD_UC.png", units="in", width=15, height=5, res=300)
par(mfrow = c(1,3))
plot(-log10(hweCD.sel$P1),-log10(hweUC.sel$P1),xlab="HWE P-value based on CD patients",ylab="HWE P-value based on UC patients")
abline(v=-log10(0.05),h=-log10(0.05),col="red")
plot(-log10(hweCD.sel$P1),-log10(hweUC.sel$P1),xlim=c(0,50),ylim=c(0,50),xlab="HWE P-value based on CD patients",ylab="HWE P-value based on UC patients")
abline(v=-log10(0.05),h=-log10(0.05),col="red")
plot(-log10(hweCD.sel$P1[frq.sel$MAF>0.001]),-log10(hweUC.sel$P1[frq.sel$MAF > 0.001]),xlim=c(0,50),ylim=c(0,50),xlab="HWE P-value based on CD patients",ylab="HWE P-value based on UC patients")
abline(v=-log10(0.05),h=-log10(0.05),col="red")
dev.off()


