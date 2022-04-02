args<- commandArgs(T)
wkdir <- as.character(args[1])
filehead <- as.character(args[2])
gwassmyfile <- as.character(args[3])

setwd(wkdir)

library(qqman)

read.table(paste0("eur/",filehead,"_eur.hwe"),header=T,stringsAsFactors=F)->hwe

read.csv(gwassmyfile,sep="\t",stringsAsFactor=F)->smy
smy$id1<- paste0("chr",smy$hm_chrom,":",smy$hm_pos,"_",smy$hm_other_allele,"_",smy$hm_effect_allele)
smy$id2<- paste0("chr",smy$hm_chrom,":",smy$hm_pos,"_",smy$hm_effect_allele,"_",smy$hm_other_allele)

intername1<- intersect(smy$id1, hwe$SNP)
intername2<- intersect(smy$id2, hwe$SNP)

smy1<- smy[smy$id1 %in% intername1,]
smy2<- smy[smy$id2 %in% intername2,]

smy1$id<- smy1$id1
smy2$id<- smy2$id2
smy.sel <- rbind(smy1,smy2)

hwe.sel <- hwe[ hwe$SNP %in%c(intername1,intername2),]

rownames(hwe.sel)<- hwe.sel$SNP
hwe.sel<- hwe.sel[smy.sel$id,]

hwe.sel$hm_pos <- smy.sel$hm_pos

hwe.sel$P1<- hwe.sel$P
hwe.sel[ hwe.sel$P1<1e-150,"P1"]<- 1e-150

smy.sel$p_value1<- smy.sel$p_value
smy.sel[ smy.sel$p_value1<1e-150,"p_value1"]<- 1e-150


png("../../figures/hwe_gwas5e-8.png", units="in", width=16, height=12, res=300)

manhattan(hwe.sel, chr = "CHR", bp = "hm_pos", p = "P1", snp = "SNP",ylim=c(2,150),cex.lab=1.5,font.lab=2,font.axis=2,cex.axis=1.6,las=2,font=2,ylab="-log10(P-value) from HWE test",genomewideline=12,suggestiveline=F,highlight = hwe.sel[smy.sel$p_value<5e-8,"SNP"])

dev.off()

#png("../figures/hwe.png", units="in", width=16, height=12, res=300)
#
#par(mfrow=c(2,1))
#par(mar=c(0,5,3,3))
#manhattan(hwe.sel, chr = "CHR", bp = "hm_pos", p = "P1", snp = "SNP",ylim=c(2,100),cex=2.2,cex.lab=1.5,font.lab=2,font.axis=2,cex.axis=1.6,las=2,font=2,ylab="-log10(P-value) from HWE test")
#par(mar=c(5,5,3,3))
#manhattan(smy.sel,chr = "hm_chrom", bp = "hm_pos", p = "p_value1", snp = "hm_rsid",ylim=c(100,2),cex=2.2,cex.lab=1.5,font.lab=2,font.axis=2,cex.axis=1.6,las=2,font=2,xlab="",xaxt="n",ylab="-log10(P-value) from GWAS")
#dev.off()

# get MHC region from https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38.p13, chr6:28,510,120-33,480,577
ancestry<-"eur"
tmp <- read.table(paste0(i,"/",filehead,"_",i,".hwe"),header=T,stringsAsFactors=F)
pos<- unlist(strsplit(tmp$SNP,split=":|_"))
pos<- pos[seq(2,length(pos),4)]
tmp$pos<- as.double(pos)
tmp<- tmp[!(tmp$CHR==6 & tmp$pos>=28510120 & tmp$pos<=33480577),]
tmp.sel<- tmp[ tmp$P<1e-12,]
hwe <- rbind(hwe,tmp.sel)

hwe.snp<- tmp.sel$SNP

write.table(hwe.snp,file="list_var_exclude_hwe",col.names=F,quote=F,row.names=F)
