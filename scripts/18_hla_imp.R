args<-commandArgs(T)
anc <- as.character(args[1])
hlagene <- as.character(args[2])
digitalNum<- as.double(args[3])

#BiocManager::install("HIBAG")
setwd("/lustre/scratch123/hgi/mdt2/projects/ibdgwas_bioresource/qz2/pre_imputation_b04b06b15b19b20")

library(HIBAG)

## load pre-fitted model
model.list <- get(load(paste0("../ref/HLA/AffyAxiomUKB-",anc,"-HLA",digitalNum,"-hg19.RData")))

geno <- hlaBED2Geno(bed.fn=paste0("HLA/",anc,"/hla_hg19.bed"), fam.fn=paste0("HLA/",anc,"/hla_hg19.fam"), bim.fn=paste0("HLA/",anc,"/hla_hg19.bim"))

#hla.ids <- c("A","B","C","DRB1","DQA1","DQB1","DPB1")

model <- hlaModelFromObj(model.list[[hlagene]])

# best-guess genotypes and all posterior probabilities
pred.guess <- hlaPredict(model, geno, type="response+prob")

saveRDS(pred.guess,file=paste0("HLA/imp/",anc,"/",digitalNum,"/",hlagene,".rds"))

write.csv(pred.guess$value, file=paste0("HLA/imp/",anc,"/",digitalNum,"/",hlagene,".csv"),quote=F,row.names=F)
