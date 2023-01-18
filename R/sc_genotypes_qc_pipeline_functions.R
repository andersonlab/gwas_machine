## Scripts for Tobi QC ##
library(data.table)
library(stringr)

#Add in the SNPs to keep
system("rm snps_to_keep")
#same SNP, just in different formats
system('echo "chr16:50729870_C_CC" >> snps_to_keep')
system('echo "chr16:50729870_CC_C" >> snps_to_keep')
system('echo "AX-96079897" >> snps_to_keep')

merge_batches = function(path_to_plate_1,path_to_plate_2){
  sprintf("/software/team152/bcftools-1.9/bcftools merge -m id %s %s --force-samples -Oz -o merged_plates.vcf.gz", path_to_plate_1, path_to_plate_2)
  sprintf("/software/team152/plink_linux_x86_64_20181202/plink --vcf merged_plates.vcf.gz --allow-extra-chr --output-chr MT --make-bed --out merged_plates")
}

# Rename VCF IDs for Tobi's manifest

# update_sample_names = function(path_to_plate_1,path_to_plate_2){
#   #Get the sample IDs for the two plates
#   system("/software/team152/bcftools-1.9/bcftools query -l /lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/genotypes/UKX16_Results/Sanger_Plate_2_Oct/VCF/plate_2_oct_joint_vcf/plate_2_oct_joint_vcf > /lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/genotypes/UKX16_Results/Sanger_Plate_2_Oct/VCF/plate_2_oct_joint_vcf/plate_2_sample_IDs")
#   system("/software/team152/bcftools-1.9/bcftools query -l /lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/genotypes/UKX16_Results/Sanger_Plate_1_Oct/VCF/plate_1_oct_joint_vcf/plate_1_oct_joint.vcf.gz >  /lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/genotypes/UKX16_Results/Sanger_Plate_1_Oct/VCF/plate_1_oct_joint_vcf/plate_1_sample_IDs")
#
#   #Read in the sample IDs:
#   sample_plate_1_ID = fread("/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/genotypes/UKX16_Results/Sanger_Plate_1_Oct/VCF/plate_1_oct_joint_vcf/plate_1_sample_IDs", header = F)
#   sample_plate_2_ID = fread("/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/genotypes/UKX16_Results/Sanger_Plate_2_Oct/VCF/plate_2_oct_joint_vcf/plate_2_sample_IDs", header = F)
#
#   #Read in the sample manifest (talk to Tobi to see how this gets auto-updated)
#   sample_manifest = fread("/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/genotypes/UKX16_Results/meta_data_nov_21_2022.csv")
#
#   #Only consider unique genotype IDs
#   sample_manifest =  sample_manifest[!duplicated(sample_manifest[,c("Genotyping_ID")]),]
#
#   #Add new columns to this chart, 'genotyping_sample_id_modified' that will be used to identify this
#
#   #First, only consider everything that is a number (which will be converted wrongly)
#   sample_manifest$genotyping_sample_id_modified = as.numeric(sample_manifest$Genotyping_ID)
#
#   #Convert the manifest IDs to the mangled scientific notation ones (e12)
#   #Note these will be characters!
#   sample_manifest$genotyping_sample_id_modified = formatC(as.numeric(sample_manifest$genotyping_sample_id_modified ), format = "e", digits = 10)
#
#   #Match on the numeric characters, since manifest has different notation and weird +'s as escape characters?
#   sample_manifest$genotyping_sample_id_modified  = str_split(sample_manifest$genotyping_sample_id_modified, "e", simplify=T)[,1]
#
#   #Split the string on the *plate ID* to get the same string lengths
#   sample_plate_1_ID_scientific_match = str_split(sample_plate_1_ID$V1, "e", simplify=T)[,1]
#
#   substr_1 = substring(sample_plate_1_ID_scientific_match,1, nchar(sample_plate_1_ID_scientific_match)-2)
#   substr_2 = substring(sample_manifest$genotyping_sample_id_modified,1,nchar(sample_manifest$genotyping_sample_id_modified))
#
#   #Perform the matching between the IDs on the plate and the manifest, merge as an additonal column!
#   sample_plate_1_rename     = matrix(data = NA, nrow = 95, ncol = 2)
#   sample_plate_1_rename[,1] = sample_plate_1_ID$V1
#   #Loop through the manifest:
#   for (i in 1:length(sample_plate_1_ID_scientific_match)) {
#     #Get the index
#     index = which(sample_manifest$genotyping_sample_id_modified %in% sample_plate_1_ID_scientific_match[i])
#     if (length(index) > 0) {
#       sample_plate_1_rename[i,2] = sample_manifest$Genotyping_ID[index]
#     }
#     else{
#       sample_plate_1_rename[i,2] = NA
#     }
#   }
#
#   #Manually refill in the missing genotype IDs
#   #3.9815768467e+012.CEL = 3981576846650
#   #3.9815768627e+012.CEL = 3981576862650
#   #3.9815768637e+012.CEL = 3981576860861X
#   #3.9815768647e+012.CEL = 3981576861875Y
#
#   sample_plate_1_rename[73,2] = "3981576846650"
#   sample_plate_1_rename[89,2] = "3981576862650"
#   sample_plate_1_rename[90,2] = "3981576863664"
#   sample_plate_1_rename[91,2] = "3981576864678"
#
#   #Write out this correspondence:
#   write.table(sample_plate_1_rename, "/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/genotypes/UKX16_Results/Sanger_Plate_1_Oct/VCF/ID_change",quote = F, col.names = F, row.names = F)
#
#   #Recode the VCF!
#   bcftools reheader -s names.txt in.vcf > out.vcf
#
#   ### Go on to using these indexes to rename the VCF on plate 1
# }

remove_non_ATCG = function(merged_plates)
{

  #Read in .bim file to get list of SNPs
  bim<-read.table("merged_plates.bim")

  #Get list of indels
  indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]

  #Get list of SNPs on non 1-22 & sex chr's (usually mitochondrial)
  otherchr <- bim[ !(bim$V1 %in% c(1:22,"X","Y")),]

  #Get joint list of SNPs to remove
  all_remove<-rbind(indels,otherchr)
  #Remove duplicate SNPs in list
  all_remove<-all_remove[!duplicated(all_remove),]

  #Generate list of SNPs to remove
  write.table(all_remove[,"V2"],"list_indel_var_exclude", col.names=F,row.names=F,quote=F,sep="\t")

  #Remove SNPs in this list from merged file:
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile merged_plates --allow-extra-chr --output-chr MT --exclude list_indel_var_exclude --make-bed --out merged_plates_ATCG")
}

align_positive_strand(){
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile merged_plates_ATCG --allow-extra-chr --output-chr MT --allow-no-sex --recode vcf --out merged_plates_ATCG")
  #Note this alignment step takes a *long* time ~ 1hr, so make yourself a coffee while it runs and don't panic if you don't see an output
  system("/software/team152/bcftools-1.9/bcftools +/software/team152/bcftools-1.9/plugins/fixref.so  merged_plates_ATCG.vcf -Oz -o merged_plates_ATCG_aligned.vcf.gz -- -f /lustre/scratch123/hgi/projects/ibdgwas/IIBDGC/resources/hg38/hg38_edited.fa -m top 2>&1 | tee alignment.log")
  #After alignment, format vcf to bed file
  system("/software/team152/plink_linux_x86_64_20181202/plink --vcf merged_plates_ATCG_aligned.vcf.gz --keep-allele-order --id-delim --make-bed --out merged_plates_ATCG_aligned")
}

update_variant_ids = function(){
  #Generate missingness across SNPs and samples
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile merged_plates_ATCG_aligned --missing --out merged_plates_ATCG_aligned")
  # Update the IDs to chr:pos:ref:alt
  system("zcat merged_plates_ATCG_aligned.vcf.gz|grep -v ^"#"|cut -f '1-5' | awk '{print $3,"chr"$1":"$2"_"$4"_"$5}' > merged_plates_ATCG_aligned")
}

remove_duplicated_variants = function(){

  #Read in bim table for SNPs
  bim<-read.table("merged_plates_ATCG_aligned.bim", sep="\t",head=F)
  #Read in IDs
  ids<-read.table("merged_plates_ATCG_aligned",header=F,stringsAsFactor=F)

  #Read in variant missingness
  varmiss<-read.table("merged_plates_ATCG_aligned.lmiss",sep="",head=T)

  #Rename column
  colnames(ids)[2]<-"ids"

  #The major allele is set to A2 by default by Plink, keep ids with real ref/alt as in vcf using the ids file
  bim.1<-cbind(bim,ids[,"ids",drop=F])

  bim.1<-merge(bim.1,varmiss[,c("SNP","F_MISS")],by.x="V2",by.y="SNP",sort=F)

  bim.1$ids<-as.character(bim.1$ids)

  # identify duplicated variants (same chr position ref and alt)
  dups<-bim.1[which(duplicated(bim.1$ids)),"ids"]
  dups <- unique(dups)

  bim.1$V2 <- as.character(bim.1$V2)
  bim.1$idx <- 1:nrow(bim.1)

  #Helper function to remove duplicated variants
  rm_dup <- function(bim,dupvar){
    tmp<-bim[bim$ids == dupvar,]
    keep<-tmp[tmp$F_MISS==min(tmp$F_MISS),]
    #if there are many duplicated variants with same missing rates, select the first one
    exclude<-tmp[tmp$V2!= keep[1,"V2"],"idx"]
    return(exclude)
  }

  bim.1.dup <- bim.1[bim.1$ids%in%dups,]
  excludes <- lapply(dups,function(x) rm_dup(bim.1.dup,x))
  excludes <- unlist(excludes)

  bim.1$ids[which(bim.1$idx %in% excludes)]<-paste(bim.1$ids[which(bim.1$idx %in% excludes)],"_rm",sep="")

  duplicated_variants<-bim.1[grep("_rm",bim.1$ids),"ids"]

  write.table(duplicated_variants,"list_duplicated_var_exclude",col.names=F,row.names=F,quote=F,sep="\t")
  write.table(bim.1[,c(2,7,3:6)],"merged_plates_ATCG_aligned_edited.bim",col.names=F,row.names=F,quote=F,sep="\t")

  #Remove duplicated variants with PLINK
  system("/software/team152/plink_linux_x86_64_20181202/plink --bed merged_plates_ATCG_aligned.bed --bim merged_plates_ATCG_aligned_edited.bim --fam merged_plates_ATCG_aligned.fam --exclude list_duplicated_var_exclude --keep-allele-order --make-bed --out merged_plates_ATCG_aligned_nodup")
}

#Compare Frequency
compare_freq_1KG = function(){
  #Note that all the reference data is located in: /lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/genotypes/qc_output/1KG_data

  # Compare MAF between 1KG EUR and study, and flip allele
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile  /lustre/scratch119/realdata/mdt2/projects/sc-eqtl-ibd/data/genotypes/qc_output/1KG_data/1000GP_EUR_b38_study_variants --freq --out /lustre/scratch119/realdata/mdt2/projects/sc-eqtl-ibd/data/genotypes/qc_output/1KG_data/1000GP_EUR_b38_study_variants")
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile merged_plates_ATCG_aligned_nodup --freq --out merged_plates_ATCG_aligned_nodup")

  #setwd("/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/genotypes/qc_output/")

  gp<-read.table("/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/genotypes/qc_output/1KG_data/1000GP_EUR_b38_study_variants.frq",head=T)
  g1<-read.table("merged_plates_ATCG_aligned_nodup.frq",head=T)

  colnames(gp)[3:6]<-paste(colnames(gp)[3:6],"_gp",sep="")
  colnames(g1)[3:6]<-paste(colnames(g1)[3:6],"_g1",sep="")

  all<-merge(g1,gp[,2:6],by="SNP",all=T)

  # variants with different minor allele were checked
  check<-all[which(all$A1_g1!=all$A1_gp),]

  # Keep only A/T C/G
  check<-check[which( (check$A1_g1=="G" & check$A2_g1=="C") | (check$A1_g1=="C" & check$A2_g1=="G") | (check$A1_g1=="A" & check$A2_g1=="T") | (check$A1_g1=="T" & check$A2_g1=="A")),]

  # LIST OF VARIANTS TO REMOVE, WE CANNOT REALLY BE SURE WHETHER THERE IS STRAND ISSUE OR NOT
  remove<-check[which(check$MAF_g1>=0.45),]

  # LIST OF VARIANTS TO FLIP:
  flip<-check[which(check$MAF_g1<0.45),]

  flip<-flip[order(flip$MAF_g1,decreasing=T),]

  remove_2<-flip[which(flip$MAF_g1>0.2 & flip$MAF_gp<0.1),]
  remove_3<-flip[which(flip$MAF_g1<0.1 & flip$MAF_gp>0.2),]

  remove<-rbind(remove,remove_2,remove_3)

  write.table(remove[,"SNP"],"list_variants_to_remove_AT_CG",col.names=F,row.names=F,quote=F,sep="\t")

  flip<-flip[which(!flip$SNP %in% remove$SNP),]

  write.table(flip[,"SNP"],"list_variants_to_flip_AT_CG",col.names=F,row.names=F,quote=F,sep="\t")

  system("grep -v -w -f snps_to_keep list_variants_to_remove_AT_CG > list_variants_to_remove_AT_CG_rmkeep")
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile merged_plates_ATCG_aligned_nodup --exclude list_variants_to_remove_AT_CG_rmkeep --flip list_variants_to_flip_AT_CG --make-bed --out merged_plates_ATCG_aligned_nodup_flip")
}

#Remove y chromosome
remove_y_chr = function(){
  system("awk '{if($1==24)print $2}' < merged_plates_ATCG_aligned_nodup_flip.bim > list_chry_variants_toremove")
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile merged_plates_ATCG_aligned_nodup_flip --allow-no-sex --exclude list_chry_variants_toremove --make-bed --out merged_plates_ATCG_aligned_nodup_flip_nochry")
}

# Remove samples and variants with low call rate
remove_low_call_rate = function(){

  # Sample Call Rate < 80%
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile merged_plates_ATCG_aligned_nodup_flip_nochry --allow-no-sex --mind 0.20 --make-bed --out merged_plates_ATCG_aligned_nodup_flip_nochry_scr0.8")

  # Variant Call Rate < 80%
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile merged_plates_ATCG_aligned_nodup_flip_nochry_scr0.8 --allow-no-sex --missing --out merged_plates_ATCG_aligned_nodup_flip_nochry_scr0.8_vcr0.8")
  system("awk '{if($5>0.2) print $2}' merged_plates_ATCG_aligned_nodup_flip_nochry_scr0.8_vcr0.8.lmiss > list_variants_to_remove_vcr0.8")
  system("grep -v -w -f snps_to_keep list_variants_to_remove_vcr0.8 > list_variants_to_remove_vcr0.8_rmkeep")
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile merged_plates_ATCG_aligned_nodup_flip_nochry_scr0.8  --allow-no-sex --exclude list_variants_to_remove_vcr0.8_rmkeep --make-bed --out merged_plates_ATCG_aligned_nodup_flip_nochry_scr0.8_vcr0.8")

  # Sample Call Rate < 95%
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile merged_plates_ATCG_aligned_nodup_flip_nochry_scr0.8_vcr0.8 --allow-no-sex --mind 0.05 --make-bed --out merged_plates_ATCG_aligned_nodup_flip_nochry_scr0.8_vcr0.8_scr0.95")

  # Variant Call Rate < 95%
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile merged_plates_ATCG_aligned_nodup_flip_nochry_scr0.8_vcr0.8_scr0.95 --allow-no-sex --missing --out merged_plates_ATCG_aligned_nodup_flip_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95")
  system("awk '{if($5>0.05) print $2}' merged_plates_ATCG_aligned_nodup_flip_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95.lmiss > list_variants_to_remove_vcr0.95")
  system("grep -v -w -f snps_to_keep list_variants_to_remove_vcr0.95 > list_variants_to_remove_vcr0.95_rmkeep")
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile merged_plates_ATCG_aligned_nodup_flip_nochry_scr0.8_vcr0.8_scr0.95 --allow-no-sex --exclude list_variants_to_remove_vcr0.95_rmkeep --make-bed --out merged_plates_ATCG_aligned_nodup_flip_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95")

  # Variant MAF<0.01 AND Call Rate <98%
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile merged_plates_ATCG_aligned_nodup_flip_nochry_scr0.8_vcr0.8_scr0.95 --allow-no-sex --missing --out merged_plates_ATCG_aligned_nodup_flip_nochry_scr0.8_vcr0.8_scr0.95")
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile merged_plates_ATCG_aligned_nodup_flip_nochry_scr0.8_vcr0.8_scr0.95 --allow-no-sex --freq --out merged_plates_ATCG_aligned_nodup_flip_nochry_scr0.8_vcr0.8_scr0.95")

  sample_miss<-read.table("merged_plates_ATCG_aligned_nodup_flip_nochry_scr0.8_vcr0.8_scr0.95.imiss",head=T)
  var_miss<-read.table("merged_plates_ATCG_aligned_nodup_flip_nochry_scr0.8_vcr0.8_scr0.95.lmiss",head=T)
  frq<-read.table("merged_plates_ATCG_aligned_nodup_flip_nochry_scr0.8_vcr0.8_scr0.95.frq",head=T)

  var<-merge(frq[,c(2:6)],var_miss,by="SNP")
  var.1<-var[which(var$MAF<0.01 & var$F_MISS>0.02),]
  write.table(var.1[,"SNP",drop=F],"list_monomorphic_vcr098maf0.01_var_exclude",col.names=F,row.names=F,quote=F,sep="\t")

  system("grep -v -w -f snps_to_keep list_monomorphic_vcr098maf0.01_var_exclude > list_monomorphic_vcr098maf0.01_var_exclude_rmkeep")
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile merged_plates_ATCG_aligned_nodup_flip_nochry_scr0.8_vcr0.8_scr0.95 --allow-no-sex --exclude list_monomorphic_vcr098maf0.01_var_exclude_rmkeep --make-bed --out merged_plates_ATCG_aligned_nodup_flip_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.01")
}

#Filter heterozygosity
filter_heterozygosity = function(){
# To check
}

#
filter_monomorphic_variants = function(){

  plink --bfile data/ancestry/eur/${filenamedup}_eur_hwe1e-12_het \
  --keep-allele-order --allow-no-sex \
  --freq \
  --out data/ancestry/eur/${filenamedup}_eur_hwe1e-12_het

  awk '{if($5==0)print $2}' < data/ancestry/eur/${filenamedup}_eur_hwe1e-12_het.frq > data/ancestry/eur/study_list_variants_toexclude_monomorphic

  plink --bfile data/ancestry/eur/${filenamedup}_eur_hwe1e-12_het \
  --allow-no-sex \
  --exclude data/ancestry/eur/study_list_variants_toexclude_monomorphic \
  --make-bed --out data/ancestry/eur/${filenamedup}_eur_hwe1e-12_het_nomonom

}





