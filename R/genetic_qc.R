#' @import stringr

#Download reference 1000G from Sanger central resources
prepare_reference_data = function()
{
  # Prepare reference data
  # #mkdir ref
  #1KG raw file can be downloaded here: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/
  #change to plink format
  #mkdir ref/1KG
  #mkdir ref/log
  ### I set the variant ID in format: chr:pos_ref_alt, only the first 20 characters were used for alleles with length longer than 20
  #for((i=1;i<23;i++))
  #do
  #	plink2 --vcf /lustre/scratch115/resources/1000g/release/20201028/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.filtered.shapeit2-duohmm-phased.vcf.gz --set-all-var-ids chr@:\#_\$r_\$a --new-id-max-allele-len 20 truncate --max-alleles 2 --rm-dup force-first --keep-allele-order --make-bed --out ref/1KG/chr${i}
  #done
  #
  ## X chromosome
  ### remove header line at first, otherwise there would be error like "Duplicate FORMAT:GT header line in --vcf file."
  #zcat /lustre/scratch115/resources/1000g/release/20201028/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.vcf.gz|grep -v "^##" > ref/1KG/chr23_tmp
  #plink2 --vcf ref/1KG/chr23_tmp --set-all-var-ids chr@:\#_\$r_\$a --new-id-max-allele-len 20 truncate --max-alleles 2 --keep-allele-order --rm-dup force-first --make-bed --out ref/1KG/chr23
  #rm ref/1KG/chr23_tmp

  ##get high LD region
  #wget https://raw.githubusercontent.com/meyer-lab-cshl/plinkQC/master/inst/extdata/high-LD-regions-hg38-GRCh38.txt -O ref/high-LD-regions-hg38-GRCh38.txt
  #sed -i 's/chr//g' ref/high-LD-regions-hg38-GRCh38.txt
  #wget https://raw.githubusercontent.com/meyer-lab-cshl/plinkQC/master/inst/extdata/high-LD-regions-hg19-GRCh37.txt -O ref/high-LD-regions-hg19-GRCh37.txt
  ##get sample ancestry information
  #wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel -O ref/integrated_call_samples_v3.20130502.ALL.panel
  #
  ## get reference data for liftover from hg38 to hg19
  #wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz' -O ref/hg38ToHg19.over.chain.gz
  #wget --timestamping 'https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz' -O ref/hg19ToHg38.over.chain.gz
  #
  ## get reference data for alignment to hg19
  #wget ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz -O ref/human_g1k_v37.fasta.gz
  #gunzip ref/human_g1k_v37.fasta.gz

  ##get reference data for HLA imputation (using HIBAG)
  #mkdir ref/HLA
  #wget https://hibag.s3.amazonaws.com/download/HLARES/AffyAxiomUKB-European-HLA2-hg19.RData -O ref/HLA/AffyAxiomUKB-European-HLA2-hg19.RData
  #wget https://hibag.s3.amazonaws.com/download/HLARES/AffyAxiomUKB-Broad-HLA2-hg19.RData -O ref/HLA/AffyAxiomUKB-Broad-HLA2-hg19.RData
  #wget https://hibag.s3.amazonaws.com/download/HLARES/AffyAxiomUKB-European-HLA4-hg19.RData -O ref/HLA/AffyAxiomUKB-European-HLA4-hg19.RData
  #wget https://hibag.s3.amazonaws.com/download/HLARES/AffyAxiomUKB-Broad-HLA4-hg19.RData -O ref/HLA/AffyAxiomUKB-Broad-HLA4-hg19.RData

  # prepare ref data for the following analysis
  ## 1KG 30x on GRCh38, description on here:https://www.internationalgenome.org/data-portal/data-collection/30x-grch38

  ## I set the variant ID in format: chr:pos_ref_alt, only the first 20 characters were used for alleles with length longer than 20
  # for (i in 1:22) {
  # system(sprintf("/software/team152/plink2 --vcf /lustre/scratch118/humgen/resources/1000g/release/20201028/CCDG_14151_B01_GRM_WGS_2020-08-05_chr%s.filtered.shapeit2-duohmm-phased.vcf.gz --set-all-var-ids chr@:\\#_\\$r_\\$a --new-id-max-allele-len 20 truncate --max-alleles 2 --rm-dup force-first --make-bed --out ref/1KG/chr%s", i,i))
  # }

  ## X chromosome
  ### remove header line at first, otherwise there would be error like "Duplicate FORMAT:GT header line in --vcf file."
  #system('zcat /lustre/scratch118/humgen/resources/1000g/release/20201028/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.vcf.gz|grep -v "^##" > ref/1KG/chr23_tmp')
  #system("/software/team152/plink2 --vcf ref/1KG/chr23_tmp --set-all-var-ids chr@:\\#_\\$r_\\$a --new-id-max-allele-len 20 truncate --max-alleles 2 --rm-dup force-first --make-bed --out ref/1KG/chr23")
  system("rm ref/1KG/chr23_tmp")

  # #get high LD region
  #system("wget https://raw.githubusercontent.com/meyer-lab-cshl/plinkQC/master/inst/extdata/high-LD-regions-hg38-GRCh38.txt -O ref/high-LD-regions-hg38-GRCh38.txt")
  #system("sed -i 's/chr//g' ref/high-LD-regions-hg38-GRCh38.txt")
}

#Convert a VCF to PLINK files
convert_vcf_to_plink = function(raw_vcf_data_directory, # This points to the full path where the raw VCFs we intend to QC are located on lustre
                                batch_numbers,          # This refers to the batch numbers from IBDBR
                                output_directory)       # Which directory (relative) do we want to store this data?
{
  #Think about what working directory we wish to be located in?
  setwd("/lustre/scratch123/hgi/projects/ibdgwas_bioresource/gwas_machine_testing")

  #Loop through the batches and convert the VCF files to PLINK format
  for (i in 1:length(batch_numbers)) {

    sprintf("Processing VCF to PLINK for batch %s",batch_numbers[i])

    #Create a subdirectory to store the batch PLINK files
    system(sprintf("mkdir -p %s", batch_numbers[i]))

    #Remove any SNPs not from chromosome 1-22,X,Y, SNPs (including MT)
    system(sprintf("mkdir -p ./%s/data", batch_numbers[i]))

    #Remove SNPs that are not in autosomal or X/Y/MT chromosomes
    system(sprintf("zcat %s/v2chip_%s_V3_calls.vcf.gz| grep -E -w ^'[1-9]|1[0-9]|2[0-2]|X|Y|MT|#CHROM' > ./%s/tmp",raw_vcf_data_directory,batch_numbers[i],batch_numbers[i]))

    #Convert from VCF to PLINK
    system(sprintf("/software/team152/plink_linux_x86_64_20181202/plink --vcf ./%s/tmp --allow-extra-chr --output-chr MT --make-bed --out ./%s/data/%s",
                                                                                                                             batch_numbers[i],batch_numbers[i],batch_numbers[i] ))
    #Remove temporary data
    system(sprintf("rm ./%s/tmp", batch_numbers[i]))
  }
}

keep_specific_snps = function(){
  system("rm ./snps_to_keep")
  system("echo 'chr16:50729870_C_CC' >> './snps_to_keep'")
  system("echo 'chr16:50729870_CC_C' >> './snps_to_keep'")
  system("echo 'AX-96079897' >> './snps_to_keep'")
}

#Keep only SNPs that are ATCG across batches
keep_atcg_snps = function(output_directory, #Where the results are being stores
                          batch_numbers)    #This refers to the batch numbers from IBDBR
{

  #Loop across batches
  for (i in 1:length(batch_numbers)) {

  #Read in bim tables from PLINK
  bim<-fread(sprintf("%s/%s/data/%s.bim", output_directory, batch_numbers[i], batch_numbers[i]),header = F)

  #Extract the indels
  indels<-bim[which(!bim$V5 %in% c("A","G","C","T") | !bim$V6 %in% c("A","G","C","T")),]

  #Extract other chr SNPs
  otherchr <- bim[ !(bim$V1 %in% c(1:22,"X","Y")),]

  #Remove indels and SNPs on different chrs
  all_remove<-rbind(indels,otherchr)
  all_remove<-all_remove[!duplicated(all_remove),]
  #Write out the indels we have excluded in each batch
  write.table(all_remove[,"V2"],paste0("./", batch_numbers[i],"/data", "/list_indel_var_exclude_",batch_numbers[i]), col.names=F,row.names=F,quote=F,sep="\t")
  system(sprintf("/software/team152/plink_linux_x86_64_20181202/plink --bfile %s/data/%s --allow-extra-chr --output-chr MT --exclude %s/data/list_indel_var_exclude_%s --make-bed --out ./%s/data/%s_ATCG",
                batch_numbers[i], batch_numbers[i], batch_numbers[i], batch_numbers[i], batch_numbers[i], batch_numbers[i]))
  }
}

#Align SNPs to the positive strand
align_positive_strand = function(batch_numbers,
                                 output_directory){

  #Loop through the batches
  for (i in 1:length(batch_numbers)) {

    #Convert PLINK files to VCF to use bcftools
    system(sprintf("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./%s/data/%s_ATCG --allow-extra-chr --output-chr MT --allow-no-sex --recode vcf --out ./%s/data/%s_ATCG",
                       batch_numbers[i], batch_numbers[i],batch_numbers[i], batch_numbers[i]))

    #Run bcftools to align the strands
    system(sprintf("/software/team152/bcftools-1.9/bcftools +/software/team152/bcftools-1.9/plugins/fixref.so ./%s/data/%s_ATCG.vcf -Oz -o ./%s/data/%s_ATCG_aligned.vcf.gz -- -f /lustre/scratch123/hgi/projects/ibdgwas/IIBDGC/resources/hg38/hg38_edited.fa -m top",
                     batch_numbers[i], batch_numbers[i],batch_numbers[i], batch_numbers[i]))

    # After alignment, convert vcf to bed file
    system(sprintf("/software/team152/plink_linux_x86_64_20181202/plink --vcf ./%s/data/%s_ATCG_aligned.vcf.gz --vcf-half-call m --id-delim --make-bed --out ./%s/data/%s_ATCG_aligned",
                   batch_numbers[i], batch_numbers[i],batch_numbers[i], batch_numbers[i]))

    # Delete intermediary file
    system(sprintf("rm ./%s/data/%s_ATCG.vcf", batch_numbers[i], batch_numbers[i]))
  }
}


# #Update the variant IDs to chrX:XXX nomenclature
update_variant_id = function(batch_numbers,
                             output_directory)
{
  for (i in 1:length(batch_numbers)) {
    system(sprintf("zcat ./%s/data/%s_ATCG_aligned.vcf.gz|grep -v ^\"#\"|cut -f '1-5' | awk '{print $3,\"chr\"$1\":\"$2\"_\"$4\"_\"$5}' > ./%s/data/%s_ATCG_aligned",batch_numbers[i],batch_numbers[i],batch_numbers[i],batch_numbers[i]))
  }
}

#Remove missing SNPs
remove_missing_snps = function(batch_numbers,
                               output_directory)
{
  for (i in 1:length(batch_numbers)){
    #Use PLINK default missingness threshold (10%)
    system(sprintf("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./%s/data/%s_ATCG_aligned --missing --out ./%s/data/%s_ATCG_aligned",batch_numbers[i],batch_numbers[i],batch_numbers[i],batch_numbers[i]))
  }
}

#Remove duplicated variants
remove_duplicated_variants = function(batch_numbers,
                                      output_directory)
{
  for (i in 1:length(batch_numbers)) {

    bim<-read.table(sprintf("./%s/data/%s_ATCG_aligned.bim",batch_numbers[i],batch_numbers[i]),sep="\t",head=F)
    ids<-read.table(sprintf("./%s/data/%s_ATCG_aligned",batch_numbers[i],batch_numbers[i]),header=F,stringsAsFactor=F) #skip line number might be differ

    varmiss<-read.table(sprintf("./%s/data/%s_ATCG_aligned.lmiss",batch_numbers[i],batch_numbers[i]),sep="",head=T)

    colnames(ids)[2]<-"ids"

    #the major allele is set to A2 by default by Plink, keep ids with real ref/alt as in vcf using the ids file
    bim.1<-cbind(bim,ids[,"ids",drop=F])

    bim.1<-merge(bim.1,varmiss[,c("SNP","F_MISS")],by.x="V2",by.y="SNP",sort=F)

    bim.1$ids<-as.character(bim.1$ids)

    # identify duplicated variants (same chr position ref and alt)
    dups<-bim.1[which(duplicated(bim.1$ids)),"ids"]
    dups <- unique(dups)

    bim.1$V2 <- as.character(bim.1$V2)
    bim.1$idx <- 1:nrow(bim.1)

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

    write.table(duplicated_variants,sprintf("%s/data/list_duplicated_var_exclude_%s",batch_numbers[i], batch_numbers[i]),col.names=F,row.names=F,quote=F,sep="\t")
    write.table(bim.1[,c(2,7,3:6)],sprintf("%s/data/%s_ATCG_aligned_edited.bim",batch_numbers[i], batch_numbers[i]),col.names=F,row.names=F,quote=F,sep="\t")

    system(sprintf("/software/team152/plink_linux_x86_64_20181202/plink --bed %s/data/%s_ATCG_aligned.bed --bim %s/data/%s_ATCG_aligned_edited.bim --fam %s/data/%s_ATCG_aligned.fam --exclude %s/data/list_duplicated_var_exclude_%s --make-bed --out %s/data/%s_ATCG_aligned_nodup",
                   batch_numbers[i],batch_numbers[i],batch_numbers[i],batch_numbers[i],batch_numbers[i],batch_numbers[i],batch_numbers[i],batch_numbers[i],batch_numbers[i],batch_numbers[i]))
  }
}

merge_batches = function(batch_numbers,
                         output_directory)
{
  #Remove old merge file for clean run
  system("rm ./merge.list")

  #Get the batches (except the most recent as a merged list for PLINK 1.9)
  for (i in 1:(length(batch_numbers)-1)) {
    system(sprintf("echo ./%s/data/%s_ATCG_aligned_nodup.bed ./%s/data/%s_ATCG_aligned_nodup.bim ./%s/data/%s_ATCG_aligned_nodup.fam >> ./merge.list",
                   batch_numbers[i],batch_numbers[i],batch_numbers[i],batch_numbers[i],batch_numbers[i],batch_numbers[i]))
  }
  #Merge batches together with most recent batch
  system(sprintf("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./%s/data/%s_ATCG_aligned_nodup --merge-list ./merge.list --make-bed --out ./AllBatch_ATCG_aligned_nodup",
        batch_numbers[length(batch_numbers)],batch_numbers[length(batch_numbers)]))
}

add_gender =  function(output_directory)
{
  read.csv("/lustre/scratch123/hgi/mdt2/projects/ibdgwas_bioresource/phenotype_data/release_20220323/raw/IBD_BioRes_phenotypes_20220323_1.csv",stringsAsFactor=F,sep="~",header=T)->data

  read.table("./AllBatch_ATCG_aligned_nodup.fam",stringsAsFactor=F)->fam

  famnew <- merge(fam,data[,c("SP0011_id_1","gender")],by.y="SP0011_id_1",by.x="V1")
  ## 1: male, 2: female
  famnew$V5 <- famnew$gender
  famnew <- famnew[,c(1:6)]

  ## for samples without gender info, set it to zero
  famnew[is.na(famnew$V5),"V5"]<- 0

  write.table(famnew,file= "./AllBatch_ATCG_aligned_nodup_gender",row.names=F,quote=F,col.names=F)

  system("/software/team152/plink_linux_x86_64_20181202/plink --bed ./AllBatch_ATCG_aligned_nodup.bed --bim ./AllBatch_ATCG_aligned_nodup.bim --fam ./AllBatch_ATCG_aligned_nodup_gender --make-bed --out ./AllBatch_ATCG_aligned_nodup_gender")
}

calculate_missigness = function(){
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./AllBatch_ATCG_aligned_nodup --missing --out ././AllBatch_ATCG_aligned_nodup")
}

identify_duplicates = function(){

  #Run KING inference
  system("/software/team152/king -b ./AllBatch_ATCG_aligned_nodup_gender.bed --related --degree 1 --prefix ./AllBatch_Duplicates")

  #### Identify the intra family/individual relatedness ###
  remove_intra_family_relatedness = c()
  kin                              = read.table("AllBatch_Duplicates.kin",header = T)

  ## Check for samples with same FID, but not large kinship, remove these since we are unsure about genotype and can't match to phenotype
  intra_family_no_duplicate = c()
  intra_family_no_duplicate = kin[which(kin$Kinship<=0.354),c("FID","ID1","ID2","Kinship")]

  #Need to check for the samples between both the 1st and 2nd sample in pair (ID1 & ID2)
  remove_intra_family_relatedness_ID1            = intra_family_no_duplicate[,c("FID","ID1")]
  remove_intra_family_relatedness_ID2            = intra_family_no_duplicate[,c("FID","ID2")]
  colnames(remove_intra_family_relatedness_ID1)  = c("FID","IID")
  colnames(remove_intra_family_relatedness_ID2)  = c("FID","IID")
  remove_intra_family_relatedness                = rbind(remove_intra_family_relatedness_ID1,remove_intra_family_relatedness_ID2)
  remove_intra_family_relatedness$reason         = "IBDBR_same_FID_low_kinship"

  #### Identify INTRA family relatedness, remove sample with lower call-rate ###
  duplicate_pairs_intra_family = kin[which(kin$Kinship>0.354),c("FID","ID1","ID2","Kinship")]
  sample_miss     = read.table("./AllBatch_ATCG_aligned_nodup.imiss",head=T)

  #Extract the unique FID/ID pairs from the inter family kinship
  duplicate_pairs_intra_family_ID1           = duplicate_pairs_intra_family[,c("FID","ID1")]
  duplicate_pairs_intra_family_ID2           = duplicate_pairs_intra_family[,c("FID","ID2")]
  colnames(duplicate_pairs_intra_family_ID1) = c("FID","IID")
  colnames(duplicate_pairs_intra_family_ID2) = c("FID","IID")
  duplicate_pairs_intra_family               = rbind(duplicate_pairs_intra_family_ID1,duplicate_pairs_intra_family_ID2)
  duplicate_pairs_intra_family               = unique(duplicate_pairs_intra_family)
  #Change the FID/ID id to be FID_ID coded
  rownames(duplicate_pairs_intra_family)     = paste0(duplicate_pairs_intra_family$FID,"_",duplicate_pairs_intra_family$IID)
  rownames(sample_miss)                      = paste0(sample_miss$FID,"_",sample_miss$IID)
  #Extract the samples that are in the duplicate pairs to the sample missingness
  sample_miss                                = sample_miss[rownames(duplicate_pairs_intra_family),]
  #Add missingness column
  duplicate_pairs_intra_family$miss          = sample_miss$F_MISS

  #Keep track of IDs that have the higher call rate
  higher_callrate_duplicate_pairs_intra_family_ID <- c()
  #Loop through the kinship results to extract the sample with lower call rate between the pair
  for(id in unique(duplicate_pairs_intra_family$FID)){
      #Extract all rows which have the same FID
      datatmp                              = duplicate_pairs_intra_family[duplicate_pairs_intra_family$FID==id,]
      #Pick the sample that has the higher call-rate
      higher_callrate_duplicate_pairs_intra_family_ID = c(higher_callrate_duplicate_pairs_intra_family_ID,rownames(datatmp[ datatmp$miss<max(datatmp$miss),]))
  }

  #Remove the higher call-rate IDs in the duplicate pairs intra family
  duplicate_pairs_intra_family                   = duplicate_pairs_intra_family[higher_callrate_duplicate_pairs_intra_family_ID,]
  duplicate_pairs_intra_family$reason            = "IBDBR_reported_est_dup"
  remove_intra_family_relatedness                = rbind(remove_intra_family_relatedness,duplicate_pairs_intra_family[,c(1,2,4)])

  #### Identify the INTER family/individual relatedness ###
  remove_inter_family_relatedness = c()
  kin                             = read.table("./AllBatch_Duplicates.kin0",head=T)
  #change FID as FID_ID, the reason to do this is sometimes different samples could have same FID, but different IID
  kin$FID1                        = paste0(kin$FID1,"_",kin$ID1)
  kin$FID2                        = paste0(kin$FID2,"_",kin$ID2)

  #exclude the samples if they have been removed from the INTRA familial relatedness calculations
  if(!is.null(remove_intra_family_relatedness)){
    kin = kin[!(kin$FID1%in%paste0(remove_intra_family_relatedness$FID,"_",remove_intra_family_relatedness$IID)|kin$FID2%in%paste0(remove_intra_family_relatedness$FID,"_",remove_intra_family_relatedness$IID)),]
  }

  #Extract the duplicates from the kinship file
  duplicate_pairs_inter_family    = kin[which(kin$Kinship>0.354),c("FID1","FID2","Kinship")]
  duplicate_pairs_inter_family_ID = c(as.character(duplicate_pairs_inter_family$FID1),as.character(duplicate_pairs_inter_family$FID2))

  #Read the fam file and add phenotype & gender information
  fam                = read.table("./AllBatch_ATCG_aligned_nodup_gender.fam",head=F)
  #change FID to FID_IID
  fam$V1             = paste0(fam$V1,"_",fam$V2)
  #Add in phenotype and sex information
  colnames(fam)[5:6] = c("sex","pheno")
  sample_miss        = read.table("./AllBatch_ATCG_aligned_nodup.imiss",head=T)
  #change FID as FID_IID
  sample_miss$FID    = paste0(sample_miss$FID,"_",sample_miss$IID)
  #The all file contains the FAM information along with the phenotype/gender information (ie. ALL information)
  all<-merge(fam[,c(1,5,6)],sample_miss[,c("FID","F_MISS")],by.x="V1",by.y="FID",all.x=T,sort=F)

  #Only include IDs which are in the full fam + phenotype file (WHY???)
  tmp                             = all[which(all$V1 %in% duplicate_pairs_inter_family_ID),]
  duplicate_pairs_inter_family_ID = duplicate_pairs_inter_family_ID[match(tmp$V1,duplicate_pairs_inter_family_ID)]
  rm(tmp)

  #Loop through duplicate IDs for the INTER family relatedness
    for (i in 1:length(duplicate_pairs_inter_family_ID)) {
      #Extract the rows in the KING duplicate table that match the FID for both FID1 and FID2
      tmp1 = duplicate_pairs_inter_family[which(duplicate_pairs_inter_family$FID1==duplicate_pairs_inter_family_ID[i]),]
      tmp2 = duplicate_pairs_inter_family[which(duplicate_pairs_inter_family$FID2==duplicate_pairs_inter_family_ID[i]),]
      colnames(tmp2)[1:2]<-colnames(tmp2)[2:1]
      #Bind together the results
      tmp<-rbind(tmp1,tmp2)
      ids_tmp<-c(as.character(tmp$FID1),as.character(tmp$FID2))
      #Remove any duplicated IDs between these rows
      ids_tmp<-ids_tmp[!duplicated(ids_tmp)]
      #keep same order as in dup_ids and in all
      ids_tmp<-ids_tmp[match(duplicate_pairs_inter_family_ID[which(duplicate_pairs_inter_family_ID %in% ids_tmp)],ids_tmp)]

      #Get number of possible combinations between all IDs that we've extracted
      #This deals with cases when we have pairs like: ID1 - ID2 pair & ID1 - ID3 pair we would expect to see ID2 -ID3 pair
      n_possible_combinations<-nrow(permutations(length(ids_tmp), 2))/2
      if (nrow(duplicate_pairs_inter_family[which((duplicate_pairs_inter_family$FID1 %in% ids_tmp) | (duplicate_pairs_inter_family$FID2 %in% ids_tmp)),])==n_possible_combinations ) {
        #For duplicated samples that have same phenotype and sex, remove the ones with lower call rate
        data<-as.data.frame(matrix(ncol=1,nrow=length(ids_tmp)))
        data$V1<-ids_tmp
        #Merge together the phenotype/gender information with these duplicated samples
        data<-merge(data,all,by="V1",all.x=T,sort=F)
        if ((dim(table(data$sex))==1 & dim(table(data$pheno))==1) ) {
          #Concordant gender and/or phenotype between duplicate pairs
          if(!exists("remove_inter_family_relatedness")) {
            #Remove sample with lower call rate
            keep_sample<-data$V1[which(data$F_MISS==min(data$F_MISS))][1]
            remove_inter_family_relatedness<-data[which(!data$V1 %in% keep_sample),]
          } else {
            #Remove sample with the lower call rate
            keep_sample<-data$V1[which(data$F_MISS==min(data$F_MISS))][1]
            remove_inter_family_relatedness<-rbind(remove_inter_family_relatedness,data[which(!data$V1 %in% keep_sample),])
          }
        } else {
          #Disconcordant gender and/or phenotype between duplicate pairs
          data<-as.data.frame(matrix(ncol=1,nrow=length(ids_tmp)))
          data$V1<-ids_tmp
          data<-merge(data,all,by="V1",all.x=T,sort=F)

          if(!exists("remove_inter_family_relatedness")) {
            remove_inter_family_relatedness<-data
          } else {
            remove_inter_family_relatedness<-rbind(remove_inter_family_relatedness,data)
          }
          if(!exists("data_inconsist")){

            jj<-1
            data_inconsist<-data
            data_inconsist$group<-jj

          } else {
            jj<-jj+1
            data$group<-jj
            data_inconsist<-rbind(data_inconsist,data)
          }
        }

      } else {
        #Not all combinations of duplicated pairs found, show list of IDs, to manually inspect issues, why would this be the case?
        print(paste("Number of expected combinations: ",n_possible_combinations,sep=""))
        print(paste("Number of observed combinations: ",nrow(duplicate_pairs_inter_family[which((duplicate_pairs_inter_family$FID1 %in% ids_tmp) | (duplicate_pairs_inter_family$FID2 %in% ids_tmp)),]),sep=""))
        print(ids_tmp)
      }
    }

  #Remove duplicate IDs from the samples we are removing based on the inter family relatedness and extract the IDs to write out
  remove_inter_family_relatedness           = remove_inter_family_relatedness[!duplicated(remove_inter_family_relatedness$V1),]
  remove_inter_family_relatedness           = remove_inter_family_relatedness[,c(1,1)]
  #Add column names
  colnames(remove_inter_family_relatedness) = c("FID","IID")
  #Replace FID_ID coding with seperate FID and ID coding
  ids                                       = unlist(strsplit(remove_inter_family_relatedness$FID,split="_"))
  remove_inter_family_relatedness$FID       = ids[seq(1,length(ids),2)]
  remove_inter_family_relatedness$IID       = ids[seq(2,length(ids),2)]
  #Add the reason for the removal of these IDs
  remove_inter_family_relatedness$reason    = "not_IBDBR_reported_est_dup"
  #Bind together both the inter and intra family duplicated sample IDs for removal
  duplicated_samples                        = rbind(remove_inter_family_relatedness,remove_intra_family_relatedness)
  write.table(duplicated_samples,"./list_duplicated_samples",col.names=F,row.names=F,quote=F,sep="\t")

  #Remove these duplicates for (some) of the subsequent analysis
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./AllBatch_ATCG_aligned_nodup_gender --remove ./list_duplicated_samples --make-bed --out ./AllBatch_ATCG_aligned_nodup_gender_nodupsample")
}

#Compare freq with EUR 1000GP to deal with A/T   C/G strand exclude A/T and C/G variants with MAF >0.45 (hard to determine if they are wrong)
compare_freq_1000G_EUR = function(){
  #Create the 1000G EUR data directory
  system("mkdir -p ./ref/1KG_EUR_study_variant")
  #Extract the list of variants from our merged batches
  system("cat ./AllBatch_ATCG_aligned_nodup_gender_nodupsample.bim | cut -f 2 > ./ref/1KG_EUR_study_variant/list_variants_posstr_nodup")
  #Extract the list of ...
  system('awk \'{if($3=="EUR")print 0,$1}\' < /lustre/scratch118/humgen/resources/1000g/release/20130502/integrated_call_samples_v3.20130502.ALL.panel > ./ref/1KG_EUR_study_variant/eur.sample')

  for (i in 1:23) {
    system(sprintf("/software/team152/plink2 --bfile ./ref/1KG/chr%s --keep ./ref/1KG_EUR_study_variant/eur.sample --rm-dup force-first --extract ./ref/1KG_EUR_study_variant/list_variants_posstr_nodup --make-bed --out ./ref/1KG_EUR_study_variant/chr%s",
           i,i))
  }

  system("rm ./merge.list")
  for (i in 1:23) {
    system(sprintf('echo "./ref/1KG_EUR_study_variant/chr%s.bed ./ref/1KG_EUR_study_variant/chr%s.bim ./ref/1KG_EUR_study_variant/chr%s.fam" >> ./merge.list',
                   i,i,i))
  }

  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./ref/1KG_EUR_study_variant/chr1 --merge-list ./merge.list --allow-no-sex --make-bed --out ./ref/1KG_EUR_study_variant/1000GP_EUR_b38_study_variants")

  #Compare MAF between 1KG EUR and study, and flip allele
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./ref/1KG_EUR_study_variant/1000GP_EUR_b38_study_variants --freq --out ./ref/1KG_EUR_study_variant/1000GP_EUR_b38_study_variants")
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./AllBatch_ATCG_aligned_nodup_gender_nodupsample --freq --out ./AllBatch_ATCG_aligned_nodup_gender_nodupsample")

  #Read in the frequencies of the GWAS study and the 1KG panel (EUR)
  gp<-read.table("./ref/1KG_EUR_study_variant/1000GP_EUR_b38_study_variants.frq",head=T)
  g1<-read.table("./AllBatch_ATCG_aligned_nodup_gender_nodupsample.frq",head=T)

  #Reconfigure column names
  colnames(gp)[3:6]<-paste(colnames(gp)[3:6],"_gp",sep="")
  colnames(g1)[3:6]<-paste(colnames(g1)[3:6],"_g1",sep="")

  #Merge the frequencies from 1000G and from the GWAS
  all<-merge(g1,gp[,2:6],by="SNP",all=T)

  #Extract the A1 allele
  check<-all[which(all$A1_g1!=all$A1_gp),]
  # keep only A/T C/G
  check<-check[which( (check$A1_g1=="G" & check$A2_g1=="C") | (check$A1_g1=="C" & check$A2_g1=="G") | (check$A1_g1=="A" & check$A2_g1=="T") | (check$A1_g1=="T" & check$A2_g1=="A")),]

  # LIST OF VARIANTS TO REMOVE, WE CANNOT REALLY BE SURE WHETHER THERE IS STRAND ISSUE OR NOT
  remove<-check[which(check$MAF_g1>=0.45),]

  # LIST OF VARIANTS TO FLIP:
  flip<-check[which(check$MAF_g1<0.45),]
  flip<-flip[order(flip$MAF_g1,decreasing=T),]

  #Check disconcordance in MAF between 1KG EUR and GWAS
  remove_2<-flip[which(flip$MAF_g1>0.2 & flip$MAF_gp<0.1),]
  remove_3<-flip[which(flip$MAF_g1<0.1 & flip$MAF_gp>0.2),]
  remove<-rbind(remove,remove_2,remove_3)
  write.table(remove[,"SNP"],"./list_variants_to_remove_AT_CG",col.names=F,row.names=F,quote=F,sep="\t")
  flip<-flip[which(!flip$SNP %in% remove$SNP),]

  write.table(flip[,"SNP"],"./list_variants_to_flip_AT_CG",col.names=F,row.names=F,quote=F,sep="\t")

  system("grep -v -w -f ./snps_to_keep ./list_variants_to_remove_AT_CG > ./list_variants_to_remove_AT_CG_rmkeep")
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./AllBatch_ATCG_aligned_nodup_gender_nodupsample --exclude ./list_variants_to_remove_AT_CG_rmkeep --flip ./list_variants_to_flip_AT_CG --make-bed --out ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip")
}

sex_discrepancies_check = function(){
  system("mkdir -p ./sex")
  system("mkdir -p ./sex/X")
  system("mkdir -p ./sex/Y")

  #Extract the sexes (M/F) of the samples in the FAM file
  system("awk '{if($5==1)print $1,$2}' < ./AllBatch_ATCG_aligned_nodup_gender_nodupsample.fam > ./sex/list_male_samples")
  system("awk '{if($5==2)print $1,$2}' < ./AllBatch_ATCG_aligned_nodup_gender_nodupsample.fam > ./sex/list_female_samples")

  ##X chromosome
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./AllBatch_ATCG_aligned_nodup_gender_nodupsample --keep ./sex/list_female_samples --chr 23 --make-bed --out ./sex/X/study_check_sex_females_only")
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./sex/X/study_check_sex_females_only --missing --hardy --freq --out ./sex/X/study_check_sex_females_only")

  ### get high quality variants
  bim      = read.table("./sex/X/study_check_sex_females_only.bim",sep="\t",head=F)
  hwe      = read.table("./sex/X/study_check_sex_females_only.hwe",head=T)
  frq      = read.table("./sex/X/study_check_sex_females_only.frq", sep="",head=T)
  var_miss = read.table("./sex/X/study_check_sex_females_only.lmiss",sep="",head=T)

  all = merge(hwe[,c("SNP","P")],frq[,c("SNP","MAF")],by="SNP",sort=F)
  all = merge(all,var_miss[,c("SNP","F_MISS")],by="SNP",sort=F)

  all = all[which(all$MAF>0.05 & all$F_MISS<0.01 & all$P>1e-12),]

  write.table(all,"./sex/X/list_good_chrX_variants",col.names=F,row.names=F,quote=F,sep="\t")

  ##Y chromosome
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./AllBatch_ATCG_aligned_nodup_gender_nodupsample --keep ./sex/list_male_samples --chr 24 --make-bed --out ./sex/Y/study_check_sex_males_only")
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./sex/Y/study_check_sex_males_only --freq --missing --out ./sex/Y/study_check_sex_males_only")
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./AllBatch_ATCG_aligned_nodup_gender_nodupsample --keep ./sex/list_female_samples --chr 24 --recode A --out ./sex/Y/study_check_sex_females_only_chr24")

  ### get high quality variants
  frq      = read.table("./sex/Y/study_check_sex_males_only.frq",head=T)
  var_miss = read.table("./sex/Y/study_check_sex_males_only.lmiss",head=T)

  all = merge(frq[,c("SNP","MAF","NCHROBS")],var_miss[,c("SNP","F_MISS")],by="SNP",sort=F)

  ## females, find variants with less % of calls:
  ped<-read.table("./sex/Y/study_check_sex_females_only_chr24.raw",head=T,check.names=F)

  dat<-matrix(nrow=nrow(all),ncol=2)
  dat<-as.data.frame(dat)
  colnames(dat)<-c("variant","percentage_NA")

  for (i in 1:nrow(dat)) {
    tmp<-ped[,6+i,drop=F]
    dat$variant[i]<-gsub("_[A-Z]{1}$","",colnames(tmp))
    dat$percentage_NA[i]<-nrow(tmp[which(is.na(tmp)),,drop=F])/nrow(tmp)
  }

  all<-all[which(all$SNP %in% dat$variant[which(dat$percentage_NA>0.95)]),]

  #Keep variants with large number of calls in males
  all<-all[which(all$NCHROBS>=max(all$NCHROBS,na.rm=T)-(max(all$NCHROBS,na.rm=T)*0.005)),]

  write.table(all,"./sex/Y/list_good_chrY_variants",col.names=F,row.names=F,quote=F,sep="\t")

  system("cat ./sex/X/list_good_chrX_variants ./sex/Y/list_good_chrY_variants > ./sex/list_good_chrXY_variants")

  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./AllBatch_ATCG_aligned_nodup_gender_nodupsample --extract ./sex/list_good_chrXY_variants --make-bed --out ./sex/study_check_sex_good_chrXY")
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./sex/study_check_sex_good_chrXY --check-sex ycount --out ./sex/study_check_sex_good_chrXY")
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./sex/study_check_sex_good_chrXY --check-sex --out ./sex/study_check_sex_good_chrX")

  system("mkdir -p ./figures")

  sex<-read.table("./sex/study_check_sex_good_chrXY.sexcheck",sep="",head=T)
  #change FID to FID_IID since sometimes samples could have same FID
  sex$FID <- paste0(sex$FID,"_",sex$IID)
  sex$PEDSEX<-as.factor(sex$PEDSEX)

  # X thresholds:
  fmin_male<-mean(sex[which(sex$PEDSEX==1),"F"])-(4*sd(sex[which(sex$PEDSEX==1),"F"]))
  fmax_female<-mean(sex[which(sex$PEDSEX==2),"F"])+(4*sd(sex[which(sex$PEDSEX==2),"F"]))

  # Y thresholds:
  ymin_male<-mean(sex[which(sex$PEDSEX==1),"YCOUNT"])-(4*sd(sex[which(sex$PEDSEX==1),"YCOUNT"]))
  ymax_female<-mean(sex[which(sex$PEDSEX==2),"YCOUNT"])+(4*sd(sex[which(sex$PEDSEX==2),"YCOUNT"]))

  p1n<-ggplot(sex[which(sex$PEDSEX==1),],aes(x=F,y=YCOUNT)) + geom_point(colour = "#00AFBB") +
    theme_bw() +
    xlim(min(sex$F),max(sex$F)) +
    ylim(0,max(sex$YCOUNT)) +
    geom_hline(yintercept=ymin_male, linetype="dashed", color = "#00AFBB")+
    geom_vline(xintercept=fmin_male, linetype="dashed", color = "#00AFBB") +
    theme(axis.text= element_text(size=14,face="bold"),axis.title=element_text(size=14,face="bold"))

  p2n<-ggplot(sex[which(sex$PEDSEX==2),],aes(x=F,y=YCOUNT)) + geom_point(colour = "#E7B800") +
    theme_bw() +
    xlim(min(sex$F),max(sex$F)) +
    ylim(0,max(sex$YCOUNT)) +
    geom_hline(yintercept=ymax_female, linetype="dashed", color = "#E7B800")+
    geom_vline(xintercept=fmax_female, linetype="dashed", color = "#E7B800") +
    theme(axis.text= element_text(size=14,face="bold"), axis.title=element_text(size=14,face="bold"))

  pn<-ggarrange(p1n,p2n,ncol=2,labels=c("Male","Female"))

  ggsave(file="./figures/F_ycount.png",width=12,height=6,dpi=300)

  sex$SNPSEX2<-0
  sex$SNPSEX2[which(sex$F<=fmax_female & sex$YCOUNT<=ymax_female)]<-2
  sex$SNPSEX2[which(sex$F>=fmin_male & sex$YCOUNT>=ymin_male)]<-1

  sex_chrx<-read.table("./sex/study_check_sex_good_chrX.sexcheck",head=T)
  #change FID to FID_IID since sometimes samples could have same FID
  sex_chrx$FID <- paste0(sex_chrx$FID,"_",sex_chrx$IID)

  colnames(sex_chrx)[4]<-"SNPSEX_chrx"
  sexm<-merge(sex,sex_chrx[,c("FID","SNPSEX_chrx")],by="FID",sort=F)

  samples_exclude<-sex[which((sex$PEDSEX==2 & sex$SNPSEX2==1) | (sex$PEDSEX==1 & sex$SNPSEX2==2)),]

  if(nrow(samples_exclude)>0){
    id <- unlist(strsplit(samples_exclude$FID,split="_"))
    samples_exclude$FID <- id[seq(1,length(id),2)]
  }

  write.table(samples_exclude[,c("FID","IID")],"./list_samples_wrong_gender",col.names=T,row.names=F,quote=F,sep="\t")
  write.table(samples_exclude[,c("FID","PEDSEX","SNPSEX2","F","YCOUNT")],"./study_list_samples_sex_discrepancy",col.names=T,row.names=F,quote=F,sep="\t")

  #Recode the rest:
  fam<-read.table("./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip.fam",sep="",head=F)
  fam$FID_IID <- paste0(fam$V1,"_",fam$V2)
  fam_ed<-merge(fam,sex[,c("FID","SNPSEX2")],by.x="FID_IID",by.y="FID",sort=F)
  fam_ed<-fam_ed[,c("V1","V2","V3","V4","SNPSEX2","V6")]

  write.table(fam_ed,"./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_edited.fam",col.names=F,row.names=F,quote=F,sep="\t")

  system("/software/team152/plink_linux_x86_64_20181202/plink --bed ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip.bed --bim ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip.bim --fam ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_edited.fam --allow-no-sex --remove ./list_samples_wrong_gender --make-bed --out ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck")
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck --allow-no-sex --set-hh-missing --make-bed --out ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh")
}

remove_Y_chr = function(){
  system("awk '{if($1==24)print $2}' < ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh.bim > ./list_chry_variants_toremove")
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh --allow-no-sex --exclude ./list_chry_variants_toremove --make-bed --out ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry")
}

remove_low_call_rate = function(){
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry --allow-no-sex --mind 0.20 --make-bed --out ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8")

  # Variant Call Rate < 80%
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8 --allow-no-sex --missing --out ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8")
  system("awk '{if($5>0.2) print $2}' ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8.lmiss > ./list_variants_to_remove_vcr0.8")
  system("grep -v -w -f ./snps_to_keep ./list_variants_to_remove_vcr0.8 > ./list_variants_to_remove_vcr0.8_rmkeep")
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8 --allow-no-sex --exclude ./list_variants_to_remove_vcr0.8_rmkeep --make-bed --out ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8")

  # Sample Call Rate < 95%
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8 --allow-no-sex --mind 0.05 --make-bed --out ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95")

  #Variant Call Rate < 95%
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95 --allow-no-sex --missing --out ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95")
  system("awk '{if($5>0.05) print $2}' ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95.lmiss > ./list_variants_to_remove_vcr0.95")
  system("grep -v -w -f ./snps_to_keep ./list_variants_to_remove_vcr0.95 > ./list_variants_to_remove_vcr0.95_rmkeep")
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95 --allow-no-sex --exclude ./list_variants_to_remove_vcr0.95_rmkeep --make-bed --out ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95")

  # Variant MAF<0.01 AND Call Rate <98%
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95 --allow-no-sex --missing --out ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95")
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95 --allow-no-sex --freq --out ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95")

  sample_miss = read.table("./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95.imiss",head=T)
  var_miss    = read.table("./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95.lmiss",head=T)
  frq         = read.table("./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95.frq",head=T)

  var   = merge(frq[,c(2:6)],var_miss,by="SNP")
  var.1 = var[which(var$MAF<0.01 & var$F_MISS>0.02),]

  write.table(var.1[,"SNP",drop=F],"./list_monomorphic_vcr098maf0.01_var_exclude",col.names=F,row.names=F,quote=F,sep="\t")

  system("grep -v -w -f ./snps_to_keep ./list_monomorphic_vcr098maf0.01_var_exclude > ./list_monomorphic_vcr098maf0.01_var_exclude_rmkeep")
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95 --allow-no-sex --exclude ./list_monomorphic_vcr098maf0.01_var_exclude_rmkeep --make-bed --out ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_vcr0.98maf0.01")
}

variant_missigness_across_batches = function(batches = c("b04","b06","b15","b19")){
  for( i in 1:length(batches)){
    system(sprintf("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95 --missing --keep ./%s/data/%s.fam --out ./tmp_%s",
    batches[i],batches[i],batches[i]))
  }
  ## three parameters: 1. work_dir; 2. batches investigated (only b04,b06,b15,b19 were considered); 3. summary of variant missingness across batches; 4. file name of list of variants with missingness > 0.1 per study/batch
  #/software/R-4.1.0/bin/Rscript scripts/11_missing_per_study.R ${wkdir}/data "b04 b06 b15 b19" "table_breakdown_number_missing_variants_perstudy.csv" "list_variants_per_study_missing_0.1"

  for (i in 1:length(batches)){
    if(i==1){
      tmp<-read.table(sprintf("./tmp_%s.lmiss",batches[i]),head=T)
      a1<-as.data.frame(table(cut(tmp$F_MISS,breaks=c(-1,0.05,0.1,0.2,0.3,0.4),labels=c("0-0.05","0.05-0.1","0.1-0.2","0.2-0.3","0.3-0.4"))) )
      colnames(a1)<-c("numbers",batches[i])
    }else{
      tmp<-read.table(sprintf("./tmp_%s.lmiss",batches[i]),head=T)
      a2<-as.data.frame(table(cut(tmp$F_MISS,breaks=c(-1,0.05,0.1,0.2,0.3,0.4),labels=c("0-0.05","0.05-0.1","0.1-0.2","0.2-0.3","0.3-0.4"))) )
      colnames(a2)<-c("numbers",batches[i])
      a1<-merge(a1,a2,by="numbers",sort=F)
    }
  }

  write.table(a1,"./table_breakdown_number_missing_variants_perstudy.csv",col.names=T,row.names=F,quote=F,sep=",")

  for (i in 1:length(batches)){
    if (i==1) {
      tmp<-read.table(sprintf("tmp_%s.lmiss",batches[i]),head=T)
      tmp$study<-batches[i]
      tmp<-tmp[,c("SNP","F_MISS","study")]
      all<-tmp
    } else {
      tmp<-read.table(sprintf("tmp_%s.lmiss",batches[i]),head=T)
      tmp$study<-batches[i]
      tmp<-tmp[,c("SNP","F_MISS","study")]
      all<-rbind(all,tmp)
    }
  }

  variants<-all[which(all$F_MISS>0.1),"SNP",drop=F]
  variants<-variants[!duplicated(variants$SNP),,drop=F]
  browser()
  write.table(variants,"./list_variants_per_study_missing_0.1",col.names=F,row.names=F,quote=F)

  system("grep -v -w -f ./snps_to_keep ./list_variants_per_study_missing_0.1 > ./list_variants_per_study_missing_0.1_rmkeep")
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95 --allow-no-sex --exclude ./list_variants_per_study_missing_0.1 --make-bed --out ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy")

  system("rm ./tmp_*")
}

PCA_Ancestry = function(){

  system("mkdir -p ./PCA")
  #Remove variants in regions with high LD
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy --exclude ./ref/high-LD-regions-hg38-GRCh38.txt --range --make-bed --out PCA/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_noHighLD")

  #Keep independent SNPs, make sure these SNPs also in 1KG
  system("awk '{print $2}' < ./ref/1KG_EUR_study_variant/1000GP_EUR_b38_study_variants.bim > ./PCA/1KG.snp")
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile PCA/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_noHighLD --indep-pairwise 50 5 0.2 --extract PCA/1KG.snp --out PCA/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_noHighLD")

  system('grep -v /"^chrX/" PCA/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_noHighLD.prune.in > PCA/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_noHighLD_nosex.prune.in')

  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile PCA/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_noHighLD --extract PCA/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_noHighLD_nosex.prune.in --make-bed --out PCA/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_noHighLD_pruned")

  #Get 1KG with study variants
  system("mkdir -p ./PCA/1KG")

  for (i in 1:22) {
    system(sprintf("/software/team152/plink2 --bfile ./ref/1KG/chr%s --extract ./PCA/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_noHighLD.prune.in --make-bed --out chr%s",i,i))
  }

  system("rm ./merge.list")

  for (i in 2:22){
  system(sprintf('echo "chr%s.bed chr%s.bim chr%s.fam" >> merge.list',i,i,i))
  }

  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile chr1 --merge-list merge.list --allow-no-sex --make-bed --out ./ref/1KG/1000GP_ALL_b38_study_variants")

  #Combine 1KG with study sample
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./ref/1KG/1000GP_ALL_b38_study_variants --bmerge ./PCA/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_noHighLD.bed ./PCA/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_noHighLD.bim ./PCA/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_noHighLD.fam --make-bed --out ./ref/1KG/1000GP_ALL_b38_plus_study")

  #PCA analysis
  system("bsub -J ./PCA/PCA_AllBatch_1KG -q normal -o ./PCA/PCA_AllBatch_1KG.log -e ./PCA/PCA_AllBatch_1KG.err -W 10:00 -m \"modern_hardware\" -M 40000 -R'select[mem>40000] rusage[mem=40000]' /software/team152/plink_linux_x86_64_20181202/plink --memory 39000 --bfile ./ref/1KG/1000GP_ALL_b38_plus_study --pca --out ./PCA/AllBatch_1KG")

  #Get study samples with different ancestry
  # ## four parameters: 1. work_dir; 2. 1KG population file; 3. prefix of input file(s); 4. file name of PCA result
  # /software/R-4.1.0/bin/Rscript scripts/13_keep_all_ancestry.R "${wkdir}/PCA" "/lustre/scratch115/resources/1000g/release/20130502/integrated_call_samples_v3.20130502.ALL.panel" "${filename}_noHighLD_pruned" "AllBatch_1KG"

  #Read in the 1000G population samples
  samples1000 = read.table("/lustre/scratch118/humgen/resources/1000g/release/20130502/integrated_call_samples_v3.20130502.ALL.panel",header=T,stringsAsFactors=F)
  fam         = read.table("./PCA/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_noHighLD_pruned.fam",header=F,stringsAsFactor=F)

  #Rename columns and add label for the population and super population
  colnames(fam)[1]<-c("sample")
  fam$super_pop<- NA
  fam$pop<- NA
  fam<- fam[,c(1,7,8)]

  samples1000<-samples1000[,c("sample","super_pop","pop")]
  fam<-rbind(fam,samples1000)

  #### CHANGE THIS #####
  #Read in PCA results (from Qians since his has finished running)
  pca<-read.table("/lustre/scratch123/hgi/mdt2/projects/ibdgwas_bioresource/qz2/pre_imputation_b04b06b15b19/PCA/AllBatch_1KG.eigenvec",head=F)
  colnames(pca)<-c("FID","IID",paste0("PC",c(1:20)))
  #Merge together the results with the fam file
  pca<-merge(pca,fam,by.x="IID",by.y="sample",all.x=T)

  #Find the variance explained by each PC
  eigenval         = read.table("/lustre/scratch123/hgi/mdt2/projects/ibdgwas_bioresource/qz2/pre_imputation_b04b06b15b19/PCA/AllBatch_1KG.eigenval",head=F)
  eigenval$var_exp = NA

  for (i in 1:nrow(eigenval)){
    eigenval$var_exp[i]<-(eigenval$V1[i] / sum(eigenval$V1))*100
  }

  #Plot the PCA clusters
  eigenval$PC<- 1:20
  eigenval$cumval<- cumsum(eigenval$var_exp)
  p1<-ggplot(data=eigenval,aes(x=PC,y=var_exp))+geom_bar(stat="identity")+theme_bw()+theme(axis.title=element_text(size=14,face="bold"),axis.text=element_text(size=14,face="bold"))+xlab("Principal Component (PC)")+ylab("Proportion of variance explained by each PC (%)")
  p2<-ggplot(data=eigenval,aes(x=PC,y=cumval))+geom_point(shape=16)+geom_line()+theme_bw()+theme(axis.title=element_text(size=14,face="bold"),axis.text=element_text(size=14,face="bold"))+xlab("Principal Component (PC)")+ylab("Cumulative proportion of variance explained by top PCs (%)")
  p3<-arrangeGrob(p1, p2,ncol=2, nrow = 1, widths = c(2.7, 3), heights = 2.3)
  ggsave("./figures/pca_var.png",dpi=300,width=13,height=6,p3)

  #Infer Populations
  pop<-c("EUR","AMR","AFR","EAS","SAS")
  pca$inferred_population<-NA
  for (i in 1:length(pop)) {
    pca$inferred_population[which(pca$PC1>=min(pca$PC1[which(pca$super_pop==pop[i])]) & pca$PC1<=max(pca$PC1[which(pca$super_pop==pop[i])]) &
                                    pca$PC2>=min(pca$PC2[which(pca$super_pop==pop[i])]) & pca$PC2<=max(pca$PC2[which(pca$super_pop==pop[i])]) &
                                    pca$PC3>=min(pca$PC3[which(pca$super_pop==pop[i])]) & pca$PC3<=max(pca$PC3[which(pca$super_pop==pop[i])]) &
                                    pca$PC4>=min(pca$PC4[which(pca$super_pop==pop[i])]) & pca$PC4<=max(pca$PC4[which(pca$super_pop==pop[i])]) &
                                    pca$PC5>=min(pca$PC5[which(pca$super_pop==pop[i])]) & pca$PC5<=max(pca$PC5[which(pca$super_pop==pop[i])]) &
                                    pca$PC6>=min(pca$PC6[which(pca$super_pop==pop[i])]) & pca$PC6<=max(pca$PC6[which(pca$super_pop==pop[i])]) &
                                    pca$PC7>=min(pca$PC7[which(pca$super_pop==pop[i])]) & pca$PC7<=max(pca$PC7[which(pca$super_pop==pop[i])]) &
                                    pca$PC8>=min(pca$PC8[which(pca$super_pop==pop[i])]) & pca$PC8<=max(pca$PC8[which(pca$super_pop==pop[i])]) &
                                    pca$PC9>=min(pca$PC9[which(pca$super_pop==pop[i])]) & pca$PC9<=max(pca$PC9[which(pca$super_pop==pop[i])]) &
                                    pca$PC10>=min(pca$PC10[which(pca$super_pop==pop[i])]) & pca$PC10<=max(pca$PC10[which(pca$super_pop==pop[i])]) )]<-pop[i]
  }

  # EXPAND LIMITS FOR EUR SUBSET OF SAMPLES
  pc1max<-max(pca$PC1[which(pca$super_pop=="EUR")])+((max(pca$PC1[which(pca$super_pop=="EUR")])-min(pca$PC1[which(pca$super_pop=="EUR")]))*0.3)
  pc1min<-min(pca$PC1[which(pca$super_pop=="EUR")])-((max(pca$PC1[which(pca$super_pop=="EUR")])-min(pca$PC1[which(pca$super_pop=="EUR")]))*0.3)
  pc2max<-max(pca$PC2[which(pca$super_pop=="EUR")])+((max(pca$PC2[which(pca$super_pop=="EUR")])-min(pca$PC2[which(pca$super_pop=="EUR")]))*0.3)
  pc2min<-min(pca$PC2[which(pca$super_pop=="EUR")])-((max(pca$PC2[which(pca$super_pop=="EUR")])-min(pca$PC2[which(pca$super_pop=="EUR")]))*0.3)
  pc3max<-max(pca$PC3[which(pca$super_pop=="EUR")])+((max(pca$PC3[which(pca$super_pop=="EUR")])-min(pca$PC3[which(pca$super_pop=="EUR")]))*0.3)
  pc3min<-min(pca$PC3[which(pca$super_pop=="EUR")])-((max(pca$PC3[which(pca$super_pop=="EUR")])-min(pca$PC3[which(pca$super_pop=="EUR")]))*0.3)

  pcaplot<- pca
  colnames(pcaplot)[23]<-"Population"
  pcaplot[is.na(pcaplot$Population),"Population"]<- "This study"

  p1<-ggplot(pcaplot,aes(x=PC1,y=PC2,col=Population))+geom_point()+theme_bw()+theme(axis.title=element_text(size=12,face="bold"),axis.text=element_text(size=12,face="bold"),legend.position="none")+geom_rect(xmin=pc1min,xmax=pc1max,ymin=pc2min,ymax=pc2max,col="red",fill = "transparent")
  p2<-ggplot(pcaplot,aes(x=PC2,y=PC3,col=Population))+geom_point()+theme_bw()+theme(axis.title=element_text(size=12,face="bold"),axis.text=element_text(size=12,face="bold"),legend.position="right",legend.title=element_text(size=12,face="bold"),legend.text=element_text(size=12,face="bold"))+geom_rect(xmin=pc2min,xmax=pc2max,ymin=pc3min,ymax=pc3max,col="red",fill = "transparent")

  grid.arrange(p1, p2,ncol=2, nrow = 1, widths = c(2.7, 3), heights = 2.7)

  p3<-arrangeGrob(p1, p2,ncol=2, nrow = 1, widths = c(2.7, 3), heights = 2.7)

  ggsave("./figures/PC123.png",dpi=300,width=13,height=6,p3)

  #Reclassify
  pca$inferred_population[which(is.na(pca$inferred_population) & pca$PC1<=pc1max & pca$PC1>=pc1min & pca$PC2<=pc2max & pca$PC2>=pc2min& pca$PC3<=pc3max & pca$PC3>=pc3min)]<-"EUR"

  # EXPAND LIMITS FOR SAS SUBSET OF SAMPLES
  pc1max_sas<-max(pca$PC1[which(pca$super_pop=="SAS")])+((max(pca$PC1[which(pca$super_pop=="SAS")])-min(pca$PC1[which(pca$super_pop=="SAS")]))*0.3)
  pc1min_sas<-min(pca$PC1[which(pca$super_pop=="SAS")])-((max(pca$PC1[which(pca$super_pop=="SAS")])-min(pca$PC1[which(pca$super_pop=="SAS")]))*0.3)
  pc2max_sas<-max(pca$PC2[which(pca$super_pop=="SAS")])+((max(pca$PC2[which(pca$super_pop=="SAS")])-min(pca$PC2[which(pca$super_pop=="SAS")]))*0.3)
  pc2min_sas<-min(pca$PC2[which(pca$super_pop=="SAS")])-((max(pca$PC2[which(pca$super_pop=="SAS")])-min(pca$PC2[which(pca$super_pop=="SAS")]))*0.3)
  pc3max_sas<-max(pca$PC3[which(pca$super_pop=="SAS")])+((max(pca$PC3[which(pca$super_pop=="SAS")])-min(pca$PC3[which(pca$super_pop=="SAS")]))*0.3)
  pc3min_sas<-min(pca$PC3[which(pca$super_pop=="SAS")])-((max(pca$PC3[which(pca$super_pop=="SAS")])-min(pca$PC3[which(pca$super_pop=="SAS")]))*0.3)

  pca$inferred_population[which(is.na(pca$inferred_population) & pca$PC1<=pc1max_sas & pca$PC1>=pc1min_sas & pca$PC2<=pc2max_sas & pca$PC2>=pc2min_sas
                                & pca$PC3<=pc3max_sas & pca$PC3>=pc3min_sas)]<-"SAS"

  # EXPAND LIMITS FOR EAS SUBSET OF SAMPLES
  pc1max_eas<-max(pca$PC1[which(pca$super_pop=="EAS")])+((max(pca$PC1[which(pca$super_pop=="EAS")])-min(pca$PC1[which(pca$super_pop=="EAS")]))*0.3)
  pc1min_eas<-min(pca$PC1[which(pca$super_pop=="EAS")])-((max(pca$PC1[which(pca$super_pop=="EAS")])-min(pca$PC1[which(pca$super_pop=="EAS")]))*0.3)
  pc2max_eas<-max(pca$PC2[which(pca$super_pop=="EAS")])+((max(pca$PC2[which(pca$super_pop=="EAS")])-min(pca$PC2[which(pca$super_pop=="EAS")]))*0.3)
  pc2min_eas<-min(pca$PC2[which(pca$super_pop=="EAS")])-((max(pca$PC2[which(pca$super_pop=="EAS")])-min(pca$PC2[which(pca$super_pop=="EAS")]))*0.3)
  pc3max_eas<-max(pca$PC3[which(pca$super_pop=="EAS")])+((max(pca$PC3[which(pca$super_pop=="EAS")])-min(pca$PC3[which(pca$super_pop=="EAS")]))*0.3)
  pc3min_eas<-min(pca$PC3[which(pca$super_pop=="EAS")])-((max(pca$PC3[which(pca$super_pop=="EAS")])-min(pca$PC3[which(pca$super_pop=="EAS")]))*0.3)

  pca$inferred_population[which(is.na(pca$inferred_population) & pca$PC1<=pc1max_eas & pca$PC1>=pc1min_eas & pca$PC2<=pc2max_eas & pca$PC2>=pc2min_eas
                                & pca$PC3<=pc3max_eas & pca$PC3>=pc3min_eas)]<-"EAS"

  # EXPAND LIMITS FOR AFR SUBSET OF SAMPLES
  pc1max_afr<-max(pca$PC1[which(pca$super_pop=="AFR")])+((max(pca$PC1[which(pca$super_pop=="AFR")])-min(pca$PC1[which(pca$super_pop=="AFR")]))*0.3)
  pc1min_afr<-min(pca$PC1[which(pca$super_pop=="AFR")])-((max(pca$PC1[which(pca$super_pop=="AFR")])-min(pca$PC1[which(pca$super_pop=="AFR")]))*0.3)
  pc2max_afr<-max(pca$PC2[which(pca$super_pop=="AFR")])+((max(pca$PC2[which(pca$super_pop=="AFR")])-min(pca$PC2[which(pca$super_pop=="AFR")]))*0.3)
  pc2min_afr<-min(pca$PC2[which(pca$super_pop=="AFR")])-((max(pca$PC2[which(pca$super_pop=="AFR")])-min(pca$PC2[which(pca$super_pop=="AFR")]))*0.3)
  pc3max_afr<-max(pca$PC3[which(pca$super_pop=="AFR")])+((max(pca$PC3[which(pca$super_pop=="AFR")])-min(pca$PC3[which(pca$super_pop=="AFR")]))*0.3)
  pc3min_afr<-min(pca$PC3[which(pca$super_pop=="AFR")])-((max(pca$PC3[which(pca$super_pop=="AFR")])-min(pca$PC3[which(pca$super_pop=="AFR")]))*0.3)

  pca$inferred_population[which(is.na(pca$inferred_population) & pca$PC1<=pc1max_afr & pca$PC1>=pc1min_afr & pca$PC2<=pc2max_afr & pca$PC2>=pc2min_afr
                                & pca$PC3<=pc3max_afr & pca$PC3>=pc3min_afr)]<-"AFR"

  # EXPAND LIMITS FOR AMR SUBSET OF SAMPLES
  pc1max_amr<-max(pca$PC1[which(pca$super_pop=="AMR")])+((max(pca$PC1[which(pca$super_pop=="AMR")])-min(pca$PC1[which(pca$super_pop=="AMR")]))*0.3)
  pc1min_amr<-min(pca$PC1[which(pca$super_pop=="AMR")])-((max(pca$PC1[which(pca$super_pop=="AMR")])-min(pca$PC1[which(pca$super_pop=="AMR")]))*0.3)
  pc2max_amr<-max(pca$PC2[which(pca$super_pop=="AMR")])+((max(pca$PC2[which(pca$super_pop=="AMR")])-min(pca$PC2[which(pca$super_pop=="AMR")]))*0.3)
  pc2min_amr<-min(pca$PC2[which(pca$super_pop=="AMR")])-((max(pca$PC2[which(pca$super_pop=="AMR")])-min(pca$PC2[which(pca$super_pop=="AMR")]))*0.3)
  pc3max_amr<-max(pca$PC3[which(pca$super_pop=="AMR")])+((max(pca$PC3[which(pca$super_pop=="AMR")])-min(pca$PC3[which(pca$super_pop=="AMR")]))*0.3)
  pc3min_amr<-min(pca$PC3[which(pca$super_pop=="AMR")])-((max(pca$PC3[which(pca$super_pop=="AMR")])-min(pca$PC3[which(pca$super_pop=="AMR")]))*0.3)

  pca$inferred_population[which(is.na(pca$inferred_population) & pca$PC1<=pc1max_amr & pca$PC1>=pc1min_amr & pca$PC2<=pc2max_amr & pca$PC2>=pc2min_amr
                                & pca$PC3<=pc3max_amr & pca$PC3>=pc3min_amr)]<-"AMR"


  pca$inferred_population<-as.factor(pca$inferred_population)
  pca$inferred_population<-factor(pca$inferred_population, levels=c("EUR","AFR","EAS","SAS","AMR"))

  fam        = read.table("./PCA/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_noHighLD_pruned.fam",header=F,stringsAsFactor=F)
  pca_nodup  = pca[which(pca$IID %in% fam$V2),]

  #Write out ancestries
  write.table(pca_nodup[which(pca_nodup$inferred_population=="EUR"),c("FID","IID")],"./PCA/list_eur_ancestry_samples_noDuplicates",col.names=F,row.names=F,quote=F,sep="\t")
  write.table(pca_nodup[which(pca_nodup$inferred_population=="AFR"),c("FID","IID")],"./PCA/list_afr_ancestry_samples_noDuplicates",col.names=F,row.names=F,quote=F,sep="\t")
  write.table(pca_nodup[which(pca_nodup$inferred_population=="EAS"),c("FID","IID")],"./PCA/list_eas_ancestry_samples_noDuplicates",col.names=F,row.names=F,quote=F,sep="\t")
  write.table(pca_nodup[which(pca_nodup$inferred_population=="SAS"),c("FID","IID")],"./PCA/list_sas_ancestry_samples_noDuplicates",col.names=F,row.names=F,quote=F,sep="\t")
  write.table(pca_nodup[which(pca_nodup$inferred_population=="AMR"),c("FID","IID")],"./PCA/list_amr_ancestry_samples_noDuplicates",col.names=F,row.names=F,quote=F,sep="\t")

  system("mkdir -p ./ancestry")
  ancestry = c("eur","amr","eas","afr","sas")
  for (i in 1:length(ancestry)) {
    system(sprintf("mkdir -p ./ancestry/%s", ancestry[i]))
    system(sprintf("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy --keep PCA/list_%s_ancestry_samples_noDuplicates --make-bed --out ./ancestry/%s/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_%s",
           ancestry[i],ancestry[i],ancestry[i]))
  }
}

HWE_filter = function(){
  #Calculate HWE (in EUR samples only)
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./ancestry/eur/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_eur --hardy --out ./ancestry/eur/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_eur")
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./ancestry/eur/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_eur --hardy --filter-females --chr 23 --out ./ancestry/eur/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_eur_females")

  #Check whether variants with significant HWE P-value also have significant GWAS P-value, and found most of such variants are in MHC region
  hwe = read.table("./ancestry/eur/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_eur.hwe",header=T,stringsAsFactors=F)

  smy = fread("/lustre/scratch123/hgi/mdt2/teams/anderson/qz2/scratch115/proj1/data_23_11_2020/lange_smy/28067908-GCST004131-EFO_0003767.h.tsv",sep="\t",stringsAsFactor=F)
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

  png("./figures/hwe_gwas5e-8.png", units="in", width=16, height=12, res=300)

  # get MHC region from https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38.p13, chr6:28,510,120-33,480,577
  ancestry<-"eur"

  ###ASK QIAN ABOUT THIS
  #tmp <- read.table(paste0(i,"/",AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_eur,"_",i,".hwe"),header=T,stringsAsFactors=F)
  #pos<- unlist(strsplit(tmp$SNP,split=":|_"))
  #pos<- pos[seq(2,length(pos),4)]
  #tmp$pos<- as.double(pos)
  #tmp<- tmp[!(tmp$CHR==6 & tmp$pos>=28510120 & tmp$pos<=33480577),]
  #tmp.sel<- tmp[ tmp$P<1e-12,]
  #hwe <- rbind(hwe,tmp.sel)
  #hwe.snp<- tmp.sel$SNP
  write.table(hwe,file="./list_var_exclude_hwe",col.names=F,quote=F,row.names=F)

  # 14.3 get SNPs not in MHC region but with HWE P-value < 1e-12 (1e-12 is used since they are cases)
  ## three parameters: 1. work_dir; 2. prefix of input file(s); 3. output file name (list of variants with extreme HWE P-value)
  #/software/R-4.1.0/bin/Rscript scripts/14_get_hwe.R ${wkdir}/data/ancestry ${filename} "list_var_exclude_hwe"

  hwe = read.table("./ancestry/eur/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_eur.hwe",header=T,stringsAsFactors=F)

  pos<- unlist(strsplit(hwe$SNP,split=":|_"))
  pos<- pos[seq(2,length(pos),4)]
  hwe$pos<- as.double(pos)
  hwe_auto<- hwe[!(hwe$CHR==6 & hwe$pos>=28510120 & hwe$pos<=33480577) & hwe$CHR!=23,]
  hwe_auto.sel<- hwe_auto[ hwe_auto$P<1e-12,]
  hwe_auto.snp<- hwe_auto.sel$SNP

  #ASK QIAN ABOUT THIS
  # hwe_sex.snp <- c()
  # filename <- "./ancestry/eur/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_eur_females.hwe"
  # if(file.exists(filename)){
  #   read.table(filename,header=T,stringsAsFactors=F)->hwe
  #   pos<- unlist(strsplit(hwe$SNP,split=":|_"))
  #   pos<- pos[seq(2,length(pos),4)]
  #   hwe$pos<- as.double(pos)
  #   hwe_sex.sel<- hwe_sex[ hwe_sex$P<1e-12,]
  #   hwe_sex.snp<- hwe_sex.sel$SNP
  # }

  hwe.snp <- c(hwe_auto.snp)
  write.table(hwe.snp,file="./ancestry/list_var_exclude_hwe",col.names=F,quote=F,row.names=F)

  #Remove SNPs did not pass HWE test
  system("grep -v -w -f ./snps_to_keep ./ancestry/list_var_exclude_hwe > ./ancestry/list_var_exclude_hwe_rmkeep")

  ancestry= c("eur","amr","eas","afr","sas")
  for (i in 1:length(ancestry)) {
    system(sprintf("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./ancestry/%s/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_%s --exclude ./ancestry/list_var_exclude_hwe_rmkeep --make-bed --out ./ancestry/%s/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_%s_hwe1e-12",
                   ancestry[i],ancestry[i],ancestry[i],ancestry[i]))
  }
}

Heterozygosity_filter = function(){

  #Loop through all ancestries
  ancestry = c("eur","afr","sas")
  for (i in 1:length(ancestry)) {
    system(sprintf("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./ancestry/%s/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_%s_hwe1e-12 --allow-no-sex --het --out ./ancestry/%s/het",
            ancestry[i],ancestry[i],ancestry[i]))

    het = read.table(sprintf("./ancestry/%s/het.het",ancestry[i]),header=T,stringsAsFactors=F)

    het$het<-(het$N.NM.-het$O.HOM)/het$N.NM.

    # #### CREATE AN HISTOGRAM

    # number of bins in histogram
    fd=function(x) {
      n=length(x)
      r=IQR(x)
        2*r/n^(1/3)
    }

    limits.lower<-mean(het$het)-4*sd(het$het)
    limits.upper<-mean(het$het)+4*sd(het$het)

    ymax <- max(table(cut(het$het,seq(min(het$het),max(het$het),dist(range(het$het))/100))))

    ggplot(data=het,aes(het))+geom_histogram(color="darkblue", fill="lightblue",binwidth = dist(range(het$het))/100)+ylab("Number of samples")+xlab("Heterozygosity rate")+theme_bw()+theme(axis.title=element_text(size=14,face="bold"),axis.text=element_text(size=14,face="bold"))+geom_vline(xintercept=limits.lower,lty=2,col="violetred")+geom_vline(xintercept=limits.upper,col="violetred",lty=2)+geom_vline(xintercept=mean(het$het),col="violetred",lty=2)+annotate(geom="text", x=limits.lower, y=ymax*0.95, label="-4SD",color="violetred",size=4,hjust=0)+annotate(geom="text", x=limits.upper, y=ymax*0.95, label="+4SD",color="violetred",size=4,hjust=0)+annotate(geom="text", x=mean(het$het), y=ymax*0.95, label="Average",color="violetred",size=4,hjust=0)
    ggsave(sprintf("./figures/%s_het.png",ancestry[i]),dpi=300,width=8,height=6)

    het.rm<- het[ (het$het > limits.upper)|(het$het < limits.lower),]
    write.table(het.rm[,c("FID","IID")],sprintf("./ancestry/%s/het_rm.list",ancestry[i]),col.names=F,quote=F,row.names=F)

    system(sprintf("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./ancestry/%s/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_%s_hwe1e-12 --allow-no-sex --remove ./ancestry/%s/het_rm.list --make-bed --out ./ancestry/%s/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_%s_hwe1e-12_het",
                   ancestry[i],ancestry[i],ancestry[i],ancestry[i],ancestry[i]))
  }
}

split_cohort = function(){
  # Split cohort in broad ancestry sets, Exclude monomorphic variants

  #Split into EUR
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./ancestry/eur/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_eur_hwe1e-12_het --keep-allele-order --allow-no-sex --freq --out ./ancestry/eur/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_eur_hwe1e-12.het")

  #Remove monomorphic variants
  system("awk '{if($5==0)print $2}' < ./ancestry/eur/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_eur_hwe1e-12.het.frq > ./ancestry/eur/study_list_variants_toexclude_monomorphic")
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./ancestry/eur/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_eur_hwe1e-12_het --allow-no-sex --exclude ./ancestry/eur/study_list_variants_toexclude_monomorphic --make-bed --out ./ancestry/eur/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_eur_hwe1e-12.het_nomonom")

  #Split into NON-EUR
  system("mkdir -p ./ancestry/non-eur")
  system("cat ./PCA/list_eur_ancestry_samples_noDuplicates ./ancestry/afr/het_rm.list ./ancestry/sas/het_rm.list > ./ancestry/non-eur/eur_sashet_afrhet_samples_noDuplicates")

  ###Generate non-EUR samples, including SAS(pass HET QC), AFR(pass HET QC), AMR, EAS and unknown
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy --allow-no-sex --exclude ./ancestry/list_var_exclude_hwe_rmkeep   --remove ./ancestry/non-eur/eur_sashet_afrhet_samples_noDuplicates --make-bed --out ./ancestry/non-eur/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_non-eur_hwe1e-12_het")

  #Remove monomorphic SNPs
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./ancestry/non-eur/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_non-eur_hwe1e-12_het --keep-allele-order --allow-no-sex --freq --out ./ancestry/non-eur/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_non-eur_hwe1e-12_het")
  system("awk '{if($5==0)print $2}' < ./ancestry/non-eur/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_non-eur_hwe1e-12_het.frq > ./ancestry/non-eur/study_list_variants_toexclude_monomorphic")
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./ancestry/non-eur/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_non-eur_hwe1e-12_het --allow-no-sex --exclude ./ancestry/non-eur/study_list_variants_toexclude_monomorphic --make-bed --out ./ancestry/non-eur/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_non-eur_hwe1e-12_het_nomonom")
}

prepare_sanger_imputation = function()
{
  # 17. Update to B37 (liftover) and double check with b37 fasta file, need to do this since sanger imputation requires b37
  system("mkdir -p ./submit_for_imputation")

  ## 17.1 EUR
  system("mkdir -p ./submit_for_imputation/eur")
  system("awk '{print \"chr\"$1,$4,$4+1,$2}' ./ancestry/eur/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_eur_hwe1e-12.het_nomonom.bim|sed 's/chr23/chrX/g' > ./ancestry/eur/study_eur_hg38_postqc.bed")
  system("/software/team152/liftover/liftOver ./ancestry/eur/study_eur_hg38_postqc.bed ./ref/hg38ToHg19.over.chain.gz ./submit_for_imputation/eur/study_hg38_postqc_lifted_hg19 ./submit_for_imputation/eur/study_hg38_postqc_no_lifted_hg19")

  system("cut -f 4 ./submit_for_imputation/eur/study_hg38_postqc_no_lifted_hg19 | sed \"/^#/d\" > ./submit_for_imputation/eur/nonlifted_hg19_variants_to_exclude_tmp.dat")

  system("grep \"alt|random\" ./submit_for_imputation/eur/study_hg38_postqc_lifted_hg19 | cut -f 4 | cat - ./submit_for_imputation/eur/nonlifted_hg19_variants_to_exclude_tmp.dat > ./submit_for_imputation/eur/nonlifted_hg19_variants_to_exclude.dat")

  ## 17.2  exclude not lifted variants and change position
  system("grep -w -v -f ./submit_for_imputation/eur/nonlifted_hg19_variants_to_exclude.dat ./submit_for_imputation/eur/study_hg38_postqc_lifted_hg19|awk '{print $4,$2}' > ./submit_for_imputation/eur/study_hg38_postqc_lifted_hg19_liftedvariants.map")

  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./ancestry/eur/AllBatch_ATCG_aligned_nodup_gender_nodupsample_flip_sexcheck_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95_perstudy_eur_hwe1e-12.het_nomonom --exclude ./submit_for_imputation/eur/nonlifted_hg19_variants_to_exclude.dat --update-map ./submit_for_imputation/eur/study_hg38_postqc_lifted_hg19_liftedvariants.map --make-bed --out ./submit_for_imputation/eur/study_postqc_lifted_hg19")

  ## 17.3 force A1 and A2 to be ref and alt alleles
  system("zcat ./b04/data/b04_ATCG_aligned.vcf.gz|cut -f '1-5'| awk '{print \"chr\"$1\":\"$2\"_\"$4\"_\"$5,$4}'|sed '/^chr##/d'  > ./submit_for_imputation/eur/list_variants_study_hg38_posstrandaligned_with_A1")

  system("sort ./submit_for_imputation/eur/list_variants_study_hg38_posstrandaligned_with_A1  | uniq > ./submit_for_imputation/eur/list_variants_study_hg38_posstrandaligned_with_A1_ed")

  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./submit_for_imputation/eur/study_postqc_lifted_hg19 --allow-no-sex --a2-allele ./submit_for_imputation/eur/list_variants_study_hg38_posstrandaligned_with_A1_ed --make-bed --out ./submit_for_imputation/eur/study_postqc_lifted_hg19_RefAlt")

  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./submit_for_imputation/eur/study_postqc_lifted_hg19_RefAlt --allow-no-sex --keep-allele-order --output-chr MT --recode vcf-fid --out ./submit_for_imputation/eur/study_postqc_lifted_hg19_RefAlt")

  ## 17.4 alignment (Long running time)
  system("/software/team152/bcftools-1.9/bcftools +/software/team152/bcftools-1.9/plugins/fixref.so  ./submit_for_imputation/eur/study_postqc_lifted_hg19_RefAlt.vcf -Oz -o ./submit_for_imputation/eur/study_postqc_lifted_hg19_RefAlt_posstrandaligned.vcf.gz -- -f ./ref/human_g1k_v37.fasta -m top 2>&1 | tee ./submit_for_imputation/eur/alignment_1.log")

  #### when use --check-ref x, it will exclude incorrect or missing REF allele is encountered, the unsolved variants from previous step would be removed
  #### it should be noted that INDEL would be normalized (left aligned) at this step, the position of indel therefore could be changed ####
  system("/software/team152/bcftools-1.9/bcftools norm --check-ref x -f ./ref/human_g1k_v37.fasta ./submit_for_imputation/eur/study_postqc_lifted_hg19_RefAlt_posstrandaligned.vcf.gz --threads 10 -Oz -o ./submit_for_imputation/eur/study_postqc_lifted_hg19_RefAlt_posstrandaligned_2.vcf.gz 2>&1 |tee ./submit_for_imputation/eur/alignment_2.log")

  ## 17.5 change variant ID based on the lifted variant position
  system("/software/team152/plink_linux_x86_64_20181202/plink --vcf ./submit_for_imputation/eur/study_postqc_lifted_hg19_RefAlt_posstrandaligned_2.vcf.gz --keep-allele-order --allow-no-sex --double-id --make-bed --out ./submit_for_imputation/eur/study_postqc_lifted_hg19_RefAlt_posstrandaligned")
  system("zcat ./submit_for_imputation/eur/study_postqc_lifted_hg19_RefAlt_posstrandaligned_2.vcf.gz|grep -v ^\"#\" | cut -f '1-5' | awk '{print $1\":\"$2\"_\"$4\"_\"$5,$3}'|sed 's/^X/23/g' > ./submit_for_imputation/eur/list_variants_study_hg19_posstrandaligned")

  system("mkdir -p ./submit_for_imputation/ready_for_imp_hg19")
  system("mkdir submit_for_imputation/ready_for_imp_hg19/eur")
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./submit_for_imputation/eur/study_postqc_lifted_hg19_RefAlt_posstrandaligned --update-name submit_for_imputation/eur/list_variants_study_hg19_posstrandaligned 1 2 --keep-allele-order --make-bed --out ./submit_for_imputation/ready_for_imp_hg19/eur/study_hg19")
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./submit_for_imputation/ready_for_imp_hg19/eur/study_hg19 --allow-no-sex --keep-allele-order --output-chr MT --recode vcf-fid --out ./submit_for_imputation/ready_for_imp_hg19/eur/study_hg19")

  ## 17.6 check if there is any error before submit for imputation
  system("/software/team152/bcftools-1.9/bcftools norm -ce -f ./ref/human_g1k_v37.fasta ./submit_for_imputation/ready_for_imp_hg19/eur/study_hg19.vcf -Ou -o /dev/null")
}
