#' @import stringr

prepare_reference_data = function()
{
  #Create directory for the reference data
  #system(sprintf("mkdir %s/ref"))
  #Download the raw 1K raw file
  #system("wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/")
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
  #system("rm ref/1KG/chr23_tmp")

  # #get high LD region
  #system("wget https://raw.githubusercontent.com/meyer-lab-cshl/plinkQC/master/inst/extdata/high-LD-regions-hg38-GRCh38.txt -O ref/high-LD-regions-hg38-GRCh38.txt")
  #system("sed -i 's/chr//g' ref/high-LD-regions-hg38-GRCh38.txt")
}

merge_batches_ibdbr = function(){
  system("mkdir -p ./data")
  system("rm -f ./data/old_batch_merge.list")
  system("rm -f ./data/new_batch_merge.list")
  system("rm -f ./data/merge.list")

  rawDataDir1 = "/lustre/scratch123/hgi/projects/ibdgwas_bioresource/genotyped_data/raw/data_transfer_20220228"
  rawDataDir2 = "/lustre/scratch123/hgi/projects/ibdgwas_bioresource/genotyped_data/raw/data_transfer_20220324"
  batch_old   = c("b04","b06","b08","b09","b10","b12","b15","b17","b18","b19")

  # Loop through all the previous batches and add their names to a list of data we wish to merge
  for (i in 1:length(batch_old)) {
    system(sprintf("echo %s/v2chip_%s_V3_calls.vcf.gz >> ./data/old_batch_merge.list", rawDataDir1, batch_old[i]))
    }

  # Add the new batches on from the new directory
  batch_new = c("b20")
  for (i in 1:length(batch_new)) {
    system(sprintf("echo %s/v2chip_%s_V3_calls.vcf.gz >> ./data/new_batch_merge.list", rawDataDir2, batch_new[i]))
  }

  #Merge the two files
  system("cat ./data/new_batch_merge.list ./data/old_batch_merge.list > ./data/merge.list")

  #Merge the VCF files together
  system("/software/team152/bcftools-1.9/bcftools merge -m id --file-list ./data/merge.list --force-samples -Oz -o ./data/AllBatch.vcf.gz")

  #Convert the VCF to PLINK format
  system("/software/team152/plink_linux_x86_64_20181202/plink --vcf ./data/AllBatch.vcf.gz --allow-extra-chr --output-chr MT --make-bed --out ./data/AllBatch_tmp")
}

update_sample_ids = function(){
  system('/software/R-4.1.0/bin/Rscript ./scripts/0_update_sample_ID.R ./data/AllBatch_tmp "sample_id_update"
  /software/team152/plink_linux_x86_64_20181202/plink --bfile ./data/AllBatch_tmp --allow-extra-chr --output-chr MT --update-ids ./data/sample_id_update --make-bed --out ./data/AllBatch')
}

keep_specific_snps = function(){
  system("rm -f ./data/snps_to_keep")
  system("echo 'chr16:50729870_C_CC' >> './data/snps_to_keep'")
  system("echo 'chr16:50729870_CC_C' >> './data/snps_to_keep'")
  system("echo 'AX-96079897' >> './data/snps_to_keep'")
}

#Keep only ATCG SNPs
keep_atcg_snps = function()
{
  system("/software/R-4.1.0/bin/Rscript scripts/1_rm_non_ATCG.R  ./data/AllBatch
  grep -v -w -f ./data/snps_to_keep data/snps_to_keep > ./data/snps_to_keep
  /software/team152/plink_linux_x86_64_20181202/plink --bfile ./data/AllBatch --allow-extra-chr --output-chr MT --exclude ./data/list_indel_var_exclude_rmkeep --make-bed --out ./data/AllBatch_ATCG")
}

#Align SNPs to the positive strand
align_positive_strand = function(){
  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./data/AllBatch_ATCG --allow-extra-chr --output-chr MT --allow-no-sex --recode vcf --out ./data/AllBatch_ATCG")
  system("/software/team152/bcftools-1.9/bcftools +/software/team152/bcftools-1.9/plugins/fixref.so ./data/AllBatch_ATCG.vcf -Oz -o ./data/AllBatch_ATCG_aligned.vcf.gz -- -f /lustre/scratch123/hgi/projects/ibdgwas/IIBDGC/resources/hg38/hg38_edited.fa -m top 2>&1 | tee ./data/alignment.log")
  # after alignment, format vcf to bed file
  system("/software/team152/plink_linux_x86_64_20181202/plink --vcf ./data/AllBatch_ATCG_aligned.vcf.gz --keep-allele-order --allow-extra-chr --id-delim --make-bed --out ./data/AllBatch_ATCG_aligned")
  system("rm ./data/AllBatch_ATCG.vcf")
}

#Update the variant IDs to chrX:XXX nomenclature
update_variant_id = function(){
  system("/software/team152/plink_linux_x86_64_20181202/plink --allow-extra-chr --bfile ./data/AllBatch_ATCG_aligned --missing --out ./data/AllBatch_ATCG_aligned")
  system('zcat ./data/AllBatch_ATCG_aligned.vcf.gz|grep -v ^"#"|cut -f \'1-5\' | awk \'{print $3,"chr"$1":"$2"_"$4"_"$5}\' > ./data/AllBatch_ATCG_aligned')

  #For duplicated variants(same chr, pos, ref, alt), remove the one with more missingness, for variants with same missingness, keep the first one
  system("/software/R-4.1.0/bin/Rscript ./scripts/5_rm_dup_var.R ./data/AllBatch_ATCG_aligned")
  system("grep -v -w -f ./data/snps_to_keep ./data/list_duplicated_var_exclude > ./data/list_duplicated_var_exclude_rmkeep")
  system("/software/team152/plink_linux_x86_64_20181202/plink --allow-extra-chr --bed ./data/AllBatch_ATCG_aligned.bed --bim ./data/AllBatch_ATCG_aligned_edited.bim --fam ./data/AllBatch_ATCG_aligned.fam --exclude ./data/list_duplicated_var_exclude_rmkeep --keep-allele-order --make-bed --out ./data/AllBatch_ATCG_aligned_nodup")
}

add_gender =  function(){
  system('/software/R-4.1.0/bin/Rscript scripts/5_add_gender.R ./data/AllBatch_ATCG_aligned_nodup "/lustre/scratch123/hgi/mdt2/projects/ibdgwas_bioresource/phenotype_data/release_20220404/raw/IBD_BioRes_phenotypes_20220404.txt"')
  system("/software/team152/plink_linux_x86_64_20181202/plink --allow-extra-chr --bed ./data/AllBatch_ATCG_aligned_nodup.bed --bim ./data/AllBatch_ATCG_aligned_nodup.bim --fam ./data/AllBatch_ATCG_aligned_nodup_gender --keep-allele-order --make-bed --out ./data/AllBatch_ATCG_aligned_nodup_gender")
}

compare_freq_1000G_EUR = function(){
  system("mkdir -p 1KG_EUR_study_variant")
  system("cat ./data/AllBatch_ATCG_aligned_nodup_gender.bim | cut -f 2 > ./1KG_EUR_study_variant/list_variants_posstr_nodup")

  system("cd 1KG_EUR_study_variant")
  system("awk \'{if($3==\"EUR\")print 0,$1}\' < ./ref/integrated_call_samples_v3.20130502.ALL.panel > ./1KG_EUR_study_variant/eur.sample")

  for (i in 1:23) {
    system(sprintf("/software/team152/plink2 --bfile ./ref/1KG/chr%s --keep ./1KG_EUR_study_variant/eur.sample --rm-dup force-first --extract ./1KG_EUR_study_variant/list_variants_posstr_nodup --make-bed --out ./1KG_EUR_study_variant/chr%s",i,i))
  }

  system("rm -f ./data/merge.list")

  for (i in 2:23) {
    system(sprintf("echo ./1KG_EUR_study_variant/chr%s.bed ./1KG_EUR_study_variant/chr%s.bim ./1KG_EUR_study_variant/chr%s.fam >> ./data/merge.list", i,i,i))
  }

  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./1KG_EUR_study_variant/chr1 --merge-list ./data/merge.list --allow-no-sex --make-bed --out ./1KG_EUR_study_variant/1000GP_EUR_b38_study_variants")

  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./1KG_EUR_study_variant/1000GP_EUR_b38_study_variants --freq --out ./1KG_EUR_study_variant/1000GP_EUR_b38_study_variants")
  system("/software/team152/plink_linux_x86_64_20181202/plink --allow-extra-chr --bfile ./data/AllBatch_ATCG_aligned_nodup_gender --freq --out ./data/AllBatch_ATCG_aligned_nodup_gender")

  system('/software/R-4.1.0/bin/Rscript ./scripts/6_maf_comp.R "./1KG_EUR_study_variant/1000GP_EUR_b38_study_variants.frq" "./data/AllBatch_ATCG_aligned_nodup_gender.frq"')
  system("grep -v -w -f ./data/snps_to_keep ./data/list_variants_to_remove_AT_CG > ./data/list_variants_to_remove_AT_CG_rmkeep")
  system("/software/team152/plink_linux_x86_64_20181202/plink  --allow-extra-chr --bfile ./data/AllBatch_ATCG_aligned_nodup_gender --exclude ./data/list_variants_to_remove_AT_CG_rmkeep --flip ./data/list_variants_to_flip_AT_CG --make-bed --out ./data/AllBatch_ATCG_aligned_nodup_gender_flip")
}

sex_discrepancies_check = function(){
  system("mkdir -p ./data/sex")
  system("mkdir -p ./data/sex/X")
  system("mkdir -p ./data/sex/Y")

  system("awk \'{if($5==1)print $1,$2}\' < ./data/AllBatch_ATCG_aligned_nodup_gender_flip.fam > ./data/sex/list_male_samples")
  system("awk \'{if($5==2)print $1,$2}\' < ./data/AllBatch_ATCG_aligned_nodup_gender_flip.fam > ./data/sex/list_female_samples")

  ## X chromosome
  system("/software/team152/plink_linux_x86_64_20181202/plink --allow-extra-chr --bfile ./data/AllBatch_ATCG_aligned_nodup_gender_flip --keep ./data/sex/list_female_samples --chr 23 --make-bed --out ./data/sex/X/study_check_sex_females_only")

  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./data/sex/X/study_check_sex_females_only --missing --hardy --freq --out ./data/sex/X/study_check_sex_females_only")

  ### get high quality variants
  ## three parameters: 1. work_dir; 2. hwe threshold; 3. prefix of file input; 4.output file name
  system('/software/R-4.1.0/bin/Rscript ./scripts/7_X_var.R 1e-12 "./data/sex/X/study_check_sex_females_only" "./data/sex/X/list_good_chrX_variants"')

  #Y chromosome
  system("/software/team152/plink_linux_x86_64_20181202/plink --allow-extra-chr --bfile ./data/AllBatch_ATCG_aligned_nodup_gender_flip --keep ./data/sex/list_male_samples --chr 24 --make-bed --out ./data/sex/Y/study_check_sex_males_only")

  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./data/sex/Y/study_check_sex_males_only --freq --missing --out ./data/sex/Y/study_check_sex_males_only")

  system("/software/team152/plink_linux_x86_64_20181202/plink --allow-extra-chr --bfile ./data/AllBatch_ATCG_aligned_nodup_gender_flip --keep ./data/sex/list_female_samples --chr 24 --recode A --out ./data/sex/Y/study_check_sex_females_only_chr24")

  ### get high quality variants
  ## three parameters: 1. work_dir; 2. prefix of file input (frq/lmiss); 3. dosage file name; 4.output file name
  system('/software/R-4.1.0/bin/Rscript ./scripts/7_Y_var.R  "./data/sex/Y/study_check_sex_males_only" "./data/sex/Y/study_check_sex_females_only_chr24.raw" "./data/sex/Y/list_good_chrY_variants"')

  system("cat ./data/sex/X/list_good_chrX_variants ./data/sex/Y/list_good_chrY_variants > ./data/sex/list_good_chrXY_variants")

  system("/software/team152/plink_linux_x86_64_20181202/plink --allow-extra-chr --bfile ./data/AllBatch_ATCG_aligned_nodup_gender_flip --extract ./data/sex/list_good_chrXY_variants --make-bed --out ./data/sex/study_check_sex_good_chrXY")

  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./data/sex/study_check_sex_good_chrXY --check-sex ycount --out ./data/sex/study_check_sex_good_chrXY")

  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./data/sex/study_check_sex_good_chrXY --check-sex --out ./data/sex/study_check_sex_good_chrX")

  system("mkdir figures")
  ## three parameters: 1. work_dir; 2. filename of good XY variants; 3. filename of good X variants; 4.prefix of file
  system('/software/R-4.1.0/bin/Rscript scripts/7_check_sex.R "./data/sex/study_check_sex_good_chrXY.sexcheck" "./data/sex/study_check_sex_good_chrX.sexcheck" ./data/AllBatch_ATCG_aligned_nodup_gender_flip')

  system("/software/team152/plink_linux_x86_64_20181202/plink --allow-extra-chr --bed ./data/AllBatch_ATCG_aligned_nodup_gender_flip.bed --bim ./data/AllBatch_ATCG_aligned_nodup_gender_flip.bim --fam ./data/AllBatch_ATCG_aligned_nodup_gender_flip_edited.fam --allow-no-sex --remove ./data/sex/list_samples_wrong_gender --make-bed --out ./data/AllBatch_ATCG_aligned_nodup_gender_flip_sexcheck")

  system("/software/team152/plink_linux_x86_64_20181202/plink --allow-extra-chr --bfile data/AllBatch_ATCG_aligned_nodup_gender_flip --allow-no-sex --set-hh-missing --make-bed --out ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh")
}

#Remove Y chromosome
remove_Y_chr = function(){
  system("awk '{if($1==24)print $2}' < ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh.bim > ./data/list_chry_variants_toremove")

  system("/software/team152/plink_linux_x86_64_20181202/plink --allow-extra-chr --bfile ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh --allow-no-sex --exclude ./data/list_chry_variants_toremove --make-bed --out ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry")
}

remove_low_call_rate = function(){

  # Sample Call Rate < 80%
  system("/software/team152/plink_linux_x86_64_20181202/plink --allow-extra-chr --bfile ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry --allow-no-sex --mind 0.20 --make-bed --out ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8")

  # Variant Call Rate < 80%
  system("/software/team152/plink_linux_x86_64_20181202/plink --allow-extra-chr --bfile ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8 --allow-no-sex --missing --out ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8")
  system("awk '{if($5>0.2) print $2}' ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8.lmiss > ./data/list_variants_to_remove_vcr0.8")
  system("grep -v -w -f ./data/snps_to_keep ./data/list_variants_to_remove_vcr0.8 > ./data/list_variants_to_remove_vcr0.8_rmkeep")
  system("/software/team152/plink_linux_x86_64_20181202/plink --allow-extra-chr --bfile ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8 --allow-no-sex --exclude ./data/list_variants_to_remove_vcr0.8_rmkeep --make-bed --out ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8")

  # Sample Call Rate < 95%
  system("/software/team152/plink_linux_x86_64_20181202/plink --allow-extra-chr --bfile ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8 --allow-no-sex --mind 0.05 --make-bed --out ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95")

  # Variant Call Rate < 95%
  system("/software/team152/plink_linux_x86_64_20181202/plink --allow-extra-chr --bfile ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95 --allow-no-sex --missing --out ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95")
  system("awk '{if($5>0.05) print $2}' ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95.lmiss > ./data/list_variants_to_remove_vcr0.95")
  system("grep -v -w -f ./data/snps_to_keep data/list_variants_to_remove_vcr0.95 > ./data/list_variants_to_remove_vcr0.95_rmkeep")
  system("/software/team152/plink_linux_x86_64_20181202/plink --allow-extra-chr --bfile ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95 --allow-no-sex --exclude ./data/list_variants_to_remove_vcr0.95_rmkeep --make-bed --out ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.95")

  # Variant MAF<0.01 AND Call Rate <98%
  system("/software/team152/plink_linux_x86_64_20181202/plink --allow-extra-chr --bfile ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95 --allow-no-sex --missing --out ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95")
  system("/software/team152/plink_linux_x86_64_20181202/plink --allow-extra-chr --bfile ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95 --allow-no-sex --freq --out ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95")

  ## three parameters: 1. work_dir; 2. prefix of input files; 3. output filename
  system('/software/R-4.1.0/bin/Rscript ./scripts/9_callrate.R ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95  "./data/list_monomorphic_vcr098maf0.01_var_exclude"')
  system("grep -v -w -f ./data/snps_to_keep ./data/list_monomorphic_vcr098maf0.01_var_exclude > ./data/list_monomorphic_vcr098maf0.01_var_exclude_rmkeep")
  system("/software/team152/plink_linux_x86_64_20181202/plink --allow-extra-chr --bfile ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95 --allow-no-sex --exclude ./data/list_monomorphic_vcr098maf0.01_var_exclude_rmkeep --make-bed --out ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.01")
}

variant_missingness = function(){
  ## I did not do it across studies, but do it across batches (only five main batches including b04,b06,b15,b19,b20 were considered)
  batch_old=c("b04", "b06", "b08","b09","b10", "b12", "b15", "b17", "b18", "b19")
  system("rm -f ./data/id_batch ./data/id_old_batch ./data/id_new_batch")

  rawDataDir1="/lustre/scratch123/hgi/projects/ibdgwas_bioresource/genotyped_data/raw/data_transfer_20220228"
  rawDataDir2="/lustre/scratch123/hgi/projects/ibdgwas_bioresource/genotyped_data/raw/data_transfer_20220324"

  for (i in 1:length(batch_old)) {
    system(sprintf('/software/team152/bcftools-1.9/bcftools query -l %s/v2chip_%s_V3_calls.vcf.gz|awk -v bnum="%s" \'{print $0,bnum}\' >> ./data/id_old_batch', rawDataDir1,batch_old[i],batch_old[i]))
  }

  batch_new=c("b20")
  for (i in 1:length(batch_new)) {
    system(sprintf('/software/team152/bcftools-1.9/bcftools query -l %s/v2chip_%s_V3_calls.vcf.gz|awk -v bnum="%s" \'{print $0,bnum}\' >> data/id_new_batch',rawDataDir2,batch_new,batch_new))
  }

  system('/software/R-4.1.0/bin/Rscript scripts/11_sample_batch.R "./data/id_old_batch" "./data/id_new_batch"')

  batch_used=c("b04","b06","b15","b19","b20")
  for (i in 1:length(batch_used)) {
    system(sprintf("/software/team152/plink_linux_x86_64_20181202/plink --allow-extra-chr --bfile ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.01 --missing --keep ./data/%s.sample --out ./data/tmp_%s",batch_used[i],batch_used[i]))
  }

  system('/software/R-4.1.0/bin/Rscript scripts/11_missing_per_study.R "b04 b06 b15 b19 b20" "./data/table_breakdown_number_missing_variants_perstudy.csv" "./data/list_variants_per_study_missing_0.1"')
  system("grep -v -w -f ./data/snps_to_keep ./data/list_variants_per_study_missing_0.1 > ./data/list_variants_per_study_missing_0.1_rmkeep")
  system("/software/team152/plink_linux_x86_64_20181202/plink --allow-extra-chr --bfile ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.01 --allow-no-sex --exclude ./data/list_variants_per_study_missing_0.1_rmkeep --make-bed --out ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.01_perstudy")

  system("rm -f ./data/tmp_*")
}

duplicate_samples = function(){

  system("/software/team152/plink_linux_x86_64_20181202/plink --allow-extra-chr --bfile ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.01_perstudy --missing --out ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.01_perstudy")

  system("/software/team152/king -b ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.01_perstudy.bed --related --cpus 1 --prefix ./data/study_king")

  #Four parameters: 1. work_dir; 2. prefix of file name of king's result (within and across family); 3. prefix of file(s); 4. output file name
  system('/software/R-4.1.0/bin/Rscript ./scripts/12_rm_dup_ind.R  "./data/study_king" ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.01_perstudy "./data/list_duplicated_samples"')

  # Extract non-dupulicated samples
  system("/software/team152/plink_linux_x86_64_20181202/plink --allow-extra-chr --bfile ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.01_perstudy --remove ./data/list_duplicated_samples --make-bed --out ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.01_perstudy_nodupsample")
}

PCA = function(){
  system("mkdir -p ./PCA")

  #Remove variants in regions with high LD (using non-duplicated samples to calculate LD)
  system("/software/team152/plink_linux_x86_64_20181202/plink --allow-extra-chr --bfile ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.01_perstudy_nodupsample --exclude /lustre/scratch123/hgi/projects/ibdgwas_bioresource/qz2/ref/high-LD-regions-hg38-GRCh38.txt --range --make-bed --out ./PCA/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.01_perstudy_nodupsample_noHighLD")

  #Keep independent SNPs, make sure these SNPs also in 1KG
  system("awk '{print $2}' < ./1KG_EUR_study_variant/1000GP_EUR_b38_study_variants.bim > ./PCA/1KG.snp")
  system("/software/team152/plink_linux_x86_64_20181202/plink --allow-extra-chr --bfile ./PCA/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.01_perstudy_nodupsample_noHighLD --indep-pairwise 50 5 0.2 --extract ./PCA/1KG.snp --out ./PCA/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.01_perstudy_nodupsample_noHighLD")

  system("grep -v \"^chrX\" ./PCA/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.01_perstudy_nodupsample_noHighLD.prune.in > ./PCA/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.01_perstudy_nodupsample_noHighLD_noHighLD_nosex.prune.in")

  system("/software/team152/plink_linux_x86_64_20181202/plink --allow-extra-chr --bfile ./PCA/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.01_perstudy_nodupsample_noHighLD --extract ./PCA/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.01_perstudy_nodupsample_noHighLD.prune.in --make-bed --out ./PCA/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.01_perstudy_nodupsample_noHighLD_pruned")

  #Get 1KG with study variants
  system("mkdir -p ./PCA/1KG")

  ###get 2504 unrelated 1KG sample IDs, additional 698 samples that are related to 2504 samples are not considered. see details here https://www.internationalgenome.org/data-portal/data-collection/30x-grch38
  system("awk '{if(NR!=1)print 0,$1}' ./ref/integrated_call_samples_v3.20130502.ALL.panel > ./PCA/1KG/1KG_urel.id")

  for (i in 1:23) {
    system(sprintf("/software/team152/plink2 --bfile ./ref/1KG/chr%s --keep ./PCA/1KG/1KG_urel.id --extract ./PCA/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.01_perstudy_nodupsample_noHighLD.prune.in --make-bed --out ./PCA/1KG/chr%s",i,i))
  }

  system("rm ./PCA/1KG/merge.list")

  for (i in 2:23) {
    system(sprintf('echo "./PCA/1KG/chr%s.bed ./PCA/1KG/chr%s.bim ./PCA/1KG/chr%s.fam" >> ./PCA/1KG/merge.list',i,i,i))
  }

  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./PCA/1KG/chr1 --merge-list ./PCA/1KG/merge.list --allow-no-sex --make-bed --out ./PCA/1KG/1000GP_ALL_b38_study_variants")

  system("/software/team152/plink_linux_x86_64_20181202/plink --allow-extra-chr --bfile ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.01_perstudy --extract ./PCA/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.01_perstudy_nodupsample_noHighLD.prune.in --make-bed --out ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.01_perstudy_noHighLD_pruned")

  system("/software/team152/plink_linux_x86_64_20181202/plink --bfile ./PCA/1KG/1000GP_ALL_b38_study_variants --bmerge ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.01_perstudy_noHighLD_pruned.bed ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.01_perstudy_noHighLD_pruned.bim ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.01_perstudy_noHighLD_pruned.fam --make-bed --out ./PCA/1KG/1000GP_ALL_b38_plus_study_dup")

  #PCA

  ## Run this on the Cluster!
  #bsub -J ./PCA/1KG/PCA_AllBatch_1KG_dup -q normal -o ./PCA/1KG/PCA_AllBatch_1KG_dup.log -e ./PCA/1KG/PCA_AllBatch_1KG_dup.err -W 10:00 -m "modern_hardware" -M 40000 -R'select[mem>40000] rusage[mem=40000]' /software/team152/plink2 --memory 39000 --bfile ./PCA/1KG/1000GP_ALL_b38_plus_study_dup --pca --out ./PCA/1KG/AllBatch_1KG_dup

  #Get study samples with different ancestry
  system('/software/R-4.1.0/bin/Rscript ./scripts/13_keep_all_ancestry.R "./ref/integrated_call_samples_v3.20130502.ALL.panel" "./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.01_perstudy_noHighLD_pruned" "./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.01_perstudy_noHighLD_pruned" "/lustre/scratch123/hgi/projects/ibdgwas_bioresource/qz2/pre_imputation_b04b06b15b19b20/PCA/AllBatch_1KG_dup"')

  system("mkdir -p ./data/ancestry")
  ancestry=c("eur","amr","eas","afr","sas")

  for (i in 1:length(ancestry)) {
    system(sprintf("mkdir -p ./data/ancestry/%s",ancestry[i]))
    system(sprintf("/software/team152/plink_linux_x86_64_20181202/plink --allow-extra-chr --bfile ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.01_perstudy --keep ./PCA/list_%s_ancestry_samples_withDuplicates --make-bed --out data/ancestry/%s/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.01_perstudy_noHighLD_pruned_%s",ancestry[i],ancestry[i],ancestry[i]))
    system(sprintf("/software/team152/plink_linux_x86_64_20181202/plink --allow-extra-chr --bfile ./data/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.01_perstudy --keep PCA/list_%s_ancestry_samples_noDuplicates --make-bed --out data/ancestry/%s/AllBatch_ATCG_aligned_nodup_gender_flip_missinghh_nochry_scr0.8_vcr0.8_scr0.95_vcr0.98maf0.01_perstudy_%s", ancestry[i], ancestry[i], ancestry[i]))
  }
}
