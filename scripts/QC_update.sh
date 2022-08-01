#!/bin/bash
#This document follows the pipeline shared by Laura, 18-01-2021: https://drive.google.com/file/d/1RmDcq90uPAWZMZaKy3l5Vy3JtG3jbwq5/view?usp=sharing
#A few changes have been made to fit this data.
#Qian Zhang, 01-03-2022

#cd /lustre/scratch123/hgi/projects/ibdgwas_bioresource/qz2/
#
## prepare ref data for the following analysis
### 1KG 30x on GRCh38, description on here:https://www.internationalgenome.org/data-portal/data-collection/30x-grch38
#mkdir ref
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

# data preparation
# 0. merge data and transfer from VCF to plink
wkdir="/lustre/scratch123/hgi/projects/ibdgwas_bioresource/qz2/pre_imputation_b04b06b15b19b20"

cd ${wkdir}
mkdir data
rawDataDir1="/lustre/scratch123/hgi/projects/ibdgwas_bioresource/genotyped_data/raw/data_transfer_20220228"
rawDataDir2="/lustre/scratch123/hgi/projects/ibdgwas_bioresource/genotyped_data/raw/data_transfer_20220324"
batch_old=(b04 b06 b08 b09 b10 b12 b15 b17 b18 b19)
rm data/merge.list
for i in ${batch_old[@]}
do
	echo ${rawDataDir1}/v2chip_${i}_V3_calls.vcf.gz >> data/merge.list
done

batch_new="b20"
for i in ${batch_new[@]}
do
	echo ${rawDataDir2}/v2chip_${i}_V3_calls.vcf.gz >> data/merge.list
done

filename="AllBatch"

bcftools merge -m id --file-list data/merge.list --force-samples -Oz -o data/${filename}.vcf.gz

plink --vcf data/${filename}.vcf.gz --allow-extra-chr --output-chr MT --make-bed --out data/${filename}_tmp

## update sample ID (this is for UKIBDBR only, for samples with same FID, change their IID as FID-1/2/3)
#three parameters: 1. work_dir; 2. prefix of file; 3. output file
/software/R-4.1.0/bin/Rscript scripts/0_update_sample_ID.R ${wkdir} ${filename}_tmp "sample_id_update"
plink --bfile data/${filename}_tmp --allow-extra-chr --output-chr MT --update-ids data/sample_id_update --make-bed --out data/${filename}

rm data/${filename}_tmp*

rm data/snps_to_keep
#same SNP, just in different formats
echo "chr16:50729870_C_CC" >> data/snps_to_keep
echo "chr16:50729870_CC_C" >> data/snps_to_keep
echo "AX-96079897" >> data/snps_to_keep

# 1. Remove Indels, Mitochondrial Variants (None ATCG SNPs)
#two parameters: 1. work_dir; 2. prefix of file
/software/R-4.1.0/bin/Rscript scripts/1_rm_non_ATCG.R ${wkdir} ${filename}
grep -v -w -f data/snps_to_keep data/list_indel_var_exclude > data/list_indel_var_exclude_rmkeep
plink --bfile data/${filename} --allow-extra-chr --output-chr MT --exclude data/list_indel_var_exclude_rmkeep --make-bed --out data/${filename}_ATCG

filename=${filename}_ATCG

# 2. remove samples with no phenotype data or samples that have withdrawn consent (did not perform here, applied after data merge)

# 3. Update to B37 (not done here since our data is hg38)

# 4. Align to (+) strand
plink --bfile data/${filename} --allow-extra-chr --output-chr MT --allow-no-sex --recode vcf --out data/${filename}
bcftools +fixref data/${filename}.vcf -Oz -o data/${filename}_aligned.vcf.gz -- -f /lustre/scratch123/hgi/projects/ibdgwas/IIBDGC/resources/hg38/hg38_edited.fa -m top 2>&1 | tee data/alignment.log
# after alignment, format vcf to bed file
plink --vcf data/${filename}_aligned.vcf.gz --keep-allele-order --id-delim --make-bed --out data/${filename}_aligned
rm data/${filename}.vcf

filename=${filename}_aligned

# 5. Update variant ids (using chr:position_ref_alt); Remove duplicated variants (same Chr pos ref alt, remove ones with more missingness)
# 5.1. update varaint id and remove duplicated variants
plink --bfile data/${filename} --missing --out data/${filename}
zcat data/${filename}.vcf.gz|grep -v ^"#"|cut -f '1-5' | awk '{print $3,"chr"$1":"$2"_"$4"_"$5}' > data/${filename}

# for duplicated variants(same chr, pos, ref, alt), remove the one with more missingness, for variants with same missingness, keep the first one
# two parameters: 1. work_dir; 2. prefix of file
/software/R-4.1.0/bin/Rscript scripts/5_rm_dup_var.R ${wkdir} ${filename}
grep -v -w -f data/snps_to_keep data/list_duplicated_var_exclude > data/list_duplicated_var_exclude_rmkeep
plink --bed data/${filename}.bed --bim data/${filename}_edited.bim --fam data/${filename}.fam --exclude data/list_duplicated_var_exclude_rmkeep --keep-allele-order --make-bed --out data/${filename}_nodup

filename=${filename}_nodup

## 5.2. add gender info
## three parameters: 1. work_dir; 2. prefix of filename; 3. phenotype file with gender information
/software/R-4.1.0/bin/Rscript scripts/5_add_gender.R ${wkdir} ${filename} "/lustre/scratch123/hgi/mdt2/projects/ibdgwas_bioresource/phenotype_data/release_20220404/raw/IBD_BioRes_phenotypes_20220404.txt"

plink --bed data/${filename}.bed --bim data/${filename}.bim --fam data/${filename}_gender --keep-allele-order --make-bed --out data/${filename}_gender

filename=${filename}_gender

# 6. Compare freq with EUR 1000GP to deal with A/T   C/G strand exclude A/T and C/G variants with MAF >0.45 (hard to determine if they are wrong)
## 6.1. get 1KG EUR ref
cd ${wkdir}
mkdir 1KG_EUR_study_variant
cat data/${filename}.bim | cut -f 2 > 1KG_EUR_study_variant/list_variants_posstr_nodup

cd 1KG_EUR_study_variant
awk '{if($3=="EUR")print 0,$1}' < /lustre/scratch123/hgi/projects/ibdgwas_bioresource/qz2/ref/integrated_call_samples_v3.20130502.ALL.panel > eur.sample

for((i=1;i<24;i++))
do
	plink2 --bfile /lustre/scratch123/hgi/projects/ibdgwas_bioresource/qz2/ref/1KG/chr${i} --keep eur.sample --rm-dup force-first --extract list_variants_posstr_nodup --make-bed --out chr${i}
done

rm merge.list
for((i=2;i<24;i++))
do
	echo "chr${i}.bed chr${i}.bim chr${i}.fam" >> merge.list
done

plink --bfile chr1 --merge-list merge.list --allow-no-sex --make-bed --out 1000GP_EUR_b38_study_variants

# 6.2. Compare MAF between 1KG EUR and study, and flip allele
cd ${wkdir}

plink --bfile 1KG_EUR_study_variant/1000GP_EUR_b38_study_variants --freq --out 1KG_EUR_study_variant/1000GP_EUR_b38_study_variants
plink --bfile data/${filename} --freq --out data/${filename}

## three parameters: 1. work_dir; 2. filename of 1KG EUR allele freq; 3. filename of study sample allele freq
/software/R-4.1.0/bin/Rscript scripts/6_maf_comp.R ${wkdir} "1KG_EUR_study_variant/1000GP_EUR_b38_study_variants.frq" "data/${filename}.frq"
grep -v -w -f data/snps_to_keep data/list_variants_to_remove_AT_CG > data/list_variants_to_remove_AT_CG_rmkeep
plink --bfile data/${filename} --exclude data/list_variants_to_remove_AT_CG_rmkeep --flip data/list_variants_to_flip_AT_CG --make-bed --out data/${filename}_flip

filename=${filename}_flip

# 7. Sex Discrepancies check
mkdir data/sex
mkdir data/sex/X
mkdir data/sex/Y

awk '{if($5==1)print $1,$2}' < data/${filename}.fam > data/sex/list_male_samples
awk '{if($5==2)print $1,$2}' < data/${filename}.fam > data/sex/list_female_samples

## 7.1. for X chromosome
plink --bfile data/${filename} --keep data/sex/list_female_samples --chr 23 --make-bed --out data/sex/X/study_check_sex_females_only

plink --bfile data/sex/X/study_check_sex_females_only --missing --hardy --freq --out data/sex/X/study_check_sex_females_only

### get high quality variants
## three parameters: 1. work_dir; 2. hwe threshold; 3. prefix of file input; 4.output file name
/software/R-4.1.0/bin/Rscript scripts/7_X_var.R ${wkdir}/data/sex/X  1e-12 "study_check_sex_females_only" "list_good_chrX_variants"

## 7.2. for Y chromosome
plink --bfile data/${filename} --keep data/sex/list_male_samples --chr 24 --make-bed --out data/sex/Y/study_check_sex_males_only

plink --bfile data/sex/Y/study_check_sex_males_only --freq --missing --out data/sex/Y/study_check_sex_males_only

plink --bfile data/${filename} --keep data/sex/list_female_samples --chr 24 --recode A --out data/sex/Y/study_check_sex_females_only_chr24

### get high quality variants
## three parameters: 1. work_dir; 2. prefix of file input (frq/lmiss); 3. dosage file name; 4.output file name
/software/R-4.1.0/bin/Rscript scripts/7_Y_var.R ${wkdir}/data/sex/Y "study_check_sex_males_only" "study_check_sex_females_only_chr24.raw" "list_good_chrY_variants"

cat data/sex/X/list_good_chrX_variants data/sex/Y/list_good_chrY_variants > data/sex/list_good_chrXY_variants

plink --bfile data/${filename} --extract data/sex/list_good_chrXY_variants --make-bed --out data/sex/study_check_sex_good_chrXY

plink --bfile data/sex/study_check_sex_good_chrXY --check-sex ycount --out data/sex/study_check_sex_good_chrXY

plink --bfile data/sex/study_check_sex_good_chrXY --check-sex --out data/sex/study_check_sex_good_chrX

mkdir figures
## three parameters: 1. work_dir; 2. filename of good XY variants; 3. filename of good X variants; 4.prefix of file
/software/R-4.1.0/bin/Rscript scripts/7_check_sex.R ${wkdir}/data/sex "study_check_sex_good_chrXY.sexcheck" "study_check_sex_good_chrX.sexcheck" ${filename}

plink --bed data/${filename}.bed --bim data/${filename}.bim --fam data/${filename}_edited.fam --allow-no-sex --remove data/sex/list_samples_wrong_gender --make-bed --out data/${filename}_sexcheck
filename=${filename}_sexcheck

plink --bfile data/${filename} --allow-no-sex --set-hh-missing --make-bed --out data/${filename}_missinghh
filename=${filename}_missinghh

# 8. remove Y chromosome
awk '{if($1==24)print $2}' < data/${filename}.bim > data/list_chry_variants_toremove

plink --bfile data/${filename} --allow-no-sex --exclude data/list_chry_variants_toremove --make-bed --out data/${filename}_nochry
filename=${filename}_nochry

# 9. remove samples and variants with low call rate
# Sample Call Rate < 80%
plink --bfile data/${filename} --allow-no-sex --mind 0.20 --make-bed --out data/${filename}_scr0.8
filename=${filename}_scr0.8

# Variant Call Rate < 80%
plink --bfile data/${filename} --allow-no-sex --missing --out data/${filename}_vcr0.8
awk '{if($5>0.2) print $2}' data/${filename}_vcr0.8.lmiss > data/list_variants_to_remove_vcr0.8
grep -v -w -f data/snps_to_keep data/list_variants_to_remove_vcr0.8 > data/list_variants_to_remove_vcr0.8_rmkeep
plink --bfile data/${filename} --allow-no-sex --exclude data/list_variants_to_remove_vcr0.8_rmkeep --make-bed --out data/${filename}_vcr0.8
filename=${filename}_vcr0.8

# Sample Call Rate < 95%
plink --bfile data/${filename} --allow-no-sex --mind 0.05 --make-bed --out data/${filename}_scr0.95
filename=${filename}_scr0.95

# Variant Call Rate < 95%
plink --bfile data/${filename} --allow-no-sex --missing --out data/${filename}_vcr0.95
awk '{if($5>0.05) print $2}' data/${filename}_vcr0.95.lmiss > data/list_variants_to_remove_vcr0.95
grep -v -w -f data/snps_to_keep data/list_variants_to_remove_vcr0.95 > data/list_variants_to_remove_vcr0.95_rmkeep
plink --bfile data/${filename} --allow-no-sex --exclude data/list_variants_to_remove_vcr0.95_rmkeep --make-bed --out data/${filename}_vcr0.95
filename=${filename}_vcr0.95

# Variant MAF<0.01 AND Call Rate <98%
plink --bfile data/${filename} --allow-no-sex --missing --out data/${filename}
plink --bfile data/${filename} --allow-no-sex --freq --out data/${filename}
## three parameters: 1. work_dir; 2. prefix of input files; 3. output filename
/software/R-4.1.0/bin/Rscript scripts/9_callrate.R ${wkdir}/data ${filename} "list_monomorphic_vcr098maf0.01_var_exclude"
grep -v -w -f data/snps_to_keep data/list_monomorphic_vcr098maf0.01_var_exclude > data/list_monomorphic_vcr098maf0.01_var_exclude_rmkeep
plink --bfile data/${filename} --allow-no-sex --exclude data/list_monomorphic_vcr098maf0.01_var_exclude_rmkeep --make-bed --out data/${filename}_vcr0.98maf0.01
filename=${filename}_vcr0.98maf0.01

# 10. Variant Missingness discrepancy between disease status (Ca/Ctr) Cases vs Ctr (Fisher exact p-value <1E-4)
# not done here since we only have case samples

##!!!! the following scripts are from Laura's pipeline and have not been tested in my analysis. A few changes have been made so it can be easy to read
mkdir data/disease
mkdir "data/disease/case"
mkdir "data/disease/control"

awk '{if($6==1)print $1,$2}' < data/${filename}.fam > data/disease/list_control_samples
awk '{if($6==2)print $1,$2}' < data/${filename}.fam > data/disease/list_case_samples

plink --bfile data/${filename} --allow-no-sex --keep disease/list_control_samples --missing --out data/${filename}_ctr

plink --bfile data/${filename} --allow-no-sex --keep disease/list_case_samples --missing --out data/${filename}_cases

## three parameters: 1. work_dir; 2. prefix of input files; 3. output file name1; 4. output file name2
/software/R-4.1.0/bin/Rscript scripts/10_dis_discrp.R ${wkdir}/data/disease "comparison_missingness_cases_ctr.lmiss" "list_variants_exclude_by_missingness_cases_ctr.lmiss"
grep -v -w -f data/snps_to_keep data/disease/list_variants_exclude_by_missingness_cases_ctr.lmiss > data/disease/list_variants_exclude_by_missingness_cases_ctr.lmiss_rmkeep
plink --bfile data/${filename} --allow-no-sex --exclude data/disease/list_variants_exclude_by_missingness_cases_ctr.lmiss_rmkeep --make-bed --out data/${filename}_nomissidiscrep
filename=${filename}_nomissidiscrep

# 11. Variant Missingness discrepancy between studies
## I did not do it across studies, but do it across batches (only five main batches including b04,b06,b15,b19,b20 were considered)
batch_old=(b04 b06 b08 b09 b10 b12 b15 b17 b18 b19)
rm data/id_batch data/id_old_batch data/id_new_batch
for i in ${batch_old[@]}
do
	bcftools query -l ${rawDataDir1}/v2chip_${i}_V3_calls.vcf.gz|awk -v bnum="${i}" '{print $0,bnum}' >> data/id_old_batch
done

batch_new="b20"
for i in ${batch_new[@]}
do
	bcftools query -l ${rawDataDir2}/v2chip_${i}_V3_calls.vcf.gz|awk -v bnum="${i}" '{print $0,bnum}' >> data/id_new_batch
done

## three parameters: 1. work_dir; 2. old batch id file; 3. new batch id file
/software/R-4.1.0/bin/Rscript scripts/11_sample_batch.R ${wkdir}/data "id_old_batch" "id_new_batch"

batch_used=(b04 b06 b15 b19 b20)
for i in ${batch_used[@]}
do
        plink --bfile data/${filename} --missing --keep data/${i}.sample --out data/tmp_${i}
done

## three parameters: 1. work_dir; 2. batches investigated (only b04,b06,b15,b19,b20 were considered); 3. summary of variant missingness across batches; 4. file name of list of variants with missingness > 0.1 per study/batch

/software/R-4.1.0/bin/Rscript scripts/11_missing_per_study.R ${wkdir}/data "b04 b06 b15 b19 b20" "table_breakdown_number_missing_variants_perstudy.csv" "list_variants_per_study_missing_0.1"
grep -v -w -f data/snps_to_keep data/list_variants_per_study_missing_0.1 > data/list_variants_per_study_missing_0.1_rmkeep
plink --bfile data/${filename} --allow-no-sex --exclude data/list_variants_per_study_missing_0.1_rmkeep --make-bed --out data/${filename}_perstudy

filename=${filename}_perstudy
rm data/tmp_*

# 12. check if there are duplicates (Kinship>=0.354), when sex are equal, keep the sample with largest call rate(Laura also check pheno, which I did not), when sex/pheno don't match exclude all. For >1 pair duplicated samples, check pair combinations are detected. 1st-degree samples (Intra- Inter- study) Do not exclude but identify and report

plink --bfile data/${filename} --missing --out data/${filename}

# 12.1 1) for samples reported to be duplicates (same FID), but not estimated to be duplicates in our analysis, remove both; 2) for samples reported to be duplicates, and estimated to be duplicates, keep the one with highest call rate; 3) for duplicates not reported by UKIBDBR, but estimated to be duplicates here, select the one with highest call rate.

#### king v2.2.5 not working (run forever), king v2.2.7 works
king -b data/${filename}.bed --related --cpus 1 --prefix data/study_king

#Four parameters: 1. work_dir; 2. prefix of file name of king's result (within and across family); 3. prefix of file(s); 4. output file name
/software/R-4.1.0/bin/Rscript scripts/12_rm_dup_ind.R ${wkdir}/data "study_king" ${filename} "list_duplicated_samples"

# 12.2 extract non-dupulicated samples
plink --bfile data/${filename} --remove data/list_duplicated_samples --make-bed --out data/${filename}_nodupsample
filenamedup=${filename}
filename=${filename}_nodupsample

# 13. PCA analysis to estimate the ancestry (based on variants from autosome chromosomes, on duplicated samples)
mkdir PCA

# 13.1. remove variants in regions with high LD (using non-duplicated samples to calculate LD)
plink --bfile data/${filename} --exclude /lustre/scratch123/hgi/projects/ibdgwas_bioresource/qz2/ref/high-LD-regions-hg38-GRCh38.txt --range --make-bed --out PCA/${filename}_noHighLD

# 13.2. keep independent SNPs, make sure these SNPs also in 1KG
awk '{print $2}' < 1KG_EUR_study_variant/1000GP_EUR_b38_study_variants.bim > PCA/1KG.snp
plink --bfile PCA/${filename}_noHighLD --indep-pairwise 50 5 0.2 --extract PCA/1KG.snp --out PCA/${filename}_noHighLD

grep -v "^chrX" PCA/${filename}_noHighLD.prune.in > PCA/${filename}_noHighLD_nosex.prune.in

plink --bfile PCA/${filename}_noHighLD --extract PCA/${filename}_noHighLD_nosex.prune.in --make-bed --out PCA/${filename}_noHighLD_pruned

# 13.3. Get 1KG with study variants
mkdir PCA/1KG
cd PCA/1KG

###get 2504 unrelated 1KG sample IDs, additional 698 samples that are related to 2504 samples are not considered. see details here https://www.internationalgenome.org/data-portal/data-collection/30x-grch38
awk '{if(NR!=1)print 0,$1}' /lustre/scratch123/hgi/projects/ibdgwas_bioresource/qz2/ref/integrated_call_samples_v3.20130502.ALL.panel > 1KG_urel.id

for((i=1;i<23;i++))
do
        plink2 --bfile /lustre/scratch123/hgi/projects/ibdgwas_bioresource/qz2/ref/1KG/chr${i} --keep 1KG_urel.id --extract ../${filename}_noHighLD_nosex.prune.in --make-bed --out chr${i}
done

rm merge.list
for((i=2;i<23;i++))
do
        echo "chr${i}.bed chr${i}.bim chr${i}.fam" >> merge.list
done

plink --bfile chr1 --merge-list merge.list --allow-no-sex --make-bed --out 1000GP_ALL_b38_study_variants

# 13.4. Combine 1KG with study sample
cd ..

plink --bfile ../data/${filenamedup} --extract ${filename}_noHighLD_nosex.prune.in --make-bed --out ${filenamedup}_noHighLD_pruned
plink --bfile 1KG/1000GP_ALL_b38_study_variants --bmerge ${filenamedup}_noHighLD_pruned.bed ${filenamedup}_noHighLD_pruned.bim ${filenamedup}_noHighLD_pruned.fam --make-bed --out 1000GP_ALL_b38_plus_study_dup

# 13.5 PCA

bsub -J PCA_AllBatch_1KG_dup -q normal -o PCA_AllBatch_1KG_dup.log -e PCA_AllBatch_1KG_dup.err -W 10:00 -m "modern_hardware" -M 40000 -R'select[mem>40000] rusage[mem=40000]' plink2 --memory 39000 --bfile 1000GP_ALL_b38_plus_study_dup --pca --out AllBatch_1KG_dup

bsub -J PCA_AllBatch_1KG_dup_proxy -q normal -o PCA_AllBatch_1KG_dup_proxy.log -e PCA_AllBatch_1KG_dup_proxy.err -W 10:00 -m "modern_hardware" -M 40000 -R'select[mem>40000] rusage[mem=40000]' plink2 --memory 39000 --bfile 1000GP_ALL_b38_plus_study_dup --pca approx --out AllBatch_1KG_dup_proxy


# 13.6. get study samples with different ancestry
cd ${wkdir}
## four parameters: 1. work_dir; 2. 1KG population file; 3. prefix of input file with duplicated samples; 4. prefix of input file with non-duplicated samples; 5. file name of PCA result
/software/R-4.1.0/bin/Rscript scripts/13_keep_all_ancestry.R "${wkdir}/PCA" "/lustre/scratch123/hgi/projects/ibdgwas_bioresource/qz2/ref/integrated_call_samples_v3.20130502.ALL.panel" "${filenamedup}_noHighLD_pruned" "${filename}_noHighLD_pruned" "AllBatch_1KG_dup"

mkdir data/ancestry
ancestry=(eur amr eas afr sas)
for i in ${ancestry[@]}
do
	mkdir data/ancestry/${i}
	plink --bfile data/${filenamedup} --keep PCA/list_${i}_ancestry_samples_withDuplicates --make-bed --out data/ancestry/${i}/${filenamedup}_${i}
	plink --bfile data/${filename} --keep PCA/list_${i}_ancestry_samples_noDuplicates --make-bed --out data/ancestry/${i}/${filename}_${i}
done

# 14. Calculation of HWE P-value (in EUR only, using non-duplicated samples), and check whether variants with small P-value are more likely to show significance in Case/Control analysis

# 14.1 Calculation of HWE P-value (in EUR only) in CD and UC seperately (this is to avoid removing SNPs associated with CD but not UC, or vise verse)
#plink --bfile data/ancestry/eur/${filenamedup}_eur --keep PCA/list_eur_ancestry_samples_noDuplicates --hardy --out data/ancestry/eur/${filename}_eur
#plink --bfile data/ancestry/eur/${filenamedup}_eur --keep PCA/list_eur_ancestry_samples_noDuplicates --hardy --filter-females --chr 23 --out data/ancestry/eur/${filename}_eur_females


awk -F"~" '{if($7==1) print $1}' /lustre/scratch123/hgi/projects/ibdgwas_bioresource/phenotype_data/release_20220404/raw/IBD_BioRes_phenotypes_20220404.txt|grep -w -f - PCA/list_eur_ancestry_samples_noDuplicates > PCA/list_eur_ancestry_samples_noDuplicates_CD

awk -F"~" '{if($7==2) print $1}' /lustre/scratch123/hgi/projects/ibdgwas_bioresource/phenotype_data/release_20220404/raw/IBD_BioRes_phenotypes_20220404.txt|grep -w -f - PCA/list_eur_ancestry_samples_noDuplicates > PCA/list_eur_ancestry_samples_noDuplicates_UC

# test from Qian, try to check whether the HWE P is associated with GWAS P, based on LD independent SNPs
## get LD independent SNPs
plink --bfile data/ancestry/eur/${filenamedup}_eur --keep PCA/list_eur_ancestry_samples_noDuplicates --indep-pairwise 1000 100 0.01 --out data/ancestry/eur/ld

plink --bfile data/ancestry/eur/${filenamedup}_eur --keep PCA/list_eur_ancestry_samples_noDuplicates --freq --out data/ancestry/eur/${filename}_eur
plink --bfile data/ancestry/eur/${filenamedup}_eur --keep PCA/list_eur_ancestry_samples_noDuplicates --hardy --out data/ancestry/eur/${filename}_eur

plink --bfile data/ancestry/eur/${filenamedup}_eur --keep PCA/list_eur_ancestry_samples_noDuplicates_CD --hardy --out data/ancestry/eur/${filename}_eur_CD
plink --bfile data/ancestry/eur/${filenamedup}_eur --keep PCA/list_eur_ancestry_samples_noDuplicates_UC --hardy --out data/ancestry/eur/${filename}_eur_UC

plink --bfile data/ancestry/eur/${filenamedup}_eur --keep PCA/list_eur_ancestry_samples_noDuplicates --freq --filter-females --chr 23 --out data/ancestry/eur/${filename}_eur_females
plink --bfile data/ancestry/eur/${filenamedup}_eur --keep PCA/list_eur_ancestry_samples_noDuplicates --hardy --filter-females --chr 23 --out data/ancestry/eur/${filename}_eur_females
plink --bfile data/ancestry/eur/${filenamedup}_eur --keep PCA/list_eur_ancestry_samples_noDuplicates_CD --hardy --filter-females --chr 23 --out data/ancestry/eur/${filename}_eur_females_CD
plink --bfile data/ancestry/eur/${filenamedup}_eur --keep PCA/list_eur_ancestry_samples_noDuplicates_UC --hardy --filter-females --chr 23 --out data/ancestry/eur/${filename}_eur_females_UC

# 14.2 Check whether variants with significant HWE P-value also have significant GWAS P-value, and found most of such variants are in MHC region
## three parameters: 1. work_dir; 2. prefix of input file(s); 3. IBD GWAS summary statistics from de Lange et al.
/software/R-4.1.0/bin/Rscript scripts/14_check_hwe.R ${wkdir}/data/ancestry ${filename} "/lustre/scratch123/hgi/mdt2/teams/anderson/qz2/scratch115/proj1/data_23_11_2020/lange_smy/28067908-GCST004131-EFO_0003767.h.tsv"

# 14.3 get SNPs not in MHC region but with HWE P-value < 1e-12 (1e-12 is used since they are cases)
## three parameters: 1. work_dir; 2. prefix of input file(s); 3. output file name (list of variants with extreme HWE P-value)
/software/R-4.1.0/bin/Rscript scripts/14_get_hwe.R ${wkdir}/data/ancestry ${filename} "list_var_exclude_hwe"

# 14.4 remove SNPs did not pass HWE test
grep -v -w -f data/snps_to_keep data/ancestry/list_var_exclude_hwe > data/ancestry/list_var_exclude_hwe_rmkeep
ancestry=(eur amr eas afr sas)
for i in ${ancestry[@]}
do
	plink --bfile data/ancestry/${i}/${filenamedup}_${i} --exclude data/ancestry/list_var_exclude_hwe_rmkeep --make-bed --out data/ancestry/${i}/${filenamedup}_${i}_hwe1e-12
done

# 15. Heterozygosity +/-  4SD (in EUR, SAS and AFR only, others have small sample size, threshold of Heterozygosity is determined based on non-duplicated samples)
ancestry=(eur afr sas)
for i in ${ancestry[@]}
do
	plink --bfile data/ancestry/${i}/${filenamedup}_${i}_hwe1e-12 --allow-no-sex --het --out data/ancestry/${i}/het
	## three parameters: 1. work_dir; 2. ancestry; 3. list of non-duplicated samples
	/software/R-4.1.0/bin/Rscript scripts/15_rm_het.R ${wkdir} ${i} PCA/list_${i}_ancestry_samples_noDuplicates
	plink --bfile data/ancestry/${i}/${filenamedup}_${i}_hwe1e-12 --allow-no-sex --remove data/ancestry/${i}/het_rm.list --make-bed --out data/ancestry/${i}/${filenamedup}_${i}_hwe1e-12_het
done

#EUR, non-dup
i="eur"
plink --bfile data/ancestry/${i}/${filename}_${i} --allow-no-sex --exclude data/ancestry/list_var_exclude_hwe_rmkeep --remove data/ancestry/${i}/het_rm.list --make-bed --out data/ancestry/${i}/${filename}_${i}_hwe1e-12_het_tmp
awk '{print $1,$2,$1,$1}' data/ancestry/${i}/${filename}_${i}_hwe1e-12_het_tmp.fam > data/ancestry/${i}/${filename}_${i}_id.map
plink --bfile data/ancestry/${i}/${filename}_${i}_hwe1e-12_het_tmp --update-ids data/ancestry/${i}/${filename}_${i}_id.map --make-bed --out data/ancestry/${i}/${filename}_${i}_hwe1e-12_het
rm data/ancestry/${i}/${filename}_${i}_hwe1e-12_het_tmp*

###################################################################################
## compare MAF with Gnomad, did not perform here since all our samples are cases ##
###################################################################################

# 16. Split cohort into broad ancestry sets, Exclude monomorphic variants
## 16.1 EUR
plink --bfile data/ancestry/eur/${filenamedup}_eur_hwe1e-12_het \
	--keep-allele-order --allow-no-sex \
	--freq \
	--out data/ancestry/eur/${filenamedup}_eur_hwe1e-12_het

awk '{if($5==0)print $2}' < data/ancestry/eur/${filenamedup}_eur_hwe1e-12_het.frq > data/ancestry/eur/study_list_variants_toexclude_monomorphic

plink --bfile data/ancestry/eur/${filenamedup}_eur_hwe1e-12_het \
        --allow-no-sex \
        --exclude data/ancestry/eur/study_list_variants_toexclude_monomorphic \
        --make-bed --out data/ancestry/eur/${filenamedup}_eur_hwe1e-12_het_nomonom


## 16.2 non-EUR
mkdir data/ancestry/non-eur
cat PCA/list_eur_ancestry_samples_withDuplicates data/ancestry/afr/het_rm.list data/ancestry/sas/het_rm.list > data/ancestry/non-eur/eur_sashet_afrhet_samples_withDuplicates

### 16.2.1 generate non-EUR samples, including SAS(pass HET QC), AFR(pass HET QC), AMR, EAS and unknown
plink --bfile data/${filenamedup} \
        --allow-no-sex \
        --exclude data/ancestry/list_var_exclude_hwe_rmkeep \
        --remove data/ancestry/non-eur/eur_sashet_afrhet_samples_withDuplicates \
	--make-bed \
        --out data/ancestry/non-eur/${filenamedup}_non-eur_hwe1e-12_het

### 16.2.2 remove monomorphic SNPs
plink --bfile data/ancestry/non-eur/${filenamedup}_non-eur_hwe1e-12_het \
        --keep-allele-order \
	--allow-no-sex \
        --freq \
        --out data/ancestry/non-eur/${filenamedup}_non-eur_hwe1e-12_het

awk '{if($5==0)print $2}' < data/ancestry/non-eur/${filenamedup}_non-eur_hwe1e-12_het.frq > data/ancestry/non-eur/study_list_variants_toexclude_monomorphic

plink --bfile data/ancestry/non-eur/${filenamedup}_non-eur_hwe1e-12_het \
        --allow-no-sex \
        --exclude data/ancestry/non-eur/study_list_variants_toexclude_monomorphic \
        --make-bed --out data/ancestry/non-eur/${filenamedup}_non-eur_hwe1e-12_het_nomonom

## 16.3 for ID mathcing exercise: MERGE EUR and non-EUR (use SNPs in EUR, if include SNPs exist in non-EUR but not in EUR, it will make that SNP has large missing rate).
mkdir data/ancestry/merge
awk '{print $2}' data/ancestry/eur/${filenamedup}_eur_hwe1e-12_het_nomonom.bim > data/ancestry/eur.snps

plink --bfile data/ancestry/eur/${filenamedup}_eur_hwe1e-12_het_nomonom --bmerge data/ancestry/non-eur/${filenamedup}_non-eur_hwe1e-12_het.bed data/ancestry/non-eur/${filenamedup}_non-eur_hwe1e-12_het.bim data/ancestry/non-eur/${filenamedup}_non-eur_hwe1e-12_het.fam --make-bed --out data/ancestry/merge/${filenamedup}_all_hwe1e-12_het_nomonom --extract data/ancestry/eur.snps


# 17. Update to B37 (liftover) and double check with b37 fasta file, need to do this since Sanger imputation requires b37
mkdir submit_for_imputation

## 17.1 EUR, ***** maybe should do on all samples ******
mkdir submit_for_imputation/eur
awk '{print "chr"$1,$4,$4+1,$2}' data/ancestry/eur/${filenamedup}_eur_hwe1e-12_het_nomonom.bim|sed 's/chr23/chrX/g' > data/ancestry/eur/study_eur_hg38_postqc.bed
liftOver data/ancestry/eur/study_eur_hg38_postqc.bed ../ref/hg38ToHg19.over.chain.gz submit_for_imputation/eur/study_hg38_postqc_lifted_hg19 submit_for_imputation/eur/study_hg38_postqc_no_lifted_hg19

cut -f 4 submit_for_imputation/eur/study_hg38_postqc_no_lifted_hg19 | sed "/^#/d" > submit_for_imputation/eur/nonlifted_hg19_variants_to_exclude_tmp.dat

grep "alt\|random" submit_for_imputation/eur/study_hg38_postqc_lifted_hg19 | cut -f 4 | cat - submit_for_imputation/eur/nonlifted_hg19_variants_to_exclude_tmp.dat > submit_for_imputation/eur/nonlifted_hg19_variants_to_exclude.dat

wc -l submit_for_imputation/eur/nonlifted_hg19_variants_to_exclude.dat

## 17.2  exclude not lifted variants and change position
grep -w -v -f submit_for_imputation/eur/nonlifted_hg19_variants_to_exclude.dat submit_for_imputation/eur/study_hg38_postqc_lifted_hg19|awk '{print $4,$2}' > submit_for_imputation/eur/study_hg38_postqc_lifted_hg19_liftedvariants.map

plink --bfile data/ancestry/eur/${filenamedup}_eur_hwe1e-12_het_nomonom --exclude submit_for_imputation/eur/nonlifted_hg19_variants_to_exclude.dat --update-map submit_for_imputation/eur/study_hg38_postqc_lifted_hg19_liftedvariants.map --make-bed --out submit_for_imputation/eur/study_postqc_lifted_hg19

## 17.3 force A1 and A2 to be ref and alt alleles
zcat data/AllBatch_ATCG_aligned.vcf.gz|cut -f '1-5'| awk '{print "chr"$1":"$2"_"$4"_"$5,$4}'|sed '/^chr##/d'  > submit_for_imputation/eur/list_variants_study_hg38_posstrandaligned_with_A2
wc -l submit_for_imputation/eur/list_variants_study_hg38_posstrandaligned_with_A2
sort submit_for_imputation/eur/list_variants_study_hg38_posstrandaligned_with_A2  | uniq > submit_for_imputation/eur/list_variants_study_hg38_posstrandaligned_with_A2_ed

plink --bfile submit_for_imputation/eur/study_postqc_lifted_hg19 --allow-no-sex --a2-allele submit_for_imputation/eur/list_variants_study_hg38_posstrandaligned_with_A2_ed --make-bed --out submit_for_imputation/eur/study_postqc_lifted_hg19_RefAlt

plink --bfile submit_for_imputation/eur/study_postqc_lifted_hg19_RefAlt --allow-no-sex --keep-allele-order --output-chr MT --recode vcf --out submit_for_imputation/eur/study_postqc_lifted_hg19_RefAlt

## 17.4 alignment
bcftools +fixref submit_for_imputation/eur/study_postqc_lifted_hg19_RefAlt.vcf -Oz -o submit_for_imputation/eur/study_postqc_lifted_hg19_RefAlt_posstrandaligned.vcf.gz -- -f ../ref/human_g1k_v37.fasta -m top 2>&1 | tee submit_for_imputation/eur/alignment_1.log

#### when use --check-ref x, it will exclude incorrect or missing REF allele is encountered, the unsolved variants from previous step would be removed
#### it should be noted that INDEL would be normalized (left aligned) at this step, the position of indel therefore could be changed ####
bcftools norm --check-ref x -f ../ref/human_g1k_v37.fasta submit_for_imputation/eur/study_postqc_lifted_hg19_RefAlt_posstrandaligned.vcf.gz --threads 10 -Oz -o submit_for_imputation/eur/study_postqc_lifted_hg19_RefAlt_posstrandaligned_2.vcf.gz 2>&1 |tee submit_for_imputation/eur/alignment_2.log

## 17.5 change variant ID based on the lifted variant position
plink --vcf submit_for_imputation/eur/study_postqc_lifted_hg19_RefAlt_posstrandaligned_2.vcf.gz --keep-allele-order --allow-no-sex --make-bed --out submit_for_imputation/eur/study_postqc_lifted_hg19_RefAlt_posstrandaligned

zcat submit_for_imputation/eur/study_postqc_lifted_hg19_RefAlt_posstrandaligned_2.vcf.gz|grep -v ^"#" | cut -f '1-5' | awk '{print $1":"$2"_"$4"_"$5,$3}'|sed 's/^X/23/g' > submit_for_imputation/eur/list_variants_study_hg19_posstrandaligned

mkdir submit_for_imputation/ready_for_imp_hg19
mkdir submit_for_imputation/ready_for_imp_hg19/eur
plink --bfile submit_for_imputation/eur/study_postqc_lifted_hg19_RefAlt_posstrandaligned --update-name submit_for_imputation/eur/list_variants_study_hg19_posstrandaligned 1 2 --keep-allele-order --make-bed --out submit_for_imputation/ready_for_imp_hg19/eur/study_hg19

plink --bfile submit_for_imputation/ready_for_imp_hg19/eur/study_hg19 --allow-no-sex --keep-allele-order --output-chr MT --recode vcf --out submit_for_imputation/ready_for_imp_hg19/eur/study_hg19

## 17.6 check if there is any error before submit for imputation
bcftools norm -ce -f ../ref/human_g1k_v37.fasta submit_for_imputation/ready_for_imp_hg19/eur/study_hg19.vcf -Ou -o /dev/null




# 18. HLA imputation
mkdir HLA

## 18.1 Extract varaints from MHC region (chr6:28,477,797-33,448,354 from https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37)
## 18.1.1 EUR ancestry
mkdir HLA/European
plink --bfile submit_for_imputation/ready_for_imp_hg19/eur/study_hg19 --chr 6 --from-bp 28477797 --to-bp 33448354 --make-bed --out HLA/European/hla_hg19

## 18.1.2 all ancestry
mkdir HLA/Broad
plink --bfile data/ancestry/merge/${filenamedup}_all_hwe1e-12_het_nomonom --chr 6 --from-bp 28510120 --to-bp 33480577 --make-bed --out HLA/Broad/hla_hg38

awk '{print "chr"$1,$4,$4+1,$2}' HLA/Broad/hla_hg38.bim > HLA/Broad/hla_hg38_postqc.bed

liftOver HLA/Broad/hla_hg38_postqc.bed ../ref/hg38ToHg19.over.chain.gz HLA/Broad/hla_hg38_postqc_lifted_hg19 HLA/Broad/hla_hg38_postqc_no_lifted_hg19
wc -l HLA/all/hla_hg38_postqc_no_lifted_hg19
# 0

awk '{print $4,$2}' HLA/Broad/hla_hg38_postqc_lifted_hg19 > HLA/Broad/hla_hg38_postqc_lifted_hg19_liftedvariants.map

plink --bfile HLA/Broad/hla_hg38 --update-map HLA/Broad/hla_hg38_postqc_lifted_hg19_liftedvariants.map --make-bed --out HLA/Broad/hla_hg38_postqc_lifted_hg19

plink --bfile HLA/Broad/hla_hg38_postqc_lifted_hg19 --allow-no-sex --a2-allele submit_for_imputation/eur/list_variants_study_hg38_posstrandaligned_with_A2_ed --make-bed --out HLA/Broad/hla_hg38_postqc_lifted_hg19_RefAlt

plink --bfile HLA/Broad/hla_hg38_postqc_lifted_hg19_RefAlt --allow-no-sex --keep-allele-order --output-chr MT --recode vcf --out HLA/Broad/hla_hg38_postqc_lifted_hg19_RefAlt

bcftools +fixref HLA/Broad/hla_hg38_postqc_lifted_hg19_RefAlt.vcf -Oz -o HLA/Broad/hla_hg38_postqc_lifted_hg19_RefAlt_posstrandaligned.vcf.gz -- -f ../ref/human_g1k_v37.fasta -m top 2>&1 | tee HLA/Broad/alignment_1.log

plink --vcf HLA/Broad/hla_hg38_postqc_lifted_hg19_RefAlt_posstrandaligned.vcf.gz --keep-allele-order --allow-no-sex --make-bed --out HLA/Broad/hla_hg38_postqc_lifted_hg19_RefAlt_posstrandaligned

zcat HLA/Broad/hla_hg38_postqc_lifted_hg19_RefAlt_posstrandaligned.vcf.gz|grep -v ^"#" | cut -f '1-5' | awk '{print $1":"$2"_"$4"_"$5,$3}' > HLA/Broad/list_variants_study_hg19_posstrandaligned

plink --bfile HLA/Broad/hla_hg38_postqc_lifted_hg19_RefAlt_posstrandaligned --update-name HLA/Broad/list_variants_study_hg19_posstrandaligned 1 2 --keep-allele-order --make-bed --out HLA/Broad/hla_hg19

## 18.2 imputation
mkdir HLA/imp
mkdir HLA/imp/log

anc="European"
mkdir HLA/imp/${anc}
hla_ids=("A" "B" "C" "DRB1" "DQA1" "DQB1" "DPB1")
for i in ${hla_ids[@]}
do
	bsub -J HLA_IMP_${anc}_${i} -q normal -o HLA/imp/log/${anc}_${i}.log -e HLA/imp/log/${anc}_${i}.err -W 10:00 -m "modern_hardware" -M 20000 -R'select[mem>20000] rusage[mem=20000]' /software/R-4.1.0/bin/Rscript scripts/18_hla_imp.R ${anc} ${i}
done

anc="Broad"
mkdir HLA/imp/${anc}
for i in ${hla_ids[@]}
do
        bsub -J HLA_IMP_${anc}_${i} -q normal -o HLA/imp/log/${anc}_${i}.log -e HLA/imp/log/${anc}_${i}.err -W 10:00 -m "modern_hardware" -M 20000 -R'select[mem>20000] rusage[mem=20000]' /software/R-4.1.0/bin/Rscript scripts/18_hla_imp.R ${anc} ${i}
done


