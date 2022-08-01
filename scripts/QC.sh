#!/bin/bash
#This document follows the pipeline shared by Laura, 18-01-2021: https://drive.google.com/file/d/1RmDcq90uPAWZMZaKy3l5Vy3JtG3jbwq5/view?usp=sharing
#A few changes have been made to fit this data.
#Qian Zhang, 01-03-2022

wkdir="/lustre/scratch123/hgi/projects/ibdgwas_bioresource/qz2/pre_imputation_b04b06b15b19b20"

cd /lustre/scratch123/hgi/projects/ibdgwas_bioresource/qz2/

## prepare ref data for the following analysis
### 1KG 30x on GRCh38, description on here:https://www.internationalgenome.org/data-portal/data-collection/30x-grch38
#mkdir ref
#mkdir ref/1KG
#mkdir ref/log
### I set the variant ID in format: chr:pos_ref_alt, only the first 20 characters were used for alleles with length longer than 20
#
#for((i=1;i<23;i++))
#do
#	plink2 --vcf /lustre/scratch115/resources/1000g/release/20201028/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${i}.filtered.shapeit2-duohmm-phased.vcf.gz --set-all-var-ids chr@:\#_\$r_\$a --new-id-max-allele-len 20 truncate --max-alleles 2 --rm-dup force-first --make-bed --out ref/1KG/chr${i}
#done
### X chromosome
#### remove header line at first, otherwise there would be error like "Duplicate FORMAT:GT header line in --vcf file."
#zcat /lustre/scratch115/resources/1000g/release/20201028/CCDG_14151_B01_GRM_WGS_2020-08-05_chrX.filtered.eagle2-phased.vcf.gz|grep -v "^##" > ref/1KG/chr23_tmp
#plink2 --vcf ref/1KG/chr23_tmp --set-all-var-ids chr@:\#_\$r_\$a --new-id-max-allele-len 20 truncate --max-alleles 2 --rm-dup force-first --make-bed --out ref/1KG/chr23
#rm ref/1KG/chr23_tmp
#
##get high LD region
#wget https://raw.githubusercontent.com/meyer-lab-cshl/plinkQC/master/inst/extdata/high-LD-regions-hg38-GRCh38.txt -O ref/high-LD-regions-hg38-GRCh38.txt
#sed -i 's/chr//g' ref/high-LD-regions-hg38-GRCh38.txt
#
#
## get reference data for liftover from hg38 to hg19
#wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz' -O ref/hg38ToHg19.over.chain.gz
#
## get reference data for alignment to hg19
#wget ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.gz -O ref/human_g1k_v37.fasta.gz
#gunzip ref/human_g1k_v37.fasta.gz


# data preparation
# 0. transfer from vcf to plink bed format
cd ${wkdir}
rawDataDir="/lustre/scratch123/hgi/projects/ibdgwas_bioresource/genotyped_data/raw/data_transfer_20220324"

mkdir data
batch="b20"
for i in ${batch[@]}
do
	mkdir data/${i}
	# extract SNPs from chromosome 1-22,X,Y,MT, SNPs from other chromosomes were removed, otherwise there would be error like "Invalid chromosome code"
	zcat ${rawDataDir}/v2chip_${i}_V3_calls.vcf.gz| grep -E -w ^'[1-9]|1[0-9]|2[0-2]|X|Y|MT|#CHROM' > data/${i}/tmp
	# change to plink format
	plink --vcf data/${i}/tmp --allow-extra-chr --output-chr MT --make-bed --out data/${i}/${i}
	rm data/${i}/tmp
done


rm data/snps_to_keep
#same SNP, just in different formats
echo "chr16:50729870_C_CC" >> data/snps_to_keep
echo "chr16:50729870_CC_C" >> data/snps_to_keep
echo "AX-96079897" >> data/snps_to_keep

# 1. Remove Indels, Mitochondrial Variants (None ATCG SNPs)
for i in ${batch[@]}
do
	#two parameters: 1. work_dir; 2.batch
	/software/R-4.1.0/bin/Rscript scripts/1_rm_non_ATCG.R ${wkdir} ${i}
	grep -v -w -f data/snps_to_keep data/${i}/list_indel_var_exclude_${i} > data/${i}/list_indel_var_exclude_rmkeep_${i}
	plink --bfile data/${i}/${i} --allow-extra-chr --output-chr MT --exclude data/${i}/list_indel_var_exclude_rmkeep_${i} --make-bed --out data/${i}/${i}_ATCG
done

# 2. remove samples with no phenotype data or samples that have withdrawn consent (did not perform here, applied after data merge)

# 3. Update to B37 (not done here since our data is hg38)

# 4. Align to (+) strand
for i in ${batch[@]}
do
	plink --bfile data/${i}/${i}_ATCG --allow-extra-chr --output-chr MT --allow-no-sex --recode vcf --out data/${i}/${i}_ATCG
	bcftools +fixref data/${i}/${i}_ATCG.vcf -Oz -o data/${i}/${i}_ATCG_aligned.vcf.gz -- -f /lustre/scratch123/hgi/projects/ibdgwas/IIBDGC/resources/hg38/hg38_edited.fa -m top 2>&1 | tee data/${i}/alignment_${i}.log
	# after alignment, format vcf to bed file
	plink --vcf data/${i}/${i}_ATCG_aligned.vcf.gz --id-delim --make-bed --out data/${i}/${i}_ATCG_aligned
	rm data/${i}/${i}_ATCG.vcf
done

# 5. Update variant ids (using chr:position_ref_alt); Remove duplicated variants (same Chr pos ref alt, remove ones with more missingness), merge batches
# 5.1. update varaint id and remove duplicated variants
for i in ${batch[@]}
do
	plink --bfile data/${i}/${i}_ATCG_aligned --missing --out data/${i}/${i}_ATCG_aligned
	zcat data/${i}/${i}_ATCG_aligned.vcf.gz|grep -v ^"#"|cut -f '1-5' | awk '{print $3,"chr"$1":"$2"_"$4"_"$5}' > data/${i}/${i}_ATCG_aligned

	# for duplicated variants(same chr, pos, ref, alt), remove the one with more missingness, for variants with same missingness, keep the first one
	# two parameters: 1. work_dir; 2.batch
	/software/R-4.1.0/bin/Rscript scripts/5_rm_dup_var.R ${wkdir} ${i}
	grep -v -w -f data/snps_to_keep data/${i}/list_duplicated_var_exclude_${i} > data/${i}/list_duplicated_var_exclude_rmkeep_${i}
	plink --bed data/${i}/${i}_ATCG_aligned.bed --bim data/${i}/${i}_ATCG_aligned_edited.bim --fam data/${i}/${i}_ATCG_aligned.fam --exclude data/${i}/list_duplicated_var_exclude_rmkeep_${i} --make-bed --out data/${i}/${i}_ATCG_aligned_nodup
done

## 5.2. merge batches together
rm data/merge.list
batch_old=(b04 b06 b08 b09 b10 b12 b15 b17 b18 b19)
for i in ${batch_old[@]}
do
        echo "/lustre/scratch123/hgi/mdt2/projects/ibdgwas_bioresource/qz2/pre_imputation_b04b06b15b19/data/${i}/${i}_ATCG_aligned_nodup.bed /lustre/scratch123/hgi/mdt2/projects/ibdgwas_bioresource/qz2/pre_imputation_b04b06b15b19/data/${i}/${i}_ATCG_aligned_nodup.bim /lustre/scratch123/hgi/mdt2/projects/ibdgwas_bioresource/qz2/pre_imputation_b04b06b15b19/data/${i}/${i}_ATCG_aligned_nodup.fam" >> data/merge.list
done

filename="AllBatch_ATCG_aligned_nodup"
plink --bfile data/b20/b20_ATCG_aligned_nodup --merge-list data/merge.list --make-bed --out data/${filename}


## 5.3. add gender info
## three parameters: 1. work_dir; 2. prefix of filename; 3. phenotype file with gender information
/software/R-4.1.0/bin/Rscript scripts/5_add_gender.R ${wkdir} ${filename} "/lustre/scratch123/hgi/mdt2/projects/ibdgwas_bioresource/phenotype_data/release_20220323/raw/IBD_BioRes_phenotypes_20220323_1.csv"

plink --bed data/${filename}.bed --bim data/${filename}.bim --fam data/${filename}_gender --make-bed --out data/${filename}_gender

filename=${filename}_gender

## 5.4. check if there are duplicates (Kinship>=0.354), when sex are equal, keep the sample with largest call rate(Laura also check pheno, which I did not), when sex/pheno don't match exclude all. For >1 pair duplicated samples, check pair combinations are detected. 1st-degree samples (Intra- Inter- study) Do not exclude but identify and report

plink --bfile data/${filename} --missing --out data/${filename}

### 5.4.1 1) for samples reported to be duplicates (same FID), but not estimated to be duplicates in our analysis, remove both; 2) for samples reported to be duplicates, and estimated to be duplicates, keep the one with highest call rate; 3) for duplicates not reported by UKIBDBR, but estimated to be duplicates here, select the one with highest call rate.

#### king v2.2.5 not working (run forever), king v2.2.7 works
king -b data/${filename}.bed --related --cpus 1 --prefix data/study_king

#Four parameters: 1. work_dir; 2. prefix of file name of king's result (within and across family); 3. prefix of file(s); 4. output file name
/software/R-4.1.0/bin/Rscript scripts/5_rm_dup_ind.R ${wkdir}/data "study_king" ${filename} "list_duplicated_samples"

### 5.4.2 extract non-dupulicated samples
plink --bfile data/${filename} --remove data/list_duplicated_samples --make-bed --out data/${filename}_nodupsample
filename=${filename}_nodupsample

# 6. Compare freq with EUR 1000GP to deal with A/T   C/G strand exclude A/T and C/G variants with MAF >0.45 (hard to determine if they are wrong)
## 6.1. get 1KG EUR ref
cd ${wkdir}
mkdir 1KG_EUR_study_variant
cat data/${filename}.bim | cut -f 2 > 1KG_EUR_study_variant/list_variants_posstr_nodup

cd 1KG_EUR_study_variant
awk '{if($3=="EUR")print 0,$1}' < /lustre/scratch115/resources/1000g/release/20130502/integrated_call_samples_v3.20130502.ALL.panel > eur.sample

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
## I did not do it across studies, but do it across batches (only four main batches including b04,b06,b15,b19 were considered)
batch=(b04 b06 b15 b19)
for i in ${batch[@]}
do
        plink --bfile data/${filename} --missing --keep data/${i}/${i}.fam --out data/tmp_${i}
done

## three parameters: 1. work_dir; 2. batches investigated (only b04,b06,b15,b19 were considered); 3. summary of variant missingness across batches; 4. file name of list of variants with missingness > 0.1 per study/batch
/software/R-4.1.0/bin/Rscript scripts/11_missing_per_study.R ${wkdir}/data "b04 b06 b15 b19" "table_breakdown_number_missing_variants_perstudy.csv" "list_variants_per_study_missing_0.1"
grep -v -w -f data/snps_to_keep data/list_variants_per_study_missing_0.1 > data/list_variants_per_study_missing_0.1_rmkeep
plink --bfile data/${filename} --allow-no-sex --exclude data/list_variants_per_study_missing_0.1 --make-bed --out data/${filename}_perstudy

filename=${filename}_perstudy
rm data/tmp_*

# 12. Remove duplidated samples
# (Kinship>=0.354, double sex (Laura also check pheno, which I did not), when equal keep the sample with largest call rate, when sex/pheno don't match exclude all) For >1 pair duplicated samples, check pair combinations are detected. 1st-degree samples (Intra- Inter- study ) Do not exclude but identify and report

king -b data/${filename}.bed --related --cpus 1 --prefix data/study_king

plink --bfile data/${filename} --allow-no-sex --missing --out data/${filename}

## four parameters: 1. work_dir; 2. file name of king's result; 3. prefix of file(s); 4. output file name
/software/R-4.1.0/bin/Rscript scripts/12_rm_dup_ind.R ${wkdir}/data "study_king.kin0" ${filename} "list_duplicated_samples"

plink --bfile data/${filename} --remove data/list_duplicated_samples --make-bed --out data/${filename}_nodupsample
filename=${filename}_nodupsample

# 13. PCA analysis and select EUR samples (autosome chromosomes)
mkdir PCA

# 13.1. remove variants in regions with high LD
plink --bfile data/${filename} --exclude /lustre/scratch123/hgi/projects/ibdgwas_bioresource/qz2/ref/high-LD-regions-hg38-GRCh38.txt --range --make-bed --out PCA/${filename}_noHighLD

# 13.2. keep independent SNPs, make sure these SNPs also in 1KG
awk '{print $2}' < 1KG_EUR_study_variant/1000GP_EUR_b38_study_variants.bim > PCA/1KG.snp
plink --bfile PCA/${filename}_noHighLD --indep-pairwise 50 5 0.2 --extract PCA/1KG.snp --out PCA/${filename}_noHighLD

grep -v "^chrX" PCA/${filename}_noHighLD.prune.in > PCA/${filename}_noHighLD_nosex.prune.in

plink --bfile PCA/${filename}_noHighLD --extract PCA/${filename}_noHighLD_nosex.prune.in --make-bed --out PCA/${filename}_noHighLD_pruned

# 13.3. Get 1KG with study variants
mkdir PCA/1KG
cd PCA/1KG

###get 2504 unrelated 1KG sample IDs, additional 698 samples are related to 2504 samples. https://www.internationalgenome.org/data-portal/data-collection/30x-grch38
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

plink --bfile 1KG/1000GP_ALL_b38_study_variants --bmerge ${filename}_noHighLD_pruned.bed ${filename}_noHighLD_pruned.bim ${filename}_noHighLD_pruned.fam --make-bed --out 1000GP_ALL_b38_plus_study

# 13.5. PCA analysis
bsub -J PCA_AllBatch_1KG -q normal -o PCA_AllBatch_1KG.log -e PCA_AllBatch_1KG.err -W 10:00 -m "modern_hardware" -M 40000 -R'select[mem>40000] rusage[mem=40000]' plink --memory 39000 --bfile 1000GP_ALL_b38_plus_study --pca --out AllBatch_1KG

# 13.6. get study samples with different ancestry
cd ${wkdir}
## four parameters: 1. work_dir; 2. 1KG population file; 3. prefix of input file(s); 4. file name of PCA result
/software/R-4.1.0/bin/Rscript scripts/13_keep_all_ancestry.R "${wkdir}/PCA" "/lustre/scratch123/hgi/projects/ibdgwas_bioresource/qz2/ref/integrated_call_samples_v3.20130502.ALL.panel" "${filename}_noHighLD_pruned" "AllBatch_1KG"

mkdir data/ancestry
ancestry=(eur amr eas afr sas)
for i in ${ancestry[@]}
do
	mkdir data/ancestry/${i}
	plink --bfile data/${filename} --keep PCA/list_${i}_ancestry_samples_noDuplicates --make-bed --out data/ancestry/${i}/${filename}_${i}
done

# 14. Calculation of HWE P-value (in EUR only), and check whether variants with small P-value are more likely to show significance in Case/Control analysis

# 14.1 Calculation of HWE P-value (in EUR only)
plink --bfile data/ancestry/eur/${filename}_eur --hardy --out data/ancestry/eur/${filename}_eur
plink --bfile data/ancestry/eur/${filename}_eur --hardy --filter-females --chr 23 --out data/ancestry/eur/${filename}_eur_females

# 14.2 Check whether variants with significant HWE P-value also have significant GWAS P-value, and found most of such variants are in MHC region
## three parameters: 1. work_dir; 2. prefix of input file(s); 3. GWAS summary statistics of IBD
/software/R-4.1.0/bin/Rscript scripts/14_check_hwe.R ${wkdir}/data/ancestry ${filename} "/lustre/scratch123/hgi/mdt2/teams/anderson/qz2/scratch115/proj1/data_23_11_2020/lange_smy/28067908-GCST004131-EFO_0003767.h.tsv"

# 14.3 get SNPs not in MHC region but with HWE P-value < 1e-12 (1e-12 is used since they are cases)
## three parameters: 1. work_dir; 2. prefix of input file(s); 3. output file name (list of variants with extreme HWE P-value)
/software/R-4.1.0/bin/Rscript scripts/14_get_hwe.R ${wkdir}/data/ancestry ${filename} "list_var_exclude_hwe"


# 14.4 remove SNPs did not pass HWE test
grep -v -w -f data/snps_to_keep data/ancestry/list_var_exclude_hwe > data/ancestry/list_var_exclude_hwe_rmkeep
ancestry=(eur amr eas afr sas)
for i in ${ancestry[@]}
do
	plink --bfile data/ancestry/${i}/${filename}_${i} --exclude data/ancestry/list_var_exclude_hwe_rmkeep --make-bed --out data/ancestry/${i}/${filename}_${i}_hwe1e-12
done

# 15. Heterozygosity +/-  4SD (in EUR, SAS and AFR only, others have small sample size)
ancestry=(eur afr sas)
for i in ${ancestry[@]}
do
	plink --bfile data/ancestry/${i}/${filename}_${i}_hwe1e-12 --allow-no-sex --het --out data/ancestry/${i}/het
	## two parameters: 1. work_dir; 2. ancestry
	/software/R-4.1.0/bin/Rscript scripts/15_rm_het.R ${wkdir} ${i}
	plink --bfile data/ancestry/${i}/${filename}_${i}_hwe1e-12 --allow-no-sex --remove data/ancestry/${i}/het_rm.list --make-bed --out data/ancestry/${i}/${filename}_${i}_hwe1e-12_het
done

# 16. Split cohort in broad ancestry sets, Exclude monomorphic variants
## 16.1 EUR
plink --bfile data/ancestry/eur/${filename}_eur_hwe1e-12_het \
	--keep-allele-order --allow-no-sex \
	--freq \
	--out data/ancestry/eur/${filename}_eur_hwe1e-12_het

awk '{if($5==0)print $2}' < data/ancestry/eur/${filename}_eur_hwe1e-12_het.frq > data/ancestry/eur/study_list_variants_toexclude_monomorphic

plink --bfile data/ancestry/eur/${filename}_eur_hwe1e-12_het \
        --allow-no-sex \
        --exclude data/ancestry/eur/study_list_variants_toexclude_monomorphic \
        --make-bed --out data/ancestry/eur/${filename}_eur_hwe1e-12_het_nomonom

## 16.2 non-EUR
mkdir data/ancestry/non-eur
cat PCA/list_eur_ancestry_samples_noDuplicates data/ancestry/afr/het_rm.list data/ancestry/sas/het_rm.list > data/ancestry/non-eur/eur_sashet_afrhet_samples_noDuplicates

### 16.2.1 generate non-EUR samples, including SAS(pass HET QC), AFR(pass HET QC), AMR, EAS and unknown
plink --bfile data/${filename} \
        --allow-no-sex \
        --exclude data/ancestry/list_var_exclude_hwe_rmkeep \
        --remove data/ancestry/non-eur/eur_sashet_afrhet_samples_noDuplicates \
	--make-bed \
        --out data/ancestry/non-eur/${filename}_non-eur_hwe1e-12_het

### 16.2.2 remove monomorphic SNPs
plink --bfile data/ancestry/non-eur/${filename}_non-eur_hwe1e-12_het \
        --keep-allele-order \
	--allow-no-sex \
        --freq \
        --out data/ancestry/non-eur/${filename}_non-eur_hwe1e-12_het

awk '{if($5==0)print $2}' < data/ancestry/non-eur/${filename}_non-eur_hwe1e-12_het.frq > data/ancestry/non-eur/study_list_variants_toexclude_monomorphic

plink --bfile data/ancestry/non-eur/${filename}_non-eur_hwe1e-12_het \
        --allow-no-sex \
        --exclude data/ancestry/non-eur/study_list_variants_toexclude_monomorphic \
        --make-bed --out data/ancestry/non-eur/${filename}_non-eur_hwe1e-12_het_nomonom

## 16.3 for ID mathcing exercise: MERGE EUR and non-EUR (use SNPs in EUR, if also include SNPs exist in non-EUR but not in EUR, that will make that SNP has large missing rate).
mkdir data/ancestry/merge
awk '{print $2}' data/ancestry/eur/${filename}_eur_hwe1e-12_het_nomonom.bim > data/ancestry/eur.snps

plink --bfile data/ancestry/eur/${filename}_eur_hwe1e-12_het_nomonom --bmerge data/ancestry/non-eur/${filename}_non-eur_hwe1e-12_het.bed data/ancestry/non-eur/${filename}_non-eur_hwe1e-12_het.bim data/ancestry/non-eur/${filename}_non-eur_hwe1e-12_het.fam --make-bed --out data/ancestry/merge/${filename}_all_hwe1e-12_het_nomonom --extract data/ancestry/eur.snps


# 17. Update to B37 (liftover) and double check with b37 fasta file, need to do this since sanger imputation requires b37
mkdir submit_for_imputation

## 17.1 EUR
mkdir submit_for_imputation/eur
awk '{print "chr"$1,$4,$4+1,$2}' data/ancestry/eur/${filename}_eur_hwe1e-12_het_nomonom.bim|sed 's/chr23/chrX/g' > data/ancestry/eur/study_eur_hg38_postqc.bed
liftOver data/ancestry/eur/study_eur_hg38_postqc.bed ../ref/hg38ToHg19.over.chain.gz submit_for_imputation/eur/study_hg38_postqc_lifted_hg19 submit_for_imputation/eur/study_hg38_postqc_no_lifted_hg19

cut -f 4 submit_for_imputation/eur/study_hg38_postqc_no_lifted_hg19 | sed "/^#/d" > submit_for_imputation/eur/nonlifted_hg19_variants_to_exclude_tmp.dat

grep "alt\|random" submit_for_imputation/eur/study_hg38_postqc_lifted_hg19 | cut -f 4 | cat - submit_for_imputation/eur/nonlifted_hg19_variants_to_exclude_tmp.dat > submit_for_imputation/eur/nonlifted_hg19_variants_to_exclude.dat

wc -l submit_for_imputation/eur/nonlifted_hg19_variants_to_exclude.dat

## 17.2  exclude not lifted variants and change position
grep -w -v -f submit_for_imputation/eur/nonlifted_hg19_variants_to_exclude.dat submit_for_imputation/eur/study_hg38_postqc_lifted_hg19|awk '{print $4,$2}' > submit_for_imputation/eur/study_hg38_postqc_lifted_hg19_liftedvariants.map

plink --bfile data/ancestry/eur/${filename}_eur_hwe1e-12_het_nomonom --exclude submit_for_imputation/eur/nonlifted_hg19_variants_to_exclude.dat --update-map submit_for_imputation/eur/study_hg38_postqc_lifted_hg19_liftedvariants.map --make-bed --out submit_for_imputation/eur/study_postqc_lifted_hg19

## 17.3 force A1 and A2 to be ref and alt alleles
zcat data/b04/b04_ATCG_aligned.vcf.gz|cut -f '1-5'| awk '{print "chr"$1":"$2"_"$4"_"$5,$4}'|sed '/^chr##/d'  > submit_for_imputation/eur/list_variants_study_hg38_posstrandaligned_with_A1
wc -l submit_for_imputation/eur/list_variants_study_hg38_posstrandaligned_with_A1
sort submit_for_imputation/eur/list_variants_study_hg38_posstrandaligned_with_A1  | uniq > submit_for_imputation/eur/list_variants_study_hg38_posstrandaligned_with_A1_ed

plink --bfile submit_for_imputation/eur/study_postqc_lifted_hg19 --allow-no-sex --a2-allele submit_for_imputation/eur/list_variants_study_hg38_posstrandaligned_with_A1_ed --make-bed --out submit_for_imputation/eur/study_postqc_lifted_hg19_RefAlt

plink --bfile submit_for_imputation/eur/study_postqc_lifted_hg19_RefAlt --allow-no-sex --keep-allele-order --output-chr MT --recode vcf-fid --out submit_for_imputation/eur/study_postqc_lifted_hg19_RefAlt

## 17.4 alignment
bcftools +fixref submit_for_imputation/eur/study_postqc_lifted_hg19_RefAlt.vcf -Oz -o submit_for_imputation/eur/study_postqc_lifted_hg19_RefAlt_posstrandaligned.vcf.gz -- -f ../ref/human_g1k_v37.fasta -m top 2>&1 | tee submit_for_imputation/eur/alignment_1.log

#### when use --check-ref x, it will exclude incorrect or missing REF allele is encountered, the unsolved variants from previous step would be removed
#### it should be noted that INDEL would be normalized (left aligned) at this step, the position of indel therefore could be changed ####
bcftools norm --check-ref x -f ../ref/human_g1k_v37.fasta submit_for_imputation/eur/study_postqc_lifted_hg19_RefAlt_posstrandaligned.vcf.gz --threads 10 -Oz -o submit_for_imputation/eur/study_postqc_lifted_hg19_RefAlt_posstrandaligned_2.vcf.gz 2>&1 |tee submit_for_imputation/eur/alignment_2.log

## 17.5 change variant ID based on the lifted variant position
plink --vcf submit_for_imputation/eur/study_postqc_lifted_hg19_RefAlt_posstrandaligned_2.vcf.gz --keep-allele-order --allow-no-sex --double-id --make-bed --out submit_for_imputation/eur/study_postqc_lifted_hg19_RefAlt_posstrandaligned

zcat submit_for_imputation/eur/study_postqc_lifted_hg19_RefAlt_posstrandaligned_2.vcf.gz|grep -v ^"#" | cut -f '1-5' | awk '{print $1":"$2"_"$4"_"$5,$3}'|sed 's/^X/23/g' > submit_for_imputation/eur/list_variants_study_hg19_posstrandaligned

mkdir submit_for_imputation/ready_for_imp_hg19
mkdir submit_for_imputation/ready_for_imp_hg19/eur
plink --bfile submit_for_imputation/eur/study_postqc_lifted_hg19_RefAlt_posstrandaligned --update-name submit_for_imputation/eur/list_variants_study_hg19_posstrandaligned 1 2 --keep-allele-order --make-bed --out submit_for_imputation/ready_for_imp_hg19/eur/study_hg19

plink --bfile submit_for_imputation/ready_for_imp_hg19/eur/study_hg19 --allow-no-sex --keep-allele-order --output-chr MT --recode vcf-fid --out submit_for_imputation/ready_for_imp_hg19/eur/study_hg19

## 17.6 check if there is any error before submit for imputation
bcftools norm -ce -f ../ref/human_g1k_v37.fasta submit_for_imputation/ready_for_imp_hg19/eur/study_hg19.vcf -Ou -o /dev/null
