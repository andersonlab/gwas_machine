## 17.4 alignment
cd /lustre/scratch123/hgi/mdt2/projects/ibdgwas_bioresource/qz2/pre_imputation_b04b06b15b19b20/

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
