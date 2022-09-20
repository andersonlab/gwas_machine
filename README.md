# gwas_machine
Platform to run GWAS machinery within the Anderson lab. This will be a 'soup-to-nuts' target identification operation and will include all ingredients necessary to transform the raw genotype data coming into the Anderson lab into high quality imputed data which will be used for GWAS and post-GWAS functional analysis to delve into potential therapeutic targets.

Please email mt27@sanger.ac.uk if you have any comments, bugs, questions or feature suggestions for any of the code included within!

## Using this package

This package is aimed to be used as a fully containerized product with singularity that includes a nextflow pipeline.

Therefore to run this...

## QC Pipeline

This pipeline will involve the transformation of raw genotyped data to pre and post imputed QC'd data used for GWAS analysis. The workflow consists of a .yaml file which is defined for a specific project (eg. IBDBR phenotype analysis or single-cell QTL analysis) but can also be designed adhoc to run a user-specified QC pipeline. 

#### Functions Defintion and Usage

Briefly, here is an overview and tutorial of all the functions (some of which are only used for internal pipeline building purposes) used in the QC scripts (located in `R/genetic_qc.R`). Some of the helper scripts for each of these functions are located in `./scripts`

* `prepare_reference_data()`: Locate and if necessary, download the data required to run the QC in some of the associated steps in this pipeline. This data is downloaded to (`/lustre/scratch123/hgi/mdt2/projects/gwas_machine/gwasmachine/ref`). The steps associated in this data downloaded and pre-processing are:
  + Downloading the raw 1000G data from:         `http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/` and converting    the file into format `chr:pos_ref_alt`
  + Extract the high LD regions
  + Get the reference data for alignment to hg19
  + Get reference data for lift-over from hg38 to hg19
  + Get reference data for HLA imputation (using HIBAG)

* `update_sample_ids`: Sample IDs are updated (this is for UKIBDBR only) -- for samples with same FID, change their IID as FID-1/2/3)

* `keep_specific_snps` : Tracks list of specific SNPs to keep that may be excluded for other reasons in the QC pipeline to bring to imputation

* `keep_atcg_snps` :Remove Indels, Mitochondrial Variants (None ATCG SNPs) -> `AllBatches_ATCG`

* `align_positive_strands`: Align to (+) strand, and change format to bed file -> `AllBatches_ATCG_Aligned`

* `update_variant_id`: Update variant ids (using chr:position_ref_alt) and remove duplicated variants (same chr:pos:ref:alt, remove ones with higher missingness) -> `AllBatch_ATCG_aligned_nodup`

*`add_gender_info`: Add gender info to the fam plink files -> `AllBatch_ATCG_aligned_nodup_gender`

*`compare_freq_1000G_EUR`: Compare freq with EUR 1000GP to deal with A/T   C/G strand exclude A/T and C/G variants with MAF >0.45 then flip allele if necessary -> `AllBatch_ATCG_aligned_nodup_gender_filp`







