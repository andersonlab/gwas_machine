# These scripts will be the ones used in the separate snake-make rules
batches = c("b04","b06","b08","b09","b10","b12","b15","b17","b18","b19")

#### Read in the YAML file with the configurations
download_1000G_reference()

convert_vcf_to_plink(raw_vcf_data_directory = "/lustre/scratch123/hgi/projects/ibdgwas_bioresource/raw/data_transfer_20220228/",
                     batch_numbers          = batches,
                     output_directory       = "/lustre/scratch123/hgi/projects/ibdgwas_bioresource/gwas_machine_testing")

keep_specific_snps()

keep_atcg_snps(output_directory = "/lustre/scratch123/hgi/projects/ibdgwas_bioresource/gwas_machine_testing",
               batch_numbers    = batches)

align_positive_strand(output_directory = "/lustre/scratch123/hgi/projects/ibdgwas_bioresource/gwas_machine_testing",
                      batch_numbers    = batches)

update_variant_id(output_directory = "/lustre/scratch123/hgi/projects/ibdgwas_bioresource/gwas_machine_testing",
                   batch_numbers    = batches)

remove_missing_snps(output_directory = "/lustre/scratch123/hgi/projects/ibdgwas_bioresource/gwas_machine_testing",
                  batch_numbers    = batches)

remove_duplicated_variants(output_directory = "/lustre/scratch123/hgi/projects/ibdgwas_bioresource/gwas_machine_testing",
                           batch_numbers    = batches)

merge_batches(output_directory = "/lustre/scratch123/hgi/projects/ibdgwas_bioresource/gwas_machine_testing",
              batch_numbers    = batches)

calculate_missigness()

add_gender()

identify_duplicates()

compare_freq_1000G_EUR()

sex_discrepancy_check()

remove_Y_chr()

remove_low_call_rate()

variant_missigness_across_batches()

PCA_Ancestry()



