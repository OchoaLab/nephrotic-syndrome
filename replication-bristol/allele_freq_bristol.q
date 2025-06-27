#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=allele_freq
#SBATCH --output=allele_freq_bristol.out
#SBATCH --mem=64G
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=tiffany.tu@duke.edu
#SBATCH --mail-type=END,FAIL

# NOTE: 32G sufficed for HO and HGDP, but increased for TGP
dir='/datacommons/ochoalab/ssns_gwas/replication/saige_results/bristol/AF'
module load Plink/2.00a2LM

# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep $dir/ssns_ctrl/bristol_ssns_euro_allage.txt --extract $dir/ukbb/ssns_ctrl/clean_control_snp_id.txt --freq 'counts' --out $dir/ukbb/ssns_ctrl/bristolAC_eur
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep $dir/ssns_ctrl/bristol_ssns_sas_allage.txt --extract $dir/ukbb/ssns_ctrl/clean_control_snp_id.txt --freq 'counts' --out $dir/ukbb/ssns_ctrl/bristolAC_sas
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep $dir/ssns_ctrl/bristol_ssns_af_allage.txt --extract $dir/ukbb/ssns_ctrl/clean_control_snp_id.txt --freq 'counts' --out $dir/ukbb/ssns_ctrl/bristolAC_afr
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep $dir/ssns_ctrl/bristol_ssns_allage.txt --make-bed --out $dir/ssns_ctrl/bristol_ssns_allage

# ssns eur
#time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep $dir/ssns_ctrl/bristol_ssns_euro_allage.txt --extract $dir/ukbb/ssns_ctrl/eur/clean_control_snp_id.txt --freq 'counts' --out $dir/ukbb/ssns_ctrl/eur/bristolAC_eur
# ssns sas
#time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep $dir/ssns_ctrl/bristol_ssns_sas_allage.txt --extract $dir/ukbb/ssns_ctrl/sas/clean_control_snp_id.txt --freq 'counts' --out $dir/ukbb/ssns_ctrl/sas/bristolAC_sas



#time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute --keep $dir/ssns_ctrl/bristol_ssns_euro_allage.txt --extract $dir/ukbb/control_snp_id_ssns_all.txt --freq 'counts' --out $dir/ukbb/bristolAC_eur_nofilter
#time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute --keep $dir/ssns_ctrl/bristol_ssns_sas_allage.txt --extract $dir/ukbb/control_snp_id_ssns_all.txt --freq 'counts' --out $dir/ukbb/bristolAC_sas_nofilter
#time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute --extract $dir/ukbb/control_snp_id_ssns_all.txt --freq 'counts' --out $dir/ukbb/bristolAC_all_nofilter
#time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute --keep $dir/ssns_ctrl/bristol_ssns_af_allage.txt --extract $dir/ukbb/control_snp_id_ssns_afr.txt --freq 'counts' --out $dir/ukbb/bristolAC_afr_nofilter
#time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute --keep $dir/ssns_ctrl/bristol_ssns_sas_allage.txt --extract $dir/ukbb/control_snp_id_ssns_sas.txt --freq 'counts' --out $dir/ukbb/bristolAC_sas_nofilter

# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep $dir/ssns_ctrl/bristol_ssns_eur.txt --extract $dir/ssns_ctrl/clean_control_snp_id.txt --freq 'counts' --out $dir/ssns_ctrl/bristolAC_eur
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep $dir/ssns_ctrl/bristol_ssns_sas.txt --extract $dir/ssns_ctrl/clean_control_snp_id.txt --freq 'counts' --out $dir/ssns_ctrl/bristolAC_sas
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep $dir/ssns_ctrl/bristol_ssns_afr.txt --extract $dir/ssns_ctrl/clean_control_snp_id.txt --freq 'counts' --out $dir/ssns_ctrl/bristolAC_afr
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep $dir/ssns_ctrl/bristol_ssns_allage.txt --make-bed --out $dir/ssns_ctrl/bristol_ssns_allage

##############
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep $dir/ns_ctrl/bristol_ns_eur_allage.txt --extract $dir/ns_ctrl/clean_control_snp_id.txt --freq 'counts' --out $dir/ns_ctrl/bristolAC_eur
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep $dir/ns_ctrl/bristol_ns_sas_allage.txt --extract $dir/ns_ctrl/clean_control_snp_id.txt --freq 'counts' --out $dir/ns_ctrl/bristolAC_sas
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep $dir/ns_ctrl/bristol_ns_afr_allage.txt --extract $dir/ns_ctrl/clean_control_snp_id.txt --freq 'counts' --out $dir/ns_ctrl/bristolAC_afr
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep $dir/ns_ctrl/bristol_ns_allage.txt --make-bed --out $dir/ns_ctrl/bristol_ns_allage

# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep $dir/ns_ctrl/bristol_ns_eur_allage.txt --extract $dir/ns_ctrl/eur/clean_control_snp_id.txt --freq 'counts' --out $dir/ns_ctrl/eur/bristolAC_eur
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep $dir/ns_ctrl/bristol_ns_eur_allage.txt --make-bed --out $dir/ns_ctrl/eur/bristol_ns_eur
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep $dir/ns_ctrl/bristol_ns_sas_allage.txt --extract $dir/ns_ctrl/sas/clean_control_snp_id.txt --freq 'counts' --out $dir/ns_ctrl/sas/bristolAC_sas
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep $dir/ns_ctrl/bristol_ns_sas_allage.txt --make-bed --out $dir/ns_ctrl/sas/bristol_ns_sas

#time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep $dir/ssns_ctrl/bristol_ssns_eur.txt --extract $dir/ssns_ctrl/eur/control_snp_id.txt --freq 'counts' --out $dir/ssns_ctrl/eur/bristolAC_eur
#time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep $dir/ssns_ctrl/bristol_ssns_eur.txt --make-bed --out $dir/ssns_ctrl/eur/bristol_ssns_eur

#time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep $dir/ssns_ctrl/bristol_ssns_sas.txt --extract $dir/ssns_ctrl/sas/control_snp_id.txt --freq 'counts' --out $dir/ssns_ctrl/sas/bristolAC_sas
#time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep $dir/ssns_ctrl/bristol_ssns_sas.txt --make-bed --out $dir/ssns_ctrl/sas/bristol_ssns_sas

time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep $dir/srns_ctrl/bristol_srns_eur.txt --extract $dir/srns_ctrl/clean_control_snp_id.txt --freq 'counts' --out $dir/srns_ctrl/bristolAC_eur
time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep $dir/srns_ctrl/bristol_srns_sas.txt --extract $dir/srns_ctrl/clean_control_snp_id.txt --freq 'counts' --out $dir/srns_ctrl/bristolAC_sas
time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep $dir/srns_ctrl/bristol_srns_afr.txt --extract $dir/srns_ctrl/clean_control_snp_id.txt --freq 'counts' --out $dir/srns_ctrl/bristolAC_afr
time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep $dir/srns_ctrl/bristol_srns_allage.txt --make-bed --out $dir/srns_ctrl/bristol_srns_allage

# ssns snps from all ancestry
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --extract $dir/ssns_ctrl/snp_id_from_table.txt --freq 'counts' --out $dir/ssns_ctrl/bristolAC_ukbb
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep $dir/ssns_ctrl/bristol_ssns_sas_allage.txt --extract $dir/ssns_ctrl/sas/snp_id_from_table.txt --freq 'counts' --out $dir/ssns_ctrl/sas/bristolAC_ukbb
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep $dir/ssns_ctrl/bristol_ssns_afr_allage.txt --extract $dir/ssns_ctrl/afr/snp_id_from_table.txt --freq 'counts' --out $dir/ssns_ctrl/afr/bristolAC_ukbb
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep $dir/ssns_ctrl/bristol_ssns_eur_allage.txt --extract $dir/ssns_ctrl/eur/snp_id_from_table.txt --freq 'counts' --out $dir/ssns_ctrl/eur/bristolAC_ukbb


module unload Plink/2.00a2LM



