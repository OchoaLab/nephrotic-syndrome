#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=allele_freq
#SBATCH --output=allele_freq_ancestry.out
#SBATCH --mem=64G
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=tiffany.tu@duke.edu
#SBATCH --mail-type=END,FAIL

# NOTE: 32G sufficed for HO and HGDP, but increased for TGP

module load Plink/2.00a2LM

time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep bristol_ssns_euro_allage.txt --extract control_snp_id.txt --freq 'counts' --out b_ssns_euro_allage
time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep bristol_ssns_sas_allage.txt --extract control_snp_id.txt --freq 'counts' --out b_ssns_sas_allage
time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep bristol_ssns_af_allage.txt --extract control_snp_id.txt --freq 'counts' --out b_ssns_af_allage
time plink2 --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep bristol_ssns_allage.txt --make-bed --out bristol_ssns_allage

module unload Plink/2.00a2LM



