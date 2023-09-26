#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=AF_casecontrol
#SBATCH --output=AF_casecontrol.out
#SBATCH --mem=64G
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=tiffany.tu@duke.edu
#SBATCH --mail-type=END,FAIL

# NOTE: 32G sufficed for HO and HGDP, but increased for TGP

module load Plink/1.90  

# ssns vs ctrl
#time plink --bfile /datacommons/ochoalab/ssns_gwas/imputed/ssns_ctrl/mac20 --pheno /datacommons/ochoalab/ssns_gwas/imputed/patient-data_plink.txt.gz --pheno-name ssns_ctrl --keep-allele-order --extract snp_id/ssns_ctrl_all_snpid.txt --freq case-control --nonfounders --allow-no-sex --out ssns_ctrl_all_casecontrol

# sas
#time plink --bfile /datacommons/ochoalab/ssns_gwas/imputed/ssns_ctrl/sas/mac20 --pheno /datacommons/ochoalab/ssns_gwas/imputed/patient-data_plink.txt.gz --pheno-name ssns_ctrl --keep-allele-order --extract snp_id/ssns_ctrl_sas_snpid.txt --freq case-control --nonfounders --allow-no-sex --out ssns_ctrl_sas_casecontrol
# afr
#time plink --bfile /datacommons/ochoalab/ssns_gwas/imputed/ssns_ctrl/afr/mac20 --pheno /datacommons/ochoalab/ssns_gwas/imputed/patient-data_plink.txt.gz --pheno-name ssns_ctrl --keep-allele-order --extract snp_id/ssns_ctrl_afr_snpid.txt --freq case-control --nonfounders --allow-no-sex --out ssns_ctrl_afr_casecontrol
# eur
#time plink --bfile /datacommons/ochoalab/ssns_gwas/imputed/ssns_ctrl/eur/mac20 --pheno /datacommons/ochoalab/ssns_gwas/imputed/patient-data_plink.txt.gz --pheno-name ssns_ctrl --keep-allele-order --extract snp_id/ssns_ctrl_eur_snpid.txt --freq case-control --nonfounders --allow-no-sex --out ssns_ctrl_eur_casecontrol

# ns vs ctrl
#time plink --bfile /datacommons/ochoalab/ssns_gwas/imputed/mac20 --pheno /datacommons/ochoalab/ssns_gwas/imputed/patient-data_plink.txt.gz --pheno-name ns_ctrl --keep-allele-order --extract snp_id/ns_ctrl_all_snpid.txt --freq case-control --nonfounders --allow-no-sex --out ns_ctrl_all_casecontrol
# eur
time plink --bfile /datacommons/ochoalab/ssns_gwas/eur/imputed/mac20 --pheno /datacommons/ochoalab/ssns_gwas/imputed/patient-data_plink.txt.gz --pheno-name ns_ctrl --keep-allele-order --extract snp_id/ns_ctrl_eur_snpid.txt --freq case-control --nonfounders --allow-no-sex --out ns_ctrl_eur_casecontrol

module unload Plink/1.90  



