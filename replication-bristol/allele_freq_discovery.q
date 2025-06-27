#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=ns.allele_freq
#SBATCH --output=AF_ns.out
#SBATCH --mem=64G
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=tiffany.tu@duke.edu
#SBATCH --mail-type=END,FAIL

# NOTE: 32G sufficed for HO and HGDP, but increased for TGP

dir='/datacommons/ochoalab/ssns_gwas/replication/saige_results/bristol/AF'
module load Plink/2.00a2LM

# 
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/imputed/ssns_ctrl/mac20 --keep $dir/ns_ctrl/ctrl_sas_id.txt --extract $dir/ssns_ctrl/ukbb_control_snp_id.txt --freq counts --out $dir/ssns_ctrl/discoveryAC_sas_ukbbsnps
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/imputed/ssns_ctrl/mac20 --keep $dir/ns_ctrl/ctrl_afr_id.txt --extract $dir/ssns_ctrl/ukbb_control_snp_id.txt --freq counts --out $dir/ssns_ctrl/discoveryAC_afr_ukbbsnps
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/imputed/ssns_ctrl/mac20 --keep $dir/ns_ctrl/ctrl_eur_id.txt --extract $dir/ssns_ctrl/ukbb_control_snp_id.txt --freq counts --out $dir/ssns_ctrl/discoveryAC_eur_ukbbsnps
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/imputed/ssns_ctrl/mac20 --keep $dir/ns_ctrl/ctrl_all_id.txt --extract $dir/ssns_ctrl/ukbb_control_snp_id.txt --freq counts --out $dir/ssns_ctrl/discoveryAC_all_ukbbsnps


# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/imputed/mac20 --keep $dir/ns_ctrl/ctrl_sas_id.txt --extract $dir/ns_ctrl/ukbb_control_snp_id.txt --freq counts --out $dir/ns_ctrl/discoveryAC_sas_ukbbsnps
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/imputed/mac20 --keep $dir/ns_ctrl/ctrl_afr_id.txt --extract $dir/ns_ctrl/ukbb_control_snp_id.txt --freq counts --out $dir/ns_ctrl/discoveryAC_afr_ukbbsnps
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/imputed/mac20 --keep $dir/ns_ctrl/ctrl_eur_id.txt --extract $dir/ns_ctrl/ukbb_control_snp_id.txt --freq counts --out $dir/ns_ctrl/discoveryAC_eur_ukbbsnps
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/imputed/mac20 --keep $dir/ns_ctrl/ctrl_all_id.txt --extract $dir/ns_ctrl/ukbb_control_snp_id.txt --freq counts --out $dir/ns_ctrl/discoveryAC_all_ukbbsnps

# ns vs control - eur
#time plink2 --bfile /datacommons/ochoalab/ssns_gwas/imputed/eur/mac20 --keep $dir/ns_ctrl/ctrl_eur_id.txt --extract $dir/ns_ctrl/eur/ukbb_control_snp_id.txt --freq counts --out $dir/ns_ctrl/eur/discoveryAC_eur_ukbbsnps

#ns vs control - sas
#time plink2 --bfile /datacommons/ochoalab/ssns_gwas/imputed/sas/mac20 --keep $dir/ns_ctrl/ctrl_sas_id.txt --extract $dir/ns_ctrl/sas/ukbb_control_snp_id.txt --freq counts --out $dir/ns_ctrl/sas/discoveryAC_sas_ukbbsnps

# ssns vs control - eur
#time plink2 --bfile /datacommons/ochoalab/ssns_gwas/imputed/ssns_ctrl/eur/mac20 --keep $dir/ns_ctrl/ctrl_eur_id.txt --extract $dir/ssns_ctrl/eur/ukbb_control_snp_id.txt --freq counts --out $dir/ssns_ctrl/eur/discoveryAC_eur_ukbbsnps

# ssns vs control - sas
time plink2 --bfile /datacommons/ochoalab/ssns_gwas/imputed/ssns_ctrl/sas/mac20 --keep $dir/ns_ctrl/ctrl_sas_id.txt --extract $dir/ssns_ctrl/sas/ukbb_control_snp_id.txt --freq counts --out $dir/ssns_ctrl/sas/discoveryAC_sas_ukbbsnps

# srns vs control match10
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/imputed/srns_ctrl/mac20 --keep $dir/ns_ctrl/ctrl_sas_id.txt --extract $dir/srns_ctrl/ukbb_control_snp_id.txt --freq counts --out $dir/srns_ctrl/discoveryAC_sas_ukbbsnps
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/imputed/srns_ctrl/mac20 --keep $dir/ns_ctrl/ctrl_afr_id.txt --extract $dir/srns_ctrl/ukbb_control_snp_id.txt --freq counts --out $dir/srns_ctrl/discoveryAC_afr_ukbbsnps
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/imputed/srns_ctrl/mac20 --keep $dir/ns_ctrl/ctrl_eur_id.txt --extract $dir/srns_ctrl/ukbb_control_snp_id.txt --freq counts --out $dir/srns_ctrl/discoveryAC_eur_ukbbsnps

module unload Plink/2.00a2LM



