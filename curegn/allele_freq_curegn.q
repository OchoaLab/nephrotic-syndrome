#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=allele_freq
#SBATCH --output=allele_freq.out
#SBATCH --mem=64G
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mail-user=tiffany.tu@duke.edu
#SBATCH --mail-type=END,FAIL

# NOTE: 32G sufficed for HO and HGDP, but increased for TGP
dir='/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF'

module load Plink/2.00a2LM

# time plink2 --bfile /datacommons/ochoalab/curegn/imputed/clean/curegn_ns_mac20_clean --keep /datacommons/ochoalab/curegn/imputed/curegn_ns_afr.txt --extract /datacommons/ochoalab/ssns_gwas/replication/curegn/AF/ns/gnomad_snp_id.txt --freq 'counts' --out /datacommons/ochoalab/ssns_gwas/replication/curegn/AF/ns/curegnAC_afr
# time plink2 --bfile /datacommons/ochoalab/curegn/imputed/clean/curegn_ns_mac20_clean --keep /datacommons/ochoalab/curegn/imputed/curegn_ns_eur.txt --extract /datacommons/ochoalab/ssns_gwas/replication/curegn/AF/ns/gnomad_snp_id.txt --freq 'counts' --out /datacommons/ochoalab/ssns_gwas/replication/curegn/AF/ns/curegnAC_eur
# time plink2 --bfile /datacommons/ochoalab/curegn/imputed/clean/curegn_ns_mac20_clean --keep /datacommons/ochoalab/curegn/imputed/curegn_ns_sas.txt --extract /datacommons/ochoalab/ssns_gwas/replication/curegn/AF/ns/gnomad_snp_id.txt --freq 'counts' --out /datacommons/ochoalab/ssns_gwas/replication/curegn/AF/ns/curegnAC_sas
# time plink2 --bfile /datacommons/ochoalab/curegn/imputed/clean/curegn_ns_mac20_clean --keep /datacommons/ochoalab/curegn/imputed/curegn_ns_eas.txt --extract /datacommons/ochoalab/ssns_gwas/replication/curegn/AF/ns/gnomad_snp_id.txt --freq 'counts' --out /datacommons/ochoalab/ssns_gwas/replication/curegn/AF/ns/curegnAC_eas

# eur subanalysis
#time plink2 --bfile /datacommons/ochoalab/curegn/imputed/clean/curegn_ns_eur_mac20_clean --keep /datacommons/ochoalab/curegn/imputed/curegn_ns_eur.txt --extract /datacommons/ochoalab/ssns_gwas/replication/curegn/AF/ns/eur/gnomad_snp_id_eur.txt --freq 'counts' --out /datacommons/ochoalab/ssns_gwas/replication/curegn/AF/ns/eur/curegnAC_eur_sub
#time plink2 --bfile /datacommons/ochoalab/curegn/imputed/curegn_ns_mac20 --keep /datacommons/ochoalab/curegn/imputed/curegn_ssns_eur.txt --extract /datacommons/ochoalab/ssns_gwas/replication/curegn/AF/ssns/eur/gnomad_snp_id_eur.txt --freq 'counts' --out /datacommons/ochoalab/ssns_gwas/replication/curegn/AF/ssns/eur/curegnAC_eur_sub

# afr subanalysis
#time plink2 --bfile /datacommons/ochoalab/curegn/imputed/clean/curegn_ns_afr_mac20_clean --keep /datacommons/ochoalab/curegn/imputed/curegn_ns_afr.txt --extract /datacommons/ochoalab/ssns_gwas/replication/curegn/AF/ns/afr/gnomad_snp_id_afr.txt --freq 'counts' --out /datacommons/ochoalab/ssns_gwas/replication/curegn/AF/ns/afr/curegnAC_afr_sub
#time plink2 --bfile /datacommons/ochoalab/curegn/imputed/curegn_ns_mac20 --keep /datacommons/ochoalab/curegn/imputed/curegn_ssns_afr.txt --extract /datacommons/ochoalab/ssns_gwas/replication/curegn/AF/ssns/afr/gnomad_snp_id_afr.txt --freq 'counts' --out /datacommons/ochoalab/ssns_gwas/replication/curegn/AF/ssns/afr/curegnAC_afr_sub

## ssns
#time plink2 --bfile /datacommons/ochoalab/curegn/imputed/curegn_ns_mac20 --keep /datacommons/ochoalab/curegn/imputed/curegn_ssns_afr.txt --extract /datacommons/ochoalab/ssns_gwas/replication/curegn/AF/ssns/gnomad_snp_id.txt --freq 'counts' --out /datacommons/ochoalab/ssns_gwas/replication/curegn/AF/ssns/curegnAC_afr
#time plink2 --bfile /datacommons/ochoalab/curegn/imputed/curegn_ns_mac20 --keep /datacommons/ochoalab/curegn/imputed/curegn_ssns_eur.txt --extract /datacommons/ochoalab/ssns_gwas/replication/curegn/AF/ssns/gnomad_snp_id.txt --freq 'counts' --out /datacommons/ochoalab/ssns_gwas/replication/curegn/AF/ssns/curegnAC_eur



time plink2 --bfile /datacommons/ochoalab/curegn/merge_tgp/curegn_tgp_merge --keep /datacommons/ochoalab/curegn/imputed/curegn_ssns_afr.txt --extract AC_selected_snp_id.txt --freq 'counts' --out AC_selected_snp
#time plink2 --bfile /datacommons/ochoalab/curegn/imputed/curegn_NS_afr --keep /datacommons/ochoalab/curegn/imputed/curegn_ns_cases_id.txt --extract AC_selected_snp_id.txt --freq 'counts' --out AC_selected_snp

module unload Plink/2.00a2LM



