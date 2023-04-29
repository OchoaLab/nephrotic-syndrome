#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=sub.allele_freq
#SBATCH --output=allele_freq_disc.out
#SBATCH --mem=64G
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=tiffany.tu@duke.edu
#SBATCH --mail-type=END,FAIL

# NOTE: 32G sufficed for HO and HGDP, but increased for TGP

module load Plink/2.00a2LM

time plink2 --bfile /datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/ssns_ctr_mac20 --keep /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/ssns_merge_black_control.txt --extract gnomad_snp_id.txt --freq counts --out ssns_ctr_afr_gnomadsnps
time plink2 --bfile /datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/ssns_ctr_mac20 --keep /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/ssns_merge_white_control.txt --extract gnomad_snp_id.txt --freq counts --out ssns_ctr_euro_gnomadsnps
time plink2 --bfile /datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/ssns_ctr_mac20 --keep /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/ssns_merge_asian_control.txt --extract gnomad_snp_id.txt --freq counts --out ssns_ctr_sas_gnomadsnps

module unload Plink/2.00a2LM



