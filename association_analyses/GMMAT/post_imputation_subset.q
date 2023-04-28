#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=split
#SBATCH --output=split.out
#SBATCH --mem=264G
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=tiffany.tu@duke.edu
#SBATCH --mail-type=END,FAIL

# NOTE: 32G sufficed for HO and HGDP, but increased for TGP

module load Plink/2.00a2LM

time plink2 --bfile /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/ssns_tgp_impute --keep /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/ssns_ctr.txt --make-bed --out ssns_ctr
time plink2 --bfile /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/ssns_tgp_impute --keep /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/srns_ctr.txt --make-bed --out srns_ctr
time plink2 --bfile /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/ssns_tgp_impute --keep /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/ssns_srns.txt --make-bed --out ssns_srns

module unload Plink/2.00a2LM



