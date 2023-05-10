#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=subset
#SBATCH --output=subset_sig.out
#SBATCH --mem=64G
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=tiffany.tu@duke.edu
#SBATCH --mail-type=END,FAIL

# NOTE: 32G sufficed for HO and HGDP, but increased for TGP

module load Plink/2.00a2LM

time plink2 --bfile /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/ssns_tgp_imp_geno --extract suggsig_uniqueID.txt --make-bed --out combined_suggsig
module unload Plink/2.00a2LM



