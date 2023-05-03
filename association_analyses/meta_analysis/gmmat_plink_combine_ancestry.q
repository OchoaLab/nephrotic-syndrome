#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=combine_race_ssns
#SBATCH --output=combine_race_ssns.out
#SBATCH --mem=264G
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=tiffany.tu@duke.edu
#SBATCH --mail-type=END,FAIL

# NOTE: 32G sufficed for HO and HGDP, but increased for TGP

module load Plink/2.00a2LM

time plink2 --bfile /datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/ssns_ctr_mac20  --keep snpid_ssns_ancestry.txt --make-bed --out ssns_ctr_mac20_ancestry
module unload Plink/2.00a2LM



