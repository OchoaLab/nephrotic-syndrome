#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=ssns_grm
#SBATCH --output=ssns_grm.out
#SBATCH --mem=220G
#SBATCH --ntasks-per-node=10
#SBATCH --mail-user=tiffany.tu@duke.edu
#SBATCH --mail-type=END,FAIL

# NOTE: 32G sufficed for HO and HGDP, but increased for TGP


time gcta --bfile /datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/ssns_ctr_mac20 --autosome --make-grm --out ssns_ctr_mac20
time gcta --grm ssns_ctr_mac20 --pca 10 --out ssns_ctr_mac20 