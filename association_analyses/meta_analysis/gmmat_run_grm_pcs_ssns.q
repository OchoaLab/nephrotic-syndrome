#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=grm_ssns
#SBATCH --output=grm_ssns.out
#SBATCH --mem=220G
#SBATCH --ntasks-per-node=10
#SBATCH --mail-user=tiffany.tu@duke.edu
#SBATCH --mail-type=END,FAIL

# NOTE: 32G sufficed for HO and HGDP, but increased for TGP



time gcta --bfile ssns_ctr_mac20_ancestry --autosome --make-grm --out ssns_ctr_mac20_ancestry
time gcta --grm ssns_ctr_mac20_ancestry --pca 10 --out ssns_ctr_mac20_ancestry 