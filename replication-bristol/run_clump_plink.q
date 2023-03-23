#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=clump3
#SBATCH --output=clump3.out
#SBATCH --mem=164G
#SBATCH --ntasks-per-node=10
#SBATCH --mail-user=tiffany.tu@duke.edu
#SBATCH --mail-type=END,FAIL

# NOTE: 32G sufficed for HO and HGDP, but increased for TGP

module load Plink/1.90

time plink --bfile /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/GMMAT/ssns_ctr_hwe_mac20 --clump /datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ssns_ctr_all/replication_test_ssns_2913_1441sig.txt --clump-r2 0.3 --clump-kb 6000 --clump-p1 1 --clump-p2 1 --out clump_bristol_all_r2_003_bp6_p1p2_1441sig.txt

module unload Plink/1.90
