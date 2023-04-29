#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=clump_ssns
#SBATCH --output=clump.out
#SBATCH --mem=164G
#SBATCH --ntasks-per-node=10
#SBATCH --mail-user=tiffany.tu@duke.edu
#SBATCH --mail-type=END,FAIL

# NOTE: 32G sufficed for HO and HGDP, but increased for TGP

module load Plink/1.90

time plink --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ssns_ctr_new/bristol_ssns_allage --clump /datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ssns_ctr_new/replication_unclump_ssns_allage.txt --clump-r2 0.3 --clump-kb 6000 --clump-p1 1 --clump-p2 1 --out /datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ssns_ctr_new/clump_bristol_6205
# time plink --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ssns_ctr_all/bristol_ssns_35 --clump /datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ssns_ctr_all/replication_forLDclump_ssns_all35.txt --clump-r2 0.3 --clump-kb 6000 --clump-p1 1 --clump-p2 1 --out /datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ssns_ctr_all/clump_bristol_2913_age35
# time plink --bfile /datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ssns_ctr_all/bristol_ssns_21 --clump /datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ssns_ctr_all/replication_forLDclump_ssns_all21.txt --clump-r2 0.3 --clump-kb 6000 --clump-p1 1 --clump-p2 1 --out /datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/ssns_ctr_all/clump_bristol_2913_age21

module unload Plink/1.90
