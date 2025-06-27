#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=clump_srns
#SBATCH --output=clump.out
#SBATCH --mem=164G
#SBATCH --ntasks-per-node=10
#SBATCH --mail-user=tiffany.tu@duke.edu
#SBATCH --mail-type=END,FAIL

# NOTE: 32G sufficed for HO and HGDP, but increased for TGP
dir='/datacommons/ochoalab/ssns_gwas/replication/saige_results/bristol/AF'
module load Plink/1.90

time plink --bfile $dir/ssns_ctrl/bristol_ssns_allage --clump $dir/ssns_ctrl/replication_unclump_bristol_OR.txt --clump-r2 0.3 --clump-kb 6000 --clump-p1 1 --clump-p2 1 --out $dir/ssns_ctrl/clump_bristol_4848_OR
#time plink --bfile $dir/ssns_ctrl/age35/bristol_ssns_35 --clump $dir/ssns_ctrl/age35/replication_unclump_bristol_age35.txt --clump-r2 0.3 --clump-kb 6000 --clump-p1 1 --clump-p2 1 --out $dir/ssns_ctrl/age35/clump_bristol_2278
#time plink --bfile $dir/ssns_ctrl/age21/bristol_ssns_21 --clump $dir/ssns_ctrl/age21/replication_unclump_bristol_age21.txt --clump-r2 0.3 --clump-kb 6000 --clump-p1 1 --clump-p2 1 --out $dir/ssns_ctrl/age21/clump_bristol_2278
#time plink --bfile $dir/ssns_ctrl/eur/bristol_ssns_eur --clump $dir/ssns_ctrl/eur/replication_unclump_bristol.txt --clump-r2 0.3 --clump-kb 6000 --clump-p1 1 --clump-p2 1 --out $dir/ssns_ctrl/eur/clump_bristol_1899
#time plink --bfile $dir/ssns_ctrl/sas/bristol_ssns_sas --clump $dir/ssns_ctrl/sas/replication_unclump_bristol.txt --clump-r2 0.3 --clump-kb 6000 --clump-p1 1 --clump-p2 1 --out $dir/ssns_ctrl/sas/clump_bristol_1054

#time plink --bfile $dir/ns_ctrl/bristol_ns_allage --clump $dir/ns_ctrl/replication_unclump_bristol.txt --clump-r2 0.3 --clump-kb 6000 --clump-p1 1 --clump-p2 1 --out $dir/ns_ctrl/clump_bristol_4290
#time plink --bfile $dir/ns_ctrl/eur/bristol_ns_eur --clump $dir/ns_ctrl/eur/replication_unclump_bristol.txt --clump-r2 0.3 --clump-kb 6000 --clump-p1 1 --clump-p2 1 --out $dir/ns_ctrl/eur/clump_bristol_1164
#time plink --bfile $dir/ns_ctrl/sas/bristol_ns_sas --clump $dir/ns_ctrl/sas/replication_unclump_bristol.txt --clump-r2 0.3 --clump-kb 6000 --clump-p1 1 --clump-p2 1 --out $dir/ns_ctrl/sas/clump_bristol_928

#time plink --bfile $dir/srns_ctrl/bristol_srns_allage --clump $dir/srns_ctrl/replication_unclump_bristol.txt --clump-r2 0.3 --clump-kb 6000 --clump-p1 1 --clump-p2 1 --out $dir/srns_ctrl/clump_bristol_81

module unload Plink/1.90
