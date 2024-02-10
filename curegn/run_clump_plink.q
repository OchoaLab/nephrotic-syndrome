#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=clump
#SBATCH --output=clump.out
#SBATCH --mem=64G
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=tiffany.tu@duke.edu
#SBATCH --mail-type=END,FAIL

# NOTE: 32G sufficed for HO and HGDP, but increased for TGP

module load Plink/1.90

#time plink --bfile /datacommons/ochoalab/curegn/imputed/clean/curegn_ns_mac20_clean --clump /datacommons/ochoalab/ssns_gwas/replication/curegn/AF/ns/replication_unclump_curegn_snps.txt --clump-r2 0.3 --clump-kb 6000 --clump-p1 1 --clump-p2 1 --out /datacommons/ochoalab/ssns_gwas/replication/curegn/AF/ns/clump_curegn_5629
#time plink --bfile /datacommons/ochoalab/curegn/imputed/clean/curegn_ns_eur_mac20_clean --clump /datacommons/ochoalab/ssns_gwas/replication/curegn/AF/ns/eur/replication_unclump_curegn_snps_eur.txt --clump-r2 0.3 --clump-kb 6000 --clump-p1 1 --clump-p2 1 --out /datacommons/ochoalab/ssns_gwas/replication/curegn/AF/ns/eur/clump_curegn_1437_eur
#time plink --bfile /datacommons/ochoalab/curegn/imputed/clean/curegn_ns_afr_mac20_clean --clump /datacommons/ochoalab/ssns_gwas/replication/curegn/AF/ns/afr/replication_unclump_curegn_snps_afr.txt --clump-r2 0.3 --clump-kb 6000 --clump-p1 1 --clump-p2 1 --out /datacommons/ochoalab/ssns_gwas/replication/curegn/AF/ns/afr/clump_curegn_1201_afr

#time plink --bfile /datacommons/ochoalab/curegn/imputed/curegn_ns_mac20 --keep /datacommons/ochoalab/curegn/imputed/curegn_ssns_cases_id.txt --clump /datacommons/ochoalab/ssns_gwas/replication/curegn/AF/ssns/replication_unclump_curegn_snps.txt --clump-r2 0.3 --clump-kb 6000 --clump-p1 1 --clump-p2 1 --out /datacommons/ochoalab/ssns_gwas/replication/curegn/AF/ssns/clump_curegn_6605
#time plink --bfile /datacommons/ochoalab/curegn/imputed/curegn_ns_mac20 --keep /datacommons/ochoalab/curegn/imputed/curegn_ssns_eur.txt --clump /datacommons/ochoalab/ssns_gwas/replication/curegn/AF/ssns/eur/replication_unclump_curegn_snps_eur.txt --clump-r2 0.3 --clump-kb 6000 --clump-p1 1 --clump-p2 1 --out /datacommons/ochoalab/ssns_gwas/replication/curegn/AF/ssns/eur/clump_curegn_2347
time plink --bfile /datacommons/ochoalab/curegn/imputed/curegn_ns_mac20 --keep /datacommons/ochoalab/curegn/imputed/curegn_ssns_afr.txt --clump /datacommons/ochoalab/ssns_gwas/replication/curegn/AF/ssns/afr/replication_unclump_curegn_snps_afr.txt --clump-r2 0.3 --clump-kb 6000 --clump-p1 1 --clump-p2 1 --out /datacommons/ochoalab/ssns_gwas/replication/curegn/AF/ssns/afr/clump_curegn_1291

module unload Plink/1.90
