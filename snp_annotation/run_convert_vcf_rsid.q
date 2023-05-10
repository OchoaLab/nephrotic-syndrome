#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=vcf_rsid
#SBATCH --output=vcf_rsid.out
#SBATCH --mem=264G
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=tiffany.tu@duke.edu
#SBATCH --mail-type=END,FAIL

# NOTE: 32G sufficed for HO and HGDP, but increased for TGP

module load Plink/2.00a2LM
module load bcftools/1.4

time plink2 --bfile combined_suggsig_recode --keep-allele-order --recode vcf --out combined_suggsig_recode
time bcftools view -O z -o combined_suggsig_recode.vcf.gz combined_suggsig_recode.vcf
time bcftools index combined_suggsig_recode.vcf.gz --tbi
time bcftools annotate -a /datacommons/ochoalab/ssns_gwas/rsid/00-All.vcf.gz -c ID combined_suggsig_recode.vcf.gz -o rsid_combined_suggsig_recode.vcf.gz

module unload Plink/2.00a2LM
module load bcftools/1.4


