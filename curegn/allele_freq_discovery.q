#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=ns.allele_freq
#SBATCH --output=AF.out
#SBATCH --mem=64G
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mail-user=tiffany.tu@duke.edu
#SBATCH --mail-type=END,FAIL

# NOTE: 32G sufficed for HO and HGDP, but increased for TGP

dir='/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF'
module load Plink/2.00a2LM

# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/imputed/ssns_ctrl/mac20 --keep $dir/ssns_ctrl/array_ctrl_sas_id.txt --extract $dir/gnomad/outfile_new/snp-id/ssns_ctrl_all_snp_id.txt --freq counts --out $dir/ssns_ctrl/array_discoveryAC_sas_gnomadsnps
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/imputed/ssns_ctrl/mac20 --keep $dir/ssns_ctrl/array_ctrl_afr_id.txt --extract $dir/gnomad/outfile_new/snp-id/ssns_ctrl_all_snp_id.txt --freq counts --out $dir/ssns_ctrl/array_discoveryAC_afr_gnomadsnps
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/imputed/ssns_ctrl/mac20 --keep $dir/ssns_ctrl/array_ctrl_eur_id.txt --extract $dir/gnomad/outfile_new/snp-id/ssns_ctrl_all_snp_id.txt --freq counts --out $dir/ssns_ctrl/array_discoveryAC_eur_gnomadsnps
# time plink2 --bfile /datacommons/ochoalab/ssns_gwas/imputed/ssns_ctrl/mac20 --keep $dir/ssns_ctrl/array_ctrl_all_id.txt --extract $dir/gnomad/outfile_new/snp-id/ssns_ctrl_all_snp_id.txt --freq counts --out $dir/ssns_ctrl/array_discoveryAC_all_gnomadsnps
# # 
# 

time plink2 --bfile /datacommons/ochoalab/ssns_gwas/imputed/mac20 --keep $dir/ns_ctrl/array_ctrl_sas_id.txt --extract $dir/gnomad/outfile_new/snp-id/ns_ctrl_all_snp_id.txt --freq counts --out $dir/ns_ctrl/array_discoveryAC_sas_gnomadsnps
time plink2 --bfile /datacommons/ochoalab/ssns_gwas/imputed/mac20 --keep $dir/ns_ctrl/array_ctrl_afr_id.txt --extract $dir/gnomad/outfile_new/snp-id/ns_ctrl_all_snp_id.txt --freq counts --out $dir/ns_ctrl/array_discoveryAC_afr_gnomadsnps
time plink2 --bfile /datacommons/ochoalab/ssns_gwas/imputed/mac20 --keep $dir/ns_ctrl/array_ctrl_eur_id.txt --extract $dir/gnomad/outfile_new/snp-id/ns_ctrl_all_snp_id.txt --freq counts --out $dir/ns_ctrl/array_discoveryAC_eur_gnomadsnps
time plink2 --bfile /datacommons/ochoalab/ssns_gwas/imputed/mac20 --keep $dir/ns_ctrl/array_ctrl_all_id.txt --extract $dir/gnomad/outfile_new/snp-id/ns_ctrl_all_snp_id.txt --freq counts --out $dir/ns_ctrl/array_discoveryAC_all_gnomadsnps

module unload Plink/2.00a2LM



