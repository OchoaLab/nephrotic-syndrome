#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=processing
#SBATCH --output=processing.out
#SBATCH --mem=64G
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=tiffany.tu@duke.edu
#SBATCH --mail-type=END,FAIL

# NOTE: 32G sufficed for HO and HGDP, but increased for TGP

module load R/4.0.3-rhel8
time Rscript 01-list-to-regions.R /datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_suggsig.txt snp-list.txt


google_gnomad=https://storage.googleapis.com/gcp-public-data--gnomad/release
module load htslib/1.11
module load bcftools/1.4
module load VCFtools/0.1.15
for chr in {1..22}; do
  time tabix -h $google_gnomad/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr$chr.vcf.bgz -R snp-list.txt > gnomad-3-genome-chr$chr.vcf
  # have to compress and index for merging
  bgzip gnomad-3-genome-chr$chr.vcf
  bcftools sort -Oz gnomad-3-genome-chr$chr.vcf.gz -o gnomad-3-genome-chr$chr.sort.vcf.gz
  vcf-sort gnomad-3-genome-chr$chr.vcf.gz > gnomad-3-genome-sort-chr$chr.vcf.gz
  bcftools view -Oz -o gnomad-3-genome-sort-out-chr$chr.vcf.gz gnomad-3-genome-sort-chr$chr.vcf.gz
  bcftools index gnomad-3-genome-sort-out-chr$chr.vcf.gz
done
# 0m9.030s, 0m7.024s, 0m4.446s, 0m5.121s, 0m5.480s, 0m12.747s, ...
# now merge them for simplicity
bcftools merge gnomad-3-genome-sort-out-chr{?,??}.vcf.gz > gnomad-genome_ssns_control.vcf
#clean up intermediates
rm gnomad-3-genome-chr*.vcf.gz{,.csi}

#module unload bcftools/1.4


#awkwardly add header line, for easy parsing in R... (bcftools doesn't add it, frustrating)
echo -e "CHROM\tPOS\tID\tREF\tALT\tAC_nfe\tAN_nfe\tAC_afr\tAN_afr\tAC_sas\tAN_sas" > gnomad-genome_ssns_control.txt
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AC_nfe\t%AN_nfe\t%AC_afr\t%AN_afr\t%AC_sas\t%AN_sas\n' gnomad-genome_ssns_control.vcf >> gnomad-genome_ssns_control.txt
time Rscript 02-cleanup.R /datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/ssns_ctr_suggsig.txt gnomad-genome_ssns_control.txt

