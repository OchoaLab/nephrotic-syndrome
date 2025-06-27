#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=ns_ctrl-gnomad
#SBATCH --output=processing_ns_ctrl_all.out
#SBATCH --mem=80G
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mail-user=tiffany.tu@duke.edu
#SBATCH --mail-type=END,FAIL

# NOTE: 32G sufficed for HO and HGDP, but increased for TGP

input_file='/datacommons/ochoalab/ssns_gwas/imputed/annotate/locuszoom/ssns_ctrl_all.txt'
output_name='/datacommons/ochoalab/ssns_gwas/replication/curegn/gnomad/outfile/gnomad-genome_ssns_ctrl'

module load R/4.0.3-rhel8
time Rscript 01-list-to-regions.R $input_file snp-list.txt

#https://storage.googleapis.com/gcp-public-data--gnomad/release/4.0/vcf/genomes/gnomad.genomes.v4.0.sites.chr1.vcf.bgz
google_gnomad=https://storage.googleapis.com/gcp-public-data--gnomad/release
module load htslib/1.11
module load bcftools/1.4
module load VCFtools/0.1.15

for chr in {1..22}; do
  echo "chr"$chr
  time tabix -h $google_gnomad/4.0/vcf/genomes/gnomad.genomes.v4.0.sites.chr$chr.vcf.bgz -R snp-list.txt > gnomad-4-genome-chr$chr.vcf
  # have to compress and index for merging
  bgzip gnomad-4-genome-chr$chr.vcf
  vcf-sort gnomad-4-genome-chr$chr.vcf.gz > gnomad-4-genome-sort-chr$chr.vcf.gz
  bcftools view -Oz -o gnomad-4-genome-sort-out-chr$chr.vcf.gz gnomad-4-genome-sort-chr$chr.vcf.gz
  bcftools index gnomad-4-genome-sort-out-chr$chr.vcf.gz
done
# 0m9.030s, 0m7.024s, 0m4.446s, 0m5.121s, 0m5.480s, 0m12.747s, ...
# now merge them for simplicity
bcftools merge gnomad-4-genome-sort-out-chr{?,??}.vcf.gz > ${output_name}.vcf
#bcftools merge gnomad-3-genome-sort-out-chr1.vcf.gz gnomad-3-genome-sort-out-chr2.vcf.gz > ${output_name}.vcf

#clean up intermediates
rm gnomad-4-genome-chr*.vcf.gz{,.csi}
rm *tbi
rm gnomad-4-genome-sort-*
  
#awkwardly add header line, for easy parsing in R... (bcftools doesn't add it, frustrating)
echo -e "CHROM\tPOS\tID\tREF\tALT\tAC_afr\tAN_afr\tAC_nfe\tAN_nfe\tAC_sas\tAN_sas\tAC_eas\tAN_eas" > ${output_name}.txt
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AC_afr\t%AN_afr\t%AC_nfe\t%AN_nfe\t%AC_sas\t%AN_sas\t%AC_eas\t%AN_eas\n' ${output_name}.vcf >> ${output_name}.txt
time Rscript 02-cleanup.R $input_file ${output_name}.txt

