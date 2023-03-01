####################
### DEMOGRAPHICS ###
####################

# basic demographics report, to see what sample sizes we're dealing with for subanalyses
# reads: pheno_Bristol.csv
# creates: age.pdf , rest are STDOUT table reports manually incorporated into a presentation
Rscript 00-bristol-demographics.R


##############
### GNOMAD ###
##############

# NOTE: ran 3 versions:
# - gnomad v3 genomes: good one with largest coverage and sample sizes
# - gnomad v2 genomes and exomes: both ignored in the end, but commands stay here for reference

# convert loci to replicate into simple chr/pos table in right format (chr prefixes, no header)
time Rscript 01-list-to-regions.R 2023-02-20_replication_gene_list.csv snp-list.txt
# 0m0.897s duke-ochoa

# get gnomad info for those SNPs!

# shared variables:
# this is the Google option, there's also Amazon and Microsoft (if there's problems in the future, can switch)
google_gnomad=https://storage.googleapis.com/gcp-public-data--gnomad/release

# v2 exomes, lifted over
time tabix -h $google_gnomad/2.1.1/liftover_grch38/vcf/exomes/gnomad.exomes.r2.1.1.sites.liftover_grch38.vcf.bgz -R snp-list.txt > gnomad-2-exome.vcf
# 0m16.388s dell-xps
# sadly only 5 hits, maybe because it's exome...

# v2 genomes, lifted over
time tabix -h $google_gnomad/2.1.1/liftover_grch38/vcf/genomes/gnomad.genomes.r2.1.1.sites.liftover_grch38.vcf.bgz -R snp-list.txt > gnomad-2-genome.vcf
# 1m0.838s dell-xps

# v3 genomes
# here we've got to query each chromosome :( no merged version
for chr in {1..22}; do
    time tabix -h $google_gnomad/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr$chr.vcf.bgz -R snp-list.txt > gnomad-3-genome-chr$chr.vcf
    # have to compress and index for merging
    bgzip gnomad-3-genome-chr$chr.vcf
    bcftools index gnomad-3-genome-chr$chr.vcf.gz
done
# 0m9.030s, 0m7.024s, 0m4.446s, 0m5.121s, 0m5.480s, 0m12.747s, ...
# now merge them for simplicity
bcftools merge gnomad-3-genome-chr{?,??}.vcf.gz > gnomad-3-genome.vcf
# clean up intermediates
rm gnomad-3-genome-chr*.vcf.gz{,.csi}

# numbers of SNPs
wc -l snp-list.txt                  # 106
grep -c -v '^#' gnomad-2-exome.vcf  # 5
grep -c -v '^#' gnomad-2-genome.vcf # 140 # contains multiple hits
grep -c -v '^#' gnomad-3-genome.vcf # 188 # even more!

# also clean up temporary indexes that were automatically downloaded by tabix (all gnomad versions covered)
rm gnomad.*.tbi
# this is no longer needed, redundant
rm snp-list.txt
# move other extras to unpub, they won't be on github
mkdir unpub
mv gnomad*.vcf unpub/

# extra hits suggests alt/ref mismatches, and potentially overlaps that I'm not interested in
# look at bcftools -T, and does tabix have something similar?  Or maybe we prefer this format?

# conclusion:
# - since exomes are practically absent, maybe just stick with genomes?  In that case, might as well stick with newest one only (v3)

# gnomad notation notes
# inferred by reading INFO tags, figuring out patterns, and checked against link below for completeness
# https://gnomad.broadinstitute.org/help/what-populations-are-represented-in-the-gnomad-data

# ancestries:
# - matching our data:              = race/ethnicity labels in Bristol
#   - nfe: Non-Finnish European     = White
#   - afr: African/African-American = Black
#   - sas: South Asian              = Asian
# - Though plausible, Tiffany's admixture analysis identified none that don't just fit in the previous three categories
#   - eas: East Asian               = Asian too?
#   - amr: Latinos                  = Mixed, Unknown
#   - oth: Other                    = Mixed, Unknown
# - unlikely/ignore:
#   - ami: Amish
#   - fin: Finnish
#   - asj: Ashkenazi Jewish
#   - mid: Middle Eastern

# interesting subsets (all self-explanatory, though "overall" is only case we really want):
# - (overall, no prefix)
# - controls_and_biobanks
# - non_cancer
# - non_neuro
# - non_topmed
# - non_v2

# lastly, all have XX and XY subsets, but again we don't need, except to potentially "control" for sex this way
# - confirmed (for NFE ancestry) that there are no cases with missing sex

# There aren't overall AC/AN in this data!

# turn it into a way nicer table
# version with all ancestries, no sex stratification
# awkwardly add header line, for easy parsing in R... (bcftools doesn't add it, frustrating)
echo -e "CHROM\tPOS\tID\tREF\tALT\tAC_nfe\tAN_nfe\tAC_afr\tAN_afr\tAC_sas\tAN_sas" > gnomad-3-genome.txt
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AC_nfe\t%AN_nfe\t%AC_afr\t%AN_afr\t%AC_sas\t%AN_sas\n' unpub/gnomad-3-genome.vcf >> gnomad-3-genome.txt

# version that included more populations, probably unnecessarily
# echo -e "CHROM\tPOS\tID\tREF\tALT\tAC_nfe\tAN_nfe\tAC_afr\tAN_afr\tAC_sas\tAN_sas\tAC_eas\tAN_eas\tAC_amr\tAN_amr\tAC_oth\tAN_oth" > gnomad-3-genome.txt
# bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AC_nfe\t%AN_nfe\t%AC_afr\t%AN_afr\t%AC_sas\t%AN_sas\t%AC_eas\t%AN_eas\t%AC_amr\t%AN_amr\t%AC_oth\t%AN_oth\n' gnomad-3-genome.vcf >> gnomad-3-genome.txt

# # testing sex stratification, NFE ancestry only
# bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AC_nfe\t%AN_nfe\t%AC_nfe_XX\t%AN_nfe_XX\t%AC_nfe_XY\t%AN_nfe_XY\n' gnomad-3-genome.vcf > gnomad-3-genome.txt

# clean up the data using this script
# most of the cleanup is removing unmatched rows in gnomad output, most of which were junk
# only one query SNP was not found:
time Rscript 02-cleanup.R 2023-02-20_replication_gene_list.csv gnomad-3-genome.txt
# 0m1.003s duke-ochoa
# Removing loci unmatched in gnomad:
#   CHR       POS rsid         REF   ALT   chrpos      
# 1 chr4  9172304 rs1319121362 T     C     chr4:9172304

###################
### GNOMAD NULL ###
###################

# minimal repeat of analysis on random SNPs, to assess calibration of test

# download a large list of random SNPs to compare to, as controls for inflation analysis (not in repo)
scp $dcc:/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF/gmmat_ssns_ctr_all_random1000.txt .

# create regions for tabix from this data
time Rscript 01-list-to-regions.R gmmat_ssns_ctr_all_random1000.txt snp-list-null.txt

# query gnomad v3 only!
for chr in {1..22}; do
    time tabix -h $google_gnomad/3.1.2/vcf/genomes/gnomad.genomes.v3.1.2.sites.chr$chr.vcf.bgz -R snp-list-null.txt > gnomad-3-genome-null-chr$chr.vcf
    # have to compress and index for merging
    bgzip gnomad-3-genome-null-chr$chr.vcf
    bcftools index gnomad-3-genome-null-chr$chr.vcf.gz
done
# 0m58.218s, 0m49.737s, 0m52.019s, 0m52.127s, 0m35.224s, ...

# now merge them for simplicity
bcftools merge gnomad-3-genome-null-chr{?,??}.vcf.gz > gnomad-3-genome-null.vcf
# clean up intermediates
rm gnomad-3-genome-null-chr*.vcf.gz{,.csi}

# numbers of SNPs
wc -l snp-list-null.txt                  # 1000
grep -c -v '^#' gnomad-3-genome-null.vcf # 1384

# also clean up temporary indexes that were automatically downloaded by tabix (all gnomad versions covered)
rm gnomad.*.tbi
# this is no longer needed, redundant
rm snp-list-null.txt
# move other extras to unpub, they won't be on github
mv gnomad-3-genome-null.vcf unpub/

# make nice table out of nasty VCF data
echo -e "CHROM\tPOS\tID\tREF\tALT\tAC_nfe\tAN_nfe\tAC_afr\tAN_afr\tAC_sas\tAN_sas" > gnomad-3-genome-null.txt
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%AC_nfe\t%AN_nfe\t%AC_afr\t%AN_afr\t%AC_sas\t%AN_sas\n' unpub/gnomad-3-genome-null.vcf >> gnomad-3-genome-null.txt

# and subset properly!
time Rscript 02-cleanup.R gmmat_ssns_ctr_all_random1000.txt gnomad-3-genome-null.txt
# Removing loci unmatched in gnomad:
#   CHR         POS REF    ALT   chrposrefalt            
# 1 chr3  118562867 A      G     chr3:118562867:A:G      
# 2 chr4    2396433 G      A     chr4:2396433:G:A        
# 3 chr12 108308483 GCCTGT G     chr12:108308483:GCCTGT:G
# 4 chr19    133886 T      A     chr19:133886:T:A        
