# steps for intersecting loci between GWAS and TGP datasets, for merging

# go where the GWAS data is
cd ~/dbs/ssns_gwas/

wc -l ssns.gwas.hg37.{bim,fam}
# 1748250 ssns.gwas.hg37.bim
#     762 ssns.gwas.hg37.fam

perl -w common_loci_bim.pl ssns.gwas.hg37.bim ~/dbs/tgp/plink2/all_phase3_filt-minimal.bim out.bim
# Reading file 1: ssns.gwas.hg37.bim
# File 1 (ssns.gwas.hg37.bim): 1677135 loci passed filters
# Reading file 2: /home/viiia/dbs/tgp/plink2/all_phase3_filt-minimal.bim
# and writing: out.bim

wc -l out.bim 
# 1407537 out.bim
gzip out.bim

################################

# script for plotting allele frequencies in cases vs controls
# here we provide name of BED data in DMPI server, but analysis will work with other BED files as long as case/control status is given in FAM pheno column as 2/1 values.
Rscript af-case-ctrl.R /dmpi/analysis/Rasheed/ssns_gwas/merge_tgp/combine

# CATHGEN

cd /dmpi/analysis/Rasheed/
ls -lh
# total 5.8M
# drwxrws--- 2 root   rasheedlab    0 Aug 30  2018 Project1
# drwxrws--- 2 root   rasheedlab  160 Feb 20 15:40 gwas_plink
# lrwxrwx--- 1 root   rasheedlab   84 Nov  2 16:00 pVCF -> /dmpi/analysis/cathgen_dmpi_analysis/exome_seq/DUKE_Freeze_One_GRCh38_pVCF/data/pVCF
# -rwxrwxr-x 1 as1052 rasheedlab 3.3M Jan 24  2020 plink2
# -rwxrwx--- 1 as1052 rasheedlab 1.6M Nov  5 12:52 plink2_linux_i686.zip
# drwxrws--- 3 ao128  rasheedlab  333 Nov 16 10:32 ssns_gwas
# drwxrwxr-x 2 ao128  rasheedlab  136 Nov  5 15:16 ssns_wgs_family
cd pVCF/
wc -l all_variants/DUKE_Freeze_One_GRCh38.GL.pVCF.biallelic.{bim,fam}
# 3533081 all_variants/DUKE_Freeze_One_GRCh38.GL.pVCF.biallelic.bim
#    8783 all_variants/DUKE_Freeze_One_GRCh38.GL.pVCF.biallelic.fam
wc -l PASS_variants/DUKE_Freeze_One_GRCh38.GL.PASS.pVCF.biallelic.{bim,fam}
# 3320159 PASS_variants/DUKE_Freeze_One_GRCh38.GL.PASS.pVCF.biallelic.bim
#    8783 PASS_variants/DUKE_Freeze_One_GRCh38.GL.PASS.pVCF.biallelic.fam

cd /dmpi/analysis/Rasheed/
cd gwas_plink/
wc -l all_fixed4.{bim,fam}
# 905955 all_fixed4.bim
#   3264 all_fixed4.fam
wc -l black_inds.txt
# 687 black_inds.txt
wc -l white_inds.txt 
# 2277 white_inds.txt

###################

# gwas hack
# data from Amika, from a run with several issues but meh
# creates several plots looking at the p-values, filtering, and manhattan plots
Rscript gwas-simple-filter.R
# converg some of the big PDFs to PNG
pdf2pngScreen ../data/pca_combine3_maf10.PHENO1.glm.logistic-manhattan-all.pdf
pdf2pngScreen ../data/pca_combine3_maf10.PHENO1.glm.logistic-manhattan-clean.pdf
