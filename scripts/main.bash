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

