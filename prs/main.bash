### GMMAT Wald/ORs ###

# SRNS vs SSNS

# get ORs for top SNPs
sbatch -a 1-100 gmmat_OR_srns_ssns_batch.q
# gather batches into a single file
Rscript gmmat_OR_srns_ssns_batch-collect.R


# SSNS vs control ... Coming soon!



### PRSice ###

# NOT DONE YET!

# load preinstalled PRSice on DCC!
module load PRSice/3.2.3 # actually 2.1.4.beta (19 October 2018), older!
module load PRSice/2.3.3 # correct, 2.3.3 (2020-08-05) 
# my local version is 2.3.5 (2021-09-20) 
# online log goes up to 2.3.3, so I don't know what the difference is

# NOTES:
# -k <prevalence> used 0.2 = proportion of NS cases that are SRNS, from Rasheed's intro

time Rscript PRSice.R \
     --prisce PRSice \
#time PRSice \
    --bp POS \
    --pvalue PVAL \
    --or \
    --stat OR \
    --binary-target T \
    --clump-kb 6000 \
    --clump-r2 0.3 \
    --keep-ambig \
    -o prsice \
    --memory 1000 \
    -n 1 \
    -k 0.2 \
    -f <pheno file for target?> \
    --cov <continuous covariate file for target?> \
    --cov-factor <categorical covariate file for target?> \
    -b /datacommons/ochoalab/ssns_gwas/GMMAT_0418/srns_ctr/srns_ssns/glmm.score.bed.srns_ssns_PC_mac20_nohwe.txt \
    -t /datacommons/ochoalab/ssns_gwas/replication/bristol_data/GMMAT/srns_ssns_hwe_mac20

# add after assessing if results depend too much on this, and runtime
    --logit-perm \

# to see all options
PRSice

# base header
# CHR SNP             cM POS     A1 A2 N   AF       SCORE     VAR     PVAL
# 1   chr1:66861:C:T  0  66861   T  C  918 0.977669 -2.99743  6.86321 0.252559
# 1   chr1:598941:G:A 0  598941  A  G  918 0.989107  0.804131 2.11309 0.580139
# 1   chr1:710225:T:A 0  710225  A  T  918 0.985294 -1.82461  4.73649 0.401815
# 1   chr1:722408:G:C 0  722408  C  G  918 0.322985 -11.8305  37.437  0.0531706
# 1   chr1:722858:C:T 0  722858  T  C  918 0.987473  2.66111  3.02411 0.125954

# PROBLEMS:
# - don't have OR in this "score" test!  Have to rerun :(
# - 
