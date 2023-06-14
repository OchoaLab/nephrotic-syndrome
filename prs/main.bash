### GMMAT Wald/ORs ###

# SRNS vs SSNS

# get ORs for top SNPs
sbatch -a 1-100 gmmat_OR_srns_ssns_batch.q
# gather batches into a single file
Rscript gmmat_OR_srns_ssns_batch-collect.R

# SSNS vs control

sbatch -a 1-100 gmmat_OR_ssns_ctrl_batch.q
Rscript gmmat_OR_ssns_ctrl_batch-collect.R



### PRSice ###

# ALL run from this location!
# /datacommons/ochoalab/ssns_gwas/replication/bristol_data/GMMAT

wc -l srns_ssns_mac20.{bim,fam}
# 7891682 srns_ssns_mac20.bim
#     490 srns_ssns_mac20.fam

# base data is really small in this case, maybe should extend it, make it bigger!
zcat /datacommons/ochoalab/ssns_gwas/GMMAT_0418/PRS/glmm.wald_srns_ssns.txt.gz|wc -l
# 2651

# creates
# - srns_ssns.phen (removed bad col names of srns_ssns_pheno.txt)
# - srns_ssns_covar_prsice.txt (merged/clean of srns_ssns_covar.txt and srns_ssns.eigenvec)
module load R/4.0.0
Rscript prsice-00-reformat.R
module unload R/4.0.0

# load preinstalled PRSice on DCC!
#module load PRSice/3.2.3 # actually 2.1.4.beta (19 October 2018), older!
module load PRSice/2.3.3 # correct, 2.3.3 (2020-08-05) 
# my local version is 2.3.5 (2021-09-20) 
# online log goes up to 2.3.3, so I don't know what the difference is

# NOTES:
# -k <prevalence> used 0.2 = proportion of NS cases that are SRNS, from Rasheed's intro

# time PRSice.R \
#      --prsice /opt/apps/rhel8/PRSice-2.3.3/PRSice \
time PRSice \
       --bp POS \
       -p PVAL \
       --beta \
       --binary-target T \
       --a1 A2 \
       --a2 A1 \
       --clump-kb 6000 \
       --clump-r2 0.3 \
       --keep-ambig \
       --memory 1000 \
       -n 1 \
       -k 0.2 \
       -b /datacommons/ochoalab/ssns_gwas/GMMAT_0418/PRS/glmm.wald_srns_ssns.txt.gz \
       -t srns_ssns_mac20 \
       -f srns_ssns.phen \
       -C srns_ssns_covar_prsice.txt \
       --cov-factor sex,race \
       -o srns_ssns_mac20_prsice
# 0m15.906s

# plotting failed on DCC, downloaded the necessary outputs (non-genetic, non-identifiable info) and ran this locally to complete plots:
time ~/bin/PRSice_linux/PRSice.R \
     --plot \
     --binary-target T \
     -t srns_ssns_mac20 \
     -f srns_ssns.phen \
     -o srns_ssns_mac20_prsice

# add after assessing if results depend too much on this, and runtime, only active if --perm <n> is used
#    --logit-perm

# to see all options
# PRSice

# base header for wald test
# CHR SNP                          cM POS      A1 A2           N   AF                 BETA                SE                  PVAL                  converged
# 1   chr1:17947774:CCAAATGTAAAT:C 0  17947774 C  CCAAATGTAAAT 918 0.8055555555555556 -0.6724062622564424 0.1723690003733841  9.580810979269296e-5  TRUE
# 1   chr1:17947810:G:A            0  17947810 A  G            918 0.8050108932461874 -0.6605024881345551 0.17195180467048415 1.2243077790403302e-4 TRUE
# 1   chr1:17948129:A:G            0  17948129 G  A            918 0.7908496732026143 -0.6442018520772689 0.16264308957822604 7.468926446101956e-5  TRUE
# 1   chr1:18888482:A:G            0  18888482 G  A            918 0.982570806100218  -1.6106168821651532 0.4346107919473915  2.106552723113269e-4  TRUE
# 1   chr1:19533450:G:A            0  19533450 A  G            918 0.9842047930283224 -1.937307174485691  0.5269969683580391  2.3680274313231243e-4 TRUE
