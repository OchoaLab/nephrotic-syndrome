### GMMAT Wald/ORs ###

# run from this location:
cd /datacommons/ochoalab/ssns_gwas/GMMAT_0418/PRS/

# SRNS vs SSNS

# get ORs for top SNPs
sbatch -a 1-100 gmmat_OR_srns_ssns_batch.q
# gather batches into a single file
Rscript gmmat_OR_srns_ssns_batch-collect.R

# SSNS vs control
sbatch -a 1-100 gmmat_OR_ssns_ctrl_batch.q
Rscript gmmat_OR_ssns_ctrl_batch-collect.R

# calculate in-sample R2 for SRNS

# merge covariates (demographics and PCs) into one with standard column names
cd /datacommons/ochoalab/ssns_gwas/GMMAT_0418/srns_ctr/srns_ssns/
Rscript prsice-00-reformat-srns-discovery.R

# back to PRS subdir
cd /datacommons/ochoalab/ssns_gwas/GMMAT_0418/PRS/
mkdir srns-discovery
cd srns-discovery
base=/datacommons/ochoalab/ssns_gwas
# for genotypes, covariates
name=$base/GMMAT_0418/srns_ctr/srns_ssns/srns_ssns_mac20
# phenotype is elsewhere, lacks header as desired!
phen=$base/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/ssns_srns_pheno.txt

pcuts=1e-40,1e-39,1e-38,1e-37,1e-36,1e-35,1e-34,1e-33,1e-32,1e-31,1e-30,1e-29,1e-28,1e-27,1e-26,1e-25,1e-24,1e-23,1e-22,1e-21,1e-20,1e-19,1e-18,1e-17,1e-16,1e-15,1e-14,1e-13,1e-12,1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,2e-5,3e-5,4e-5,5e-5,6e-5,7e-5,8e-5,9e-5,1e-4,2e-4,3e-4,4e-4,5e-4,6e-4,7e-4,8e-4,9e-4

# NOTE: this case appears backwards unexpectedly
module load PRSice/2.3.3
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
       --fastscore \
       --bar-levels $pcuts \
       -b $base/GMMAT_0418/PRS/glmm.wald_srns_ssns.txt.gz \
       -t $name \
       -f $phen \
       -C ${name}_covar_prsice.txt \
       --cov-factor sex,race \
       -o prsice
# 0m51.982s
module purge

# plotting failed on DCC, downloaded the necessary outputs (non-genetic, non-identifiable info) and ran this locally to complete plots:
time ~/bin/PRSice_linux/PRSice.R \
     --plot \
     --binary-target T \
     --bar-levels $pcuts \
     -t srns_ssns_mac20 \
     -f ssns_srns_pheno.txt \
     -o prsice


# test SSNS-ctrl on SRNS-SSNS (both discovery but only SSNS are shared)
# use same $base, $name, $phen as above; all is the same except
# - base data is ssns-ctrl
# - expect direction to be reversed (srns is 1 in srns-ssns, but ssns is 1 in ssns-ctrl), so I reversed A1/A2 to compensate (back to default; that is, usually gmmat is reversed but here we want that)
cd ..
mkdir srns-discovery-ssns
cd srns-discovery-ssns

module load PRSice/2.3.3
time PRSice \
       --bp POS \
       -p PVAL \
       --beta \
       --binary-target T \
       --clump-kb 6000 \
       --clump-r2 0.3 \
       --keep-ambig \
       --memory 1000 \
       -n 1 \
       -k 0.2 \
       --fastscore \
       --bar-levels $pcuts \
       -b $base/GMMAT_0418/PRS/glmm.wald_ssns_ctrl.txt.gz \
       -t $name \
       -f $phen \
       -C ${name}_covar_prsice.txt \
       --cov-factor sex,race \
       -o prsice
# 0m51.982s
module purge

# plotting failed on DCC, downloaded the necessary outputs (non-genetic, non-identifiable info) and ran this locally to complete plots:
time ~/bin/PRSice_linux/PRSice.R \
     --plot \
     --binary-target T \
     --bar-levels $pcuts \
     -t srns_ssns_mac20 \
     -f ssns_srns_pheno.txt \
     -o prsice


### BRISTOL VALIDATIONS ###

# ALL run from this location!
cd /datacommons/ochoalab/ssns_gwas/replication/bristol_data/GMMAT

wc -l srns_ssns_mac20.{bim,fam}
# 7891682 srns_ssns_mac20.bim
#     490 srns_ssns_mac20.fam

# base data is really small in this case, maybe should extend it, make it bigger!
zcat /datacommons/ochoalab/ssns_gwas/GMMAT_0418/PRS/glmm.wald_srns_ssns.txt.gz|wc -l
# 2651
# 13972

# creates
# - srns_ssns.phen (removed bad col names of srns_ssns_pheno.txt)
# - srns_ssns_covar_prsice.txt (merged/clean of srns_ssns_covar.txt and srns_ssns.eigenvec)
module load R/4.0.0
Rscript prsice-00-reformat-srns-bristol.R
module unload R/4.0.0

# NOTE: this other module appears newer but isn't!
# module load PRSice/3.2.3 # actually 2.1.4.beta (19 October 2018), older!

# NOTES:
# -k <prevalence> used 0.2 = proportion of NS cases that are SRNS, from Rasheed's intro

# put the output in a separate dir
mkdir srns-ssns-bristol-srns-ssns-discovery
cd srns-ssns-bristol-srns-ssns-discovery

# NOTE: this case appears backwards unexpectedly
module load PRSice/2.3.3
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
     --fastscore \
     --bar-levels $pcuts \
     -b $base/GMMAT_0418/PRS/glmm.wald_srns_ssns.txt.gz \
     -t ../srns_ssns_mac20 \
     -f ../srns_ssns.phen \
     -C ../srns_ssns_covar_prsice.txt \
     --cov-factor sex,race \
     -o prsice
module purge
# 0m15.906s slice1
# 0m27.905s, 0m9.434s slice2

# plotting failed on DCC, downloaded the necessary outputs (non-genetic, non-identifiable info) and ran this locally to complete plots:
time ~/bin/PRSice_linux/PRSice.R \
     --plot \
     --binary-target T \
     --bar-levels $pcuts \
     -t srns_ssns_mac20 \
     -f srns_ssns.phen \
     -o prsice

# version with ssns gwas!
cd ..
mkdir srns-ssns-bristol-ssns-ctrl-discovery
cd srns-ssns-bristol-ssns-ctrl-discovery
module load PRSice/2.3.3
time PRSice \
       --bp POS \
       -p PVAL \
       --beta \
       --binary-target T \
       --clump-kb 6000 \
       --clump-r2 0.3 \
       --keep-ambig \
       --memory 1000 \
       -n 1 \
       -k 0.2 \
       --fastscore \
       --bar-levels $pcuts \
       -b $base/GMMAT_0418/PRS/glmm.wald_ssns_ctrl.txt.gz \
       -t ../srns_ssns_mac20 \
       -f ../srns_ssns.phen \
       -C ../srns_ssns_covar_prsice.txt \
       --cov-factor sex,race \
       -o prsice
module purge

# plot locally
time ~/bin/PRSice_linux/PRSice.R \
     --plot \
     --binary-target T \
     --bar-levels $pcuts \
     -t srns_ssns_mac20 \
     -f srns_ssns.phen \
     -o prsice

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


###############
### LDPRED2 ###
###############

# go where the data is
cd /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp
# make working dir for this project, to keep things separate, since many temp files are created by ldpred
mkdir prs
cd prs
ln -s ../bristol_impute_mac20.bed data.bed
ln -s ../bristol_impute_mac20.bim data.bim
ln -s ../bristol_impute_mac20.fam data.fam

# for reruns
cd /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/prs
