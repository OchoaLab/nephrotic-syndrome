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

# based on this tutorial
# https://privefl.github.io/bigsnpr/articles/LDpred2.html

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

# start interactive shell if needed
srun --mem 16G -p ochoalab --account ochoalab --pty bash -i
module load R/4.1.1-rhel8 

# add phenotype to fam, which simplifies things later
time Rscript ldpred-00-fix-pheno-fam.R

# add posg to bim, helps set accurate/dynamic window sizes for LD calculations
time Rscript bim-add-posg.R data 38
# 1m20.104s DCC

# create RDS versions of train+test data (here merged)
time Rscript prs-new-04-make-rds.R data
# 2m12.065s DCC

# clean summary stats (convert scores to betas, subset to array)
time Rscript prs-new-05-sumstats-clean.R ssns_ctrl
# Array has these many variants: 761366
# 20,838,869 variants to be matched.
# 82,356 ambiguous SNPs have been removed.
# 672,362 variants have been matched; 0 were flipped and 0 were reversed.
# 2m3.935s DCC
time Rscript prs-new-05-sumstats-clean.R ssns_srns
# Array has these many variants: 761366
# 12,274,557 variants to be matched.
# 76,584 ambiguous SNPs have been removed.
# 635,789 variants have been matched; 0 were flipped and 0 were reversed.
# 1m34.423s DCC

# TODO: here vignette suggests QC on sumstats, but no code is provided (there's equations and a massive repo, may have to revisit)

# subset again to match to training data
time Rscript prs-new-06-sumstats-match.R ssns_ctrl
# 672,362 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 508,082 variants have been matched; 0 were flipped and 0 were reversed.
time Rscript prs-new-06-sumstats-match.R ssns_srns
# 635,789 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 505,231 variants have been matched; 0 were flipped and 0 were reversed.

# calculate LD matrices for training data
# since they have to match SNP sets and these vary by type, run for both base gwas types
type=ssns_ctrl; sbatch -J ld-$type -o ld-$type.out --export=type=$type prs-new-07-ld-matched-snps.q
# 9m35.203s/163m1.993s DCC
type=ssns_srns; sbatch -J ld-$type -o ld-$type.out --export=type=$type prs-new-07-ld-matched-snps.q
# 7m16.877s/120m55.839s DCC

# randomly select the subset of Bristol individuals used for testing vs training
time Rscript ldpred-000-draw-test-train-inds.R

# then run ldpred-inf version, test a grid of heritabilities to determine quickly what is more promising
type=ssns_ctrl; sbatch -J ldpred-01-inf-$type -o ldpred-01-inf-$type.out --export=type=$type ldpred-01-inf.q
# 2m38.207s DCC
type=ssns_srns; sbatch -J ldpred-01-inf-$type -o ldpred-01-inf-$type.out --export=type=$type ldpred-01-inf.q
# 2m5.644s DCC

# fit parameters using training data (inf version)
time Rscript ldpred-01-inf-fit.R ssns_ctrl
# 1m20.346s DCC
time Rscript ldpred-01-inf-fit.R ssns_srns
# 1m20.708s DCC

# run grid version, which is more computationally intensive
type=ssns_ctrl; sbatch -J ldpred-03-grid-$type -o ldpred-03-grid-$type.out --export=type=$type ldpred-03-grid.q
# 31m56.981s DCC
type=ssns_srns; sbatch -J ldpred-03-grid-$type -o ldpred-03-grid-$type.out --export=type=$type ldpred-03-grid.q
# 31m54.374s DCC

# fit parameters using training data (grid version)
time Rscript ldpred-04-grid-fit.R ssns_ctrl
# 3m13.340s DCC
time Rscript ldpred-04-grid-fit.R ssns_srns
# 3m19.394s DCC

# run auto version
type=ssns_ctrl; sbatch -J ldpred-05-auto-$type -o ldpred-05-auto-$type.out --export=type=$type ldpred-05-auto.q 
type=ssns_srns; sbatch -J ldpred-05-auto-$type -o ldpred-05-auto-$type.out --export=type=$type ldpred-05-auto.q 
# 3m10.090s DCC OLD

# run lassosum version
type=ssns_ctrl; sbatch -J ldpred-06-lassosum-$type -o ldpred-06-lassosum-$type.out --export=type=$type ldpred-06-lassosum.q
# untimed because ran interactively debugging
type=ssns_srns; sbatch -J ldpred-06-lassosum-$type -o ldpred-06-lassosum-$type.out --export=type=$type ldpred-06-lassosum.q
# 3m30.440s DCC

# fit parameters using training data (lassosum version)
type=ssns_ctrl; sbatch -J ldpred-07-lassosum-fit-$type -o ldpred-07-lassosum-fit-$type.out --export=type=$type ldpred-07-lassosum-fit.q
# 2m23.840s DCC
type=ssns_srns; sbatch -J ldpred-07-lassosum-fit-$type -o ldpred-07-lassosum-fit-$type.out --export=type=$type ldpred-07-lassosum-fit.q
# 2m29.938s DCC

# get correlation values that actually reveal which value was best
type=ssns_ctrl; sbatch -J ldpred-02-score-$type -o ldpred-02-score-$type.out --export=type=$type ldpred-02-score.q
# 0m52.456s DCC
type=ssns_srns; sbatch -J ldpred-02-score-$type -o ldpred-02-score-$type.out --export=type=$type ldpred-02-score.q
# 0m52.168s DCC

# TODO: include PCs in all evaluations!
# '/datacommons/ochoalab/ssns_gwas/replication/bristol_data/GMMAT/ssns_srns/ssns_srns_mac20.eigenvec'


#######################
### LDPRED2 PRS-NEW ###
#######################

# reboot redoing discovery, separating into two subsets for improved validation

module load R/4.1.1-rhel8 
module load Plink/2.00a3LM

# define subsets (lists of individuals in each set, */ids.txt below)
time Rscript prs-new-00-create-subsets.R

# actually create data
cd /datacommons/ochoalab/ssns_gwas/imputed/prs-new
time plink2 --bfile ../mac20 --keep base/ids.txt --mac 20 --make-bed --out base/mac20
# 2m11.571s DCC
time plink2 --bfile ../mac20 --keep train/ids.txt --mac 20 --make-bed --out train/mac20
# 1m14.038s DCC
# onl testing data comes from Bristol, which is in a totally different, awkward path
time plink2 --bfile ../../replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep test/ids.txt --mac 20 --make-bed --out test/mac20
# 0m13.507s DCC

# confirm dimensions
wc -l */ids.txt
# 4085 base/ids.txt
#  519 test/ids.txt
#  386 train/ids.txt

wc -l */*.{bim,fam}
# 20511795 base/mac20.bim
#  8043559 test/mac20.bim
#  9549985 train/mac20.bim
#     4085 base/mac20.fam
#      514 test/mac20.fam
#      386 train/mac20.fam

# NOTES:
# - only testing (Bristol) lost individuals (5 total), is it because my rawer list hadn't been de-duped yet?  This is fine.

# add phenotype and sex to fam files (otherwise blank), which will simplify PRS trait handling later
time Rscript prs-new-01-fix-pheno-fam.R
# 0m10.434s DCC

# calculate GRMs and PCs for every dataset
# (PCs can be used to adjust correlations)
dir=base;  sbatch -J grm-$dir -o grm-$dir.out --export=dir=$dir prs-new-02-grm.q
# 88m22.320s DCC
dir=train; sbatch -J grm-$dir -o grm-$dir.out --export=dir=$dir prs-new-02-grm.q
# 0m46.235s DCC
dir=test;  sbatch -J grm-$dir -o grm-$dir.out --export=dir=$dir prs-new-02-grm.q
# 0m56.275s DCC

# runs GMMAT on base data only
sbatch prs-new-03-gmmat.q
# 1721m16.303s/1699m25.061s DCC glmmkin and GDS parts
# 953m41.492s/18874m47.095s DCC glmm.score part
# compress output
gzip base/mac20-glmm-score.txt

# add posg to bim, helps set accurate/dynamic window sizes for LD calculations
# (really only required for training data, but could add to all of them for completeness)
time Rscript bim-add-posg.R base/mac20 38
# 4m34.560s DCC
time Rscript bim-add-posg.R train/mac20 38
# 1m53.121s DCC
time Rscript bim-add-posg.R test/mac20 38
# 1m35.606s DCC

# create RDS versions of train and test (needed for both)
time Rscript prs-new-04-make-rds.R train/mac20
# 1m42.520s DCC
time Rscript prs-new-04-make-rds.R test/mac20
# 1m38.067s DCC

# clean summary stats (convert scores to betas, subset to array)
time Rscript prs-new-05-sumstats-clean.R
# Array has these many variants: 761366
# 20,511,795 variants to be matched.
# 82,354 ambiguous SNPs have been removed.
# 672,360 variants have been matched; 0 were flipped and 0 were reversed.
# 2m10.932s DCC

# subset again to match to training data
time Rscript prs-new-06-sumstats-match.R
# 672,360 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 528,964 variants have been matched; 0 were flipped and 0 were reversed.

# calculate LD matrix for training data (check that it works if $type is left undefined, as it should be here!)
sbatch -J ld -o ld.out prs-new-07-ld-matched-snps.q
# 6m17.402s/86m49.116s DCC




