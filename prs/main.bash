###############
### LDPRED2 ###
###############

# based on this tutorial
# https://privefl.github.io/bigsnpr/articles/LDpred2.html

# start interactive shell if needed
srun --mem 16G -p ochoalab --account ochoalab --pty bash -i
module load R/4.1.1-rhel8 
module load Plink/2.00a3LM

# location of data and scripts
cd /datacommons/ochoalab/ssns_gwas/imputed/prs-new

# define subsets to split Discovery mainly, but also cleans up Bristol minimally
# (lists of individuals in each set, */ids.txt below)
time Rscript prs-new-00-create-subsets.R

# actually create data
time plink2 --bfile ../mac20 --keep base/ids.txt --mac 20 --make-bed --out base/mac20
# 2m11.571s DCC
time plink2 --bfile ../mac20 --keep train/ids.txt --mac 20 --make-bed --out train/mac20
# 1m14.038s DCC
# only testing data comes from Bristol, which is in a totally different, awkward path
time plink2 --bfile ../../replication/bristol_data/imputation/post_imp/bristol_impute_mac20 --keep test/ids.txt --mac 20 --make-bed --out test/mac20
# 0m13.507s DCC

# confirm dimensions
wc -l */ids.txt
# 4085 base/ids.txt
#  522 test/ids.txt
#  386 train/ids.txt
 
wc -l */*.{bim,fam}
# 20511795 base/mac20.bim
#  8053402 test/mac20.bim
#  9549985 train/mac20.bim
#     4085 base/mac20.fam
#      517 test/mac20.fam
#      386 train/mac20.fam

# NOTES:
# - only testing (Bristol) lost individuals (5 total), because my rawer list had duplicates and high-missingness individuals.  This is fine.

# alternatively, make working copies of the original base data (not split in this weird way)
# make a copy because we edit positions, calculate LD, create other temp files
for name in ssns_ctrl ssns_srns; do
    # make new base name
    mkdir base-$name
    cd base-$name
    # link genotype data, to estimate LD from
    ln -s {../../$name/,}mac20.bed
    ln -s {../../$name/,}mac20.bim
    ln -s {../../$name/,}mac20.fam
    # don't need GRM and PCs, but we could link them too otherwise
    # link summary stats too
    ln -s {../../../saige/$name/,}saige_output.txt # SAIGE
    # go back down
    cd ..
done

# the same for curegn as a training dataset, and later also testing in a non-overlapping case
mkdir train-curegn
cd train-curegn
# link genotype data and PCs
curegn=../../../../curegn/subanalysis/curegn_ssns_srns_
ln -s {${curegn},}mac20.bed
ln -s {${curegn},}mac20.bim
ln -s {${curegn},}mac20.fam
ln -s {${curegn},}mac20.eigenvec
# # assess missingness (only CureGN has any because it wasn't imputed)
# plink2 --bfile mac20 --missing --out mac20
# # look at AFs, to consider dumb imputation schemes
# plink2 --bfile mac20 --freq --out mac20

# recently source files got moved or disappeared for weird reasons; bed/bim/fam are fine (were processed and overwritten), just eigenvec file needed regenerating, do that here!
time plink2 --bfile mac20 --pca --out mac20
# 2m28.241s DCC
cd ..

# add phenotype and sex to fam files (otherwise blank), which will simplify PRS trait handling later
time Rscript prs-new-01-fix-pheno-fam.R
# 0m10.434s DCC

# calculate GRMs and PCs for every new dataset
# (PCs can be used to adjust correlations)
dir=base;  sbatch -J grm-$dir -o grm-$dir.out --export=dir=$dir prs-new-02-grm.q
# 88m22.320s DCC
dir=train; sbatch -J grm-$dir -o grm-$dir.out --export=dir=$dir prs-new-02-grm.q
# 0m46.235s DCC
dir=test;  sbatch -J grm-$dir -o grm-$dir.out --export=dir=$dir prs-new-02-grm.q
# 0m56.275s DCC

# more recently I needed PCs again and had to normalize the format, so redid eigenvec with plink2 instead of using GCTA as above
cd test
# keep the original files just in case
mv mac20.eigenval mac20-gcta.eigenval
mv mac20.eigenvec mac20-gcta.eigenvec
# make new version as above
time plink2 --bfile mac20 --pca --out mac20
# 2m58.211s DCC
cd ..

# Tiffany ran SAIGE on "base" using this code:
sbatch saige.q


# add posg to bim, helps set accurate/dynamic window sizes for LD calculations
# (really only required for training data, but could add to all of them for completeness)
time Rscript bim-add-posg.R base/mac20 38
# 4m34.560s DCC
time Rscript bim-add-posg.R train/mac20 38
# 1m53.121s DCC
time Rscript bim-add-posg.R test/mac20 38
# 1m35.606s DCC
time Rscript bim-add-posg.R base-ssns_ctrl/mac20 38
# 3m11.185s DCC
time Rscript bim-add-posg.R base-ssns_srns/mac20 38
# 1m53.496s DCC
time Rscript bim-add-posg.R train-curegn/mac20 38
# 1m14.851s DCC

# since cureGN has low amounts of missingness, as a hack, apply some semi-random imputation (based on global AFs, but unlikely to matter)
time Rscript impute-dumb.R train-curegn/mac20
# 7m42.698s DCC
# NOTE: we will use original PCs estimated from data with missingnes, which is better; nothing else is affected in this case (we're not calculating LD for CureGN)

# create RDS versions of all data
time Rscript prs-new-04-make-rds.R train/mac20
# 1m42.520s DCC
time Rscript prs-new-04-make-rds.R test/mac20
# 1m38.067s DCC
time Rscript prs-new-04-make-rds.R base/mac20
# 16m13.505s DCC
time Rscript prs-new-04-make-rds.R base-ssns_ctrl/mac20
# 16m49.657s DCC
time Rscript prs-new-04-make-rds.R base-ssns_srns/mac20
# 3m15.834s DCC
time Rscript prs-new-04-make-rds.R train-curegn/mac20
# 1m33.345s DCC

# now is a good time to link (first) curegn test and train versions
# NOTES: do after imputing genotypes and creating RDS!  Also inherits posg, though that is unimportant
mkdir test-curegn
cd test-curegn
ln -s {../train-curegn/,}mac20.bed
ln -s {../train-curegn/,}mac20.bim
ln -s {../train-curegn/,}mac20.fam
ln -s {../train-curegn/,}mac20.eigenvec
ln -s {../train-curegn/,}mac20.bk
ln -s {../train-curegn/,}mac20.rds
cd ..

# clean summary stats (subset to array)
time Rscript prs-new-05-sumstats-clean.R base
# Array has these many variants: 761366
# 20,511,795 variants to be matched.
# 82,354 ambiguous SNPs have been removed.
# 672,360 variants have been matched; 0 were flipped and 0 were reversed.
# 2m2.685s DCC
time Rscript prs-new-05-sumstats-clean.R base-ssns_ctrl
# Array has these many variants: 761366
# 20,838,869 variants to be matched.
# 82,356 ambiguous SNPs have been removed.
# 672,362 variants have been matched; 0 were flipped and 0 were reversed.
# 2m5.272s DCC
time Rscript prs-new-05-sumstats-clean.R base-ssns_srns
# Array has these many variants: 761366
# 12,274,557 variants to be matched.
# 76,584 ambiguous SNPs have been removed.
# 635,789 variants have been matched; 0 were flipped and 0 were reversed.
# 1m26.980s DCC

# subset again to match to training data, and to self-base in all cases to get LD in base
time Rscript prs-new-06-sumstats-match.R base train
# 672,360 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 528,964 variants have been matched; 0 were flipped and 0 were reversed.
# 1m36.754s DCC
time Rscript prs-new-06-sumstats-match.R base base
# 672,360 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 672,360 variants have been matched; 0 were flipped and 0 were reversed.
# 2m44.363s DCC
time Rscript prs-new-06-sumstats-match.R base-ssns_ctrl base-ssns_ctrl
# 672,362 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 672,362 variants have been matched; 0 were flipped and 0 were reversed.
# 2m40.015s DCC
time Rscript prs-new-06-sumstats-match.R base-ssns_srns base-ssns_srns
# 635,789 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 635,789 variants have been matched; 0 were flipped and 0 were reversed.
# 1m48.589s DCC
# even in new setup, need these to align to training data
time Rscript prs-new-06-sumstats-match.R base-ssns_ctrl train-curegn
# 672,362 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 511,666 variants have been matched; 307 were flipped and 881 were reversed.
# 1m25.818s DCC
time Rscript prs-new-06-sumstats-match.R base-ssns_srns train-curegn
# 635,789 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 510,987 variants have been matched; 306 were flipped and 881 were reversed.
# 1m25.505s DCC

# calculate LD matrix
base=base; sbatch -J ld-$base -o ld-$base.out --export=base=$base prs-new-07-ld-matched-snps.q
# 80m57.891s/1240m14.603s DCC
base=base-ssns_ctrl; sbatch -J ld-$base -o ld-$base.out --export=base=$base prs-new-07-ld-matched-snps.q
# 149m20.374s/1265m22.475s DCC
base=base-ssns_srns; sbatch -J ld-$base -o ld-$base.out --export=base=$base prs-new-07-ld-matched-snps.q
# 15m27.211s/233m21.564s DCC

# before scoring, align training and testing data
time Rscript prs-new-08-match-train-test.R base train test
# 528,964 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 466,359 variants have been matched; 0 were flipped and 0 were reversed.
# 1m32.168s DCC
time Rscript prs-new-08-match-train-test.R base train test-curegn
# 528,964 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 491,349 variants have been matched; 292 were flipped and 872 were reversed.
# 1m41.491s DCC
time Rscript prs-new-08-match-train-test.R base-ssns_ctrl train-curegn test
# 511,666 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 470,468 variants have been matched; 213 were flipped and 853 were reversed.
# 1m32.108s DCC
time Rscript prs-new-08-match-train-test.R base-ssns_srns train-curegn test
# 510,987 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 469,912 variants have been matched; 212 were flipped and 853 were reversed.
# 1m27.549s DCC

# then run ldpred-inf version, test a grid of heritabilities to determine quickly what is more promising
base=base; sbatch -J ldpred-01-inf-$base -o ldpred-01-inf-$base.out --export=base=$base ldpred-01-inf.q
# 3m58.094s DCC
base=base-ssns_ctrl; sbatch -J ldpred-01-inf-$base -o ldpred-01-inf-$base.out --export=base=$base ldpred-01-inf.q
# 4m6.373s DCC
base=base-ssns_srns; sbatch -J ldpred-01-inf-$base -o ldpred-01-inf-$base.out --export=base=$base ldpred-01-inf.q
# 3m14.548s DCC

# fit parameters using training data (inf version)
time Rscript ldpred-01-inf-fit.R base train
# 1m22.593s DCC
time Rscript ldpred-01-inf-fit.R base-ssns_ctrl train-curegn
# 1m14.105s DCC
time Rscript ldpred-01-inf-fit.R base-ssns_srns train-curegn
# 1m18.862s DCC
# plot that combines results from all three fits, for paper
time Rscript ldpred-01-inf-fit-all.R
# 0m4.767s DCC

# run grid version, which is more computationally intensive
base=base; sbatch -J ldpred-03-grid-$base -o ldpred-03-grid-$base.out --export=base=$base ldpred-03-grid.q
# 44m53.794s DCC
base=base-ssns_ctrl; sbatch -J ldpred-03-grid-$base -o ldpred-03-grid-$base.out --export=base=$base ldpred-03-grid.q
# 46m2.997s DCC
base=base-ssns_srns; sbatch -J ldpred-03-grid-$base -o ldpred-03-grid-$base.out --export=base=$base ldpred-03-grid.q
# 44m3.382s DCC

# fit parameters using training data (grid version)
time Rscript ldpred-04-grid-fit.R base train
# 3m6.011s DCC
time Rscript ldpred-04-grid-fit.R base-ssns_ctrl train-curegn
# 2m58.704s DCC
time Rscript ldpred-04-grid-fit.R base-ssns_srns train-curegn
# 2m52.887s DCC
# plot that combines results from all three fits, for paper
time Rscript ldpred-04-grid-fit-all.R

# run auto version
base=base; sbatch -J ldpred-05-auto-$base -o ldpred-05-auto-$base.out --export=base=$base ldpred-05-auto.q
# 15m28.283s DCC
base=base-ssns_ctrl; sbatch -J ldpred-05-auto-$base -o ldpred-05-auto-$base.out --export=base=$base ldpred-05-auto.q
# 15m57.995s DCC
base=base-ssns_srns; sbatch -J ldpred-05-auto-$base -o ldpred-05-auto-$base.out --export=base=$base ldpred-05-auto.q
# 8m49.049s DCC

# base version needs one more step to map data to training SNPs
time Rscript ldpred-05-auto-map-base-to-train.R base train
# 0m12.843s DCC
time Rscript ldpred-05-auto-map-base-to-train.R base-ssns_ctrl train-curegn
# 0m12.028s DCC
time Rscript ldpred-05-auto-map-base-to-train.R base-ssns_srns train-curegn
# 0m11.622s DCC

# run lassosum version
base=base; sbatch -J ldpred-06-lassosum-$base -o ldpred-06-lassosum-$base.out --export=base=$base ldpred-06-lassosum.q
# 5m49.648s DCC
base=base-ssns_ctrl; sbatch -J ldpred-06-lassosum-$base -o ldpred-06-lassosum-$base.out --export=base=$base ldpred-06-lassosum.q
# 5m1.334s DCC
base=base-ssns_srns; sbatch -J ldpred-06-lassosum-$base -o ldpred-06-lassosum-$base.out --export=base=$base ldpred-06-lassosum.q
# 6m46.425s DCC

# fit parameters using training data (lassosum version)
base=base; train=train; sbatch -J ldpred-07-lassosum-fit-$base-$train -o ldpred-07-lassosum-fit-$base-$train.out --export=base=$base,train=$train ldpred-07-lassosum-fit.q
# 2m13.012s DCC
base=base-ssns_ctrl; train=train-curegn; sbatch -J ldpred-07-lassosum-fit-$base-$train -o ldpred-07-lassosum-fit-$base-$train.out --export=base=$base,train=$train ldpred-07-lassosum-fit.q
# 1m55.754s DCC
base=base-ssns_srns; train=train-curegn; sbatch -J ldpred-07-lassosum-fit-$base-$train -o ldpred-07-lassosum-fit-$base-$train.out --export=base=$base,train=$train ldpred-07-lassosum-fit.q
# 2m2.664s DCC
# make combined figure
time Rscript ldpred-07-lassosum-fit-all.R

# run regular and "stacked" CT (clump and threshold)
# clumping step is base-only (uses LD information)
base=base; sbatch -J ldpred-10-clump-$base -o ldpred-10-clump-$base.out --export=base=$base ldpred-10-clump.q
# 162m37.282s/1164m7.533s DCC
base=base-ssns_ctrl; sbatch -J ldpred-10-clump-$base -o ldpred-10-clump-$base.out --export=base=$base ldpred-10-clump.q
# 171m58.583s/1213m16.160s DCC
base=base-ssns_srns; sbatch -J ldpred-10-clump-$base -o ldpred-10-clump-$base.out --export=base=$base ldpred-10-clump.q
# 65m48.709s/336m46.577s DCC

# actual training happens now
base=base; train=train; sbatch -J ldpred-11-ct-fit-$base-$train -o ldpred-11-ct-fit-$base-$train.out --export=base=$base,train=$train ldpred-11-ct-fit.q
# 14m13.541s DCC
base=base-ssns_ctrl; train=train-curegn; sbatch -J ldpred-11-ct-fit-$base-$train -o ldpred-11-ct-fit-$base-$train.out --export=base=$base,train=$train ldpred-11-ct-fit.q
# 14m6.112s DCC
base=base-ssns_srns; train=train-curegn; sbatch -J ldpred-11-ct-fit-$base-$train -o ldpred-11-ct-fit-$base-$train.out --export=base=$base,train=$train ldpred-11-ct-fit.q
# 14m11.926s DCC
# make combined figure
time Rscript ldpred-11-ct-fit-all.R

# construct CT-specific report of best nonzero betas
# allows us to count predictor SNPs and determine where they're located
time Rscript ldpred-12-ct-fit-info.R base train
# 0m7.177s DCC
time Rscript ldpred-12-ct-fit-info.R base-ssns_ctrl train-curegn
# 0m9.549s DCC
time Rscript ldpred-12-ct-fit-info.R base-ssns_srns train-curegn
# 0m9.371s DCC

# another quick way to get number of non-zero betas for CT outputs:
zgrep -c -v '^0$' train/betas-base-ldpred2-ct-best.txt.gz # 8
zgrep -c -v '^0$' train-curegn/betas-base-ssns_ctrl-ldpred2-ct-best.txt.gz # 6
zgrep -c -v '^0$' train-curegn/betas-base-ssns_srns-ldpred2-ct-best.txt.gz # 12

# get PRSs and correlation values that actually reveal which model was best
base=base; train=train; test=test; sbatch -J ldpred-02-score-$base-$train-$test -o ldpred-02-score-$base-$train-$test.out --export=base=$base,train=$train,test=$test ldpred-02-score.q
# 1m19.138s DCC
base=base; train=train; test=test-curegn; sbatch -J ldpred-02-score-$base-$train-$test -o ldpred-02-score-$base-$train-$test.out --export=base=$base,train=$train,test=$test ldpred-02-score.q
# 1m40.356s DCC
base=base-ssns_ctrl; train=train-curegn; test=test; sbatch -J ldpred-02-score-$base-$train-$test -o ldpred-02-score-$base-$train-$test.out --export=base=$base,train=$train,test=$test ldpred-02-score.q
# 1m11.484s DCC
base=base-ssns_srns; train=train-curegn; test=test; sbatch -J ldpred-02-score-$base-$train-$test -o ldpred-02-score-$base-$train-$test.out --export=base=$base,train=$train,test=$test ldpred-02-score.q
# 1m10.048s DCC

# make nice plot that summarize testing results
# (these are less polished in many ways though)
time Rscript ldpred-08-test-plot.R base train test
# 0m8.598s DCC
time Rscript ldpred-08-test-plot.R base train test-curegn
# 0m12.047s DCC
time Rscript ldpred-08-test-plot.R base-ssns_ctrl train-curegn test
# 0m14.383s DCC
time Rscript ldpred-08-test-plot.R base-ssns_srns train-curegn test
# 0m6.913s DCC

# combine correlation data across both test datasets, which involves some non-trivial calculations
mkdir test-bristol-curegn
time Rscript ldpred-17-combine-pcor.R
# 0m3.821s DCC

# makes a single plot combining the data from all the previous cases
time Rscript ldpred-09-test-plot-combined.R
# 0m5.039s DCC

# perform ancestry subanalyses too!
time Rscript ldpred-02-score-anc.R base train test
# 1m8.761s DCC
time Rscript ldpred-02-score-anc.R base-ssns_ctrl train-curegn test
# 1m16.066s DCC
time Rscript ldpred-02-score-anc.R base-ssns_srns train-curegn test
# 1m9.333s DCC
time Rscript ldpred-02-score-anc.R base train test-curegn
# untimed

# plot for ancestry subanalysis, for a single method
time Rscript ldpred-13-test-plot-combined-anc.R grid-h0.1-best
# 0m7.819s DCC
time Rscript ldpred-13-test-plot-combined-anc.R ct-best
# 0m7.013s DCC

# look in more detail at scores by themselves, and as predictors (binned into quantiles)
time Rscript ldpred-14-prs-or-quantiles.R test base-train-ldpred2-grid-h0.1-best
# not timed
time Rscript ldpred-14-prs-or-quantiles.R test base-ssns_ctrl-train-curegn-ldpred2-grid-h0.1-best
# 0m21.026s DCC
time Rscript ldpred-14-prs-or-quantiles.R test base-ssns_srns-train-curegn-ldpred2-grid-h0.1-best
# 0m21.497s DCC
time Rscript ldpred-14-prs-or-quantiles.R test-curegn base-train-ldpred2-grid-h0.1-best
# 0m21.301s DCC

# do some for the CT method, hopefully it is similar qualitatively
# do only these two, which are the most relevant cases
# scale is different but results otherwise sort of similar
time Rscript ldpred-14-prs-or-quantiles.R test base-train-ldpred2-ct-best
time Rscript ldpred-14-prs-or-quantiles.R test-curegn base-train-ldpred2-ct-best

# version that combines test (bristol) and test-curegn, for larger sample sizes, now that we know that they were both similar when looked at separately
time Rscript ldpred-15-prs-or-quantiles-combined.R base-train-ldpred2-grid-h0.1-best
time Rscript ldpred-15-prs-or-quantiles-combined.R base-train-ldpred2-ct-best

# paper figures, with multiple panels, combining best of earlier scripts (14-15)
time Rscript ldpred-16-prs-or-quantiles-paper-panels.R base-train-ldpred2-grid-h0.1-best
# AUC values
#   Bristol         CureGN Bristol+CureGN 
# 0.6264282      0.6496568      0.6382799 
# 0m5.551s DCC
time Rscript ldpred-16-prs-or-quantiles-paper-panels.R base-train-ldpred2-ct-best
# AUC values
#   Bristol         CureGN Bristol+CureGN 
# 0.6265045      0.6302012      0.6282716 
# 0m4.619s DCC

# create "score" files to share, easy to use with plink2
time Rscript ldpred-20-prs-betas-plink-score.R base train
# confirm that scores using plink2 match (scaled versions anyway) of the ldpred numbers
method=ct; time plink2 --bfile test/mac20 --score train/betas-base-ldpred2-$method-best-plink-score.txt header --out test/prs-base-train-ldpred2-$method-best-plink
# 0m6.421s DCC
method=grid-h0.1; time plink2 --bfile test/mac20 --score train/betas-base-ldpred2-$method-best-plink-score.txt header --out test/prs-base-train-ldpred2-$method-best-plink
# 0m20.160s DCC
# confirms that plink and ldpred scores match!  (they are not negatively correlated, for example)
time Rscript ldpred-21-prs-betas-plink-score-validate.R base train test

# make a big table with info for Debo to evaluate his diagnistic test
time Rscript ldpred-22-mk-table-for-debo.R


#######################
### LDPRED2 CLEANUP ###
#######################

# remove bulky file-backed matrices and things of that sort that are redundant with rawer data and easy to regenerate if needed
# these are huge and satisfies the above:
rm */mac20.{rds,bk}
# remove less complete copies of bim/fam
# some of these are softlinks (that didn't take up any space), but others were big files, especially BIM
rm */mac20.bim~
rm */mac20_ORIG.fam
# these are only used as precursors of PCs, not used at all afterwards
rm */mac20.grm.*
# this we want to keep, but we'll compress for now (only this one isn't a softlink).  Will have to fix path if we want to use it later again in the same script
gzip base/saige_output.txt
# boring pca logs
rm */mac20.log
# LD stuff too, it can be slow to regenerate but not extremely so, and it is quite huge
rm */ld*.{sbk,RData}


#######################
### PROPOSAL PRELIM ###
#######################

### U01 2024-05

# this shows age is predictive of SSNS vs SRNS
# looked at pediatric cases only
time Rscript age-vs-steroid-response.R

# creates clean XLSX for Rasheed, with all discovery and Bristol samples that should have ages of onset (most do but some are too high), sorted by age
time Rscript age-02-make-xlsx.R

# look at dominance at top loci
# creates figures for each locus and a table with p-values and coefficients (with all loci together)
time Rscript dominance.R

### 2025-10 PRS grant

# manually ran this analysis, converted 8 C+T loci (picked under additive model) to dom or rec, but R2 are worse for those resulting PRS (refitting coeffs on Bristol)
prelim-00-rec-ct.R
# repeat for 4 top loci in paper; sadly though the results are better (higher R2), additive > dom > rec again.
prelim-01-rec-top.R
# repeated for top HLA haplotype; more of the same, though this model is worst overall
prelim-02-hla-top-haplotype.R
# repeat with HLA types (singles); still, the same
prelim-03-hla-types.R

# limited test of epistasis, using top loci again
prelim-10-epistasis-top.R

# repeat non-prs(trained models) tests on discovery data! (Bristol is much too small)
# same patterns as in Bristol (additive > dom > rec)
prelim-20-disc-rec-top.R
# VICTORY!  For top haplotype, dom > add > rec!  Difference is small though, but significant
prelim-21-disc-hla-top-haplotype.R
# DITTO!  For all 5 types, dom > add > rec, and significantly so for 3/5
prelim-22-disc-hla-types.R

# epistasis on discovery data
# using top loci, we see significant interactions!
prelim-30-disc-epistasis-top.R
# using top hla types, also see some significant interactions!
prelim-31-disc-epistasis-hla-types.R

# apply domrec to the old base data
cd base
time domrec mac20 mac20-rec rec
# 23m53.931s DCC
ln -s mac20.fam mac20-rec.fam
ln -s mac20.bim mac20-rec.bim
time domrec mac20 mac20-dom dom
# 23m58.953s DCC
ln -s mac20.fam mac20-dom.fam
ln -s mac20.bim mac20-dom.bim

# to run SAIGE, need to link covariates file
# (use same additive PCs for dom/rec, it'd be confusing to change them)
ln -s ../../../saige/ssns_ctrl/covar_ssns_ctrl.txt .

# done at first submission, not totally sure if needed or not; subsequent submissions were ok without this
# source /hpc/group/ochoalab/tt207/miniconda3/etc/profile.d/conda.sh
# conda activate RSAIGE

# submit job!
name=mac20-rec; sbatch -J saige-$name -o saige-$name.out --export=name=$name saige-alex.q
# 364m24.425s/8222m19.559s DCC step 1
# 247m25.098s DCC step 2
name=mac20-dom; sbatch -J saige-$name -o saige-$name.out --export=name=$name saige-alex.q
########## RUNNING

# are top p-values better?  are the top loci different?
# add vs rec: add is better practically always, but there are 3 loci outside chr6 where rec is much better
prelim-40-saige-domrec.R


#########################
### PRS for discovery ###
#########################

# apply trained PRS to discovery individuals, for testing out some ideas

# first use LDPred2 to calculate them, reusing our previous pipeline as much as possible!

# the existing sets don't work if we want scores for absolutely everybody, so let's start all over here
mkdir discovery
cd discovery
# link the most complete discovery ns_ctrl data, still excludes Bristol but meh
ln -s ../../mac20.bed .
ln -s ../../mac20.bim .
ln -s ../../mac20.fam .
cd ..

# match trained scores to discovery
# this is exactly what we want!  Just treat discovery as a testing set
time Rscript prs-new-08-match-train-test.R base train discovery
# 528,964 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 528,964 variants have been matched; 0 were flipped and 0 were reversed.
# 2m30.370s DCC

# make RDS just before scoring
time Rscript prs-new-04-make-rds.R discovery/mac20
# 17m9.991s DCC

# finally, get PRS ("scores") but skip correlations (so also don't need true traits or PCs), only here that case makes sense
base=base; train=train; test=discovery; just_score=1; sbatch -J ldpred-02-score-$base-$train-$test -o ldpred-02-score-$base-$train-$test.out --export=base=$base,train=$train,test=$test,just_score=$just_score ldpred-02-score.q
# 38m45.241s/2m7.509s DCC (there's something weird about how slow this was, but oh well)

# now merge data of interest into simple table for myself and collaborators!
time Rscript ldpred-18-report-prs.R

# cleanup
cd discovery
rm mac20.{bk,rds}

# locally validate by contrasting PRS to diagnosis labels
time Rscript ldpred-19-report-prs-validate.R
