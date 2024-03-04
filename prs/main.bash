###############
### LDPRED2 ###
###############

# based on this tutorial
# https://privefl.github.io/bigsnpr/articles/LDpred2.html

# start interactive shell if needed
srun --mem 16G -p ochoalab --account ochoalab --pty bash -i
module load R/4.1.1-rhel8 
module load Plink/2.00a3LM

# define subsets to split Discovery mainly, but also cleans up Bristol minimally
# (lists of individuals in each set, */ids.txt below)
time Rscript prs-new-00-create-subsets.R

# actually create data
cd /datacommons/ochoalab/ssns_gwas/imputed/prs-new
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
    ln -s {../../$name/,}mac20-glmm-score.txt
    # go back down
    cd ..
done

# the same for curegn as a training dataset
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
cd ..

# since cureGN has low amounts of missingness, as a hack, apply some semi-random imputation (based on global AFs, but unlikely to matter)
time Rscript impute-dumb.R train-curegn/mac20
# 7m42.698s DCC
# NOTE: we will use original PCs estimated from data with missingnes, which is better; nothing else is affected in this case (we're not calculating LD for CureGN)

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

# runs GMMAT on new "base" only
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
time Rscript bim-add-posg.R base-ssns_ctrl/mac20 38
# 3m11.185s DCC
time Rscript bim-add-posg.R base-ssns_srns/mac20 38
# 1m53.496s DCC
time Rscript bim-add-posg.R train-curegn/mac20 38
# 1m14.851s DCC

# create RDS versions of train and test (needed for both)
time Rscript prs-new-04-make-rds.R train/mac20
# 1m42.520s DCC
time Rscript prs-new-04-make-rds.R test/mac20
# 1m38.067s DCC
# do base too, to calculate LD in it
time Rscript prs-new-04-make-rds.R base/mac20
# 16m13.505s DCC
time Rscript prs-new-04-make-rds.R base-ssns_ctrl/mac20
# 16m49.657s DCC
time Rscript prs-new-04-make-rds.R base-ssns_srns/mac20
# 3m15.834s DCC
time Rscript prs-new-04-make-rds.R train-curegn/mac20
# 1m33.345s DCC

# clean summary stats (convert scores to betas, subset to array)
time Rscript prs-new-05-sumstats-clean.R base
# Array has these many variants: 761366
# 20,511,795 variants to be matched.
# 82,354 ambiguous SNPs have been removed.
# 672,360 variants have been matched; 0 were flipped and 0 were reversed.
# 2m10.932s DCC
time Rscript prs-new-05-sumstats-clean.R base-ssns_ctrl
# Array has these many variants: 761366
# 20,838,869 variants to be matched.
# 82,356 ambiguous SNPs have been removed.
# 672,362 variants have been matched; 0 were flipped and 0 were reversed.
# 2m28.642s DCC
time Rscript prs-new-05-sumstats-clean.R base-ssns_srns
# Array has these many variants: 761366
# 12,274,557 variants to be matched.
# 76,584 ambiguous SNPs have been removed.
# 635,789 variants have been matched; 0 were flipped and 0 were reversed.
# 1m25.307s DCC

# subset again to match to training data, and to self-base in all cases to get LD in base
time Rscript prs-new-06-sumstats-match.R base train
# 672,360 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 528,964 variants have been matched; 0 were flipped and 0 were reversed.
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
# 512,349 variants have been matched; 308 were flipped and 881 were reversed.
# 1m25.818s DCC
time Rscript prs-new-06-sumstats-match.R base-ssns_srns train-curegn
# 635,789 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 511,666 variants have been matched; 307 were flipped and 881 were reversed.
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
# 466,208 variants have been matched; 0 were flipped and 0 were reversed.
# 1m32.168s DCC
time Rscript prs-new-08-match-train-test.R base-ssns_ctrl train-curegn test
# 512,349 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 470,424 variants have been matched; 213 were flipped and 853 were reversed.
# 1m32.108s DCC
time Rscript prs-new-08-match-train-test.R base-ssns_srns train-curegn test
# 511,666 variants to be matched.
# 0 ambiguous SNPs have been removed.
# 469,868 variants have been matched; 212 were flipped and 853 were reversed.
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

# run grid version, which is more computationally intensive
base=base; sbatch -J ldpred-03-grid-$base -o ldpred-03-grid-$base.out --export=base=$base ldpred-03-grid.q
# 34m41.436s DCC
base=base-ssns_ctrl; sbatch -J ldpred-03-grid-$base -o ldpred-03-grid-$base.out --export=base=$base ldpred-03-grid.q
# 35m23.125s DCC
base=base-ssns_srns; sbatch -J ldpred-03-grid-$base -o ldpred-03-grid-$base.out --export=base=$base ldpred-03-grid.q
# 31m46.973s DCC

# fit parameters using training data (grid version)
time Rscript ldpred-04-grid-fit.R base train
# 3m6.011s DCC
time Rscript ldpred-04-grid-fit.R base-ssns_ctrl train-curegn
# 2m58.704s DCC
time Rscript ldpred-04-grid-fit.R base-ssns_srns train-curegn
# 2m52.887s DCC

# run auto version
base=base; sbatch -J ldpred-05-auto-$base -o ldpred-05-auto-$base.out --export=base=$base ldpred-05-auto.q
# 8m8.392s DCC
base=base-ssns_ctrl; sbatch -J ldpred-05-auto-$base -o ldpred-05-auto-$base.out --export=base=$base ldpred-05-auto.q
# 7m58.486s DCC
base=base-ssns_srns; sbatch -J ldpred-05-auto-$base -o ldpred-05-auto-$base.out --export=base=$base ldpred-05-auto.q
# 6m53.674s DCC

# base version needs one more step to map data to training SNOs
time Rscript ldpred-05-auto-map-base-to-train.R base train
# 0m12.843s DCC
time Rscript ldpred-05-auto-map-base-to-train.R base-ssns_ctrl train-curegn
# 0m12.028s DCC
time Rscript ldpred-05-auto-map-base-to-train.R base-ssns_srns train-curegn
# 0m11.622s DCC

# run lassosum version
base=base; sbatch -J ldpred-06-lassosum-$base -o ldpred-06-lassosum-$base.out --export=base=$base ldpred-06-lassosum.q
# 4m3.944s DCC
base=base-ssns_ctrl; sbatch -J ldpred-06-lassosum-$base -o ldpred-06-lassosum-$base.out --export=base=$base ldpred-06-lassosum.q
# 3m41.389s DCC
base=base-ssns_srns; sbatch -J ldpred-06-lassosum-$base -o ldpred-06-lassosum-$base.out --export=base=$base ldpred-06-lassosum.q
# 5m3.267s DCC

# fit parameters using training data (lassosum version)
base=base; train=train; sbatch -J ldpred-07-lassosum-fit-$base-$train -o ldpred-07-lassosum-fit-$base-$train.out --export=base=$base,train=$train ldpred-07-lassosum-fit.q
# 2m13.012s DCC
base=base-ssns_ctrl; train=train-curegn; sbatch -J ldpred-07-lassosum-fit-$base-$train -o ldpred-07-lassosum-fit-$base-$train.out --export=base=$base,train=$train ldpred-07-lassosum-fit.q
# 1m55.754s DCC
base=base-ssns_srns; train=train-curegn; sbatch -J ldpred-07-lassosum-fit-$base-$train -o ldpred-07-lassosum-fit-$base-$train.out --export=base=$base,train=$train ldpred-07-lassosum-fit.q
# 2m2.664s DCC

# get correlation values that actually reveal which value was best
base=base; train=train; test=test; sbatch -J ldpred-02-score-$base-$train-$test -o ldpred-02-score-$base-$train-$test.out --export=base=$base,train=$train,test=$test ldpred-02-score.q
# 1m1.389s DCC
base=base-ssns_ctrl; train=train-curegn; test=test; sbatch -J ldpred-02-score-$base-$train-$test -o ldpred-02-score-$base-$train-$test.out --export=base=$base,train=$train,test=$test ldpred-02-score.q
# 1m28.522s DCC
base=base-ssns_srns; train=train-curegn; test=test; sbatch -J ldpred-02-score-$base-$train-$test -o ldpred-02-score-$base-$train-$test.out --export=base=$base,train=$train,test=$test ldpred-02-score.q
# 1m2.556s DCC

# make nice plot that summarize testing results
time Rscript ldpred-08-test-plot.R base train test
# 0m8.598s DCC
time Rscript ldpred-08-test-plot.R base-ssns_ctrl train-curegn test
# 0m14.383s DCC
time Rscript ldpred-08-test-plot.R base-ssns_srns train-curegn test
# 0m6.913s DCC

# makes a single plot combining the data from the previous three
time Rscript ldpred-09-test-plot-combined.R
