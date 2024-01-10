# start interactive shell
# cd /datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/prs
# srun --mem 16G -p ochoalab --account ochoalab --pty bash -i
# module load R/4.1.1-rhel8
# R

# based on this tutorial
# https://privefl.github.io/bigsnpr/articles/LDpred2.html

library(bigsnpr)
library(genio)
library(ochoalabtools)
library(ggplot2)
library(dplyr)
library(readr)

# run from here
setwd('/datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/prs')
# name of data to predict on
name <- 'data'
#file_phen <- '../../../../bristol_covar.csv'
file_phen <- '/datacommons/ochoalab/ssns_gwas/array/patient-data.txt.gz' 
# base summary stats
#file_sumstats <- '/datacommons/ochoalab/ssns_gwas/imputed/ssns_srns/mac20-glmm-score.txt'
file_sumstats <- '/datacommons/ochoalab/ssns_gwas/imputed/ssns_ctrl/mac20-glmm-score.txt'
# file to subset to clean array SNPs (which passed QC already)
name_array <- '/datacommons/ochoalab/ssns_gwas/array/ssns_tgp_merge_clean'

# Tiffany said this PCs file is updated
# '/datacommons/ochoalab/ssns_gwas/replication/bristol_data/GMMAT/ssns_srns/ssns_srns_mac20.eigenvec'

# load validation dataset
# generate this if it doesn't exist
rds <- paste0( name, '.rds' )
if ( !file.exists( rds ) )
    # this generates .bk and .rds files
    snp_readBed( paste0( name, '.bed' ) )
# Attach the "bigSNP" object in R session
obj.bigSNP <- snp_attach( rds )

# read phen, only file with actual trait
# load patient data, which excludes TGP info
data <- read_tsv( file_phen, show_col_types = FALSE )
#phen <- read_csv( file_phen )
#phen <- read_phen( file_phen )

stopifnot( all( obj.bigSNP$fam$sample.ID %in% data$id ) )

## # these two are perfectly aligned!
## stopifnot( all( phen$id == obj.bigSNP$fam$sample.ID ) )

# subset and reorder
indexes <- match( obj.bigSNP$fam$sample.ID, data$id )
# reorder phen
data <- data[ indexes, ]
stopifnot( all( obj.bigSNP$fam$sample.ID == data$id ) )

# make numerical version
data$pheno <- NA
data$pheno[ data$diagnosis == 'SRNS' ] <- 1
data$pheno[ data$diagnosis == 'SSNS' ] <- 0

## table( phen$pheno ) # oldest stuff, which is obsolete; it appears to have some bug, all other versions agree with each other but not this one
##   0   1 
## 327 163
## table( phen$diagnosis ) # file Tiffany actually wanted me to use most recently
## NS UNCLASSIFIED             SNS            SRNS            SSNS 
##              70               1             149             364 
## table( data$diagnosis ) # my re-cleaned up data, after subsetting
## NS UNCLASSIFIED            SRNS            SSNS 
##              70             149             365 
## table( data$pheno ) # numerical version of above
##   0   1 
## 365 149 

# transfer to this other object
obj.bigSNP$fam$affection <- data$pheno

# Get aliases for useful slots
G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
y   <- obj.bigSNP$fam$affection
#NCORES <- nb_cores() # on DCC this is 1 with default settings, which is probably fine because ours is a small dataset
NCORES <- 1 # was needed locally because of automatic blas multithreading leading to two levels of multithreading, something like that
#NCORES <- 20 # for DCC run for precalculating LD matrix
message( NCORES )

# rename columns to what snp_match wants
map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))

# load summary statistics (big file because it's from imputed data)
sumstats <- bigreadr::fread2( file_sumstats )
# apply Debo's transformation
sumstats$beta <- sumstats$SCORE / sumstats$VAR
sumstats$beta_se <- 1 / sqrt( sumstats$VAR )
# OLD: exclude "cM", "converged"
#sumstats <- setNames( sumstats[ -c(3, 12) ], c('chr', 'rsid', 'pos', 'a1', 'a0', 'n_eff', 'af', 'beta', 'beta_se', 'p') )
# exclude "MISSRATE" (there's no missingness here), SCORE, VAR
stopifnot( all( sumstats$MISSRATE == 0 ) )
# NOTE: a1, a0 could be reversed!
#sumstats <- setNames( sumstats[ -c(7, 9, 10) ], c('rsid', 'chr', 'pos', 'a1', 'a0', 'n_eff', 'af', 'p', 'beta', 'beta_se') )
sumstats <- setNames( sumstats[ -c(7, 9, 10) ], c('rsid', 'chr', 'pos', 'a0', 'a1', 'n_eff', 'af', 'p', 'beta', 'beta_se') )

# if using ssns-ctrl, reverse signs!  (ssns-srns was fine though)
sumstats$beta <- -sumstats$beta

# load array SNP set
bim <- read_bim( name_array )
# change names to match snp_match
#names( bim ) <- c('chr', 'rsid', 'a0', 'a1')
names( bim )[ names( bim ) == 'id' ] <- 'rsid'
names( bim )[ names( bim ) == 'alt' ] <- 'a1'
names( bim )[ names( bim ) == 'ref' ] <- 'a0'
# remove unneeded column
bim[ names( bim ) == 'posg' ] <- NULL
# convert chr to integer
bim$chr <- as.integer( bim$chr )
# find matching subset
sumstats2 <- snp_match( sumstats, bim )
### ssns-srns
## 12,274,557 variants to be matched.
## 76,584 ambiguous SNPs have been removed.
## 635,789 variants have been matched; 0 were flipped and 0 were reversed.
### ssns-ctrl OLD and NEWEST
## 20,838,869 variants to be matched.
## 82,356 ambiguous SNPs have been removed.
## 672,362 variants have been matched; 0 were flipped and 0 were reversed.
nrow( bim ) # 761,366 # for reference

# before next round, correct some funny business with snp_match's output
sumstats2[ names( sumstats2 ) == 'rsid' ] <- NULL # IDs from BIM, don't want
names( sumstats2 )[ names( sumstats2 ) == 'rsid.ss' ] <- 'rsid' # restore original IDs
sumstats2[ names( sumstats2 ) == '_NUM_ID_.ss' ] <- NULL # new column added, also don't want
sumstats2[ names( sumstats2 ) == '_NUM_ID_' ] <- NULL # new column added, also don't want

# map assuming things are largely aligned already
df_beta <- snp_match( sumstats2, map )
### ssns-srns
## 635,789 variants to be matched.
## 0 ambiguous SNPs have been removed.
## 492,982 variants have been matched; 0 were flipped and 0 were reversed.
### ssns-ctrl OLD and NEWEST
## 672,362 variants to be matched.
## 0 ambiguous SNPs have been removed.
## 495,160 variants have been matched; 0 were flipped and 0 were reversed. # OLD
## 508,082 variants have been matched; 0 were flipped and 0 were reversed. # NEWEST

# here processing suggests QC on sumstats, but no code is provided (there's equations and a massive repo, may have to revisit)

# corr backing file, to reload if needed
corr_file <- paste0( name, '-corr' )
# NOTE: code blindly adds another .sbk extension if I specify one (leading to double '.sbk.sbk"), so don't include it
# these ones have added extension for use internally
corr_file_rdata <- paste0( corr_file, '.RData' )

if ( file.exists( corr_file_rdata ) ) {
    # load from a different session
    load( corr_file_rdata )
} else {
    # convert basepairs to genetic position this way
    POS2 <- snp_asGeneticPos(CHR, POS, ncores = NCORES)
    
    # compute correlations (LD matrices) on the fly
    for (chr in 1:22) {
        message(chr)
        ## indices in 'df_beta'
        ind.chr <- which(df_beta$chr == chr)
        ## indices in 'G'
        ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
        # calculate sparse symmetric correlation matrix
        corr0 <- snp_cor(G, ind.col = ind.chr2, size = 3 / 1000,
                         infos.pos = POS2[ind.chr2], ncores = NCORES)
        # calculate LD scores
        ld0 <- Matrix::colSums(corr0^2)
        # combine data across chromosomes
        if (chr == 1) {
            ld <- ld0
            # set up backing file, first time only
            corr <- as_SFBM(corr0, corr_file, compact = TRUE)
        } else {
            ld <- c(ld, ld0)
            corr$add_columns(corr0, nrow(corr))
        }
    }
    # save object, keep backing file, so we can load it in another session
    save( corr, ld, file = corr_file_rdata )
}

# size of LD data, in GB
file.size(corr$sbk) / 1024^3  # file size in GB
#[1] 4.118546 # ssns-srns
#[1] 4.153004 # ssns-ctrl OLD
#[1] 4.373985 # ssns-ctrl NEWEST

# decide how many individuals are in each subset... example had 503 indivs total, our data has slightly fewer, 490
# only used for LDpred2-grid and lassosum2
set.seed(1)
ind.val <- sample( nrow(G), nrow(G) * 0.7 )   # example 350 (70%)
ind.test <- setdiff( rows_along(G), ind.val ) # example 153 (30%)
# NEWEST:
length(ind.val)  # [1] 408
length(ind.test) # [1] 176


### ldpred2-inf

# Estimate of h2 from LD Score regression
ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2, sample_size = n_eff, blocks = NULL))
h2_est <- ldsc[["h2"]]
h2_est
# [1] 1.171371 # ssns-ctrl OLD
# [1] 1.245419 # ssns-ctrl NEWEST
# NEWEST only, this is BS, let's set heritability to something on the high end but not above 1 when this happens
if ( h2_est > 1 )
    h2_est <- 0.8
beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
pred_inf <- big_prodVec(G, beta_inf, ind.row = ind.test, ind.col = df_beta[["_NUM_ID_"]])
pcor(pred_inf, y[ind.test], NULL)
# [1]  0.04728706 -0.11549036  0.20759114 # ssns-srns
# [1]  0.19839389  0.03771186  0.34907536 # ssns-ctrl OLD
# [1]  0.10088376 -0.05563708  0.25256512 # ssns-ctrl NEWEST
#cor(pred_inf, y[ind.test])

### ldpred2-grid

#h2_seq <- round(h2_est * c(0.3, 0.7, 1, 1.4), 4) # OLD
h2_seq <- c(0.1, 0.3, 0.5, 0.7, 0.9)
p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2)
params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE))
beta_grid <- snp_ldpred2_grid(corr, df_beta, params, ncores = NCORES)
pred_grid <- big_prodMat(G, beta_grid, ind.col = df_beta[["_NUM_ID_"]])
params$score <- apply(pred_grid[ind.val, ], 2, function(x) {
  if (all(is.na(x))) return(NA)
  summary(lm(y[ind.val] ~ x))$coef["x", 3]
  # summary(glm(y[ind.val] ~ x, family = "binomial"))$coef["x", 3]
})

fig_start( 'ldpred2-grid', width = 6 )
ggplot(params, aes(x = p, y = score, color = as.factor(h2))) +
    theme_bigstatsr() +
    geom_point() +
    geom_line() +
    scale_x_log10(breaks = 10^(-5:0), minor_breaks = params$p) +
    facet_wrap(~ sparse, labeller = label_both) +
    labs(y = "GLM Z-Score", color = "h2") +
    theme(legend.position = "top", panel.spacing = unit(1, "lines"))
fig_end()

params %>%
    mutate(sparsity = colMeans(beta_grid == 0), id = row_number()) %>%
    arrange(desc(score)) %>%
    mutate_at(c("score", "sparsity"), round, digits = 3) %>%
    slice(1:10)
### ssns-srns
##          p     h2 sparse score sparsity  id
## 1  1.0e-05 0.3391  FALSE 5.059    0.000   1
## 2  5.6e-05 0.7913   TRUE 4.930    0.921 109
## 3  3.2e-05 1.5826   TRUE 4.619    0.962 150
## 4  1.0e-04 0.3391  FALSE 4.322    0.000   5
## 5  1.8e-04 1.5826  FALSE 4.280    0.000  69
## 6  5.6e-05 0.3391  FALSE 4.243    0.000   4
## 7  3.2e-05 0.3391   TRUE 4.185    0.910  87
## 8  5.6e-05 1.1305   TRUE 4.176    0.932 130
## 9  3.2e-05 1.5826  FALSE 4.153    0.000  66
## 10 3.2e-05 1.1305  FALSE 4.150    0.000  45
### ssns-ctrl OLD
##         p     h2 sparse score sparsity id
## 1  0.0056 0.3514   TRUE 3.536    0.562 96
## 2  0.0056 0.3514  FALSE 3.480    0.000 12
## 3  0.0100 0.3514   TRUE 3.431    0.554 97
## 4  0.0100 0.3514  FALSE 3.365    0.000 13
## 5  0.0180 0.3514   TRUE 3.339    0.552 98
## 6  0.0032 0.3514  FALSE 3.309    0.000 11
## 7  0.0032 0.3514   TRUE 3.237    0.582 95
## 8  0.0320 0.3514   TRUE 3.233    0.553 99
## 9  0.0180 0.3514  FALSE 3.219    0.000 14
## 10 0.0320 0.3514  FALSE 3.047    0.000 15
### ssns-ctrl NEWEST
##          p  h2 sparse score sparsity  id
## 1  0.00100 0.1   TRUE 4.036    0.645 114
## 2  0.00100 0.1  FALSE 4.017    0.000   9
## 3  0.00180 0.1   TRUE 3.976    0.625 115
## 4  0.00180 0.1  FALSE 3.973    0.000  10
## 5  0.00320 0.1  FALSE 3.892    0.000  11
## 6  0.00320 0.1   TRUE 3.868    0.613 116
## 7  0.00560 0.1   TRUE 3.800    0.608 117
## 8  0.00560 0.1  FALSE 3.758    0.000  12
## 9  0.00056 0.1  FALSE 3.756    0.000   8
## 10 0.00056 0.1   TRUE 3.685    0.676 113

best_beta_grid <- params %>%
    mutate(id = row_number()) %>%
    # filter(sparse) %>% 
    arrange(desc(score)) %>%
    slice(1) %>%
    print() %>% 
    pull(id) %>% 
    beta_grid[, .]
##       p     h2 sparse    score id
## 1 1e-05 0.3391  FALSE 5.059323  1 # ssns-srns
##        p     h2 sparse    score id
## 1 0.0056 0.3514   TRUE 3.536261 96 # ssns-ctrl OLD
##       p  h2 sparse    score  id
## 1 0.001 0.1   TRUE 4.035512 114 # ssns-ctrl NEWEST

pred <- big_prodVec(G, best_beta_grid, ind.row = ind.test,
                    ind.col = df_beta[["_NUM_ID_"]])
pcor(pred, y[ind.test], NULL)
# [1]  0.12688295 -0.03574459  0.28296373 # ssns-srns
# [1] 0.21203091 0.05191985 0.36151471    # ssns-ctrl OLD
# [1] 0.17061100 0.01537199 0.31781850    # ssns-ctrl NEWEST
#cor(pred, y[ind.test])



### ldpred2-auto

# NOTE ssns-ctrl: not sure what happened, but this whole approach fails (all multi_auto values are NA, everything else dies subsequently).  Ditto OLD and NEWEST

coef_shrink <- 0.95  # reduce this up to 0.4 if you have some (large) mismatch with the LD ref
# takes less than 2 min with 4 cores
multi_auto <- snp_ldpred2_auto(
  corr, df_beta, h2_init = h2_est,
  vec_p_init = seq_log(1e-4, 0.2, length.out = 30), ncores = NCORES,
  # use_MLE = FALSE,  # uncomment if you have convergence issues or when power is low (need v1.11.9)
  allow_jump_sign = FALSE, shrink_corr = coef_shrink)

# skipped chain visualization
# apply filters though
range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE)))
beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
pred_auto <- big_prodVec(G, beta_auto, ind.row = ind.test, ind.col = df_beta[["_NUM_ID_"]])
pcor(pred_auto, y[ind.test], NULL)
# [1]  0.07659266 -0.08637162  0.23556498 # ssns-srns
#cor(pred_auto, y[ind.test])



### lassosum2

beta_lassosum2 <- snp_lassosum2(corr, df_beta, ncores = NCORES)
params2 <- attr(beta_lassosum2, "grid_param")
pred_grid2 <- big_prodMat(G, beta_lassosum2, ind.col = df_beta[["_NUM_ID_"]])
params2$score <- apply(pred_grid2[ind.val, ], 2, function(x) {
    if (all(is.na(x))) return(NA)
    summary(lm(y[ind.val] ~ x))$coef["x", 3]
    # summary(glm(y[ind.val] ~ x, family = "binomial"))$coef["x", 3]
})

# this plot is very different than example, it just looks bad
fig_start( 'lassosum2', width = 6 )
ggplot(params2, aes(x = lambda, y = score, color = as.factor(delta))) +
    theme_bigstatsr() +
    geom_point() +
    geom_line() +
    scale_x_log10(breaks = 10^(-5:0)) +
    labs(y = "GLM Z-Score", color = "delta") +
    theme(legend.position = "top") +
    guides(colour = guide_legend(nrow = 1))
fig_end()

best_grid_lassosum2 <- params2 %>%
    mutate(id = row_number()) %>%
    arrange(desc(score)) %>%
    print() %>% 
    slice(1) %>%
    pull(id) %>% 
    beta_lassosum2[, .]
best_grid_overall <- 
  `if`(max(params2$score, na.rm = TRUE) > max(params$score, na.rm = TRUE),
       best_grid_lassosum2, best_beta_grid)

pred <- big_prodVec(G, best_grid_lassosum2, ind.row = ind.test,
                    ind.col = df_beta[["_NUM_ID_"]])
pcor(pred, y[ind.test], NULL)
# [1]  0.10867889 -0.05416744  0.26589394 # ssns-srns
# [1] 0.17228611 0.01069101 0.32511140    # ssns-ctrl OLD
# [1] 0.17373439 0.01858988 0.32070922    # ssns-ctrl NEWEST
#cor(pred, y[ind.test])

pred <- big_prodVec(G, best_grid_overall, ind.row = ind.test,
                    ind.col = df_beta[["_NUM_ID_"]])
pcor(pred, y[ind.test], NULL)
# [1]  0.12688295 -0.03574459  0.28296373 # ssns-srns
# [1] 0.17228611 0.01069101 0.32511140    # ssns-ctrl OLD
# [1] 0.17061100 0.01537199 0.31781850    # ssns-ctrl NEWEST
#cor(pred, y[ind.test])

save.image()


### ldpred2-auto, pt 2

# again, none of this worked for ssns-ctrl sumstats; nothing else was run for this case! (but it was saved)

# reuses up to `keep` above
all_h2 <- sapply(multi_auto[keep], function(auto) tail(auto$path_h2_est, 500))
quantile(all_h2, c(0.5, 0.025, 0.975))
##       50%      2.5%     97.5% 
## 0.7597156 0.1502886 1.3498973 # ssns-srns
all_p <- sapply(multi_auto[keep], function(auto) tail(auto$path_p_est, 500))
quantile(all_p, c(0.5, 0.025, 0.975))
##         50%        2.5%       97.5% 
## 0.005180223 0.000555023 0.014349786 # ssns-srns
all_alpha <- sapply(multi_auto[keep], function(auto) tail(auto$path_alpha_est, 500))
quantile(all_alpha, c(0.5, 0.025, 0.975))
##        50%       2.5%      97.5% 
## -1.2123923 -1.5000000 -0.5682036 # ssns-srns
## # this fails, unclear why, but resulting bsamp is list with "0-column" matrices
## bsamp <- lapply(multi_auto[keep], function(auto) auto$sample_beta)
## all_r2 <- do.call("cbind", lapply(seq_along(bsamp), function(ic) {
##     b1 <- bsamp[[ic]]
##     Rb1 <- apply(b1, 2, function(x)
##         coef_shrink * bigsparser::sp_prodVec(corr, x) + (1 - coef_shrink) * x)
##     b2 <- do.call("cbind", bsamp[-ic])
##     b2Rb1 <- as.matrix(Matrix::crossprod(b2, Rb1))
## }))
## quantile(all_r2, c(0.5, 0.025, 0.975))

beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
pred_auto <- big_prodVec(G, beta_auto, ind.col = df_beta[["_NUM_ID_"]])
pcor(pred_auto, y, NULL)^2
## [1] 0.017651880 0.002007581 0.047902414 # ssns-srns
## cor(pred_auto, y)^2

# skipped bits about fine-mapping, revisit when data is improved!

