# start interactive shell
# srun --mem 16G -p ochoalab --account ochoalab --pty bash -i
# module load R/4.1.1-rhel8

# based on this tutorial
# https://privefl.github.io/bigsnpr/articles/LDpred2.html

library(bigsnpr)
library(genio)
library(ochoalabtools)
library(ggplot2)
library(dplyr)

# run from here
setwd('/datacommons/ochoalab/ssns_gwas/replication/bristol_data/GMMAT')
# name of data to predict on
name <- 'srns_ssns_mac20'
# base summary stats
file_sumstats <- "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/PRS/glmm.wald_srns_ssns.txt.gz"
file_phen <- 'srns_ssns.phen'

# load validation dataset
# generate this if it doesn't exist
rds <- paste0( name, '.rds' )
if ( !file.exists( rds ) )
    # this generates .bk and .rds files
    snp_readBed( paste0( name, '.bed' ) )
# Attach the "bigSNP" object in R session
obj.bigSNP <- snp_attach( rds )

# read phen, only file with actual trait
phen <- read_phen( file_phen )
## table( phen$pheno )
##   0   1 
## 327 163 

# these match in number, not order
indexes <- match( obj.bigSNP$fam$sample.ID, phen$id )
# reorder phen
phen <- phen[ indexes, ]
stopifnot( all( obj.bigSNP$fam$sample.ID == phen$id ) )
obj.bigSNP$fam$affection <- phen$pheno

# Get aliases for useful slots
G   <- obj.bigSNP$genotypes
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos
y   <- obj.bigSNP$fam$affection
NCORES <- nb_cores() # on DCC this is 1 with default settings, which is probably fine because ours is a small dataset
#NCORES <- 1 # was needed locally because of automatic blas multithreading leading to two levels of multithreading, something like that

# decide how many individuals are in each subset... example had 503 indivs total, our data has slightly fewer, 490
# only used for LDpred2-grid and lassosum2
ind.val <- sample( nrow(G), nrow(G) * 0.7 )   # example 350 (70%)
ind.test <- setdiff( rows_along(G), ind.val ) # example 153 (30%)

# rename columns to what snp_match wants
map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
sumstats <- bigreadr::fread2( file_sumstats )
sumstats <- setNames( sumstats[ -c(3, 12) ], c('chr', 'rsid', 'pos', 'a1', 'a0', 'n_eff', 'af', 'beta', 'beta_se', 'p') )

# map assuming things are largely aligned already
df_beta <- snp_match(sumstats, map)
## 13,971 variants to be matched.
## 1,410 ambiguous SNPs have been removed.
## 8,167 variants have been matched; 0 were flipped and 0 were reversed.

# here processing suggests QC on sumstats, but no code is provided (there's equations and a massive repo, may have to revisit)

# convert basepairs to genetic position this way
POS2 <- snp_asGeneticPos(CHR, POS, ncores = NCORES)

# compute correlations (LD matrices) on the fly
for (chr in 1:22) {
    message(chr)
    ## indices in 'df_beta'
    ind.chr <- which(df_beta$chr == chr)
    ## indices in 'G'
    ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
    corr0 <- snp_cor(G, ind.col = ind.chr2, size = 3 / 1000,
                     infos.pos = POS2[ind.chr2], ncores = NCORES)
    if (chr == 1) {
        ld <- Matrix::colSums(corr0^2)
        corr <- as_SFBM(corr0, compact = TRUE)
    } else {
        ld <- c(ld, Matrix::colSums(corr0^2))
        corr$add_columns(corr0, nrow(corr))
    }
}
# yey this is small, probably because our summary stats are so sparse
file.size(corr$sbk) / 1024^3  # file size in GB
#[1] 0.02784605

### ldpred2-inf

# Estimate of h2 from LD Score regression
ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2, sample_size = n_eff, blocks = NULL))
# ok, something is obviously wrong with this absurdly large h2, we need more summary stats!
##          int           h2 
## 1.163067e+01 4.105297e+14 
#h2_est <- ldsc[["h2"]]
h2_est <- 0.5 # choose something nice in the middle, the above fails totally!
beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
## Warning message: (for original h2_est)
## In sp_solve_sym_eigen(A, b, add_to_diag, tol, maxiter) :
##   Estimated error: 0.0119041.
pred_inf <- big_prodVec(G, beta_inf, ind.row = ind.test, ind.col = df_beta[["_NUM_ID_"]])
pcor(pred_inf, y[ind.test], NULL)
# [1] NA NA NA # both ways!
cor(pred_inf, y[ind.test])
# [1] 0.06213717
# [1] 0.115047 # second run!

### ldpred2-grid

h2_seq <- round(h2_est * c(0.3, 0.7, 1, 1.4), 4)
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
##          p   h2 sparse score sparsity  id
## 1  0.00180 0.70   TRUE 1.351    0.011 157
## 2  0.00100 0.70   TRUE 1.351    0.004 156
## 3  0.00100 0.35   TRUE 1.351    0.042 114
## 4  0.00100 0.50   TRUE 1.351    0.087 135
## 5  0.00056 0.35   TRUE 1.351    0.133 113
## 6  0.00056 0.70   TRUE 1.351    0.102 155
## 7  0.00056 0.35  FALSE 1.351    0.000  29
## 8  0.00056 0.50   TRUE 1.351    0.020 134
## 9  0.00032 0.50  FALSE 1.350    0.000  49
## 10 0.00180 0.70  FALSE 1.350    0.000  73
best_beta_grid <- params %>%
  mutate(id = row_number()) %>%
  # filter(sparse) %>% 
  arrange(desc(score)) %>%
  slice(1) %>%
  print() %>% 
  pull(id) %>% 
    beta_grid[, .]
##        p  h2 sparse   score  id
## 1 0.0018 0.7   TRUE 1.35105 157
pred <- big_prodVec(G, best_beta_grid, ind.row = ind.test,
                    ind.col = df_beta[["_NUM_ID_"]])
pcor(pred, y[ind.test], NULL)
# [1] NA NA NA
cor(pred, y[ind.test])
# [1] 0.0666571
# [1] 0.1165977 # second run!

### ldpred2-auto

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
# [1] NA NA NA
cor(pred_auto, y[ind.test])
# [1] 0.1156314

### lassosum2

beta_lassosum2 <- snp_lassosum2(corr, df_beta, ncores = NCORES)
params2 <- attr(beta_lassosum2, "grid_param")

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
# [1] NA NA NA
cor(pred, y[ind.test])
# [1] 0.1162437 # second run only

pred <- big_prodVec(G, best_grid_overall, ind.row = ind.test,
                    ind.col = df_beta[["_NUM_ID_"]])
pcor(pred, y[ind.test], NULL)
# [1] NA NA NA
cor(pred, y[ind.test])
# [1] 0.1165977 # second run only


### ldpred2-auto, pt 2

# reuses up to `keep` above
all_h2 <- sapply(multi_auto[keep], function(auto) tail(auto$path_h2_est, 500))
quantile(all_h2, c(0.5, 0.025, 0.975))
# all of these are way higher for me than they are in example
##      50%     2.5%    97.5% 
## 24.78544 23.91517 25.62840 
all_p <- sapply(multi_auto[keep], function(auto) tail(auto$path_p_est, 500))
quantile(all_p, c(0.5, 0.025, 0.975))
# higher too, I think expected because my data is enriched for most significant cases
##       50%      2.5%     97.5% 
## 0.5365942 0.3914681 0.6757716 
all_alpha <- sapply(multi_auto[keep], function(auto) tail(auto$path_alpha_est, 500))
quantile(all_alpha, c(0.5, 0.025, 0.975))
##       50%      2.5%     97.5% 
## -1.199174 -1.367569 -1.098398
## # this fails, unclear why, but resulting bsamp is list with "0-column" matrices
## bsamp <- lapply(multi_auto[keep], function(auto) auto$sample_beta)
## all_r2 <- do.call("cbind", lapply(seq_along(bsamp), function(ic) {
##   b1 <- bsamp[[ic]]
##   Rb1 <- apply(b1, 2, function(x)
##     coef_shrink * bigsparser::sp_prodVec(corr, x) + (1 - coef_shrink) * x)
##   b2 <- do.call("cbind", bsamp[-ic])
##   b2Rb1 <- as.matrix(Matrix::crossprod(b2, Rb1))
## }))
## quantile(all_r2, c(0.5, 0.025, 0.975))
beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
pred_auto <- big_prodVec(G, beta_auto, ind.col = df_beta[["_NUM_ID_"]])
pcor(pred_auto, y, NULL)^2
# [1] NA NA NA
cor(pred_auto, y)^2
# [1] 0.004926419

# skipped bits about fine-mapping, revisit when data is improved!

