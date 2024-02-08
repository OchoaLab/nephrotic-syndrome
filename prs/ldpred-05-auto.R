library(bigsnpr)
library(readr)
library(ochoalabtools)

# constants
# best value from ldpred2-inf scan, and also other previous ldpred2-grid tests.  this isn't single value tested, but tests are centered around it
h2_est <- 0.1
NCORES <- 10
# reduce this up to 0.4 if you have some (large) mismatch with the LD ref
coef_shrink <- 0.95
p_seq <- seq_log(1e-4, 0.2, length.out = 30)
# support old data for now
types_old <- c('ssns_ctrl', 'ssns_srns')

# determine which type to run
type <- args_cli()[1]
if ( is.na( type ) )
    stop( 'Usage: <type>' )

# handle old and new cases!
if ( type %in% types_old ) {
    # add a dash to separate parts of path as needed
    type_in <- paste0( '-', type )
} else {
    # all processing happens in subdirectory
    setwd( type )
    # new setup, type isn't used, this will interpolate fine in all cases below
    type_in <- ''
}

# load filtered sumstats `df_beta`!
df_beta <- read_tsv( paste0( 'betas', type_in, '-clean-matched.txt.gz' ), show_col_types = FALSE )

# load `ld` data/backing file, matching SNPs in betas
load( paste0( 'ld', type_in, '.RData' ) )

# this is random, make it reproducible
set.seed(1)

# actual run
multi_auto <- snp_ldpred2_auto(
    ld,
    df_beta,
    h2_init = h2_est,
    vec_p_init = p_seq,
    ncores = NCORES,
    # use_MLE = FALSE,  # uncomment if you have convergence issues or when power is low (need v1.11.9)
    allow_jump_sign = FALSE,
    shrink_corr = coef_shrink
)

# output `multi_auto` is a complicated list, this cleans it up to return a simple vector of betas:
# skipped chain visualization
# apply filters though

# find range of correlation estimates
range_corr <- sapply( multi_auto, function( auto ) diff( range( auto$corr_est ) ) )
# identify cases to exclude due to extreme correlation ranges
indexes_keep <- which( range_corr > ( 0.95 * quantile( range_corr, 0.95, na.rm = TRUE ) ) )
# apply filter, extract betas, average across remaining cases to produce final beta vector
betas <- rowMeans( sapply( multi_auto[ indexes_keep ], function( auto ) auto$beta_est ) )

# store results
file_out <- paste0( 'betas', type, '-ldpred2-auto-h', h2_est, '.txt.gz' )
write_lines( betas, file_out )



##########################################

## ### ldpred2-auto, pt 2

## # again, none of this worked for ssns-ctrl sumstats; nothing else was run for this case! (but it was saved)

## # reuses up to `keep` above
## all_h2 <- sapply(multi_auto[keep], function(auto) tail(auto$path_h2_est, 500))
## quantile(all_h2, c(0.5, 0.025, 0.975))
## ##       50%      2.5%     97.5% 
## ## 0.7597156 0.1502886 1.3498973 # ssns-srns
## all_p <- sapply(multi_auto[keep], function(auto) tail(auto$path_p_est, 500))
## quantile(all_p, c(0.5, 0.025, 0.975))
## ##         50%        2.5%       97.5% 
## ## 0.005180223 0.000555023 0.014349786 # ssns-srns
## all_alpha <- sapply(multi_auto[keep], function(auto) tail(auto$path_alpha_est, 500))
## quantile(all_alpha, c(0.5, 0.025, 0.975))
## ##        50%       2.5%      97.5% 
## ## -1.2123923 -1.5000000 -0.5682036 # ssns-srns
## ## # this fails, unclear why, but resulting bsamp is list with "0-column" matrices
## ## bsamp <- lapply(multi_auto[keep], function(auto) auto$sample_beta)
## ## all_r2 <- do.call("cbind", lapply(seq_along(bsamp), function(ic) {
## ##     b1 <- bsamp[[ic]]
## ##     Rb1 <- apply(b1, 2, function(x)
## ##         coef_shrink * bigsparser::sp_prodVec(ld, x) + (1 - coef_shrink) * x)
## ##     b2 <- do.call("cbind", bsamp[-ic])
## ##     b2Rb1 <- as.matrix(Matrix::crossprod(b2, Rb1))
## ## }))
## ## quantile(all_r2, c(0.5, 0.025, 0.975))

## beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
## pred_auto <- big_prodVec(G, beta_auto, ind.col = df_beta[["_NUM_ID_"]])
## pcor(pred_auto, y, NULL)^2
## ## [1] 0.017651880 0.002007581 0.047902414 # ssns-srns
## ## cor(pred_auto, y)^2

## # skipped bits about fine-mapping, revisit when data is improved!

