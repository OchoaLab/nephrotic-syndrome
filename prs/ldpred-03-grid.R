library(bigsnpr)
library(readr)
library(genio)
library(ochoalabtools)

# constants
# best value from ldpred2-inf scan, and also other previous ldpred2-grid tests.  this isn't single value tested, but tests are centered around it
h2_est <- 0.1
NCORES <- 10

# determine which type to run
type <- args_cli()[1]

# handle old and new cases!
if ( is.na( type ) ) {
    # new setup, type isn't used, this will interpolate fine in all cases below
    type <- ''
    # all processing happens in subdirectory
    setwd( 'train' )
} else {
    # add a dash to separate parts of path as needed
    type <- paste0( '-', type )
}

# load filtered sumstats `df_beta`!
df_beta <- read_tsv( paste0( 'betas', type, '-clean-matched.txt.gz' ), show_col_types = FALSE )

# load `ld` data/backing file, matching SNPs in betas
load( paste0( 'ld', type, '.RData' ) )

# set up full grid of parameters
h2_seq <- round(h2_est * c(0.3, 0.7, 1, 1.4), 4)
p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2)
params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE))

# not sure if this is random, just in case
set.seed(1)

# actual run
betas_grid <- snp_ldpred2_grid( ld, df_beta, params, ncores = NCORES )

# store results
write_matrix( paste0( 'betas', type, '-ldpred2-grid-h', h2_est ), betas_grid, ext = 'txt.gz' )
# params needed to interpret correctly
write_tsv( params, paste0( 'params', type, '-ldpred2-grid-h', h2_est, '.txt.gz' ) )
