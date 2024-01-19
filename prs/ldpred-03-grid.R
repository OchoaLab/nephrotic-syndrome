library(bigsnpr)
library(readr)
library(genio)

# constants
# best value from ldpred2-inf scan, and also other previous ldpred2-grid tests.  this isn't single value tested, but tests are centered around it
h2_est <- 0.1
NCORES <- 10

# load filtered sumstats `df_beta`!
df_beta <- read_tsv( 'betas-ssns_ctrl-array.txt.gz', show_col_types = FALSE )

# load `corr` data/backing file, matching SNPs in betas
load( 'data-corr.RData' )

# set up full grid of parameters
h2_seq <- round(h2_est * c(0.3, 0.7, 1, 1.4), 4)
p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2)
params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE))

# not sure if this is random, just in case
set.seed(1)

# actual run
betas_grid <- snp_ldpred2_grid( corr, df_beta, params, ncores = NCORES )

# store results
write_matrix( paste0( 'betas-ldpred2-grid-h', h2_est ), betas_grid, ext = 'txt.gz' )
# params needed to interpret correctly
write_tsv( params, paste0( 'params-ldpred2-grid-h', h2_est, '.txt.gz' ) )
