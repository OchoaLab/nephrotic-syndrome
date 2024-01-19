library(bigsnpr)
library(readr)
library(genio)

# constants
NCORES <- 10

# load filtered sumstats `df_beta`!
df_beta <- read_tsv( 'betas-ssns_ctrl-array.txt.gz', show_col_types = FALSE )

# load `corr` data/backing file, matching SNPs in betas
load( 'data-corr.RData' )

# not sure if this is random, just in case
set.seed(1)

# actual run
betas_grid <- snp_lassosum2( corr, df_beta, ncores = NCORES )
# save this data frame of parameters
params <- attr( betas_grid, "grid_param" )

# store results
write_matrix( 'betas-ldpred2-lassosum', betas_grid, ext = 'txt.gz' )
# params needed to interpret correctly
write_tsv( params, 'params-ldpred2-lassosum.txt.gz' )
