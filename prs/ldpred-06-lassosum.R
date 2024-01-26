library(bigsnpr)
library(readr)
library(genio)
library(ochoalabtools)

# constants
NCORES <- 10

# determine which type to run
type <- args_cli()[1]
if ( is.na( type ) )
    stop( 'Usage: <type: ssns_ctrl or ssns_srns>' )

# load filtered sumstats `df_beta`!
df_beta <- read_tsv( paste0( 'betas-', type, '-clean-matched.txt.gz' ), show_col_types = FALSE )

# load `ld` data/backing file, matching SNPs in betas
load( paste0( 'ld-', type, '.RData' ) )

# not sure if this is random, just in case
set.seed(1)

# actual run
betas_grid <- snp_lassosum2( ld, df_beta, ncores = NCORES )
# save this data frame of parameters
params <- attr( betas_grid, "grid_param" )

# store results
write_matrix( paste0( 'betas-', type, '-ldpred2-lassosum' ), betas_grid, ext = 'txt.gz' )
# params needed to interpret correctly
write_tsv( params, paste0( 'params-', type, '-ldpred2-lassosum.txt.gz' ) )
