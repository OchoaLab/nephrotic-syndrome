library(bigsnpr)
library(readr)
library(genio)
library(ochoalabtools)

# constants
# best value from ldpred2-inf scan, and also other previous ldpred2-grid tests.  this isn't single value tested, but tests are centered around it
h2_est <- 0.1
NCORES <- 10

args <- args_cli()
base <- args[1]
if ( is.na( base ) )
    stop( 'Usage: <base>' )

# all processing happens in subdirectory
setwd( base )

# paths
file_in <- paste0( 'betas-', base, '-clean-matched.txt.gz' )
file_ld <- paste0( 'ld-', base, '.RData' )
name <- paste0( base, '-ldpred2-grid-h', h2_est )
name_out <- paste0( 'betas-', name )
file_params <- paste0( 'params-', name, '.txt.gz' )

# load filtered sumstats `df_beta`!
df_beta <- read_tsv( file_in, show_col_types = FALSE )

# load `ld` data/backing file, matching SNPs in betas
load( file_ld )

# set up full grid of parameters
h2_seq <- round(h2_est * c(0.3, 0.7, 1, 1.4), 4)
p_seq <- signif(seq_log(1e-5, 1, length.out = 21), 2)
params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE))

# not sure if this is random, just in case
set.seed(1)

# actual run
betas_grid <- snp_ldpred2_grid( ld, df_beta, params, ncores = NCORES )

# store results
write_matrix( name_out, betas_grid, ext = 'txt.gz' )
# params needed to interpret correctly
write_tsv( params, file_params )
