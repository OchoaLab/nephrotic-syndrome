library(bigsnpr)
library(readr)
library(genio)
library(ochoalabtools)

# constants
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
name <- paste0( base, '-ldpred2-lassosum' )
name_out <- paste0( 'betas-', name )
file_params <- paste0( 'params-', name, '.txt.gz' )

# load filtered sumstats `df_beta`!
df_beta <- read_tsv( file_in, show_col_types = FALSE )

# load `ld` data/backing file, matching SNPs in betas
load( file_ld )

# not sure if this is random, just in case
set.seed(1)

# actual run
betas_grid <- snp_lassosum2( ld, df_beta, ncores = NCORES )
# save this data frame of parameters
params <- attr( betas_grid, "grid_param" )

# store results
write_matrix( name_out, betas_grid, ext = 'txt.gz' )
# params needed to interpret correctly
write_tsv( params, file_params )
