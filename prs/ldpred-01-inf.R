library(bigsnpr)
library(readr)
library(genio)
library(ochoalabtools)

# constants
# sequence of heritabilities to consider
# run a finer grid for smaller values after seeing that the perform best
# (doing it now is best because this method is much faster than the following ones)
herits <- c( (1:9)/100, (1:9)/10 )
# support old data for now
types_old <- c('ssns_ctrl', 'ssns_srns')

# support old data for now, expect ssns_ctrl or ssns_srns
args <- args_cli()
type_base <- args[1]
type_train <- args[2]
if ( is.na( type_base ) )
    stop( 'Usage: <type>' )

# extra step for new cases only
if ( ! type_base %in% types_old ) {
    if ( is.na( type_train ) )
        stop( 'Usage: <type_base> <type_train>' )
    # all processing happens in subdirectory
    setwd( type_train )
}

# paths
file_in <- paste0( 'betas-', type_base, '-clean-matched.txt.gz' )
file_ld <- paste0( 'ld-', type_base, '.RData' )
name_out <- paste0( 'betas-', type_base, '-ldpred2-inf' )

# load filtered sumstats `df_beta`!
df_beta <- read_tsv( file_in, show_col_types = FALSE )

# load `ld` data/backing file, matching SNPs in betas
load( file_ld )

# let's gather output into a single matrix for several heritabilities
betas_grid <- matrix( NA, nrow = nrow( df_beta ), ncol = length( herits ) )
    
# this is fast, just scan a grid of heritability values to pick a decent one!
for ( i in 1 : length( herits ) ) {
    # not sure if this is random, just in case
    set.seed(1)

    # key calculation!
    betas_grid[ , i ] <- snp_ldpred2_inf( ld, df_beta, h2 = herits[ i ] )
}

# store results
write_matrix( name_out, betas_grid, ext = 'txt.gz' )
# in this case the heritability values are so trivial we don't save them
