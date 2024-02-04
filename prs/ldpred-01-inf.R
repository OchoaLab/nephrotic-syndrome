library(bigsnpr)
library(readr)
library(genio)
library(ochoalabtools)

# constants
# sequence of heritabilities to consider
# run a finer grid for smaller values after seeing that the perform best
# (doing it now is best because this method is much faster than the following ones)
herits <- c( (1:9)/100, (1:9)/10 )

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
write_matrix( paste0( 'betas', type, '-ldpred2-inf' ), betas_grid, ext = 'txt.gz' )
# in this case the heritability values are so trivial we don't save them
