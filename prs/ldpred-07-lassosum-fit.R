library(bigsnpr)
library(genio)
library(ochoalabtools)
library(tidyverse)

# constants
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
    name_data <- 'data'
    # load indexes of training individuals
    ind_train <- as.numeric( read_lines( 'ind-train.txt.gz' ) )
} else {
    # start in this subdirectory
    setwd( type )
    # new setup, type isn't used, this will interpolate fine in all cases below
    type_in <- ''
    name_data <- 'mac20'
}

# start by loading several previous calculations and other needed data

# suffix shared by several in/out files
name_in <- paste0( type_in, '-ldpred2-lassosum' )

# load previously calculated results
# (could be in train or base)
betas_grid <- read_matrix( paste0( 'betas', name_in ), ext = 'txt.gz' )
params <- read_tsv( paste0( 'params', name_in, '.txt.gz' ), show_col_types = FALSE )

# rest always happens in training set
# (if type == 'train', then we're already there, and for old types we don't move either)
if ( type == 'base' )
    setwd( '../train' )

# load filtered sumstats `df_beta`!
df_beta <- read_tsv( paste0( 'betas', type_in, '-clean-matched.txt.gz' ), show_col_types = FALSE )

if ( type == 'base' ) {
    # if betas came from base, we need to subset them to what we actually have in training data
    # subset betas using precalculated map of SNPs from 'base' into 'train'!
    betas_grid <- betas_grid[ df_beta[["_NUM_ID_.ss"]], ]
    # also prepend to "name_in" to denote this alternative origin of data, which from here on is only used in outputs!
    name_in <- paste0( '-', type, name_in )
}

# load training dataset
# Attach the "bigSNP" object in R session
obj.bigSNP <- snp_attach( paste0( name_data, '.rds' ) )
G <- obj.bigSNP$genotypes
y <- obj.bigSNP$fam$affection
y[ y == 0 ] <- NA # in plink format, zeros are missing, translate appropriately here!
if ( type_in == '' ) {
    # in new setup, we haven't defined ind_train, do it now that we know the number of individuals (use all)
    ind_train <- 1L : length( y )
}
# load PCs, to condition on when scoring with R2
PCs <- read_eigenvec( name_data )$eigenvec

# calculate PRS for training individuals only
pred_grid <- big_prodMat( G, betas_grid, ind.row = ind_train, ind.col = df_beta[["_NUM_ID_"]] )

# use training individuals to score grid values using correlation, adjusting for PCs, determine which is best
params$cor <- apply( pred_grid, 2, function( x ) pcor( x, y[ ind_train ], PCs[ ind_train, ] )[1] )

# add more info, sort by correlation
params <- params %>%
    mutate( id = row_number() ) %>%
    arrange( desc( cor ) )

# save this table of results!
name_out <- paste0( 'eval', name_in )
write_tsv( params, paste0( name_out, '.txt.gz' ) )

# pick out the best set of parameters, to use and score out of sample later!
betas <- betas_grid[, params$id[1] ]
# save betas!
# this is a simple vector
file_out <- paste0( 'betas', name_in, '-best.txt.gz' )
write_lines( betas, file_out )

# plot results
fig_start( name_out, width = 6 )
ggplot( params, aes( x = lambda, y = cor, color = as.factor( delta ) ) ) +
    theme_classic() +
    geom_point() +
    geom_line() +
    scale_x_log10() +
    labs( x = 'LASSO lambda', y = "Correlation to trait", color = "delta" )
fig_end()
