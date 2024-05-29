library(bigsnpr)
library(genio)
library(ochoalabtools)
library(tidyverse)

# constants
# previously-chosen h2 estimate around which grid search was performed
h2_est <- 0.1
name_data <- 'mac20'

# determine which type to run
args <- args_cli()
base <- args[1]
train <- args[2]
if ( is.na( train ) )
    stop( 'Usage: <base> <train>' )

# start at base
setwd( base )

# start by loading several previous calculations and other needed data

# suffix shared by several in/out files
name_in <- paste0( base, '-ldpred2-grid-h', h2_est )

# load previously calculated results
betas_grid <- read_matrix( paste0( 'betas-', name_in ), ext = 'txt.gz' )
params <- read_tsv( paste0( 'params-', name_in, '.txt.gz' ), show_col_types = FALSE )

# rest always happens in training set
setwd( paste0( '../', train ) )

# load filtered sumstats `df_beta`!
df_beta <- read_tsv( paste0( 'betas-', base, '-clean-matched.txt.gz' ), show_col_types = FALSE )

# since betas came from base, we need to subset them to what we actually have in training data
# subset betas using precalculated map of SNPs from 'base' into 'train'!
betas_grid <- betas_grid[ df_beta[["_NUM_ID_.ss"]], ]

# load training dataset
# Attach the "bigSNP" object in R session
obj.bigSNP <- snp_attach( paste0( name_data, '.rds' ) )
G <- obj.bigSNP$genotypes
y <- obj.bigSNP$fam$affection
y[ y == 0 ] <- NA # in plink format, zeros are missing, translate appropriately here!
# load PCs, to condition on when scoring with R2
PCs <- read_eigenvec( name_data )$eigenvec

# calculate PRS for training individuals only
pred_grid <- big_prodMat( G, betas_grid, ind.col = df_beta[["_NUM_ID_"]] )

# use training individuals to score grid values using correlation, adjusting for PCs, determine which is best
data <- t( apply( pred_grid, 2, function( x ) pcor( x, y, PCs ) ) )
# incorporate back into tibble
colnames( data ) <- c('cor', 'cor_lower', 'cor_upper')
params <- bind_cols( params, as_tibble( data ) )

# add more info, sort by correlation
params <- params %>%
    mutate( sparsity = colMeans( betas_grid == 0 ), id = row_number() ) %>%
    arrange( desc( cor ) )

# save this table of results!
name_out <- paste0( 'eval-', name_in )
write_tsv( params, paste0( name_out, '.txt.gz' ) )

# pick out the best set of parameters, to use and score out of sample later!
betas <- betas_grid[, params$id[1] ]
# save betas!
# this is a simple vector
file_out <- paste0( 'betas-', name_in, '-best.txt.gz' )
write_lines( betas, file_out )

# plot results
fig_start( name_out, width = 6 )
ggplot( params, aes( x = p, y = cor, color = as.factor( h2 ) ) ) +
    theme_classic() +
    geom_point() +
    geom_line() +
    geom_errorbar( aes( ymin = cor_lower, ymax = cor_upper ), width = .1 ) +
    expand_limits( y = 0 ) +
    scale_x_log10() +
    facet_wrap( ~ sparse, labeller = label_both ) +
    labs( x = 'Proportion of causal variants', y = expression(R^2 * " to trait"), color = "Heritability" )
fig_end()
