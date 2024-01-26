library(bigsnpr)
library(genio)
library(ochoalabtools)
library(tidyverse)

# constants
# sequence of heritabilities to consider
herits <- c( (1:9)/100, (1:9)/10 )

# determine which type to run
type <- args_cli()[1]
if ( is.na( type ) )
    stop( 'Usage: <type: ssns_ctrl or ssns_srns>' )

# start by loading several previous calculations and other needed data

# suffix shared by several in/out files
name_in <- paste0( type, '-ldpred2-inf' )

# load training dataset
# Attach the "bigSNP" object in R session
obj.bigSNP <- snp_attach( 'data.rds' )
G <- obj.bigSNP$genotypes
y <- obj.bigSNP$fam$affection
y[ y == 0 ] <- NA # in plink format, zeros are missing, translate appropriately here!
# load indexes of training individuals
ind_train <- as.numeric( read_lines( 'ind-train.txt.gz' ) )
# load filtered sumstats `df_beta`!
df_beta <- read_tsv( paste0( 'betas-', type, '-clean-matched.txt.gz' ), show_col_types = FALSE )
# load previously calculated results
betas_grid <- read_matrix( paste0( 'betas-', name_in ), ext = 'txt.gz' )

# calculate PRS for training individuals only
pred_grid <- big_prodMat( G, betas_grid, ind.row = ind_train, ind.col = df_beta[["_NUM_ID_"]] )

# create a trivial params table
params <- tibble( h2 = herits )
# use training individuals to score grid values using correlation, determine which is best
params$cor <- apply( pred_grid, 2, function(x) pcor(x, y[ ind_train ], NULL)[1] )

# add more info, sort by correlation
params <- params %>%
    mutate( id = row_number() ) %>%
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
ggplot( params, aes( x = h2, y = cor ) ) +
    theme_classic() +
    geom_point() +
    geom_line() +
    labs( x = 'Heritability', y = "Correlation to trait" )
fig_end()
