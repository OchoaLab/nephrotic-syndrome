library(bigsnpr)
library(genio)
library(ochoalabtools)
library(tidyverse)

# start by loading several previous calculations and other needed data

# load training dataset
# Attach the "bigSNP" object in R session
obj.bigSNP <- snp_attach( 'data.rds' )
G <- obj.bigSNP$genotypes
y <- as.numeric( read_lines( 'pheno.txt.gz' ) )
# load indexes of training individuals
ind.val <- as.numeric( read_lines( 'ind-training.txt.gz' ) )
# load filtered sumstats `df_beta`!
df_beta <- read_tsv( 'betas-ssns_ctrl-array.txt.gz', show_col_types = FALSE )
# load previously calculated results
betas_grid <- read_matrix( 'betas-ldpred2-lassosum', ext = 'txt.gz' )
params <- read_tsv( 'params-ldpred2-lassosum.txt.gz', show_col_types = FALSE )

# calculate PRS for training individuals only
pred_grid <- big_prodMat( G, betas_grid, ind.row = ind.val, ind.col = df_beta[["_NUM_ID_"]] )

# use training individuals to score grid values using correlation, determine which is best
params$cor <- apply( pred_grid, 2, function(x) pcor(x, y[ind.val], NULL)[1] )

# add more info, sort by correlation
params <- params %>%
    mutate( id = row_number() ) %>%
    arrange( desc( cor ) )

# save this table of results!
name_out <- 'eval-ldpred2-lassosum'
write_tsv( params, paste0( name_out, '.txt.gz' ) )

# pick out the best set of parameters, to use and score out of sample later!
betas <- betas_grid[, params$id[1] ]
# save betas!
# this is a simple vector
write_lines( betas, 'betas-ldpred2-lassosum-best.txt.gz' )

# plot results
fig_start( name_out, width = 6 )
ggplot( params, aes( x = lambda, y = cor, color = as.factor( delta ) ) ) +
    theme_classic() +
    geom_point() +
    geom_line() +
    scale_x_log10() +
    labs( x = 'LASSO lambda', y = "Correlation to trait", color = "delta" )
fig_end()
