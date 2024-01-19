library(bigsnpr)
library(genio)
library(ochoalabtools)
library(tidyverse)

# constants
# previously-chosen h2 estimate around which grid search was performed
h2_est <- 0.1

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
betas_grid <- read_matrix( paste0( 'betas-ldpred2-grid-h', h2_est ), ext = 'txt.gz' )
params <- read_tsv( paste0( 'params-ldpred2-grid-h', h2_est, '.txt.gz' ), show_col_types = FALSE )

# calculate PRS for training individuals only
pred_grid <- big_prodMat(G, betas_grid, ind.row = ind.val, ind.col = df_beta[["_NUM_ID_"]])

# use training individuals to score grid values using correlation, determine which is best
params$cor <- apply( pred_grid, 2, function(x) pcor(x, y[ind.val], NULL)[1] )

# add more info, sort by correlation
params <- params %>%
    mutate( sparsity = colMeans( betas_grid == 0 ), id = row_number() ) %>%
    arrange( desc( cor ) )

# save this table of results!
name_out <- paste0( 'eval-ldpred2-grid-h', h2_est )
write_tsv( params, paste0( name_out, '.txt.gz' ) )

# pick out the best set of parameters, to use and score out of sample later!
betas <- betas_grid[, params$id[1] ]
# save betas!
# this is a simple vector
file_out <- paste0( 'betas-ldpred2-grid-h', h2_est, '-best.txt.gz' )
write_lines( betas, file_out )

# plot results
fig_start( name_out, width = 6 )
ggplot( params, aes( x = p, y = cor, color = as.factor( h2 ) ) ) +
    theme_classic() +
    geom_point() +
    geom_line() +
    scale_x_log10() +
    facet_wrap( ~ sparse, labeller = label_both ) +
    labs( x = 'Proportion of causal variants', y = "Correlation to trait", color = "Heritability" )
fig_end()
