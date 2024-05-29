library(tidyverse)
library(ochoalabtools)

# constants
# dataset print names, in order
datasets <- c(
    'base: SSNS-SRNS Discov.; train: SSNS-SRNS CureGN',
    'base: SSNS-Ctrl Discov.; train: SSNS-SRNS CureGN',
    'base: SSNS-Ctrl Discov.; train: SSNS-SRNS Discov.'
)
# previously-chosen h2 estimate around which grid search was performed
h2_est <- 0.1

# load precalculated data for all cases

# first the older scenarios where curegn was used to train data
setwd( 'train-curegn' )

# save this table of results!
params <- read_tsv( paste0( 'eval-base-ssns_srns-ldpred2-grid-h', h2_est, '.txt.gz' ), show_col_types = FALSE )
params$Dataset <- datasets[1]

params2 <- read_tsv( paste0( 'eval-base-ssns_ctrl-ldpred2-grid-h', h2_est, '.txt.gz' ), show_col_types = FALSE )
params2$Dataset <- datasets[2]
params <- bind_rows( params, params2 )

# now load final scenario
# this will also be figure's destination
setwd( '../train' )

params2 <- read_tsv( paste0( 'eval-base-ldpred2-grid-h', h2_est, '.txt.gz' ), show_col_types = FALSE )
params2$Dataset <- datasets[3]
params <- bind_rows( params, params2 )

# make sure Dataset is ordered as desired
params$Dataset <- factor( params$Dataset, levels = datasets )

# make sparsity more verbose for strip text to be more human-readable
params <- params %>% mutate( sparse = ifelse( sparse, 'Sparse', 'Dense' ) )

# plot results
wh <- fig_scale( 3/2 )
fig_start( paste0( 'eval-ALL-ldpred2-grid-h', h2_est ), width = wh[1], height = wh[2] )
ggplot( params, aes( x = p, y = cor, color = as.factor( h2 ) ) ) +
    theme_classic() +
    geom_point() +
    geom_line() +
    geom_errorbar( aes( ymin = cor_lower, ymax = cor_upper ), width = .03 ) +
    expand_limits( y = 0 ) +
    scale_x_log10() +
    facet_grid( sparse ~ Dataset ) +
    theme( strip.text = element_text( size = 5 ) ) +
    labs( x = 'Proportion of causal variants', y = expression(R^2 * " to trait"), color = "Heritability" )
fig_end()
