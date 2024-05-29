library(tidyverse)
library(ochoalabtools)

# constants
# dataset print names, in order
datasets <- c(
    'base: SSNS-SRNS Discov.; train: SSNS-SRNS CureGN',
    'base: SSNS-Ctrl Discov.; train: SSNS-SRNS CureGN',
    'base: SSNS-Ctrl Discov.; train: SSNS-SRNS Discov.'
)

# load precalculated data for all cases

# first the older scenarios where curegn was used to train data
setwd( 'train-curegn' )

# save this table of results!
params <- read_tsv( 'eval-base-ssns_srns-ldpred2-inf.txt.gz', show_col_types = FALSE )
params$Dataset <- datasets[1]

params2 <- read_tsv( 'eval-base-ssns_ctrl-ldpred2-inf.txt.gz', show_col_types = FALSE )
params2$Dataset <- datasets[2]
params <- bind_rows( params, params2 )

# now load final scenario
# this will also be figure's destination
setwd( '../train' )

params2 <- read_tsv( 'eval-base-ldpred2-inf.txt.gz', show_col_types = FALSE )
params2$Dataset <- datasets[3]
params <- bind_rows( params, params2 )

# make sure Dataset is ordered as desired
params$Dataset <- factor( params$Dataset, levels = datasets )

# plot results
wh <- fig_scale( 3 )
fig_start( 'eval-ALL-ldpred2-inf', width = wh[1], height = wh[2] )
ggplot( params, aes( x = h2, y = cor ) ) +
    theme_classic() +
    geom_point() +
    geom_line() +
    geom_errorbar( aes( ymin = cor_lower, ymax = cor_upper ), width = .03 ) +
    expand_limits( y = 0 ) + 
    labs( x = 'Heritability', y = expression(R^2 * " to trait") ) +
    facet_wrap( vars( Dataset ) ) +
    theme( strip.text = element_text( size = 6 ) )
fig_end()
