library(readr)
library(ggplot2)
library(ochoalabtools)

data <- read_tsv( 'ldpred-data-prelim.txt' )

# manually order models:
data$model <- factor( data$model, c('ldpred2-grid', 'ldpred2-inf', 'lassosum2', 'ldpred2-auto') )

# http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
# but mostly copied earlier code from /home/viiia/docs/ochoalab/popdiff/scripts/sim-05-boxplots-all-pi0.R

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge(0.5) # move them .05 to the left and right

fig_start('ldpred-data-prelim', width = 6)
ggplot( data, aes( x = model, y = R2, col = sumstats ) ) + 
    geom_errorbar( aes( ymin = lower, ymax = upper ), width = .5, position = pd ) +
#    geom_line( position = pd ) +
    geom_point( position = pd ) +
    theme_classic()
fig_end()
