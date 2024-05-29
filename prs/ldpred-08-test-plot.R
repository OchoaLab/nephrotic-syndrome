library(tidyverse)
library(ochoalabtools)

# constants

# names of cases to score in testing data
# also plot order
names_short <- c( 'ct-best', 'ct-stacked', 'inf-best', 'grid-h0.1-best', 'auto-h0.1', 'lassosum-best' )
# for plot, and combines with names_short to make names_long
name_out <- 'ldpred2'

args <- args_cli()
base <- args[1]
train <- args[2]
test <- args[3]
if ( is.na( test ) )
    stop( 'Usage: <base> <train> <test>' )

# all processing happens in subdirectory
setwd( test )

# in all cases, output name will match inputs in saying it's about correlations
# also prepend to "names" to denote this alternative origin of data, which from here on is only used in outputs!
name_out <- paste0( 'cor-', base, '-', train, '-', name_out )

# output tibble
data <- NULL

# process preexisting results
for ( name_short in names_short ) {
    name_long <- paste0( name_out, '-', name_short )
    # load PRS calculated previously
    file_in <- paste0( name_long, '.txt.gz' )
    
    # skip silently if input is missing
    if ( !file.exists( file_in ) ) next
    # report what is being processed right now
    message( name_long )

    # load vector of numbers
    cor <- as.numeric( read_lines( file_in ) )

    # reorganize as needed
    data_name <- tibble(
        model = name_short,
        cor = cor[1],
        cor_lower = cor[2],
        cor_upper = cor[3]
    )

    # append to tibble
    data <- bind_rows( data, data_name )
}

# manually order models for consistency:
data$model <- factor( data$model, names_short )

# http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
# but mostly copied earlier code from /home/viiia/docs/ochoalab/popdiff/scripts/sim-05-boxplots-all-pi0.R

fig_start( name_out, width = 6 )
ggplot( data, aes( x = model, y = cor ) ) + 
    geom_errorbar( aes( ymin = cor_lower, ymax = cor_upper ), width = .5 ) +
    geom_point( ) +
    expand_limits( y = 0 ) + 
    theme_classic() +
    labs( x = 'Model', y = expression(R^2 * " to trait") )
fig_end()
