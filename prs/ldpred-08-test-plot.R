library(tidyverse)
library(ochoalabtools)

# constants

# names of cases to score in testing data
# also plot order
# OLD: c('ldpred2-grid', 'ldpred2-inf', 'lassosum2', 'ldpred2-auto')
names_short <- c( 'inf-best', 'grid-h0.1-best', 'auto-h0.1', 'lassosum-best' )
# for plot, and combines with names_short to make names_long
name_out <- 'ldpred2'
# support old data for now
types_old <- c('ssns_ctrl', 'ssns_srns')

# determine which type to run
type <- args_cli()[1]
if ( is.na( type ) )
    stop( 'Usage: <type>' )

# handle old and new cases!
if ( type %in% types_old ) {
    name_out <- paste0( type, '-', name_out )
} else {
    # all processing happens in subdirectory
    setwd( 'test' )
    # also prepend to "names" to denote this alternative origin of data, which from here on is only used in outputs!
    if ( type == 'base' )
        name_out <- paste0( type, '-', name_out )
}

# in all cases, output name will match inputs in saying it's about correlations
name_out <- paste0( 'cor-', name_out )

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
        R2 = cor[1],
        lower = cor[2],
        upper = cor[3]
    )

    # append to tibble
    data <- bind_rows( data, data_name )
}

# manually order models for consistency:
data$model <- factor( data$model, names_short )

# http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
# but mostly copied earlier code from /home/viiia/docs/ochoalab/popdiff/scripts/sim-05-boxplots-all-pi0.R

fig_start( name_out, width = 6 )
ggplot( data, aes( x = model, y = R2 ) ) + 
    geom_errorbar( aes( ymin = lower, ymax = upper ), width = .5 ) +
    geom_point( ) +
    expand_limits( y = 0 ) + 
    theme_classic()+
    labs( x = 'Model', y = "Correlation to trait" )
fig_end()
