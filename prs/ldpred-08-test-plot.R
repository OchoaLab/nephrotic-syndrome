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

# support old data for now, expect ssns_ctrl or ssns_srns
args <- args_cli()
type_base <- args[1]
type_train <- args[2]
type_test <- args[3]
if ( is.na( type_base ) )
    stop( 'Usage: <type>' )

# handle old and new cases!
if ( type_base %in% types_old ) {
    # base is only parameter here
    type_in <- type_base
} else {
    if ( is.na( type_test ) )
        stop( 'Usage: <type_base> <type_train> <type_test>' )
    # all processing happens in subdirectory
    setwd( type_test )
    # combine base and train in new setup
    type_in <- paste0( type_base, '-', type_train )
}

# in all cases, output name will match inputs in saying it's about correlations
# also prepend to "names" to denote this alternative origin of data, which from here on is only used in outputs!
name_out <- paste0( 'cor-', type_in, '-', name_out )

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
