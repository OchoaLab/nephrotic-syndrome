library(tidyverse)
library(ochoalabtools)

# constants
# names of cases to score in testing data
# also plot order
names_short <- c( 'ct-best', 'ct-stacked', 'inf-best', 'grid-h0.1-best', 'auto-h0.1', 'lassosum-best' )
names_plot <- c( 'CT', 'CT stacked', 'Inf', 'Grid', 'Auto', 'LASSO' )
# hardcode values to process
# only use main test data (Bristol)
test <- 'test'
# consider these combinations of base and training datasets
types_in <- c('base-train', 'base-ssns_ctrl-train-curegn', 'base-ssns_srns-train-curegn')
types_short <- c('SC-DDB', 'SC-DCB', 'SR-DCB')
# in all cases, output name will match inputs in saying it's about correlations
name_out <- 'cor-ALL-ldpred2'

# all processing happens in subdirectory
setwd( test )

# output tibble
data <- NULL

# process preexisting results
for ( name_short in names_short ) {
    name_plot <- names_plot[ names_short == name_short ]
    for ( type_in in types_in ) {
        name_long <- paste0( 'cor-', type_in, '-ldpred2-', name_short )
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
            model = name_plot,
            type = types_short[ type_in == types_in ],
            R2 = cor[1],
            lower = cor[2],
            upper = cor[3]
        )

        # append to tibble
        data <- bind_rows( data, data_name )
    }
}

# manually order models for consistency:
data$model <- factor( data$model, names_plot )
data$type <- factor( data$type, types_short )

# http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
# but mostly copied earlier code from /home/viiia/docs/ochoalab/popdiff/scripts/sim-05-boxplots-all-pi0.R

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge( 0.5 ) # move them .05 to the left and right

fig_start( name_out, width = 6 )
ggplot( data, aes( x = model, y = R2, col = type ) ) + 
    geom_errorbar( aes( ymin = lower, ymax = upper ), width = .5, position = pd ) +
    geom_point( position = pd ) +
    expand_limits( y = 0 ) + 
    theme_classic() +
    labs( x = 'Model', y = "Correlation to trait" )
fig_end()
