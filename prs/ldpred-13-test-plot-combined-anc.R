library(tidyverse)
library(ochoalabtools)

# due to the added complexity of ancestry subanalyses, need to focus on a single method here
name_short <- args_cli()[1]
if ( is.na( name_short ) )
    stop( 'Usage: <name>' )

# constants
# ancestries, in desired order, same codes as they appear on the files
# these are the same for all types! (but CureGN lacks SAS)
ancs <- c('AFR', 'EUR', 'SAS', 'Other')
# names of cases to score in testing data
names_short <- c( 'ct-best', 'ct-stacked', 'inf-best', 'grid-h0.1-best', 'auto-h0.1', 'lassosum-best' )
names_plot <- c( 'CT', 'CT stacked', 'Inf', 'Grid', 'Auto', 'LASSO' )
# since we're only using one, determine which it is now, share for all cases
name_plot <- names_plot[ names_short == name_short ]
# hardcode values to process
# only use main test data (Bristol)
test <- 'test'
# actually no, use two variants of curegn as testing dataset
tests_curegn <- c('test-curegn', 'test-curegn2')
# consider these combinations of base and training datasets
types_in <- c('base-train', 'base-ssns_ctrl-train-curegn', 'base-ssns_srns-train-curegn')
types_short <- c('SC-DDB', 'SC-DCB', 'SR-DCB')
# in all cases, output name will match inputs in saying it's about correlations
name_out <- paste0( 'cor-ALL-ldpred2-ancs-', name_short )

# processing function for all methods, for each type
load_all_ancs <- function( data, type_in, test ) {
    for ( anc in ancs ) {
        name_long <- paste0( 'cor-', type_in, '-ldpred2-', name_short, '-', anc )
        # load PRS calculated previously
        file_in <- paste0( name_long, '.txt.gz' )
        
        # skip silently if input is missing
        if ( !file.exists( file_in ) ) next
        # report what is being processed right now
        message( name_long )

        # load vector of numbers
        cors <- as.numeric( read_lines( file_in ) )

        # reorganize as needed
        data_name <- tibble(
            model = name_plot,
            type = types_short[ type_in == types_in ],
            test = test,
            anc = anc,
            cor = cors[1],
            cor_lower = cors[2],
            cor_upper = cors[3]
        )

        # append to tibble
        data <- bind_rows( data, data_name )
    }
    # return updated data
    return( data )
}

# output tibble
data <- NULL

# start with hacks where curegn is used as test, but we can only properly use it if it isn't also used to train (only first type)
for ( test_curegn in tests_curegn ) {
    setwd( test_curegn )
    data <- load_all_ancs( data, types_in[1], test_curegn )
    setwd( '..' )
}

# all processing happens in subdirectory
setwd( test )

# process preexisting results
for ( type_in in types_in ) {
    data <- load_all_ancs( data, type_in, test )
}

# continue test-curegn hack
for ( test_curegn in tests_curegn ) {
    # so far the data shouldn't overlap because testing dataset is distinguished, but it won't plot correctly as-is because type does overlap with test-bristol's, for now
    indexes <- data$test == test_curegn
    stopifnot( all( data$type[ indexes ] == 'SC-DDB' ) )
    # label this test=C instead, which partly solves our problems (and mark second version with a "2" only)
    type_new <- if ( test_curegn == tests_curegn[1] ) 'SC-DDC' else 'SC-DDC2'
    data$type[ indexes ] <- type_new
    # decide which order to plot the new type.  Do last for now
    types_short <- c( types_short, type_new )
}

# manually order models for consistency:
data$anc <- factor( data$anc, ancs )
data$type <- factor( data$type, types_short )

# http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
# but mostly copied earlier code from /home/viiia/docs/ochoalab/popdiff/scripts/sim-05-boxplots-all-pi0.R

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge( 0.5 ) # move them .05 to the left and right

fig_start( name_out, width = 6 )
ggplot( data, aes( x = type, y = cor, col = anc ) ) + 
    geom_errorbar( aes( ymin = cor_lower, ymax = cor_upper ), width = .5, position = pd ) +
    geom_point( position = pd ) +
    expand_limits( y = 0 ) + 
    theme_classic() +
    labs( x = 'Evaluation', y = expression(R^2 * " to trait") ) +
    # since error bars are way bigger here, mark zero
    geom_hline( yintercept = 0, col = "gray" )
fig_end()
