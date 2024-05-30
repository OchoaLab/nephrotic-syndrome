library(tidyverse)
library(ochoalabtools)
library(ggpubr)

# constants
# names of cases to score in testing data
# also plot order
names_short <- c( 'grid-h0.1-best', 'lassosum-best', 'ct-best', 'ct-stacked', 'inf-best', 'auto-h0.1' )
names_plot <- c( 'Grid', 'LASSOsum', 'CT', 'CT stacked', 'Inf', 'Auto' )
# hardcode values to process
# only use main test data (Bristol)
test <- 'test'
# actually no, use curegn as testing too, and combined
test_curegn <- 'test-curegn'
test_combined <- 'test-bristol-curegn'
# consider these combinations of base and training datasets
types_in <- c('base-ssns_srns-train-curegn', 'base-ssns_ctrl-train-curegn', 'base-train')
# first two distinguish base, last one (and test-curegn, test-bristol-curegn additions later) distinguishes test instead, so this is organized a bit weirdly
# last version doesn't actually show these labels, but they're used to separate the data later
panels <- c(
    'train: SSNS-SRNS Discov.; test: SSNS-SRNS Bristol',
    'train: SSNS-SRNS Discov.; test: SSNS-SRNS Bristol',
    'base: SSNS-Ctrl Discov.; train: SSNS-SRNS Discov.',
    'base: SSNS-Ctrl Discov.; train: SSNS-SRNS Discov.',
    'base: SSNS-Ctrl Discov.; train: SSNS-SRNS Discov.'
)
vars <- c(
    'SSNS-SRNS',
    'SSNS-Ctrl',
    'Bristol',
    'CureGN',
    'Bristol + CureGN'
)
# in all cases, output name will match inputs in saying it's about correlations
name_out <- 'cor-ALL-ldpred2'

# processing function for all methods, for each type
load_all_names <- function( data, type_in, test, panel, var ) {
    for ( name_short in names_short ) {
        name_plot <- names_plot[ names_short == name_short ]
        name_long <- paste0( 'cor-', type_in, '-ldpred2-', name_short )
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
            test = test,
            panel = panel,
            var = var,
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

# start with hacks where curegn and combined are used as test, but we can only properly use it if it isn't also used to train (only 3rd type)
# this corresponds to the fourth case in terms of panels/vars
setwd( test_curegn )
data <- load_all_names( data, types_in[3], test_curegn, panels[4], vars[4] )
setwd( '..' )
setwd( test_combined )
data <- load_all_names( data, types_in[3], test_combined, panels[5], vars[5] )
setwd( '..' )

# all processing happens in subdirectory
setwd( test )

# process preexisting results
for ( type_in in types_in ) {
    i <- which( types_in == type_in )
    data <- load_all_names( data, type_in, test, panels[i], vars[i] )
}

# manually order models for consistency:
data$model <- factor( data$model, names_plot )
data$panel <- factor( data$panel, unique( panels ) )
data$var <- factor( data$var, unique( vars ) )

# separate the tibbles for two panels because it makes more sense to show the legends separately
data1 <- data %>% filter( panel == panels[1] )
data2 <- data %>% filter( panel != panels[1] )

# http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
# but mostly copied earlier code from /home/viiia/docs/ochoalab/popdiff/scripts/sim-05-boxplots-all-pi0.R

# The errorbars overlapped, so use position_dodge to move them horizontally
pd <- position_dodge( 0.5 ) # move them .05 to the left and right

fig_start( name_out, width = 6, height = 6 )
p1 <- ggplot( data1, aes( x = model, y = cor, col = var ) ) + 
    geom_errorbar( aes( ymin = cor_lower, ymax = cor_upper ), width = .5, position = pd ) +
    geom_point( position = pd ) +
    expand_limits( y = 0 ) + 
    theme_classic() +
    labs( x = 'Model', y = expression(R^2 * " to trait"), col = 'Base dataset' )
p2 <- ggplot( data2, aes( x = model, y = cor, col = var ) ) + 
    geom_errorbar( aes( ymin = cor_lower, ymax = cor_upper ), width = .5, position = pd ) +
    geom_point( position = pd ) +
    expand_limits( y = 0 ) + 
    theme_classic() +
    labs( x = 'Model', y = expression(R^2 * " to trait"), col = 'Test dataset' )
ggarrange( p1, p2, ncol = 1, labels = c('A', 'B'), align = "v" )
fig_end()
