library(tidyverse)
library(ochoalabtools)

# constants

# names of cases to score in testing data
# also plot order
names_in <- paste0( 'cor-', c( 'top4', 'base-train-ldpred2-grid-h0.1-best', 'base-train-ldpred2-ct-best' ) )
names_nice <- c( 'Top SNPs', 'LDpred2', 'C+T' )

args <- args_cli()
test <- args[1]
if ( is.na( test ) )
    stop( 'Usage: <test>' )

# all processing happens in subdirectory
setwd( test )

# output tibble
data <- NULL

# process preexisting results
for ( i in 1 : length( names_in ) ) {
    # load PRS calculated previously
    cors <- as.numeric( read_lines( paste0( names_in[i], '.txt.gz' ) ) )

    # reorganize as needed
    data_name <- tibble(
        model = names_nice[i],
        cor = cors[1],
        cor_lower = cors[2],
        cor_upper = cors[3]
    )

    # append to tibble
    data <- bind_rows( data, data_name )
}

# manually order models for consistency:
data$model <- factor( data$model, names_nice )

# before squaring, confirm all correlations are positive
stopifnot( all( data$cor > 0 ) )
stopifnot( all( data$cor_lower > 0 ) )

fig_start( names_in[1], width = 3 )
ggplot( data, aes( x = model, y = cor^2 ) ) + 
    geom_errorbar( aes( ymin = cor_lower^2, ymax = cor_upper^2 ), width = .5 ) +
    geom_point( ) +
    expand_limits( y = 0 ) + 
    theme_classic() +
    labs( x = 'Model', y = expression(R^2 * " to trait") )
fig_end()
