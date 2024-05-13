library(tidyverse)
library(genio)
library(ochoalabtools)
library(epitools)
source('prs_quant.R')

# combines test (bristol) and test-curegn for larger sample sizes, since results were similar in the separate datasets
# here we don't use PCs (can't do that without properly merging genotype datasets, which we won't do)

# constants
name_data <- 'mac20'

args <- args_cli()
name <- args[1]
if ( is.na( name ) )
    stop( 'Usage: <name>' )

# main input
file_prs <- paste0( 'prs-', name, '.txt.gz' )

# all processing happens in subdirectory
setwd( 'test-curegn' )

# load true phenotype, turn into zeroes and ones
y <- read_fam( name_data )$pheno - 1
# load PRS per individual
prs <- as.numeric( read_lines( file_prs ) )

# now go back to test, load and merge, and outputs will be there too
setwd( '../test' )
# load as before, but concatenate
y <- c( y, read_fam( name_data )$pheno - 1 )
prs <- c( prs, as.numeric( read_lines( file_prs ) ) )

# gather into a tibble, with string versions for ease of plotting
data <- tibble(
    Type = ifelse( y == 1, 'SSNS', 'SRNS' ),
    PRS = prs
)

# look at sort of rawest form, as histograms
fig_start( paste0( 'prs-ALL-', name, '-hist' ), width = 4 )
ggplot( data, aes( x = PRS, y= after_stat( density ), fill = Type ) ) +
    geom_histogram( position = 'dodge' ) +
    theme_classic()
fig_end()

# Debo strongly prefers kernely density estimates, make those
fig_start( paste0( 'prs-ALL-', name, '-density' ), width = 4 )
ggplot( data, aes( x = PRS, y = after_stat( density ), color = Type ) ) +
    geom_density() +
    theme_classic()
fig_end()

# this function does the magic to calculate ORs for every PRS quartile
data_quant <- prs_quant( data )

# finally, plot!
fig_start( paste0( 'prs-ALL-', name, '-or' ), width = 4 )
ggplot( data_quant, aes( x = Quantile, y = OR ) ) + 
    geom_errorbar( aes( ymin = lower, ymax = upper ), width = .5 ) + # , position = pd
    geom_point() + # position = pd
    scale_y_log10() +
    geom_hline( yintercept = 1, linetype = "dashed", color = "gray" ) +
#    expand_limits( y = 0 ) + 
    theme_classic() +
    labs( x = 'PRS Quantiles', y = "OR for SSNS vs SRNS" )
fig_end()
