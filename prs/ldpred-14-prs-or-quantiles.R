library(tidyverse)
library(genio)
library(ochoalabtools)
library(epitools)
source('prs_quant.R')

# constants
name_data <- 'mac20'

args <- args_cli()
test <- args[1]
name <- args[2]
if ( is.na( name ) )
    stop( 'Usage: <test> <name>' )

# main input
file_prs <- paste0( 'prs-', name, '.txt.gz' )

# all processing happens in subdirectory
setwd( test )

# load true phenotype
fam <- read_fam( name_data )
# turn into zeroes and ones
y <- fam$pheno - 1
# load PCs, to visualize PRS stratification
PCs <- read_eigenvec( name_data )$eigenvec
# load PRS per individual
prs <- as.numeric( read_lines( file_prs ) )

# gather into a tibble, with string versions for ease of plotting
data <- tibble(
    Type = ifelse( y == 1, 'SSNS', 'SRNS' ),
    PRS = prs,
    PC1 = PCs[ , 1 ],
    PC2 = PCs[ , 2 ]
)

# look at sort of rawest form, as histograms
fig_start( paste0( 'prs-', name, '-hist' ), width = 4 )
ggplot( data, aes( x = PRS, y= after_stat( density ), fill = Type ) ) +
    geom_histogram( position = 'dodge' ) +
    theme_classic()
fig_end()

# Debo strongly prefers kernely density estimates, make those
fig_start( paste0( 'prs-', name, '-density' ), width = 4 )
ggplot( data, aes( x = PRS, y = after_stat( density ), color = Type ) ) +
    geom_density() +
    theme_classic()
fig_end()

# now overlay true trait and PRS over PCs
# (for trait, name isn't relevant here, it's the same for all methods and sources of base and train data)
fig_start( 'prs-pca-trait', width = 4 )
ggplot( data, aes( x = PC1, y = PC2, color = Type ) ) +
    geom_point( size = 0.5 ) +
    theme_classic()
fig_end()

fig_start( paste0( 'prs-', name, '-pca-prs' ), width = 4 )
ggplot( data, aes( x = PC1, y = PC2, color = PRS ) ) +
    geom_point( size = 0.5 ) +
    theme_classic()
fig_end()

# this function does the magic to calculate ORs for every PRS quartile
data_quant <- prs_quant( data )

# finally, plot!
fig_start( paste0( 'prs-', name, '-or' ), width = 4 )
ggplot( data_quant, aes( x = Quantile, y = OR ) ) + 
    geom_errorbar( aes( ymin = lower, ymax = upper ), width = .5 ) + # , position = pd
    geom_point() + # position = pd
    scale_y_log10() +
    geom_hline( yintercept = 1, linetype = "dashed", color = "gray" ) +
#    expand_limits( y = 0 ) + 
    theme_classic() +
    labs( x = 'PRS Quantiles', y = "OR for SSNS vs SRNS" )
fig_end()
