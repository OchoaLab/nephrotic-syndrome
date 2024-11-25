library(tidyverse)
library(ggpubr)
library(genio)
library(ochoalabtools)
library(epitools)
library(PRROC)
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
# remember which dataset this came from, organize
data <- tibble(
    y = y,
    PRS = prs,
    Dataset = 'CureGN'
)

# now go back to test, load and merge, and outputs will be there too
setwd( '../test' )
# load as before
y <- read_fam( name_data )$pheno - 1
prs <- as.numeric( read_lines( file_prs ) )
data2 <- tibble(
    y = y,
    PRS = prs,
    Dataset = 'Bristol'
)
# concatenate, putting Bristol first
data <- bind_rows( data2, data )
# make trait a string, for ease of plotting
data <- data %>% mutate( Type = ifelse( y == 1, 'SSNS', 'SRNS' ) )

# create merged version by duplicating data, but relabeling dataset
data2 <- data
data2$Dataset <- 'Bristol+CureGN'
# concatenate again, putting merged last
data <- bind_rows( data, data2 )

# make sure Dataset is ordered as desired
data$Dataset <- factor( data$Dataset, levels = c('Bristol', 'CureGN', 'Bristol+CureGN') )

### ORs ###

# this function does the magic to calculate ORs for every PRS quartile
# apply to each dataset, then re-concatenate
data_quant <- NULL
datasets <- levels( data$Dataset )
for ( dataset in datasets ) {
    data_quant_i <- prs_quant( data %>% filter( Dataset == dataset ) )
    data_quant_i$Dataset <- dataset
    data_quant <- bind_rows( data_quant, data_quant_i )
}
data_quant$Dataset <- factor( data_quant$Dataset, levels = datasets )
# NOTE: levels for quantiles are a bit awkward because they're all different, but they work for this plot!
write_tsv( data_quant, paste0( 'or-quarts-ALL-panels-', name, '.txt.gz' ) )

### PLOT ###

# for plot, makes a text version combining a lot of this data
data_quant <- data_quant %>% mutate( lab = paste0( Quantile, "\nOR ", round( OR, 2 ), "\n[", round( lower, 2 ), '-', round( upper, 2 ), ']' ) )

# finally, combined plot!
wh <- fig_scale( 4/3 )
fig_start( paste0( 'prs-ALL-panels-', name ), width = wh[1], height = wh[2] )
p1 <- ggplot( data, aes( x = PRS, y = after_stat( density ), color = Type ) ) +
    geom_density() +
    theme_classic() +
    labs( y = "Density" ) +
    facet_wrap( vars( Dataset ) )
p2 <- ggplot( data_quant, aes( x = lab, y = OR ) ) + 
    geom_errorbar( aes( ymin = lower, ymax = upper ), width = .5 ) +
    geom_point() +
    scale_y_log10() +
    geom_hline( yintercept = 1, linetype = "dashed", color = "gray" ) +
    theme_classic() +
    labs( x = 'PRS Quantiles', y = "OR for SSNS vs SRNS" ) +
    facet_wrap( vars( Dataset ), scales = "free_x" ) +
    theme( axis.text.x = element_text( size = rel( 0.9 ) ) )
ggarrange( p1, p2, ncol = 1, align = "v", labels = c('A', 'B'), heights = c(1, 1.1) )
fig_end()

### ROC ###

# calculate and gather data
auc_roc <- vector( 'numeric', length( datasets ) )
names( auc_roc ) <- datasets
data_roc <- NULL
for ( dataset in datasets ) {
    data_i <- data %>% filter( Dataset == dataset )
    # positive, negative classes
    # or scores.class0 are all scores, weights.class0 are class labels
    roc <- roc.curve( data_i$PRS, weights.class0 = data_i$y, curve = TRUE )
    # store AUC value to report later
    auc_roc[ datasets == dataset ] <- roc$auc
    # store curve info in a tibble
    data_roc_i <- roc$curve
    colnames( data_roc_i ) <- c( 'FPR', 'TPR', 'PRS' )
    data_roc_i <- as_tibble( data_roc_i )
    data_roc_i$Dataset <- dataset
    data_roc <- bind_rows( data_roc, data_roc_i )
}
data_roc$Dataset <- factor( data_roc$Dataset, levels = datasets )

# make a nice ggplot
wh <- fig_scale( 3 )
fig_start( paste0( 'prs-roc-ALL-panels-', name ), width = wh[1], height = wh[2] )
ggplot( data_roc, aes( x = FPR, y = TPR ) ) +
    geom_abline( slope = 1, intercept = 0, linetype = 'dashed', color = 'gray' ) +
    geom_line() +
    theme_classic() +
    facet_wrap( ~Dataset )
fig_end()

# just report AUCs to store elsewhere
print( auc_roc )
