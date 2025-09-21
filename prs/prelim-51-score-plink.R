library(bigsnpr) # for pcor only
library(readr)
library(genio)
library(ochoalabtools)

# assumes scores have already been calculated with plink, just load them!

# constants
name_data <- 'mac20'

args <- args_cli()
test <- args[1]
if ( is.na( test ) )
    stop( 'Usage: <test>' )

# all processing happens in subdirectory
setwd( test )

# load PCs, to condition on when scoring with R2
obj <- read_eigenvec( name_data )
PCs <- obj$eigenvec
fam <- obj$fam

# load precalculated PRS for test individuals now
# this file inherits the trait too!
data <- read_tsv( 'prs-top4-plink.sscore', show_col_types = FALSE )
# confirm alignment
stopifnot( all( fam$id == data$IID ) )
# copy down values
y <- data$PHENO1
preds <- data$SCORE1_AVG

# calculate and save only correlation coefficient to truth, adjusting for PCs
cor <- pcor( preds, y, PCs )
write_lines( cor, 'cor-top4.txt.gz' )

