library(ochoalabtools)
library(readr)
library(testthat)

# constants
# the two methods to consider, will make separate files for each
methods <- c('grid-h0.1', 'ct')

# determine which type to run
args <- args_cli()
base <- args[1]
train <- args[2]
test <- args[3]
if ( is.na( test ) )
    stop( 'Usage: <base> <train> <test>' )

# all happens in training set
setwd( test )

for ( method in methods ) {
    # load both versions
    name_in <- paste0( 'prs-', base, '-', train, '-ldpred2-', method, '-best' )
    scores_ldpred <- as.numeric( read_lines( paste0( name_in, '.txt.gz' ) ) )
    data <- read_tsv( paste0( name_in, '-plink.sscore' ) )
    # this confirms that both are perfectly correlated!
    expect_equal( cor( data$SCORE1_AVG, scores_ldpred ), 1 )
    # this more strict tests confirms that one score is just scaled version of the other, and the factor is the number of alleles
    # (here tolerance has to be reduced because data has different significant digits)
    expect_equal( data$SCORE1_AVG * data$ALLELE_CT, scores_ldpred, tolerance = 1e-5 )
}
