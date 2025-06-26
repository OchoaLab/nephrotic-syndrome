library(ochoalabtools)
library(tidyverse)

# constants
# the two methods to consider, will make separate files for each
methods <- c('grid-h0.1', 'ct')

# determine which type to run
args <- args_cli()
base <- args[1]
train <- args[2]
if ( is.na( train ) )
    stop( 'Usage: <base> <train>' )

# all happens in training set
setwd( train )

# load filtered sumstats `df_beta`!
df_beta <- read_tsv( paste0( 'betas-', base, '-clean-matched.txt.gz' ), show_col_types = FALSE )
# extract columns of interest
# let's make it easy to use with plink2 --score: https://www.cog-genomics.org/plink/2.0/score
df_beta <- df_beta %>% select( id = rsid, a1 ) # chr, pos, a0, af, p, beta, 

for ( method in methods ) {
    # make a copy to edit further
    data <- df_beta
    
    # load and add best betas (simple vector) as new column!
    data$prs <- as.numeric( read_lines( paste0( 'betas-', base, '-ldpred2-', method, '-best.txt.gz' ) ) )

    # NOTES:
    # 1. There are betas for the full grid on the base data, but the best one was determined using the train data used here.  Makes most sense to restrict to the SNPs in the train data (in a sense, the other SNPs weren't evaluated)
    # 2. Train's BIM has many more SNPs (9,549,985) than the betas vector (528,964).  This is because we only considered array SNPs for PRS I think.  Anyway, df_beta is matched, use info from there instead.

    # keep only non-zero coefficients, which saves a lot of space and later computational runtime for sparse models such as CT
    # (for ldpred-grid nothing is removed, all coefficients are non-zero)
    data <- data %>% filter( prs != 0 )
    
    # save to output
    write_tsv( data, paste0( 'betas-', base, '-ldpred2-', method, '-best-plink-score.txt' ) )
}
