library(bigsnpr)
library(genio)
library(ochoalabtools)
library(readr)

# script takes summary stats that have been cleaned already, and intersects with SNPs in testing or training sets

# constants
# support old data for now
types_old <- c('ssns_ctrl', 'ssns_srns')

# support old data for now, expect ssns_ctrl or ssns_srns
type <- args_cli()[1]
if ( is.na( type ) )
    stop( 'Usage: <type>' )

# load precalculated data
# either way assume script is run from correct local path
if ( type %in% types_old ) {
    name <- 'data'
    file_sumstats_clean <- paste0( 'betas-', type, '-clean.txt.gz' )
    file_out <- paste0( 'betas-', type, '-clean-matched.txt.gz' )
} else {
    # work in desired subdirectory (usually base, train, or test)
    setwd( type )
    name <- 'mac20'
    file_out <- 'betas-clean-matched.txt.gz'
    # this one is always in base only
    file_sumstats_clean <- '../base/mac20-glmm-score-clean.txt.gz'
}

# load SNP set for training dataset
bim <- read_bim( name )
# change names to match snp_match
names( bim )[ names( bim ) == 'id' ] <- 'rsid'
names( bim )[ names( bim ) == 'alt' ] <- 'a1'
names( bim )[ names( bim ) == 'ref' ] <- 'a0'
# remove unneeded column
bim[ names( bim ) == 'posg' ] <- NULL
# convert chr to integer
bim$chr <- as.integer( bim$chr )

message( 'Reading: ', file_sumstats_clean )
sumstats2 <- read_tsv( file_sumstats_clean, show_col_types = FALSE )

message( 'Intersecting with SNPs in ', type )
# find matching subset
df_beta <- snp_match( sumstats2, bim )

# save df_beta for later use
message( 'Writing: ', file_out )
write_tsv( df_beta, file_out )

