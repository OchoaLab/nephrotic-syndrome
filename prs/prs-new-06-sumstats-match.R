library(bigsnpr)
library(genio)
library(ochoalabtools)
library(readr)

# script takes summary stats that have been cleaned already, and intersects with SNPs in testing or training sets

# constants
# support old data for now
types_old <- c('ssns_ctrl', 'ssns_srns')

# support old data for now, expect ssns_ctrl or ssns_srns
args <- args_cli()
type_base <- args[1]
type_train <- args[2]
if ( is.na( type_base ) )
    stop( 'Usage: <type>' )

# load precalculated data
# either way assume script is run from correct local path
if ( type_base %in% types_old ) {
    name_train <- 'data'
    file_sumstats_clean <- paste0( 'betas-', type_base, '-clean.txt.gz' )
} else {
    if ( is.na( type_train ) )
        stop( 'Usage: <type_base> <type_train>' )
    # this one is always in base only
    file_sumstats_clean <- paste0( '../', type_base, '/mac20-glmm-score-clean.txt.gz' )
    # work in desired training subdirectory
    setwd( type_train )
    name_train <- 'mac20'
}

# output should always specify base data being aligned
file_out <- paste0( 'betas-', type_base, '-clean-matched.txt.gz' )

# load SNP set for training dataset
bim <- read_bim( name_train )
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

message( 'Intersecting with SNPs in ', type_train )
# find matching subset
df_beta <- snp_match( sumstats2, bim )

# save df_beta for later use
message( 'Writing: ', file_out )
write_tsv( df_beta, file_out )

