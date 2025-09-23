library(bigsnpr)
library(genio)
library(ochoalabtools)
library(readr)

# script takes summary stats that have been cleaned already, and intersects with SNPs in testing or training sets

# constants
name_train <- 'mac20'

args <- args_cli()
base <- args[1]
train <- args[2]
domrec <- args[3]
if ( is.na( domrec ) )
    stop( 'Usage: <base> <train> <domrec>' )

# load precalculated data
# this one is always in base only
file_sumstats_clean <- paste0( '../', base, '/mac20-', domrec, '_saige_clean.txt.gz' )

# work in desired training subdirectory
setwd( train )

# output should always specify base data being aligned
file_out <- paste0( 'betas-', domrec, '-', base, '-clean-matched.txt.gz' )

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

message( 'Intersecting with SNPs in ', train )
# find matching subset
df_beta <- snp_match( sumstats2, bim )

# save df_beta for later use
message( 'Writing: ', file_out )
write_tsv( df_beta, file_out )

