library(bigsnpr)
library(genio)
library(readr)
library(ochoalabtools)

# loads SNPs from training set, matches them with test set (because of MAC filters and differing populations, they are not identical)

args <- args_cli()
type_base <- args[1]
type_train <- args[2]
type_test <- args[3]
if ( is.na( type_test ) )
    stop( 'Usage: <base> <train> <test>' )

# data to align to
name_test <- paste0( type_test, '/mac20' )
file_betas_in <- paste0( type_train, '/betas-', type_base, '-clean-matched.txt.gz' )
# output should always specify base and training data being aligned
file_betas_out <- paste0( type_test, '/betas-', type_base, '-', type_train, '-clean-matched.txt.gz' )

# load SNP set for test dataset
bim <- read_bim( name_test )
# change names to match snp_match
names( bim )[ names( bim ) == 'id' ] <- 'rsid'
names( bim )[ names( bim ) == 'alt' ] <- 'a1'
names( bim )[ names( bim ) == 'ref' ] <- 'a0'
# remove unneeded column
bim[ names( bim ) == 'posg' ] <- NULL
# convert chr to integer
bim$chr <- as.integer( bim$chr )

# reload precalculated df_beta
message( 'Reading: ', file_betas_in )
df_beta_in <- read_tsv( file_betas_in, show_col_types = FALSE )

# before next round, correct some funny business with snp_match's output
df_beta_in[ names( df_beta_in ) == 'rsid' ] <- NULL # IDs from BIM, don't want
names( df_beta_in )[ names( df_beta_in ) == 'rsid.ss' ] <- 'rsid' # restore original IDs
df_beta_in[ names( df_beta_in ) == '_NUM_ID_.ss' ] <- NULL # new column added, also don't want
df_beta_in[ names( df_beta_in ) == '_NUM_ID_' ] <- NULL # new column added, also don't want

message( 'Intersecting with SNPs in testing set' )
# find matching subset
df_beta_out <- snp_match( df_beta_in, bim )

# save df_beta for later use
message( 'Writing: ', file_betas_out )
write_tsv( df_beta_out, file_betas_out )

