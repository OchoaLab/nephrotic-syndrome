library(bigsnpr)
library(genio)
library(readr)

# loads SNPs from training set, matches them with test set (because of MAC filters and differing populations, they are not identical)
# NOTE: this is for new pipeline only!

# data to align to
name <- 'test/mac20'
file_betas_shared <- '/betas-clean-matched.txt.gz'
file_betas_in <- paste0( 'train', file_betas_shared ) # data to align
file_betas_out <- paste0( 'test', file_betas_shared )

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

