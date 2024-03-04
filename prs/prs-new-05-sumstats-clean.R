library(bigsnpr)
library(genio)
library(ochoalabtools)
library(readr)

# script does several things:
# - transforms score statistics into coefficient estimates
# - replaces sample size with harmonic mean estimate (ought to be more appropriate for case/control data)
# - subsets to array SNPs (to reduce total numbers)
# - renames columns

# constants
# each base has only one base dataset:
file_sumstats <- 'mac20-glmm-score.txt.gz'
file_out <- 'mac20-glmm-score-clean.txt.gz'

base <- args_cli()[1]
if ( is.na( base ) )
    stop( 'Usage: <base>' )

# work in desired subdirectory (usually base, train, or test)
setwd( base )

# one tiny hiccup is some of these cases are uncompressed...
if ( !file.exists( file_sumstats ) )
    file_sumstats <- 'mac20-glmm-score.txt'

# get sample size a different way
fam <- read_fam( 'mac20' )
counts <- table( fam$pheno - 1 ) # counts zeros and ones

# load summary statistics (big file because it's from imputed data)
message( 'Reading: ', file_sumstats )
sumstats <- bigreadr::fread2( file_sumstats )

message( 'Updating beta, beta_se, n_eff' )

# apply Debo's transformation
sumstats$beta <- sumstats$SCORE / sumstats$VAR
sumstats$beta_se <- 1 / sqrt( sumstats$VAR )
# exclude "MISSRATE" (there's no missingness here), SCORE, VAR
stopifnot( all( sumstats$MISSRATE == 0 ) )
sumstats <- setNames( sumstats[ -c(7, 9, 10) ], c('rsid', 'chr', 'pos', 'a0', 'a1', 'n_eff', 'af', 'p', 'beta', 'beta_se') )

# calculate better effective sample size, overwrite old value
sumstats$n_eff <- 4 / ( 1 / counts[[ '0' ]] + 1 / counts[[ '1' ]] )

# if using ssns-srns, reverse signs!  (ssns-ctrl was fine though)
if ( grepl( 'ssns_srns', base ) )
    sumstats$beta <- -sumstats$beta

message( 'Intersecting with array SNPs' )

# load array SNP set, to subset (which passed QC already)
bim <- read_bim( '/datacommons/ochoalab/ssns_gwas/array/ssns_tgp_merge_clean' )
message( 'Array has these many variants: ', nrow( bim ) )
# change names to match snp_match
names( bim )[ names( bim ) == 'id' ] <- 'rsid'
names( bim )[ names( bim ) == 'alt' ] <- 'a1'
names( bim )[ names( bim ) == 'ref' ] <- 'a0'
# remove unneeded column
bim[ names( bim ) == 'posg' ] <- NULL
# convert chr to integer
bim$chr <- as.integer( bim$chr )
# find matching subset
sumstats2 <- snp_match( sumstats, bim )

# before next round, correct some funny business with snp_match's output
sumstats2[ names( sumstats2 ) == 'rsid' ] <- NULL # IDs from BIM, don't want
names( sumstats2 )[ names( sumstats2 ) == 'rsid.ss' ] <- 'rsid' # restore original IDs
sumstats2[ names( sumstats2 ) == '_NUM_ID_.ss' ] <- NULL # new column added, also don't want
sumstats2[ names( sumstats2 ) == '_NUM_ID_' ] <- NULL # new column added, also don't want

# save df_beta for later use
message( 'Writing: ', file_out )
write_tsv( sumstats2, file_out )

