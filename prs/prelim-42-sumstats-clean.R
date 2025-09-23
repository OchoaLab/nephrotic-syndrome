library(bigsnpr)
library(genio)
library(ochoalabtools)
library(tidyverse)

# script does several things:
# - replaces sample size with harmonic mean estimate (ought to be more appropriate for case/control data)
# - subsets to array SNPs (to reduce total numbers)
# - renames columns

args <- args_cli()
base <- args[1]
domrec <- args[2]
if ( is.na( domrec ) )
    stop( 'Usage: <base> <domrec>' )

# work in desired subdirectory
setwd( base )

# each base has only one base dataset:
file_sumstats <- paste0( 'mac20-', domrec, '_saige.txt.gz' )
file_out <- paste0( 'mac20-', domrec, '_saige_clean.txt.gz' )

# load summary statistics (big file because it's from imputed data)
message( 'Reading: ', file_sumstats )
sumstats <- bigreadr::fread2( file_sumstats )

message( 'Updating n_eff' )

# calculate better effective sample size
sumstats <- sumstats %>% mutate( n_eff = 4 / ( 1 / N_case + 1 / N_ctrl ) )

message( 'Selecting columns' )

# select and rename a smaller set of desired columns
sumstats <- sumstats %>% select(
                             rsid = MarkerID,
                             chr = CHR,
                             pos = POS,
                             a0 = Allele1,
                             a1 = Allele2,
                             af = AF_Allele2,
                             p = p.value,
                             beta = BETA,
                             beta_se = SE,
                             n_eff
                         )

## # if using ssns-srns, reverse signs!  (ssns-ctrl was fine though)
## if ( grepl( 'ssns_srns', base ) )
##     sumstats$beta <- -sumstats$beta

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

