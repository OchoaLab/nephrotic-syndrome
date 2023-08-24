# creates phenotypes and covariates of merged array and TGP data

library(tidyverse)
library(genio)

# this is where our data is
setwd( '/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/' )

# load fam file of merged data (can be pre-imputation, they're the same individuals as post-imputation)
fam <- read_fam( 'ssns_tgp_merge_clean' )

# load patient data, which excludes TGP info
data <- read_tsv( "patient-data.txt.gz", show_col_types = FALSE )

# also load TGP assignments to superpopulations, which correspond to "race" in patient data
tgp <- read_tsv( '/datacommons/ochoalab/tgp/pops-annot.txt', comment = '#', show_col_types = FALSE )

# remove samples only in patient data, these should be Bristol and perhaps other removed individuals
# from now on we're only keeping patients present in main genetic cohort data
indexes <- data$id %in% fam$id
data <- data[ indexes, ]
# confirm that no "bristol" samples remain
stopifnot( !any( data$bristol ) )
# remove this trivial column now
data <- select( data, -bristol )
# as covariates, the smaller categories "Mixed, Other, Unknown" should be merged
data$race[ data$race == 'Mixed' ] <- 'Other'
data$race[ data$race == 'Unknown' ] <- 'Other'

# doesn't edit data, but just make sure we understand it
# confirm that no controls have age (of onset)
stopifnot( all( is.na( data$age[ data$diagnosis == 'Control' ] ) ) )
# there are NAs among non-controls, but the are rare (about 7%)
mean( is.na( data$age[ data$diagnosis != 'Control' ] ) )
# [1] 0.07403433

# identify rows of FAM that are missing from patient data
# the ones present should not have any useful information within fam, safe to exclude
indexes <- !fam$id %in% data$id
fam <- fam[ indexes, ]
# confirmed all left are TGP only (will further confirm when mapping to race)
# only these columns have any useful info for TGP, delete the rest:
fam <- select( fam, fam, id, sex )

# transform sex from numeric to character
fam$sex <- sex_to_char( fam$sex )
# and further transform to format array data follows
fam$sex[ fam$sex == 'M' ] <- 'male'
fam$sex[ fam$sex == 'F' ] <- 'female'
# (there were no other cases in TGP)

# add trivial values
fam$diagnosis <- 'Control'
fam$age <- NA

# only thing left is mapping to race
# first get superpopulation from TGP metadata
indexes <- match( fam$fam, tgp$pop )
fam$race <- tgp$superpop[ indexes ]
# manually rename cases to races as they appear in patient data
# note it folds South and East Asians together!
fam$race[ fam$race == 'AFR' ] <- 'Black'
fam$race[ fam$race == 'EUR' ] <- 'White'
fam$race[ fam$race == 'SAS' ] <- 'Asian'
fam$race[ fam$race == 'EAS' ] <- 'Asian'
fam$race[ fam$race == 'AMR' ] <- 'Hispanic'
# can now throw away original "fam" column
fam <- select( fam, -fam )

# merge tables!
data <- bind_rows( data, fam )

# these are the final distributions for these covariates
table( data$sex )
## female    male unknown 
##   2064    2418       3 
table( data$race )
## Asian    Black Hispanic    Other    White 
##  1671     1356      406       30     1022 
table( data$diagnosis )
## Control NS UNCLASSIFIED            SRNS            SSNS 
##    3553              14             193             725 

# create binarized traits (all cases)
data$ns_ctrl <- ifelse( data$diagnosis == 'Control', 0, 1 )
data$ssns_ctrl <- ifelse( data$diagnosis == 'Control', 0, ifelse( data$diagnosis == 'SSNS', 1, NA ) )
data$srns_ctrl <- ifelse( data$diagnosis == 'Control', 0, ifelse( data$diagnosis == 'SRNS', 1, NA ) )
data$ssns_srns <- ifelse( data$diagnosis == 'SSNS', 0, ifelse( data$diagnosis == 'SRNS', 1, NA ) )

# confirm encoding and get distributions
table( data$diagnosis, data$ns_ctrl )
##                    0    1
## Control         3553    0
## NS UNCLASSIFIED    0   14
## SRNS               0  193
## SSNS               0  725
table( data$diagnosis, data$ssns_ctrl )
##                    0    1
## Control         3553    0
## NS UNCLASSIFIED    0    0
## SRNS               0    0
## SSNS               0  725
table( data$diagnosis, data$srns_ctrl )
##                    0    1
## Control         3553    0
## NS UNCLASSIFIED    0    0
## SRNS               0  193
## SSNS               0    0
table( data$diagnosis, data$ssns_srns )
##                   0   1
## Control           0   0
## NS UNCLASSIFIED   0   0
## SRNS              0 193
## SSNS            725   0

# save!
write_tsv( data, 'patient-data-merged-tgp.txt.gz' )


