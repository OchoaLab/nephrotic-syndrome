library(genio)
library(readr)
library(dplyr)

# load files assumed to be in local path

# start with phenotype file
# this is not standard type, has non-standard header
phen <- read_table( 'srns_ssns_pheno.txt' )
# want to write back without header
# first transform to our usual setup
colnames( phen ) <- c('fam', 'id', 'pheno')
# now write back (to new file really)
write_phen( 'srns_ssns', phen )

## table( phen$pheno )
##   0   1 
## 327 163 
# SSNS is likely 0 (most common), SRNS 1 (least common)

# now load and combine covariates
# PCs
eigenvec <- read_eigenvec( 'srns_ssns' )
# and demographics
covar <- read_table( 'srns_ssns_covar.txt' )

# all dimensions agree!
stopifnot( nrow( phen ) == nrow( eigenvec$fam ) )
stopifnot( nrow( phen ) == nrow( covar ) )

# rename covar columns
covar <- rename( covar, fam = V1, id = AcquisitionNumber, sex = Sex, race = RACE)
# make sure fam is character (join fails otherwise)
covar$fam <- as.character( covar$fam )

# flatten eigenvec
PCs <- as_tibble( eigenvec$eigenvec )
colnames( PCs ) <- paste0( 'pc', colnames( PCs ) )
PCs <- bind_cols( eigenvec$fam, PCs )
# merge, matching by fam and id
covar <- full_join( covar, PCs, by = c( 'fam', 'id' ) )
# rename columns again to match what PRSice expects
covar <- rename( covar, FID = fam, IID = id )

# write it out!
write_tsv( covar, 'srns_ssns_covar_prsice.txt' )
