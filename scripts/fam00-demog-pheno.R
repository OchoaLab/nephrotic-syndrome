# gathers data from two excel files, harmonizes and reorders/subsets to match actual genetic data available

library(readxl)
library(genio)
library(readr)

# go where the data is
setwd( '~/dbs/ssns_wgs_family' )

# read the more standard formats first
# genetic data with those IDs
fam <- read_fam( 'data' )
# other IDs needed for mapping
fam$id2 <- read_lines( 'tmp.txt' )

# read the XLSX files
ssns <- read_excel( 'SSNS Family Study Manifest Final.xlsx' )
srns <- read_excel( 'SRNS FSGS Manifest Final.xlsx' )

# sizes
## nrow( srns ) # 218
## nrow( ssns ) # 102
## nrow( fam )  # 102

# simplify main ID
names( ssns )[ names( ssns ) == 'Sample Name' ] <- 'id'
names( srns )[ names( srns ) == 'Sample Name' ] <- 'id'

# confirm that files don't overlap
stopifnot( !any( ssns$id %in% srns$id ) )

# subset each file to samples we actually have
ssns <- ssns[ ssns$id %in% fam$id2, ]
srns <- srns[ srns$id %in% fam$id2, ]

# sizes now
## nrow( ssns ) # 94
## nrow( srns ) # 8

# confirm that no samples are missing
stopifnot( all( fam$id2 %in% c(ssns$id, srns$id) ) )

# diagnosis is a bit different between tables, normalize
srns$Diagnosis[ srns$Diagnosis == 'Unaffected' ] <- 'Unaff'

# Race is also odd for SRNS, only cases are WH (White, or White Hispanic?), will treat as White (can validate later as there aren't that many White total)
srns$Race[ srns$Race == 'WH' ] <- 'WHITE'
# merge this variant
ssns$Race[ ssns$Race == 'CAUCASIAN' ] <- 'WHITE'

# things are normalized now, merge tables!
ns <- rbind( ssns, srns )

# reorder to match fam table
ns <- ns[ match( fam$id2, ns$id ), ]
# check agreement
stopifnot( all( ns$id == fam$id2 ) )

# most of `ns` belongs in fam, let's work that in
fam$fam <- ns$Family
fam$sex <- sex_to_int( ns$Sex )
# phenotype treat both SSNS and SRNS as disease
fam$pheno[ ns$Diagnosis %in% c('SSNS', 'SRNS') ] <- 2
fam$pheno[ ns$Diagnosis == 'Unaff' ] <- 1
# only two cases that stay as -9 (missing, what input FAM already has) are "Unaff?" in `ns`, leave as-is

# only thing that doesn't fit is race and detailed diagnosis (no clear use for other IDs)

# save edited fam file
write_fam( 'data2', fam )
# and cleaned up demographics file
write_tsv( ns, 'info.txt' )
