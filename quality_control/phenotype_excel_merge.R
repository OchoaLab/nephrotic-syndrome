# Merge patient phenotype and demographic data from various spreadsheets

library(readxl)
library(tidyverse)
library(genio)

# raw input files are here
setwd( '/datacommons/ochoalab/ssns_gwas/raw/' )
# read both sheets of first file
data1 <- read_xlsx("2022-02-25_2PhenotypeDatafirstDataSets.xlsx")
data2 <- read_xlsx("2022-02-25_2PhenotypeDatafirstDataSets.xlsx", sheet = 2)
# load other file focused on newer samples
data3 <- read_xlsx("PHENOTYPEDATASECONDSETOFPLATES2021_2022.xlsx")
# Bristol had diagnosis updated by Rasheed, use those updates!
data4 <- read_xlsx('2023-02-08_pheno_Bristol_updated.xlsx')
# load genotyped data FAM table
fam <- read_fam( 'nephrotic.syndrome.gwas.all.b38.n2656.R2' )

# subset and rename columns
data1 <- data1 %>% select(id = `Acquisition Number`, sex = Sex, race = RACE, diagnosis = DIAGNOSIS, age = `AGE DIAGNOSIS`)
data2 <- data2 %>% select(id = `Acquisition Number`, sex = SEX_CDE, race = `Broad Race`, diagnosis = Diagnosis, age = `Age at Diagnosis`)
data3 <- data3 %>% select(id = `Acquisition Number`, sex = SEX_CDE, race = `Broad Race`, diagnosis = Diagnosis, age = `Age at Diagnosis`)
data4 <- data4 %>% select(id = AcquisitionNumber, sex = Sex, race = RACE, diagnosis = DIAGNOSIS, age = AGE)

# fix minor issue that some ages have "YRS" suffix in first table
data1$age <- sub( 'YRS$', '', data1$age )
# first two tables have some exact repeats, let `unique` handle it
data1 <- unique( data1 )
data2 <- unique( data2 )

# now all tables have unique IDs (confirm!)
stopifnot( nrow( data1 ) == length( unique( data1$id ) ) )
stopifnot( nrow( data2 ) == length( unique( data2$id ) ) )
stopifnot( nrow( data3 ) == length( unique( data3$id ) ) )
stopifnot( nrow( data4 ) == length( unique( data4$id ) ) )

# only data4 has two intended duplicates added to the table, and one disagrees in diagnosis, let's fix that here!
# after harmonizing, remove the intended duplicates (they will be added back later, along with duplicates for all other tables)
# perform manual edit that makes duplicates agree perfectly
data4$diagnosis[ data4$id == 'B5319' ] <- 'SRNS'
# confirm agreement in all duplicates now (rest of steps)
ids_dupes <- grep( 'd$', data4$id, value = TRUE )
ids_orig <- sub( 'd$', '', ids_dupes )
# get rows, make sure they are ordered by ID, then remove IDs since they won't match and we already know them anyway
rows1 <- data4 %>% filter( id %in% ids_orig ) %>% arrange( id ) %>% select( -id )
rows2 <- data4 %>% filter( id %in% ids_dupes ) %>% arrange( id ) %>% select( -id )
# confirm perfect agreement now
stopifnot( all( rows1 == rows2 ) )
# delete all duplicates now
data4 <- data4 %>% filter( ! id %in% ids_dupes )

# rest of processing is to identify intersections across tables
# identify the small number of IDs in first two tables
ids <- intersect( data1$id, data2$id )
# we manually found second table has one more age (NA in first), and updated diagnosis (all SSNS to SRNS)
# so remove all overlapping individuals from first copy only
data1 <- data1 %>% filter(! id %in% ids )

# data4 (Bristol update) only overlaps with data3
# all of data4 is in data3
indexes <- match( data4$id, data3$id )
data3b <- data3[ indexes, ] # temporary Bristol subset of data3 ordered to match data4
# confirm that agreement is complete except for diagnosis
rows1 <- data3b %>% select( -diagnosis )
rows2 <- data4 %>% select( -diagnosis )
stopifnot( all( rows1 == rows2 ) )
## # just for my own info, identify cases with different diagnosis
## indexes <- data3b$diagnosis != data4$diagnosis
## sum( indexes ) # 22
## table( data3b$diagnosis[ indexes ], data4$diagnosis[ indexes ] )
## ##                 SRNS SSNS
## ## NS UNCLASSIFIED    3    0
## ## SNS                0    1
## ## SSNS              18    0
## # so among changes, all unclassifieds and all SSNS became SRNS, none became SSNS or unclassified (but one typo was corrected)
# armed with this info, it's safe to replace all old bristol data with the new data, same as for data1 vs data2 above
ids <- intersect( data3$id, data4$id )
data3 <- data3 %>% filter(! id %in% ids )

# now combine tables, overwrite first one, don't use rest anymore
# there is overlap between data3 and the first two, but rows are identical so `unique` handles it
data <- bind_rows( data1, data2, data3, data4 ) %>% unique()
# confirm uniqueness
stopifnot( nrow( data ) == length( unique( data$id ) ) )

# compared to genetic data, it has more than our table because of duplicates!
nrow( data ) # 2636
nrow( fam )  # 2656

# we don't know ahead of time which samples were duplicated, best approach is to virtually duplicate them all, using the same "d" suffix used in the genetic data, then filter to match

# create entry for duplicate
data2 <- data
data2$id <- paste0( data2$id, "d" )
# combine again
data <- bind_rows( data, data2 )
# now intersect to samples present in genetic data
data <- data[ data$id %in% fam$id, ]

nrow( data ) # 2655
nrow( fam )  # 2656

# this is the only individual only present in genetic data 
fam[ !fam$id %in% data$id, ]
##   fam      id       pat   mat     sex pheno
## 1 A1014572 A1014572 0     0         0    -9

# spoiler: later we find A1014572 is a duplicate of both 19960560 and 19960560d (very high concordance), and it ends up getting removed because a different duplicate has much lower missingness, so this weird problem went away later unexpectedly!

# age isn't used as a covariate, but let's attempt to clean it anyway
# a few very common cases should be NA though they aren't properly yet:
data$age[ data$age == 'null' ] <- NA
data$age[ data$age == 'Not in File' ] <- NA
data$age[ data$age == 'Unknown' ] <- NA
# there's only one remaining case that is precise, turn that one to decimal
data$age[ data$age == '4year 2 month' ] <- 4 + 2/12

# get data that isn't already NA directly
age1 <- data$age[ !is.na( data$age ) ]
# try to convert to numeric, failures turn into NA
age2 <- as.numeric( age1 )
age_bad <- age1[ is.na( age2 ) ]
# there are the remaining cases, we have no choice but to turn all to NA
# [1] "?8YRS, ?16" "<32"        "<16"        "<43"        "<70"
# numeric or bust
data$age <- as.numeric( data$age )

# lowercase sex (all were already good except for one Unknown)
data$sex <- tolower( data$sex )

# for some quick stats, add a "bristol" column
data$bristol <- grepl( '^B', data$id )

# save the final result!
setwd( '../array/' )
write_tsv( data, 'patient-data.txt.gz' )

# generate list of Bristol IDs too, for filtering later
x <- data$id[ data$bristol ]
write_lines( paste0( x, ' ', x ), 'ids-bristol.txt' )
