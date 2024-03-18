library(readxl)
library(readr)

# read the XLSX files
# NOTE: file has two sheets, this only loads first one but that's the only one we want
data <- read_excel( 'CureGN_RNASeq_WGS_ID_Mapping_20240318.xlsx' )
nrow( data ) # 2341

# exclude withdrawn individuals
# that column only has two cases: withdrawn or NA (keep)
stopifnot( all( data$Withdrawal[ !is.na( data$Withdrawal ) ] == 'withdrawn' ) )
# this keeps individuals with NA values
data <- data[ is.na( data$Withdrawal ), ]
nrow( data ) # 2296

# many remaining rows don't have WGS IDs, exclude those too
# here non-NAs are kept
data <- data[ !is.na( data$WGS_ID ), ]
nrow( data ) # 2012 # removed close to 10% of data!

# last exclusion: two samples that are potentially swapped
data <- data[ -grep( 'excluded', data$Note ), ]
nrow( data ) # 2010

# write list of individual IDs to keep
# sorted seems nicer, they appear sorted in the VCF files too
write_lines( sort( data$WGS_ID ), 'indiv-keep.txt' )
