library(tidyverse)
library(genio)
library(janitor) # optional, for `remove_constant`

# run from here on DCC
setwd( '/datacommons/ochoalab/curegn' )

# read original fam file of actual genetic data (already filtered), to try to match to it
# only non-trivial column is `id`, rest are blank
fam <- read_fam( 'curegn-autosomes-snps' )
ids <- fam$id

# read big covariates table now
data <- read_csv( 'raw/Gbadegesin_CureGN_PATIENTS_08AUG2023.csv' )

### ID ###

# 3 individuals from `fam` are missing from `data`
# confirmed with manual text searches that these IDs are nowhere in `data` (not in other columns)
indexes <- ids %in% data$WGS_ID
write_lines( ids[ !indexes ], 'indiv-rm-round2.txt' )
# [1] "KO0353" "KO0634" "KO2033"

# filter ids to prepare the rest of the fam file
ids <- ids[ indexes ]
data <- data %>% filter( WGS_ID %in% ids )
# order `data` to match remaining `ids`
data <- data %>% arrange( factor( WGS_ID, levels = ids ) )
# confirm
stopifnot( all( data$WGS_ID == ids ) )
# confirm IDs are unique
stopifnot( nrow( data ) == length( unique( data$WGS_ID ) ) )
# these study IDs are unique too (no repeated individuals?)
stopifnot( nrow( data ) == length( unique( data$CureGNSTudyID ) ) )

# at this point it's worth removing constant columns to reduce visual clutter while exploring data
data <- remove_constant( data, quiet = FALSE )
# Removing 9 constant columns of 116 columns total (Removed: EligibleCalc, FirstBiopsyRptsSlidesAvail, ConsentCalc, RaceMiss, LevamisoleEver, GutSteroidEver, WGS_flag, FirstGutSteroidDaysEnroll, FirstLevamisoleDaysEnroll).

### SEX ###

# only one column encodes sex (based on regexp)
# only two standard values
table( data$PAT_Sex, useNA = 'i' )
## 1:Male 2:Female 
##   1047      805 
# simplify to just the standard numerical encoding
data$PAT_Sex <- sapply( strsplit( data$PAT_Sex, ':' ), function(x) x[1] )

### DIAGNOSIS ###

## table( data$LocalDiagnosis, useNA = 'i' )
## 1: MCD (includes IgM Nephropathy)                           2: FSGS 
##                               402                               485 
##                             3: MN            4: IgAN (includes HSP) 
##                               442                               507 
##                            5: C1Q                          6: Other 
##                                14                                 2 
## table( data$DiagnosisPath, useNA = 'i' )
## 1: MCD (includes IgM Nephropathy)                           2: FSGS 
##                               415                               478 
##                             3: MN            4: IgAN (includes HSP) 
##                               442                               517 

## table( data$LocalDiagnosis, data$DiagnosisPath, useNA = 'i' )
##                                   1: MCD (includes IgM Nephropathy) 2: FSGS 3: MN 4: IgAN (includes HSP)
## 1: MCD (includes IgM Nephropathy)                               379      12     3                      8
## 2: FSGS                                                          22     455     2                      6
## 3: MN                                                             4       0   437                      1
## 4: IgAN (includes HSP)                                            4       1     0                    502
## 5: C1Q                                                            5       9     0                      0
## 6: Other                                                          1       1     0                      0

# data dictionary isn't super useful here:
# LocalDiagnosis: B B2:Local Diagnosis:
# DiagnosisPath: Diagnosis Cohort: assigned using hierarchy of (1) Core Scoring CRF final disease (F17), Pathology CRF cohort (H3), Screening Form local diagnosis (B2)
# ("path" may be short for pathology)

# my impression is DiagnosisPath is the good one, let's clean it a bit
# first remove numbers and colons
data$DiagnosisPath <- sapply( strsplit( data$DiagnosisPath, ': ' ), function(x) x[2] )
# shorten long names, to remove spaces
data$DiagnosisPath[ data$DiagnosisPath == 'IgAN (includes HSP)' ] <- 'IgAN'
data$DiagnosisPath[ data$DiagnosisPath == 'MCD (includes IgM Nephropathy)' ] <- 'MCD'

# encode NS (MCD and FSGS) vs non-NS (rest)
data$ns <- grepl( 'FSGS|MCD', data$DiagnosisPath )
# confirm it worked:
## table( data$DiagnosisPath, data$ns ) 
##                                FALSE TRUE
## FSGS                               0  478
## IgAN (includes HSP)              517    0
## MCD (includes IgM Nephropathy)     0  415
## MN                               442    0

# final counts
table( data$ns, useNA = 'i' )
## FALSE  TRUE 
##   959   893

### AGE ###

## columns_age <- grep( 'age', names( data ), ignore.case = TRUE, value=TRUE )
## ## [1] "AgeAtBiopYr"     "AgeAtBiopMo"     "AgeAtScreenYr"   "AgeAtScreenMo"  
## ## [5] "AgeAtConsentYr"  "AgeAtConsentMo"  "AgeAtEnrollment"

## # by eye it appears that the first one has the lowest values, confirm...
## # it's not absolutely the smallest, but nearly always!
## mean( data$AgeAtBiopYr <= data$AgeAtScreenYr, na.rm = TRUE )   # 0.9989101
## mean( data$AgeAtScreenYr <= data$AgeAtBiopYr, na.rm = TRUE )   # 0.3814714 # reversed, confirms better to use first one
## mean( data$AgeAtBiopYr <= data$AgeAtConsentYr, na.rm = TRUE )  # 1
## mean( data$AgeAtBiopYr <= data$AgeAtEnrollment, na.rm = TRUE ) # 1

### RACE ###

# NOTE: there's no ethnicity column in this data (maybe we didn't ask for it?), and race column doesn't include Hispanics

# this is the one we want, start cleaning up
# first remove numbers and colons
data$RaceCat <- sapply( strsplit( data$RaceCat, ': ' ), function(x) x[2] )
# shorten long codes, get rid of spaces
data$RaceCat[ data$RaceCat == 'Black/African American' ] <- 'Black'
data$RaceCat[ data$RaceCat == 'White/caucasian' ] <- 'White'
data$RaceCat[ data$RaceCat == 'Pacific Islander' ] <- 'Pacific'
data$RaceCat[ data$RaceCat == 'Native American' ] <- 'NatAmr'

# these binary versions clarify mixed ancestries (somewhat)!
# only RaceUnknown maps perfectly to RaceCat==Unknown (not shown)
## table( data$RaceCat, data$RaceAsian )
##                  0: No 1: Yes
## Asian                0    154
## Black              278      0
## Multiracial         40      8
## Native American     11      0
## Pacific Islander     7      0
## Unknown             59      0
## White             1295      0
## table( data$RaceCat, data$RaceBlack )
##                    0: No 1: Yes
##   Asian              154      0
##   Black                0    278
##   Multiracial         23     25
##   Native American     11      0
##   Pacific Islander     7      0
##   Unknown             59      0
##   White             1295      0
## table( data$RaceCat, data$RaceNativAm )
##                    0: No 1: Yes
##   Asian              154      0
##   Black              278      0
##   Multiracial         26     22
##   Native American      0     11
##   Pacific Islander     7      0
##   Unknown             59      0
##   White             1295      0
## table( data$RaceCat, data$RacePacific )
##                    0: No 1: Yes
##   Asian              154      0
##   Black              278      0
##   Multiracial         46      2
##   Native American     11      0
##   Pacific Islander     0      7
##   Unknown             59      0
##   White             1295      0
## table( data$RaceCat, data$RaceWhite )
##                    0: No 1: Yes
##   Asian              154      0
##   Black              278      0
##   Multiracial          7     41
##   Native American     11      0
##   Pacific Islander     7      0
##   Unknown             59      0
##   White                0   1295

# select colums of interest now only
data2 <- select( data, id = WGS_ID, sex = PAT_Sex, race = RaceCat, age = AgeAtBiopYr, diagnosis = DiagnosisPath, ns )

# write covariates table
write_tsv( data2, 'patient-data.txt.gz' )
