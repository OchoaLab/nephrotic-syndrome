# a basic covariate analysis motivating aim3 (probabilistic prediction of SSNS vs SRNS using age, sex, race, etc)

library(tidyverse)
library(genio)

# go where the main data is
setwd( '/datacommons/ochoalab/ssns_gwas/array' )

# this one excludes tgp (fine) but includes both discovery and bristol!
data <- read_tsv( 'patient-data.txt.gz', show_col_types = FALSE )
# CureGN is only separate one
data2 <- read_tsv( '../../curegn/patient-data.txt.gz', show_col_types = FALSE )

# start harmonizing

### for discovery/bristol

# add a column denoting dataset
# (the "bristol" boolean no longer needed after that)
data <- data %>% mutate( dataset = ifelse( bristol, 'Bristol', 'Discovery' ) ) %>% select( -bristol )
# toss controls and unclassifieds
data <- data %>% filter( diagnosis %in% c('SSNS', 'SRNS') )
# copy diagnosis to diagnosis2 (for coherence with curegn)
data$diagnosis2 <- data$diagnosis
# harmonize race values, grouping smaller sets in particular
data$race[ data$race == 'Mixed' ] <- 'Other'
data$race[ data$race == 'Unknown' ] <- 'Other'

### for curegn

# add dataset
data2$dataset <- 'CureGN'
# keep NS cases only (don't need that column afterwards, toss another redundant one)
data2 <- data2 %>% filter( ns ) %>% select( -ns, -ssns_srns )
# treat MCD and FSGS as SSNS, SRNS in one version
data2 <- data2 %>% mutate( diagnosis = ifelse( diagnosis == 'MCD', 'SSNS', 'SRNS' ) )
# rename diagnoses columns (*2 is less restricted, treat as secondary)
data2 <- data2 %>% rename( diagnosis2 = diagnosis, diagnosis = diagnosis_rasheed )
# reencode sex to match data
data2 <- data2 %>% mutate( sex = sex_to_char( sex ) )
stopifnot( !any( data2$sex == 'U' ) ) # next step assumes none were unknown (in curegn)
data2 <- data2 %>% mutate( sex = ifelse( sex == 'M', 'male', 'female' ) )
# harmonize race values, grouping smaller sets in particular
data2$race[ data2$race == 'Multiracial' ] <- 'Other'
data2$race[ data2$race == 'NatAmr' ] <- 'Other'
data2$race[ data2$race == 'Pacific' ] <- 'Other'
data2$race[ data2$race == 'Unknown' ] <- 'Other'

# merge!
data <- bind_rows( data, data2 )

# not sure where to store this, put at base
setwd( '..' )
write_tsv( data, 'covars-ssns-srns-discov-bristol-curegn.txt.gz' )
