library(tidyverse)
library(ochoalabtools)

# start by splitting discovery data
setwd( '/datacommons/ochoalab/ssns_gwas/imputed/prs-new/' )

# load data
# this is discovery and bristol together
data1 <- read_tsv( '../../array/patient-data.txt.gz', show_col_types = FALSE )
# curegn is also separate!
data2 <- read_tsv( '../../../curegn/patient-data.txt.gz', show_col_types = FALSE )

# only use individuals with ages!  This automatically tosses all controls
data1 <- data1 %>% filter( !is.na( age ) )
data2 <- data2 %>% filter( !is.na( age ) )

# now filter to keep SSNS or SRNS only
data1 <- data1 %>% filter( diagnosis %in% c('SSNS', 'SRNS') )
# for CureGN, the encoding is not the same, but let's subset first, then add the same encodings
# this throws away non-pediatric cases!
data2 <- data2 %>% filter( !is.na( ssns_srns ) )
# recalculate diagnosis (currently MCD vs FSGS, move that to new slot)
data2 <- data2 %>% rename( diagnosis_orig = diagnosis, diagnosis = diagnosis_rasheed )

# plot is clearer if we filter for pediatric cases in main dataset too
data1 <- data1 %>% filter( age <= 21 )

table( data1$diagnosis )
## SRNS SSNS 
##  265  931 
table( data2$diagnosis )
## SRNS SSNS 
##  169  250

# mark dataset and combine for faceted plots
data1$dataset <- 'Duke'
data2$dataset <- 'CureGN'
# subset to few columns of interest
data1 <- data1 %>% select( id, sex, race, diagnosis, age, dataset )
data2 <- data2 %>% select( id, sex, race, diagnosis, age, dataset )
# only CureGN sex is in numeric format, change
data2 <- data2 %>% mutate( sex = ifelse( sex == 1, 'male', ifelse( sex == 2, 'female', 'unknown' ) ) )
# merge!
data <- bind_rows( data1, data2 )
# to match my other figures, rename some variables
data <- data %>% rename( Type = diagnosis, Dataset = dataset )

# make sure Dataset is ordered as desired
data$Dataset <- factor( data$Dataset, levels = c('Duke', 'CureGN') )

# start making a plot of our data
fig_start( 'age-vs-steroid-response', width = 6 )
ggplot( data, aes( x = age, y = after_stat( density ), color = Type ) ) +
    geom_density() +
    theme_classic() +
    labs( y = "Density" ) +
    facet_wrap( vars( Dataset ) )
fig_end()
