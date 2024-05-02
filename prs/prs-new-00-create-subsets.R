library(tidyverse)

# start by splitting discovery data
setwd( '/datacommons/ochoalab/ssns_gwas/imputed' )

# load data
data_dis <- read_tsv( 'patient-data.txt.gz', show_col_types = FALSE )
# bristol data is separate, load it now too
data_bri <- read_tsv( '../array/patient-data.txt.gz', show_col_types = FALSE )

# we're creating 3 subsets with these characteristics:
# - base: ssns-ctrl, from Discovery, use all discovery controls, generates GWAS data, must be the largest!
# - train: ssns-srns, from Discovery, use all discovery SRNS, and balanced SSNS
# - test: ssns-srns, entirely from Bristol, use all SSNS and SRNS

# so primary challenge is splitting discovery SSNS into two groups

# start with Bristol (easiest), immediately filter to its cases only (exclude unclassified NS)
data_bri <- data_bri %>% filter( bristol, grepl( 'S[SR]NS', diagnosis ) )
# in Bristol, p=0.3 is SRNS, 0.7 is SSNS
# mean( data_bri$diagnosis == 'SRNS' ) # [1] 0.3314176
# this is the list of IDs we want for the Bristol test (same as Bristol otherwise, just explicitly exclude unclassified NS)
ids_test <- data_bri$id

# counts of each diagnosis (SSNS, SRNS, "NS UNCLASSIFIED", and Control), by race
x <- table( data_dis$diagnosis, data_dis$race )
# rowSums( x )
## Control NS UNCLASSIFIED            SRNS            SSNS 
##    3553              14             193             725 

# if train becomes perfectly balanced set, then we want 193 SSNS samples in it, leaving 532 SSNS for base (0.15 ssns/control; original ratio is 0.20 ssns/control for full discovery set)
# numbers don't seem bad to me, let's go with perfectly balanced then!


# NOTE: here we don't know who's SAS vs EAS (among Asian cases), or whether there's an imbalance, let's just hope for the best?

# separate srns and control cases now, those are easy...
ids_base <- data_dis %>% filter( diagnosis == 'Control' ) %>% pull( id )
ids_train <- data_dis %>% filter( diagnosis == 'SRNS' ) %>% pull( id )
## length( ids_base ) # [1] 3553
## length( ids_train ) # [1] 193

# now we want a subset of data_dis that is SSNS only
data_dis_ssns <- data_dis %>% filter( diagnosis == 'SSNS' )

# sample individuals that are controls from each of these ancestries, in these schemes to maintain balance as much as possible
# navigate races
for ( racex in colnames( x ) ) {
    # get IDs for controls of racex
    ids_racex <- data_dis_ssns %>% filter( race == racex ) %>% pull( id )
    # these go to training set, 1:1 balance with SRNS (easy decision)
    ids_racex_train <- sample( ids_racex, x[ 'SRNS', racex ] )
    # rest goes to base set
    ids_racex_base <- setdiff( ids_racex, ids_racex_train )
    # add to respective big sets
    ids_base <- c( ids_base, ids_racex_base )
    ids_train <- c( ids_train, ids_racex_train )
}
## length( ids_base ) # [1] 4085
## length( ids_train ) # [1] 386

# save lists in subfolders
dir.create( 'prs-new' )
setwd( 'prs-new' )
dir.create( 'base' )
dir.create( 'train' )
dir.create( 'test' )
write_lines( ids_base , 'base/ids.txt' )
write_lines( ids_train , 'train/ids.txt' )
write_lines( ids_test , 'test/ids.txt' )
