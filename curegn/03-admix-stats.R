library(tidyverse)

# location of required admixture analysis
setwd( '/datacommons/ochoalab/curegn/merge_tgp/admixture' )

# this is the table of interest!
data <- read_tsv( 'curegn_covar_diseasesubtype.txt' )

# read main annotations to figure out one discrepancy
data2 <- read_tsv( '../../patient-data.txt.gz' )

# subset to diagnosis of interest
data <- data %>% filter( !is.na( ssns ) | !is.na( srns ) )
data2 <- data2 %>% filter( !is.na( ssns_srns ) )

nrow( data ) # [1] 420
nrow( data2 ) # [1] 419

# the smaller/newer set is a proper subset of the bigger/older one
stopifnot( all( data2$id %in% data$id ) )

# I think the problem is the presence of an individual that dropped consent
# let's just toss them in the other dataset
data <- data %>% filter( id %in% data2$id )

# now they agree!
stopifnot( nrow( data ) == nrow( data2 ) )

# raw counts
table( data$ancestry )
##   AFR   AFR_admix Asian_admix         EAS         EUR   EUR_admix 
##    64          34           2           9         218          38 
## Other         SAS 
##    44          10 
table( data$race )
## Asian       Black Multiracial      NatAmr     Pacific     Unknown 
##    21          98          21           6           2          15 
## White 
##   256 

# grouping ancestry as we report it
data$ancestry[ ! data$ancestry %in% c('AFR', 'EUR', 'EAS', 'SAS') ] <- 'Other'
x <- table( data$ancestry )
x
## AFR   EAS   EUR Other   SAS 
##  64     9   218   118    10 
round( x/sum( x ) * 100, 1 )
##  AFR   EAS   EUR Other   SAS 
## 15.3   2.1  52.0  28.2   2.4 
