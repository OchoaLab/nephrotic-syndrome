# get stats for the paper

library(readr)
data <- read_tsv( 'patient-data.txt.gz', show_col_types = FALSE )

# immediately subset to cases we classified as SSNS or SRNS, which are automatically pediatric (this was done in a previous script)
data <- data[ !is.na( data$diagnosis_rasheed ), ]

nrow( data ) # 419

# get diagnosis breakdown
x <- table( data$diagnosis )
x
## FSGS  MCD 
##  169  250 
round( x / sum( x ) * 100, 1 )
## FSGS  MCD 
## 40.3 59.7 

# and race breakdown
# first group minor ancestries
data$race[ ! data$race %in% c('White', 'Black', 'Asian') ] <- 'Other'
x <- table( data$race )
x
## Asian Black Other White 
##    21    98    44   256 
round( x / sum( x ) * 100, 1 )
## Asian Black Other White 
##   5.0  23.4  10.5  61.1 
