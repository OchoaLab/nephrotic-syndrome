library(readr)
library(tibble)
library(genio)
library(ochoalabtools)

# constants
age_cut <- 22

# go where the data is on DCC
setwd('/datacommons/ochoalab/ssns_gwas')

# load data!
data <- read_tsv( 'array/patient-data.txt.gz' )
# subset to Bristol according to ID pattern (precalculated)
data <- data[ data$bristol, ]
nrow( data ) # [1] 590

# location of Bristol data, to subset further in case some individuals were removed (dupes, and didn't pass QC)
fam <- read_fam( 'replication/bristol_data/imputation/post_imp/bristol_impute' )
nrow( fam ) # [1] 584
# intersect
data <- data[ data$id %in% fam$id, ]

# no NAs to worry about! (except for age)
stopifnot( !anyNA( data[ , names( data ) != 'age' ] ) )

table_nice <- function( x ) {
    # rearrange for nicer table
    x <- tibble( value = names( x ), count = as.integer( x ) )
    # add percentages
    x$pc <- round( x$count / sum( x$count ) * 100, 1 )
    # rank so top ancestry is shown first
    x <- x[ order( x$count, decreasing = TRUE ), ]
    # return, just show it if not assigned
    return( x )
}

# basic tables, results in presentations/
table_nice( table( data$sex ) )
table_nice( table( data$diagnosis ) )
table_nice( table( data$race ) )

# look at race by diagnosis
x <- table( data$race, data$diagnosis )
table_nice( x[ ,'SSNS' ] )
table_nice( x[ ,'SRNS' ] )

# age is most important analysis here, as it is hard to visualize as table, needs to be fig
sum( is.na( data$age ) ) # [1] 46

fig_start( 'age', width = 6 )
hist( data$age, breaks = 100, xlab = 'Age', main = '' )
fig_end()

# find out how many are children
n <- sum( data$age < age_cut, na.rm = TRUE )
pc <- round( n/nrow( data ) * 100, 1 )
message( 'children: ', n, ', ', pc, '%' )

# at this point, let's apply age filter, redo demographics there
data <- data[ !is.na( data$age ) & data$age < age_cut, ]
table_nice( table( data$sex ) )
table_nice( table( data$diagnosis ) )
table_nice( table( data$race ) )
x <- table( data$race, data$diagnosis )
table_nice( x[ ,'SSNS' ] )
table_nice( x[ ,'SRNS' ] )
