library(readr)
library(tibble)
library(ochoalabtools)

# constants
age_cut <- 18

# load data!
data <- read_tsv( 'pheno_Bristol.csv' )

# no NAs to worry about!
stopifnot( !anyNA( data ) )

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
table_nice( table( data$Sex ) )
table_nice( table( data$DIAGNOSIS ) )
table_nice( table( data$RACE ) )

# look at race by diagnosis
x <- table( data$RACE, data$DIAGNOSIS )
table_nice( x[ ,'SSNS' ] )
table_nice( x[ ,'SRNS' ] )

# age is most important analysis here, as it is hard to visualize as table, needs to be fig
# NOTE: some ages are actually missing!
data$AGE <- as.numeric( data$AGE ) # coerces "null" to NA
sum( is.na( data$AGE ) ) # [1] 47

fig_start( 'age', width = 6 )
hist( data$AGE, breaks = 100, xlab = 'Age', main = '' )
fig_end()

# find out how many are children
n <- sum( data$AGE < age_cut, na.rm = TRUE )
pc <- round( n/nrow( data ) * 100, 1 )
message( 'children: ', n, ', ', pc, '%' )

# at this point, let's apply age filter, redo demographics there
data <- data[ !is.na( data$AGE ) & data$AGE < age_cut, ]
table_nice( table( data$Sex ) )
table_nice( table( data$DIAGNOSIS ) )
table_nice( table( data$RACE ) )
x <- table( data$RACE, data$DIAGNOSIS )
table_nice( x[ ,'SSNS' ] )
table_nice( x[ ,'SRNS' ] )
