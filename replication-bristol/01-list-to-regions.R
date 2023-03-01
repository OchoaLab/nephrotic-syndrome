# creates ranges file of query loci, for input to tabix

library(readr)
library(dplyr)
library(ochoalabtools)

# get arguments from command line
args <- args_cli()
# must have input and output
file_in <- args[1]
file_out <- args[2]
# complain if they are missing!
if ( is.na( file_in ) || is.na( file_out ) )
    stop( 'Usage: <file_in> <file_out>' )

# load raw data
data <- read_tsv( file_in, show_col_types = FALSE )
# subset to the only columns we care about, has to be exactly these two only!
data <- select( data, CHR, POS )
# there are many duplicates, subset to unique
data <- distinct( data )
# resort by CHROM/POS, not strictly needed but it's just cleaner
data <- data[ order( data$CHR, data$POS ), ]
# add chr prefix for chromosomes, otherwise they're not matched!!!
data$CHR <- paste0( 'chr', data$CHR )
# save to output
# format requires no header!
write_tsv( data, file_out, col_names = FALSE )
