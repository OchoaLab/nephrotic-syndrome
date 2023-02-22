# creates ranges file of query loci, for input to tabix

library(readr)
library(dplyr)

# load raw data
data <- read_tsv( '2023-02-20_replication_gene_list.csv' )
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
write_tsv( data, 'snp-list.txt', col_names = FALSE )
