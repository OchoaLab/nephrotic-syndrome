# combines odds ratios using a crude form of meta analysis

library(genio)
library(readr)
library(testthat)
source('pcor_meta.R')

# constants
# needed to replicate CI calculations
# includes intercept (1) plus 10 PCs
df_null <- 11
# combine the first two, into the last one
tests <- c( 'test', 'test-curegn', 'test-bristol-curegn' )

# sample sizes are essentially fixed in each of the two datasets, just look at the fam tables!
df1 <- count_lines( paste0( tests[1], '/mac20.fam' ) ) - df_null
df2 <- count_lines( paste0( tests[2], '/mac20.fam' ) ) - df_null

# read the two datasets
# the part shared across datasets
file <- '/cor-top4.txt.gz'
cors1 <- as.numeric( read_lines( paste0( tests[1], file ) ) )
cors2 <- as.numeric( read_lines( paste0( tests[2], file ) ) )
# ready for our meta-analysis!
cors3 <- pcor_meta( cors1[1], df1, cors2[1], df2 )
write_lines( cors3, paste0( tests[3], file ) )
