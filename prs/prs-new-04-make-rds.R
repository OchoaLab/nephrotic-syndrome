library(bigsnpr)
library(ochoalabtools)

# name of data to predict on
name <- args_cli()[1]

# this generates .bk and .rds files
snp_readBed( paste0( name, '.bed' ) )
