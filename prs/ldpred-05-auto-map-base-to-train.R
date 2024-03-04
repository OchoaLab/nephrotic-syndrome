library(readr)
library(ochoalabtools)

# this script is for new setup only, simply maps betas from base onto train, subsetting SNPs as it is done for other trained cases

# determine which type to run
args <- args_cli()
base <- args[1]
train <- args[2]
if ( is.na( train ) )
    stop( 'Usage: <base> <train>' )

# constants
# previously-chosen h2 estimate around which auto was performed
h2_est <- 0.1

# start in this subdirectory
setwd( base )

# start by loading several previous calculations and other needed data

# suffix shared by several in/out files
name_in <- paste0( base, '-ldpred2-auto-h', h2_est )

# load previously calculated results
betas <- read_lines( paste0( 'betas-', name_in, '.txt.gz' ) )

# rest happens in training set
setwd( paste0( '../', train ) )

# load filtered sumstats `df_beta`!
df_beta <- read_tsv( paste0( 'betas-', base, '-clean-matched.txt.gz' ), show_col_types = FALSE )

# since betas came from base, we need to subset them to what we actually have in training data
# subset betas using precalculated map of SNPs from 'base' into 'train'!
betas <- betas[ df_beta[["_NUM_ID_.ss"]] ]

# save mapped betas!
# this is a simple vector
write_lines( betas, paste0( 'betas-', name_in, '.txt.gz' ) )
