library(readr)

# this script is for "base" only, simply maps betas from base onto train, subsetting SNPs as it is done for other trained cases

# constants
type <- 'base'
# previously-chosen h2 estimate around which auto was performed
h2_est <- 0.1

# start in this subdirectory
setwd( type )

# start by loading several previous calculations and other needed data

# suffix shared by several in/out files
name_in <- paste0( '-ldpred2-auto-h', h2_est )

# load previously calculated results
betas <- read_lines( paste0( 'betas', name_in, '.txt.gz' ) )

# rest happens in training set
setwd( '../train' )

# load filtered sumstats `df_beta`!
df_beta <- read_tsv( 'betas-clean-matched.txt.gz', show_col_types = FALSE )

# since betas came from base, we need to subset them to what we actually have in training data
# subset betas using precalculated map of SNPs from 'base' into 'train'!
betas <- betas[ df_beta[["_NUM_ID_.ss"]] ]

# now prepend "base" to "name_in" to denote this alternative origin of data, which from here on is only used in outputs!
name_in <- paste0( '-', type, name_in )

# save mapped betas!
# this is a simple vector
write_lines( betas, paste0( 'betas', name_in, '.txt.gz' ) )
