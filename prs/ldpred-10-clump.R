library(bigsnpr)
library(ochoalabtools)
library(readr)

# constants
NCORES <- 10

args <- args_cli()
base <- args[1]
if ( is.na( base ) )
    stop( 'Usage: <base>' )

# all processing happens in subdirectory
setwd( base )

# paths
file_in <- paste0( 'betas-', base, '-clean-matched.txt.gz' )
name <- paste0( base, '-ldpred2-ct' )
file_params <- paste0( 'params-', name, '.txt.gz' )

# load filtered sumstats `df_beta`!
df_beta <- read_tsv( file_in, show_col_types = FALSE )

# load LD dataset
# Attach the "bigSNP" object in R session
rds <- 'mac20.rds'
message( 'Reading: ', rds )
obj.bigSNP <- snp_attach( rds )
# Get aliases for useful slots
G    <- obj.bigSNP$genotypes
CHR  <- obj.bigSNP$map$chromosome
POS  <- obj.bigSNP$map$physical.pos

# beta and lpval need to have the same length as ncol(G), CHR and POS
# -> one solution is to use missing values and use the 'exclude' parameter
lpvals <- rep( NA, ncol( G ) )
lpvals[ df_beta$`_NUM_ID_` ] <- -log10(df_beta$p)

# actual run
all_keep <- snp_grid_clumping(
    G,
    CHR,
    POS,
    lpS = lpvals,
    exclude = which(is.na(lpvals)),
    ncores = NCORES
)

# this is a messy object, a list with elements for each chromosome, each of which is a list with elements for every param, each of which is a vector of integers (indexes?) of varying lengths
# so there's no way to save this as a tibble, keep it instead as an R object
save( all_keep, file = paste0( 'all-keep-', name, '.RData' ) )

## # save this data frame of parameters
## params <- attr( all_keep, "grid" )

## # params needed to interpret correctly
## write_tsv( params, file_params )
