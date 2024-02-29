library(bigsnpr)
library(ochoalabtools)
library(readr)

# script calculates LD matrix only at positions we're going to use!

# constants
# set to 1 for local runs, but DCC accepts >1
NCORES <- 30
# window size suggested in ldpred vignette, not sure if there's more typical values to consider!
size_cM <- 3 / 1000
# support old data for now
types_old <- c('ssns_ctrl', 'ssns_srns')

# support old data for now, expect ssns_ctrl or ssns_srns
args <- args_cli()
type_base <- args[1]
type_train <- args[2]
if ( is.na( type_base ) )
    stop( 'Usage: <type>' )

# load precalculated data
# either way assume script is run from correct local path
if ( type_base %in% types_old ) {
    name_train <- 'data'
} else {
    if ( is.na( type_train ) )
        stop( 'Usage: <type_base> <type_train>' )
    # work in desired training subdirectory
    setwd( type_train )
    name_train <- 'mac20'
}

# specify base used in inputs and outputs
file_betas_matched <- paste0( 'betas-', type_base, '-clean-matched.txt.gz' )
# LD backing file base (sbk extension gets added automatically)
file_ld <- paste0( 'ld-', type_base )
# to save LD object too
file_ld_rdata <- paste0( file_ld, '.RData' )


### LOAD DATA ###

# reload precalculated df_beta
message( 'Reading: ', file_betas_matched )
df_beta <- read_tsv( file_betas_matched, show_col_types = FALSE )

# load training dataset
# Attach the "bigSNP" object in R session
rds <- paste0( name_train, '.rds' )
message( 'Reading: ', rds )
obj.bigSNP <- snp_attach( rds )
# Get aliases for useful slots
G    <- obj.bigSNP$genotypes
CHR  <- obj.bigSNP$map$chromosome
POS  <- obj.bigSNP$map$physical.pos
POS2 <- obj.bigSNP$map$genetic.dist

### LD CALCULATIONS ###

# compute LD matrices
for (chr in 1:22) {
    message(chr)
    ## indices in 'df_beta'
    ind.chr <- which( df_beta$chr == chr )
    ## indices in 'G'
    ind.chr2 <- df_beta$`_NUM_ID_`[ ind.chr ]
    # calculate sparse symmetric correlation matrix
    ld0 <- snp_cor(
        G,
        ind.col = ind.chr2,
        size = size_cM,
        infos.pos = POS2[ind.chr2],
        ncores = NCORES
    )
    # combine data across chromosomes
    if ( chr == 1 ) {
        # set up backing file, first time only
        ld <- as_SFBM( ld0, file_ld, compact = TRUE )
    } else {
        ld$add_columns( ld0, nrow( ld ) )
    }
}

# save object, keep backing file, so we can load it in another session
save( ld, file = file_ld_rdata )
