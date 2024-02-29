library(bigsnpr)
library(readr)
library(genio)
library(ochoalabtools)

# constants
# support old data for now
types_old <- c('ssns_ctrl', 'ssns_srns')

# support old data for now, expect ssns_ctrl or ssns_srns
args <- args_cli()
type_base <- args[1]
type_train <- args[2]
type_test <- args[3]
if ( is.na( type_base ) )
    stop( 'Usage: <type>' )

# handle old and new cases!
if ( type_base %in% types_old ) {
    name_data <- 'data'
    # load indexes of testing individuals
    ind_test <- as.numeric( read_lines( 'ind-test.txt.gz' ) )
    # location of PRSs is local
    dir_in <- ''
    # base is only parameter here
    type_in <- type_base
} else {
    if ( is.na( type_test ) )
        stop( 'Usage: <type_base> <type_train> <type_test>' )
    # all processing happens in subdirectory
    setwd( type_test )
    name_data <- 'mac20'
    # location of PRSs is training!
    dir_in <- paste0( '../', type_train, '/' )
    # combine base and train in new setup
    type_in <- paste0( type_base, '-', type_train )
}

# add a dash to separate parts of path as needed

message( 'Loading testing dataset' )

# need to map SNPs from 'train' into 'test' here!
# load filtered sumstats `df_beta`!
file_df_beta <- paste0( 'betas-', type_in, '-clean-matched.txt.gz' )
df_beta <- read_tsv( file_df_beta, show_col_types = FALSE )

# load testing dataset
# Attach the "bigSNP" object in R session
obj.bigSNP <- snp_attach( paste0( name_data, '.rds' ) )
G <- obj.bigSNP$genotypes
y <- obj.bigSNP$fam$affection
y[ y == 0 ] <- NA # in plink format, zeros are missing, translate appropriately here!
if ( ! type_base %in% types_old )
    # in new setup, we haven't defined ind_test, do it now that we know the number of individuals (use all)
    ind_test <- 1L : length( y )
# load PCs, to condition on when scoring with R2
PCs <- read_eigenvec( name_data )$eigenvec

# names of cases to score in testing data
names <- paste0( '-ldpred2-', c( 'inf-best', 'grid-h0.1-best', 'auto-h0.1', 'lassosum-best' ) )

# process preexisting results
for ( name in names ) {
    # load PRS calculated previously
    file_in <- paste0( dir_in, 'betas-', type_base, name, '.txt.gz' )
    file_out <- paste0( 'cor-', type_in, name, '.txt.gz' )
    
    # skip costly calculations if output already exists!
    if ( file.exists( file_out ) ) next
    # skip silently if input is missing
    if ( !file.exists( file_in ) ) next
    # report what is being processed right now
    message( name )

    # load input
    betas <- as.numeric( read_lines( file_in ) )

    if ( ! type_base %in% types_old )
        # for new pipeline only, subset betas using precalculated map of SNPs from 'train' into 'test'!
        betas <- betas[ df_beta[["_NUM_ID_.ss"]] ]
    
    # calculate PRS for test individuals now
    preds <- big_prodVec( G, betas, ind.row = ind_test, ind.col = df_beta[["_NUM_ID_"]] )
    # calculate and save only correlation coefficient to truth, adjusting for PCs
    cor <- pcor( preds, y[ ind_test ], PCs[ ind_test, ] )
    write_lines( cor, file_out )
}
