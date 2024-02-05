library(bigsnpr)
library(readr)
library(genio)
library(ochoalabtools)

# determine which type to run
type <- args_cli()[1]

# handle old and new cases!
if ( is.na( type ) ) {
    # new setup, type isn't used, this will interpolate fine in all cases below
    type <- ''
    # all processing happens in subdirectory
    setwd( 'test' )
    name_data <- 'mac20'
    # location of PRSs is training!
    dir_in <- '../train/'
} else {
    # add a dash to separate parts of path as needed
    type <- paste0( '-', type )
    name_data <- 'data'
    # load indexes of testing individuals
    ind_test <- as.numeric( read_lines( 'ind-test.txt.gz' ) )
    # location of PRSs is local
    dir_in <- ''
}

message( 'Loading testing dataset' )

# TODO: need to map SNPs from 'train' into 'test' here!
# load filtered sumstats `df_beta`!
df_beta <- read_tsv( paste0( 'betas', type, '-clean-matched.txt.gz' ), show_col_types = FALSE )

# load testing dataset
# Attach the "bigSNP" object in R session
obj.bigSNP <- snp_attach( paste0( name_data, '.rds' ) )
G <- obj.bigSNP$genotypes
y <- obj.bigSNP$fam$affection
y[ y == 0 ] <- NA # in plink format, zeros are missing, translate appropriately here!
if ( type == '' ) {
    # in new setup, we haven't defined ind_test, do it now that we know the number of individuals (use all)
    ind_test <- 1L : length( y )
}
# load PCs, to condition on when scoring with R2
PCs <- read_eigenvec( name_data )$eigenvec

# names of cases to score in testing data
names <- paste0( type, '-ldpred2-', c( 'inf-best', 'grid-h0.1-best', 'auto-h0.1', 'lassosum-best' ) )

# process preexisting results
for ( name in names ) {
    # load PRS calculated previously
    file_in <- paste0( dir_in, 'betas', name, '.txt.gz' )
    file_out <- paste0( 'cor', name, '.txt.gz' )
    
    # skip costly calculations if output already exists!
    if ( file.exists( file_out ) ) next
    # skip silently if input is missing
    if ( !file.exists( file_in ) ) next
    # report what is being processed right now
    message( name )

    # load input
    betas <- as.numeric( read_lines( file_in ) )

    if ( type == '' ) {
        # for new pipeline only, subset betas using precalculated map of SNPs from 'train' into 'test'!
        betas <- betas[ df_beta[["_NUM_ID_.ss"]] ]
    }
    
    # calculate PRS for test individuals now
    preds <- big_prodVec( G, betas, ind.row = ind_test, ind.col = df_beta[["_NUM_ID_"]] )
    # calculate and save only correlation coefficient to truth, adjusting for PCs
    cor <- pcor( preds, y[ ind_test ], PCs[ ind_test, ] )
    write_lines( cor, file_out )
}
