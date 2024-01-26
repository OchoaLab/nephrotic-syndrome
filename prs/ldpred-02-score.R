library(bigsnpr)
library(readr)
library(ochoalabtools)

# determine which type to run
type <- args_cli()[1]
if ( is.na( type ) )
    stop( 'Usage: <type: ssns_ctrl or ssns_srns>' )

message( 'Loading testing dataset' )

# load filtered sumstats `df_beta`!
df_beta <- read_tsv( paste0( 'betas-', type, '-clean-matched.txt.gz' ), show_col_types = FALSE )

# load testing dataset
# Attach the "bigSNP" object in R session
obj.bigSNP <- snp_attach( 'data.rds' )
G <- obj.bigSNP$genotypes
y <- obj.bigSNP$fam$affection
y[ y == 0 ] <- NA # in plink format, zeros are missing, translate appropriately here!
# and indexes that determine testing subset, for consistency across tests
ind_test <- as.numeric( read_lines( 'ind-test.txt.gz' ) )

# add more special cases
names <- paste0( type, '-ldpred2-', c( 'inf-best', 'grid-h0.1-best', 'auto-h0.1', 'lassosum-best' ) )

# process preexisting results
for ( name in names ) {
    # load PRS calculated previously
    file_in <- paste0( 'betas-', name, '.txt.gz' )
    file_out <- paste0( 'cor-', name, '.txt.gz' )
    
    # skip costly calculations if output already exists!
    if ( file.exists( file_out ) ) next
    # skip silently if input is missing
    if ( !file.exists( file_in ) ) next
    # report what is being processed right now
    message( name )
    
    # load input
    betas <- as.numeric( read_lines( file_in ) )

    # calculate PRS for test individuals now
    preds <- big_prodVec( G, betas, ind.row = ind_test, ind.col = df_beta[["_NUM_ID_"]] )
    # calculate and save only correlation coefficient to truth
    cor <- pcor( preds, y[ ind_test ], NULL )
    write_lines( cor, file_out )
}
