library(bigsnpr)
library(readr)

# constants
# sequence of heritabilities to consider
herits <- (1:9)/10

message( 'Loading testing dataset' )

# load filtered sumstats `df_beta`!
df_beta <- read_tsv( 'betas-ssns_ctrl-array.txt.gz', show_col_types = FALSE )

# load testing dataset
# Attach the "bigSNP" object in R session
obj.bigSNP <- snp_attach( 'data.rds' )
G <- obj.bigSNP$genotypes
y <- as.numeric( read_lines( 'pheno.txt.gz' ) )
# and indexes that determine testing subset, for consistency across tests
ind.test <- as.numeric( read_lines( 'ind-testing.txt.gz' ) )

# turn heritabilities into outputs to process
names <- paste0( 'ldpred2-inf-h', herits )
# add more special cases
names <- c( names, 'ldpred2-grid-h0.1-best', 'ldpred2-auto-h0.1', 'ldpred2-lassosum-best' )

# process preexisting results
for ( name in names ) {
    # load PRS calculated previously
    file_in <- paste0( 'betas-', name, '.txt.gz' )
    file_out <- paste0( 'cor-', name, '.txt.gz' )
    
    # skip costly calculations if output already exists!
    if ( file.exists( file_out ) ) next
    # report what is being processed right now
    message( name )
    
    # load input
    betas <- as.numeric( read_lines( file_in ) )

    # calculate PRS for test individuals now
    preds <- big_prodVec( G, betas, ind.row = ind.test, ind.col = df_beta[["_NUM_ID_"]] )
    # calculate and save only correlation coefficient to truth
    cor <- pcor( preds, y[ ind.test ], NULL )
    write_lines( cor, file_out )
}
