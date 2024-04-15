library(bigsnpr)
library(readr)
library(genio)
library(ochoalabtools)

# constants
name_data <- 'mac20'

args <- args_cli()
base <- args[1]
train <- args[2]
test <- args[3]
if ( is.na( test ) )
    stop( 'Usage: <base> <train> <test>' )

# all processing happens in subdirectory
setwd( test )
# combine base and train in new setup
base_train <- paste0( base, '-', train )

# add a dash to separate parts of path as needed

message( 'Loading testing dataset' )

# need to map SNPs from 'train' into 'test' here!
# load filtered sumstats `df_beta`!
file_df_beta <- paste0( 'betas-', base_train, '-clean-matched.txt.gz' )
df_beta <- read_tsv( file_df_beta, show_col_types = FALSE )

# load testing dataset
# Attach the "bigSNP" object in R session
obj.bigSNP <- snp_attach( paste0( name_data, '.rds' ) )
G <- obj.bigSNP$genotypes
y <- obj.bigSNP$fam$affection
y[ y == 0 ] <- NA # in plink format, zeros are missing, translate appropriately here!
# load PCs, to condition on when scoring with R2
PCs <- read_eigenvec( name_data )$eigenvec

# names of cases to score in testing data
names <- paste0( '-ldpred2-', c( 'inf-best', 'grid-h0.1-best', 'auto-h0.1', 'lassosum-best', 'ct-best', 'ct-stacked' ) )

# process preexisting results
for ( name in names ) {
    # load PRS calculated previously from training
    # location of PRSs is training!
    file_in <- paste0( '../', train, '/betas-', base, name, '.txt.gz' )
    file_out <- paste0( 'cor-', base_train, name, '.txt.gz' )
    file_prs <- paste0( 'prs-', base_train, name, '.txt.gz' )
    
    # skip costly calculations if both outputs already exist!
    if ( file.exists( file_out ) && file.exists( file_prs ) ) next
    # skip silently if input is missing
    if ( !file.exists( file_in ) ) next
    # report what is being processed right now
    message( name )

    # load input
    betas <- as.numeric( read_lines( file_in ) )

    # subset betas using precalculated map of SNPs from 'train' into 'test'!
    betas <- betas[ df_beta[["_NUM_ID_.ss"]] ]
    
    # calculate PRS for test individuals now
    preds <- big_prodVec( G, betas, ind.col = df_beta[["_NUM_ID_"]] )

    # save PRS to analyze further outside this script
    write_lines( preds, file_prs )
    
    # calculate and save only correlation coefficient to truth, adjusting for PCs
    cor <- pcor( preds, y, PCs )
    write_lines( cor, file_out )
}
