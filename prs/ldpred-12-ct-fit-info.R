library(ochoalabtools)
library(readr)

# counts number of SNPs in best model, and reports their locations

# determine which type to run
args <- args_cli()
base <- args[1]
train <- args[2]
if ( is.na( train ) )
    stop( 'Usage: <base> <train>' )

# suffix shared by several in/out files
name_in <- paste0( base, '-ldpred2-ct' )

# all happens in training set only
setwd( train )

# load filtered sumstats `df_beta_train`!
df_beta <- read_tsv( paste0( 'betas-', base, '-clean-matched.txt.gz' ), show_col_types = FALSE )

# load best betas!
# this is a simple vector
betas <- read_lines( paste0( 'betas-', name_in, '-best.txt.gz' ) )

# simply overwrite the betas in the input table
df_beta$beta <- betas

# now subset to non-zero rows
df_beta <- df_beta[ df_beta$beta != 0, ]
# output this subtable
write_tsv( df_beta, paste0( 'report-nonzero-betas-', name_in, '-best.txt.gz' ) )
