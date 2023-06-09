# original downloaded from here:
# /datacommons/ochoalab/ssns_gwas/GMMAT_0418/PRS/

library(readr)
library(dplyr)

# hardcoded parameters
n_batches <- 100

# tibble to grow
data <- NULL

base <- '/datacommons/ochoalab/ssns_gwas/GMMAT_0418/PRS/glmm.wald_ssns_ctrl'

for ( batch in 1 : n_batches ) {
    file_in <- paste0( base, "_batch-", batch, '-', n_batches, ".txt" )
    data_batch <- read_tsv( file_in, show_col_types = FALSE )
    data <- bind_rows( data, data_batch )
}

file_out <- paste0( base, ".txt.gz" )
write_tsv( data, file_out )
