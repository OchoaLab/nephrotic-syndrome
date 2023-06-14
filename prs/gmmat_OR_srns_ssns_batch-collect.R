# original downloaded from here:
# /datacommons/ochoalab/ssns_gwas/GMMAT_0418/PRS/

library(readr)
library(dplyr)

# hardcoded parameters
n_batches <- 100

# tibble to grow
data <- NULL
# boolean to know if we have to sort again or not
do_sort <- FALSE

# shared base
base <- '/datacommons/ochoalab/ssns_gwas/GMMAT_0418/PRS/glmm.wald_srns_ssns'
# main output
file_out <- paste0( base, ".txt.gz" )

# if a collected output already exists, load it and merge it with the new data (should be separate slice)
if ( file.exists( file_out ) ) {
    data <- read_tsv( file_out, show_col_types = FALSE )
    # will have to sort new data before writing
    do_sort <- TRUE
    # create new output
    file_out <- paste0( base, "_merged.txt.gz" )
}

# load and merge batches
for ( batch in 1 : n_batches ) {
    file_in <- paste0( base, "_batch-", batch, '-', n_batches, ".txt" )
    data_batch <- read_tsv( file_in, show_col_types = FALSE )
    data <- bind_rows( data, data_batch )
}

if ( do_sort )
    data <- arrange( data, CHR, POS )

# write it all out
write_tsv( data, file_out )
