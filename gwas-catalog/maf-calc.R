library(tidyverse)

# need to report effective MAF filters

# two ways:
# 1) report from MAC (fixed) and sample sizes
# 2) actual minimum MAF in data

# minor allele count was fixed in all cases
mac <- 20

# scan files thoroughly, for completeness
files <- list.files( pattern='.tsv.gz$' )

# output tibble
out <- NULL

for ( file in files ) {
    data <- read_tsv( file, show_col_types = FALSE )
    
    # empirical MAF threshold
    af_range <- range( data$effect_allele_frequency )
    maf_empir <- min( af_range[1], 1 - af_range[2] )
    
    # extract sample sizes.  in theory they vary per SNP because of missingness, but because our data is imputed this should not happen
    # because there are odd numbers, I know this is individuals and not alleles (i.e. 2n)
    n_case <- max( data$N_case )
    n_ctrl <- max( data$N_ctrl )
    n <- n_case + n_ctrl
    maf_cut <- mac / ( 2 * n )

    # append to table
    out <- bind_rows(
        out,
        tibble(
            file = file,
            maf_empir = maf_empir,
            maf_cut = maf_cut,
            n = n,
            n_case = n_case,
            n_ctrl = n_ctrl
        )
    )
}

# save table when done
write_tsv( out, 'maf-cuts.txt' )
