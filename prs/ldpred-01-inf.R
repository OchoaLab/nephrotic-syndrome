library(bigsnpr)
library(readr)

# constants
# sequence of heritabilities to consider
herits <- (1:9)/10

# load filtered sumstats `df_beta`!
df_beta <- read_tsv( 'betas-ssns_ctrl-array.txt.gz', show_col_types = FALSE )

# load `corr` data/backing file, matching SNPs in betas
load( 'data-corr.RData' )

# this is fast, just scan a grid of heritability values to pick a decent one!
for ( h2 in herits ) {
    message( h2 )
    
    # not sure if this is random, just in case
    set.seed(1)

    # key calculation!
    betas <- snp_ldpred2_inf( corr, df_beta, h2 = h2 )

    # save betas!
    # this is a simple vector
    file_out <- paste0( 'betas-ldpred2-inf-h', h2, '.txt.gz' )
    write_lines( betas, file_out )
}
