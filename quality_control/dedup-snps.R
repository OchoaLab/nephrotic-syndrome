library(genio)
library(tidyverse)

# NOTE: this requires `--mem 16G` on DCC!

# identifies SNPs with duplicate IDs, which means duplicate chr:pos (format of ID)
# then looks at concordance:
# - if it's high then we keep first SNP out of the set of duplicates (because they're identical or nearly so, then practically it doesn't matter which one we keep)
# - if it's low, we throw away both (we have no way of knowing which one is the right/better one)

# analysis occurs here
setwd( '/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/' )

# read latest data, need genotypes (X) and bim table
data <- read_plink( 'ssns_gwas_maf' )
bim <- data$bim
X <- data$X

# get duplicate IDs
dup_ids <- bim %>% group_by(id) %>% filter(n()>1) %>% pull(id) %>% unique()

# determine whether to keep one or none for each set of duplicates, based on concordance
ind_remove <- c()
for ( id in dup_ids ) {
    # get indexes
    ind <- grep(id, rownames(X))
    # confirm there's only two of each duplicate (code assumes that)
    stopifnot( length( ind ) == 2L )
    # get subset
    x_temp <- X[ind,]
    # calculate concordance
    x_conc <- mean(x_temp[1,] == x_temp[2,], na.rm = TRUE)
    # keep first if they're concordant, otherwise remove both!
    if ( x_conc > 0.99 ) {
        ind_remove <- c(ind_remove, ind[2])
    } else {
        ind_remove <- c(ind_remove, ind)
    }
}

# remove as desired
X <- X[ -ind_remove, ]
bim <- bim[ -ind_remove, ]

# save new version!
write_plink( "ssns_gwas_maf_dedup", X, bim, data$fam )
