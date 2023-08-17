library(genio)
library(readr)
library(dplyr)

# pheno is already good!  Only need to process eigenvectors and other covariates
# assume local path is that of eigenvectors
# /datacommons/ochoalab/ssns_gwas/GMMAT_0418/srns_ctr/srns_ssns/
# covariates are elsewhere, will use full path for that

# load and combine covariates
# PCs
eigenvec <- read_eigenvec( 'srns_ssns_mac20' )
# and demographics, which doesn't have a header!
covar <- read_table( '/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/ssns_srns_covar.txt', col_names = c('fam', 'id', 'sex', 'race'), col_types = 'cccc' )

# all dimensions agree!
stopifnot( nrow( eigenvec$fam ) == nrow( covar ) )

# flatten eigenvec
PCs <- as_tibble( eigenvec$eigenvec )
colnames( PCs ) <- paste0( 'pc', colnames( PCs ) )
PCs <- bind_cols( eigenvec$fam, PCs )
# merge, matching by fam and id
covar <- full_join( covar, PCs, by = c( 'fam', 'id' ) )
# rename columns to match what PRSice expects
covar <- rename( covar, FID = fam, IID = id )

# write it out!
write_tsv( covar, 'srns_ssns_mac20_covar_prsice.txt' )
