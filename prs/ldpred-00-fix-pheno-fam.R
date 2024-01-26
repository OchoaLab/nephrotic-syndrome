library(genio)
library(readr)

# run from here
setwd('/datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/prs')
# name of data to predict on
name <- 'data'
file_phen <- '/datacommons/ochoalab/ssns_gwas/array/patient-data.txt.gz'

# load original fam
fam <- read_fam( name )

# read phen, only file with actual trait
data <- read_tsv( file_phen, show_col_types = FALSE )
stopifnot( all( fam$id %in% data$id ) )

# subset and reorder
indexes <- match( fam$id, data$id )
# reorder phen
data <- data[ indexes, ]
stopifnot( all( fam$id == data$id ) )

# make numerical version, in plink format
fam$pheno <- 0 # missing
fam$pheno[ data$diagnosis == 'SRNS' ] <- 1 # controls
fam$pheno[ data$diagnosis == 'SSNS' ] <- 2 # cases (to be consistent with ssns_ctrl base gwas)

# rewrite fam (here original was softlink, so whatever)
# unlink original softlink first
file.remove( paste0( name, '.fam' ) )
write_fam( name, fam )
