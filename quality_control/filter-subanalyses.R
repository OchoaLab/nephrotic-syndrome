library(readr)
library(dplyr)

# this is where our data is
setwd( '/datacommons/ochoalab/ssns_gwas/imputed/' )

# load phenotype and fixed covariates file
data <- read_tsv( 'patient-data.txt.gz', col_types = 'cccccdiiii' )

# do diagnosis subtypes first
for ( type in c('ssns_ctrl', 'srns_ctrl', 'ssns_srns') ) {
    dir.create( type )
    ids <- data$id[ !is.na( data[ , type ] ) ]
    write_lines( ids, paste0( type, '/ids.txt' ) )
}

# also do ancestry subsets
# actual filtering applies both filters in series, so here it suffices to define IDs separately
for ( my_ancestry in c('afr', 'eur', 'sas') ) {
    dir.create( my_ancestry )
    # base R filtering was failing because of comparison to NAs, this tidyverse version works though
    ids <- filter( data, ancestry == my_ancestry ) %>% pull( id )
    write_lines( ids, paste0( my_ancestry, '/ids.txt' ) )
}
