library(readr)

# this is where our data is
setwd( '/datacommons/ochoalab/ssns_gwas/imputed/' )

# load phenotype and fixed covariates file
data <- read_tsv( 'patient-data.txt.gz', col_types = 'ccccdiiii' )

# do diagnosis subtypes first
for ( type in c('ssns_ctrl', 'srns_ctrl', 'ssns_srns') ) {
    dir.create( type )
    ids <- data$id[ !is.na( data[ , type ] ) ]
    write_lines( ids, paste0( type, '/ids.txt' ) )
}

