library(tidyverse)

# this script applies additional hacks for CureGN to treat MCD/FSGS as SSNS/SRNS (without age filters, just curious)

# start by splitting discovery data
setwd( '/datacommons/ochoalab/ssns_gwas/curegn' )

# load data
data <- read_tsv( 'patient-data.txt.gz', show_col_types = FALSE )

# keep the two types of interest only, for now we only want the list of IDs and not the detailed diagnosis
ids <- data$id[ data$diagnosis %in% c('MCD', 'FSGS') ]

# outputs to elsewhere
setwd( '/datacommons/ochoalab/ssns_gwas/imputed/prs-new' )

# save lists in subfolders
dir.create( 'test-curegn2' )
write_lines( ids , 'test-curegn2/ids.txt' )
