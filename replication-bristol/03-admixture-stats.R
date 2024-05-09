library(tidyverse)
library(genio)

# read main demographics table
setwd( '/datacommons/ochoalab/ssns_gwas/array' )
data <- read_tsv( 'patient-data.txt.gz' )
# subset to Bristol
data <- data %>% filter( bristol )

# load this other table because some individuals didn't pass QC
setwd( '/datacommons/ochoalab/ssns_gwas/replication/bristol_data/imputation/post_imp/' )
fam <- read_fam( 'bristol_impute_mac20' )
# subset again
data <- data %>% filter( id %in% fam$id )

# finally, look at admixture analysis
setwd( '/datacommons/ochoalab/ssns_gwas/replication/bristol_data/AF' )
afr <- read_table( 'labels_af.txt', col_names = c('fam', 'id') )
eur <- read_table( 'labels_euro.txt', col_names = c('fam', 'id') )
sas <- read_table( 'labels_sas.txt', col_names = c('fam', 'id') )

# counts at this point
nrow( afr ) #  16 
nrow( eur ) # 428
nrow( sas ) # 111 

# assign to our bigger table
data$ancestry <- 'Other' # default
data$ancestry[ data$id %in% afr$id ] <- 'AFR'
data$ancestry[ data$id %in% eur$id ] <- 'EUR'
data$ancestry[ data$id %in% sas$id ] <- 'SAS'

# final counts!
table( data$ancestry )
## AFR   EUR Other   SAS 
##  16   428    29   111 
