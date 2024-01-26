library(tidyverse)
library(genio)

# constants
name <- 'mac20'
datasets <- c('base', 'train', 'test')

# start by splitting discovery data
setwd( '/datacommons/ochoalab/ssns_gwas/imputed/prs-new/' )

# load data
data_dis <- read_tsv( '../patient-data.txt.gz', show_col_types = FALSE )
# bristol data is separate, load it now too
data_bri <- read_tsv( '../../array/patient-data.txt.gz', show_col_types = FALSE )

# the bristol data wasn't encoded with the ssns_srns column, construct it here same as for discovery
data_bri$ssns_srns <- ifelse( data_bri$diagnosis == 'SSNS', 0, ifelse( data_bri$diagnosis == 'SRNS', 1, NA ) )

# process each case
for ( dataset in datasets ) {
    # decide which data to use
    # only test uses bristol, the rest are based on discovery
    data <- if ( dataset == 'test' ) data_bri else data_dis
    
    # processing happens in subdirectory of the same name as dataset
    setwd( dataset )
    # read each of the fam files
    fam <- read_fam( name )
    # reorder local copy of data to match current fam
    data <- data[ match( fam$id, data$id ), ] 
    # confirm alignment
    stopifnot( all( data$id == fam$id ) )
    # data only has IDs, nothing else, populate a bit more just cause it's nice (sex may not be actually used)
    fam$sex <- data$sex
    # translate sex a bit further until it's standardized
    fam <- fam %>% mutate( sex = ifelse( sex == 'male', 'M', ifelse( sex == 'female', 'F', 'U' ) ) )
    fam$sex <- sex_to_int( fam$sex )
    # finally, phenotype, actually the key thing we wanted!
    # make sure in all cases SSNS is coded as 1, rest are zero (existing ssns_srns is backwards, gets fixed below)
    # use plink format, so all get shifted by 1 (1=control, 2=case) https://www.cog-genomics.org/plink/1.9/formats#fam
    fam$pheno <- 1 + if ( dataset == 'base' ) data$ssns_ctrl else 1 - data$ssns_srns
    # if we did this right, all values are 1 or 2 (no missing cases)
    stopifnot( all( fam$pheno %in% 1:2 ) )
    
    # move original file, to preserve it in case something goes wrong
    file.rename( paste0( name, '.fam' ), paste0( name, '_ORIG.fam' ) )
    # write new version!
    write_fam( name, fam )

    setwd( '..' )
}
