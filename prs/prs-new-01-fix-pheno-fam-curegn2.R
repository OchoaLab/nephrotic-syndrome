library(tidyverse)
library(genio)

# constants
name <- 'mac20'
# datasets to process
dataset <- 'test-curegn2'

# start by splitting discovery data
setwd( '/datacommons/ochoalab/ssns_gwas/imputed/prs-new/' )

# load data
data <- read_tsv( '../../../curegn/patient-data.txt.gz', show_col_types = FALSE )

# processing happens in subdirectory of the same name as dataset
setwd( dataset )
# read each of the fam files
fam <- read_fam( name )
# reorder local copy of data to match current fam
data <- data[ match( fam$id, data$id ), ] 
# confirm alignment
stopifnot( all( data$id == fam$id ) )
# data only has IDs, nothing else, populate a bit more just cause it's nice (sex may not be actually used)
# curegn already has required numerical format for sex
fam$sex <- data$sex
# finally, phenotype, actually the key thing we wanted!
# make sure in all cases SSNS is coded as 1, rest are zero
# use plink format, so all get shifted by 1 (1=control, 2=case) https://www.cog-genomics.org/plink/1.9/formats#fam
fam$pheno <- 1 + ( data$diagnosis == 'MCD' )
# if we did this right, all values are 1 or 2 (no missing cases)
stopifnot( all( fam$pheno %in% 1:2 ) )

# move original file, to preserve it in case something goes wrong
file.rename( paste0( name, '.fam' ), paste0( name, '_ORIG.fam' ) )
# write new version!
write_fam( name, fam )
