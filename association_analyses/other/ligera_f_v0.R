library(tidyverse)
library(dplyr)
library(BEDMatrix)
library(genio)
library(qqman)
library(vroom)
library(MASS)
library(ligera)
library(popkin)

name <- '/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/ssns_gwas_maf'
# create a BEDMatrix object, which allows random access to the genotypes
X <- BEDMatrix( name )
# read the full BIM and FAM tables too
bim <- read_bim( name )
fam <- read_fam( name )
phen <- read_phen('/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/gwas/pheno_file.txt')


name_phen <- '/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/gwas/pheno_file.txt'
if ( !is.na( name_phen ) ) {
  # if the phenotype is in a separate PHEN file, load that
  phen <- read_phen( name_phen )
  # reorder phen to match fam, if needed
  phen <- phen[ match( fam$id, phen$id ), ]
  # sanity check
  # when phen has fewer individuals than fam, some are NAs after match above, so this takes care of those cases
  stopifnot( all(phen$id == fam$id, na.rm = TRUE) )
  # save trait now
  trait <- phen$pheno
} else {
  # if the trait of interest is in the fam file, extract it now
  trait <- fam$pheno
}



kinship <- popkin( X)
kinship_inv <- ginv( kinship )
print('kinship finished')
tib_f <- ligera_f( X, trait, kinship, kinship_inv, V = 0)
tib_f <- cbind( bim, tib_f )

file_out <- paste0( 'ligera_f_mafV0', '.txt' )
write_tsv( tib_f, file_out )