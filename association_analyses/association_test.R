## GCTA
### make grm
gcta --bfile /datacommons/ochoalab/ssns_gwas/imputation/merge/ssns_tgp --autosome --make-grm --out /datacommons/ochoalab/tiffany_data/DATA_NS_MERGE/ssns_tgp

gcta --mlma --bfile /datacommons/ochoalab/ssns_gwas/imputation/merge/ssns_tgp --grm /datacommons/ochoalab/tiffany_data/DATA_NS_MERGE/ssns_tgp --pheno /hpc/home/tt207/NS_data/pheno_filtered.txt --out /datacommons/ochoalab/tiffany_data/ssns_tgp_gcta --thread-num 10

## LIGERA
ligera_f.R script

## GMMAT
library(tidyverse)
library(GMMAT)
library(BEDMatrix)
.libPaths("/datacommons/ochoalab/tiffany_data")
library(genio)