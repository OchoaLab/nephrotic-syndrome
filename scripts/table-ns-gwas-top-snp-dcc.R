## # log into DCC
## ssh $dcc
## # create an interactive session in a node, with enough memory to load this data
## # default (I think 1G) is not enough, even 4G wasn't enough
## srun -p ochoalab --account ochoalab --mem=8G --pty bash -i
## # load R
## module load R/4.0.0
## R

### R code

library(genio)

# go where the data is
setwd('/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205')

# load genetic data
data <- read_plink('ssns_gwas_maf_dedup')
X <- data$X
bim <- data$bim
fam <- data$fam

# load phenotype
phen <- read_phen('gwas/pheno_file.txt')

# these are already aligned, validate assumption here
stopifnot( all( fam$id == phen$id ) )

# want this SNP only!
top_snp_id <- '11:60078020'

i <- which( bim$id == top_snp_id )
# 551307

xi <- X[ i, , drop = FALSE ]
bimi <- bim[ i, ]

xi_char <- geno_to_char( xi, bimi )

# AF agrees with GCTA's calculation, minor allele is rare
mean( xi, na.rm = TRUE ) / 2
#[1] 0.05570152

# relatively high missingness!
mean( is.na( xi ) )
# [1] 0.03937405


# this is what I really wanted!
table(xi_char, phen$pheno, useNA='i')
## xi_char   1   2
##    C/A   72 140
##    C/C  949 742
##    <NA>  28  50
