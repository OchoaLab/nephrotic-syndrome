# uses top 4 loci from main/conditional analysis to make crude PRS based on it
# an earlier analysis suggested this was better than ldpred2's!

# train fully using discovery data (base+train in PRS setup, but not using its separation of sets for gwas vs train; here since we use only top 4 loci, training is minimal/negligible; since all loci are assumed ssns loci, controls and srns are combined!)
# then fix those coeffs, use them in testing sets and calculate R2's!

# need to calculate R2 as in the paper!  Need PCs, etc
# based on ldpred-02-score.R

library(BEDMatrix)
library(genio)
library(readr)
library(tibble)
library(dplyr)

# constants
name_data <- 'mac20'

# loci to look at, from paper table 1
# only "all" populations (joint analysis), skipped the chr10 locus that was ancestry-specific
ids <- c(
    paste0( 'chr6:', c('32652506:C:T', '32689478:C:T', '31361670:A:G') ),
    'chr16:11077745:A:G'
    )

# we need the data at the base of this dir
setwd( '/datacommons/ochoalab/ssns_gwas/imputed' )

# load genotypes
X <- BEDMatrix( name_data, simple_names = TRUE )
# and PCs
obj <- read_eigenvec( name_data )
fam <- obj$fam
PCs <- obj$eigenvec
# and phenotypes
data <- read_tsv( 'patient-data.txt.gz', show_col_types = FALSE )
# data and fam are not aligned, fix!
indexes <- match( fam$id, data$id )
data <- data[ indexes, ]
stopifnot( all( data$id == fam$id ) ) # confirm alignment
# trait here is SSNS or not (combines controls and SRNS)
# treat 14 unclassifieds as NA though
y <- data$diagnosis == 'SSNS'
y[ data$diagnosis == 'NS UNCLASSIFIED' ] <- NA

# get their genotypes from the X matrix
X <- X[ , match( ids, colnames(X) ) ]

# calculate fit in linear model
mod0 <- lm( y ~ PCs )
r2_0 <- summary( mod0 )$r.squared
# [1] 0.08393301

# try version that refits coefficients.
mod2 <- lm( y ~ X + PCs )
summary( mod2 )
##                      Estimate Std. Error t value Pr(>|t|)    
## Xchr6:32652506:C:T   0.108561   0.007232  15.012  < 2e-16 ***
## Xchr6:32689478:C:T   0.078311   0.008731   8.970  < 2e-16 ***
## Xchr6:31361670:A:G  -0.120782   0.017970  -6.721 2.03e-11 ***
## Xchr16:11077745:A:G -0.042019   0.007703  -5.455 5.17e-08 ***
r2_2 <- summary( mod2 )$r.squared
# [1] 0.1604021
( r2_2 - r2_0 ) / ( 1 - r2_0 )
# [1] 0.08347547

# extract coefficients, save properly
# NOTE: these coefficients are on the linear scale, though the original analysis was logistic.  Does that make a big difference?  Is it unfair to ldpred2?
betas <- coef( mod2 )[ paste0( 'X', ids ) ]
# make plink-style table for simplicity!
data <- tibble(
    id = ids,
    a1 = sapply( strsplit( ids, ':' ), function(x) x[4] ),
    prs = betas
)

# for consistency, put under prs/train like the others
# use plink format for simplicity!
setwd( 'prs-new/train' )
write_tsv( data, 'top4-plink-score.txt' )

# CureGN has different IDs, let's see if this works
data <- data %>% mutate( id = sub( 'chr', '', id ) )
data$id <- sapply( strsplit( data$id, ':' ), function(x) paste0( x[1:2], collapse=':' ) )
write_tsv( data, 'top4-plink-score-ids-simpler.txt' )

# NOTE: top4 are all missing in the set of array SNPs!  ldpred2 fails to use them because of that!   We'll have to hack our way ahead without ldpred2
## # load this table, PRS coefficients are relative to it
## df_beta <- read_tsv( 'betas-base-clean-matched.txt.gz', show_col_types = FALSE)
## indexes <- match( ids, df_beta$rsid )
