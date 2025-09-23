# uses top 4 loci from main/conditional analysis to demonstrate that recessive model is better, and by how much?
# this version uses the discovery data for greater powr

# need to calculate R2 as in the paper!  Need PCs, etc
# based on ldpred-02-score.R

library(BEDMatrix)
library(genio)
library(readr)

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
obj <- anova( mod0, mod2 )
obj$`Pr(>F)`
# [1]           NA 8.488762e-83

# now try the same with recessive encodings!
# for simplicity re-encode to have minor alleles only (had one case)
flip <- colMeans( X )/2 > 0.5
Xm <- X
if ( any( flip ) )
    Xm[ , flip ] <- 2 - Xm[ , flip ]
# now encode recessive model
Xr <- Xm
Xr[ Xr == 1 ] <- 0
# for completeness do dominant one as well, confirm that it is worse
Xd <- Xm
Xd[ Xd == 1 ] <- 2

# recessive first.  It's a worse model than additive :(
mod2r <- lm( y ~ Xr + PCs )
summary( mod2r )
##                       Estimate Std. Error t value Pr(>|t|)    
## Xrchr6:32652506:C:T  -0.060064   0.005993 -10.022  < 2e-16 ***
## Xrchr6:32689478:C:T   0.064710   0.011353   5.700 1.28e-08 ***
## Xrchr6:31361670:A:G  -0.117536   0.046956  -2.503   0.0123 *  
## Xrchr16:11077745:A:G -0.036785   0.007216  -5.098 3.57e-07 ***
r2_2r <- summary( mod2r )$r.squared
# [1] 0.1173503
( r2_2r - r2_0 ) / ( 1 - r2_0 )
# [1] 0.03647912
obj <- anova( mod0, mod2r )
obj$`Pr(>F)`
# [1]          NA 9.07548e-35

# dominant, better than recessive but worse than the original additive one
mod2d <- lm( y ~ Xd + PCs )
summary( mod2d )
##                       Estimate Std. Error t value Pr(>|t|)    
## Xdchr6:32652506:C:T  -0.085261   0.005895 -14.463  < 2e-16 ***
## Xdchr6:32689478:C:T   0.047498   0.005306   8.952  < 2e-16 ***
## Xdchr6:31361670:A:G  -0.057863   0.009471  -6.110 1.08e-09 ***
## Xdchr16:11077745:A:G -0.021526   0.005569  -3.865 0.000112 ***
r2_2d <- summary( mod2d )$r.squared
# [1] 0.1521916
( r2_2d - r2_0 ) / ( 1 - r2_0 )
# [1] 0.0745127
obj <- anova( mod0, mod2d )
obj$`Pr(>F)`
# [1]           NA 1.977647e-73

# test nested add+dom, then do conditionals
mod2ad <- lm( y ~ X + Xd + PCs )
summary( mod2ad )
##                       Estimate Std. Error t value Pr(>|t|)    
## Xchr6:32652506:C:T    0.068607   0.012517   5.481 4.46e-08 ***
## Xchr6:32689478:C:T    0.076431   0.022968   3.328 0.000883 ***
## Xchr6:31361670:A:G   -0.176827   0.092954  -1.902 0.057196 .  
## Xchr16:11077745:A:G  -0.057837   0.014673  -3.942 8.22e-05 ***
## Xdchr6:32652506:C:T  -0.039705   0.010155  -3.910 9.38e-05 ***
## Xdchr6:32689478:C:T   0.001925   0.013944   0.138 0.890220    
## Xdchr6:31361670:A:G   0.030508   0.049047   0.622 0.533967    
## Xdchr16:11077745:A:G  0.013650   0.010560   1.293 0.196225    
r2_2ad <- summary( mod2ad )$r.squared
# [1] 0.1636301
( r2_2ad - r2_0 ) / ( 1 - r2_0 )
# 0.08699921
# dom|add
( r2_2ad - r2_2 ) / ( 1 - r2_2 )
# [1] 0.003844679
obj <- anova( mod2, mod2ad )
obj$`Pr(>F)`
# [1]          NA 0.001804729
# add|dom
( r2_2ad - r2_2d ) / ( 1 - r2_2d )
# [1] 0.01349183
obj <- anova( mod2d, mod2ad )
obj$`Pr(>F)`
# [1]           NA 2.290622e-12
