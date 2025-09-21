# uses top 4 loci from main/conditional analysis to demonstrate that recessive model is better, and by how much?

# need to calculate R2 as in the paper!  Need PCs, etc
# based on ldpred-02-score.R

# PACE is so bad right now, don't have readr/genio, so have to rely on other packages
library(BEDMatrix)

# constants
name_data <- 'mac20'
test <- 'test'

# loci to look at, from paper table 1
# only "all" populations (joint analysis)
ids <- c(
    paste0( 'chr6:', c('32652506:C:T', '32689478:C:T', '31361670:A:G') ),
#    'chr10:28810849:C:T', # missing in test data! (fails grep even without alleles)
    'chr16:11077745:A:G'
    )

# all processing happens in subdirectory
setwd( test )

# load genotypes
X <- BEDMatrix( name_data, simple_names = TRUE )
# and phenotype
fam <- read.table( paste0( name_data, '.fam' ) )
colnames( fam ) <- c('fam', 'id', 'pat', 'mat', 'sex', 'pheno')
y <- fam$pheno
y[ y == 0 ] <- NA # in plink format, zeros are missing, translate appropriately here!
# and PCs
PCs <- read.table( paste0( name_data, '.eigenvec' ) )
PCs <- PCs[ , -(1:2) ]
PCs <- as.matrix( PCs )
colnames( PCs ) <- NULL

# get their genotypes from the X matrix
X <- X[ , match( ids, colnames(X) ) ]

# calculate fit in linear model
mod0 <- lm( y ~ PCs )
r2_0 <- summary( mod0 )$r.squared
# [1] 0.03942191

# try version that refits coefficients.
mod2 <- lm( y ~ X + PCs )
# wow, this is actually a way better model than C+T earlier!
summary( mod2 )
##                      Estimate Std. Error t value Pr(>|t|)    
## Xchr6:32652506:C:T   0.127407   0.031143   4.091    5e-05 ***
## Xchr6:32689478:C:T   0.136277   0.035968   3.789 0.000170 ***
## Xchr6:31361670:A:G   0.035876   0.080488   0.446 0.655989    
## Xchr16:11077745:A:G -0.105922   0.031109  -3.405 0.000715 ***
r2_2 <- summary( mod2 )$r.squared
# [1] 0.1201819
( r2_2 - r2_0 ) / ( 1 - r2_0 )
# [1] 0.08407438

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
##                      Estimate Std. Error t value Pr(>|t|)    
## Xrchr6:32652506:C:T  -0.08891    0.03573  -2.488   0.0132 *  
## Xrchr6:32689478:C:T   0.12053    0.05295   2.276   0.0233 *  
## Xrchr6:31361670:A:G        NA         NA      NA       NA    
## Xrchr16:11077745:A:G -0.06562    0.03678  -1.784   0.0750 .  
r2_2r <- summary( mod2r )$r.squared
# [1] 0.06588368
( r2_2r - r2_0 ) / ( 1 - r2_0 )
# [1] 0.02754775

# dominant, better than recessive but worse than the original additive one
mod2d <- lm( y ~ Xd + PCs )
summary( mod2d )
##                      Estimate Std. Error t value Pr(>|t|)    
## Xdchr6:32652506:C:T  -0.07828    0.02048  -3.822 0.000149 ***
## Xdchr6:32689478:C:T   0.06578    0.02090   3.148 0.001744 ** 
## Xdchr6:31361670:A:G   0.02235    0.04036   0.554 0.580025    
## Xdchr16:11077745:A:G -0.06423    0.02028  -3.166 0.001637 ** 
r2_2d <- summary( mod2d )$r.squared
# [1] 0.1102326
( r2_2d - r2_0 ) / ( 1 - r2_0 )
# [1] 0.07371677
