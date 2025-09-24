# uses top 4 loci from main/conditional analysis to test for interaction models (epistasis)
# this version uses the discovery data for greater power

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

# try version without interactions
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

# try version with interactions
# better to name each part
hla1 <- X[,1]
hla2 <- X[,2]
hla3 <- X[,3]
clec16a <- X[,4]
# here R2 always increases with interactions, so that's not meaningful unless it were a huge difference.  The per-interaction p-values are more useful, some are significant!
mod2i <- lm( y ~ (hla1 + hla2 + hla3 + clec16a)^2 + PCs )
summary( mod2i )
##               Estimate Std. Error t value Pr(>|t|)    
## hla1          0.071621   0.012005   5.966 2.62e-09 ***
## hla2          0.025869   0.018334   1.411  0.15831    
## hla3         -0.035792   0.048720  -0.735  0.46260    
## clec16a      -0.032464   0.013222  -2.455  0.01411 *  
## hla1:hla2     0.072059   0.011848   6.082 1.29e-09 ***
## hla1:hla3    -0.009794   0.027175  -0.360  0.71856    
## hla1:clec16a  0.004647   0.009842   0.472  0.63683    
## hla2:hla3    -0.089469   0.033662  -2.658  0.00789 ** 
## hla2:clec16a -0.020392   0.011788  -1.730  0.08373 .  
## hla3:clec16a -0.041557   0.025302  -1.642  0.10058    
r2_2i <- summary( mod2i )$r.squared
# [1] 0.1695736
( r2_2i - r2_0 ) / ( 1 - r2_0 )
# [1] 0.09348727
obj <- anova( mod0, mod2i )
obj$`Pr(>F)`
# [1]           NA 1.144489e-87

# also calculate conditional version: inter|add
( r2_2i - r2_2 ) / ( 1 - r2_2 )
# [1] 0.01092366
obj <- anova( mod2, mod2i )
obj$`Pr(>F)`
# [1]           NA 7.808689e-09


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
hla1r <- Xr[,1]
hla2r <- Xr[,2]
hla3r <- Xr[,3]
clec16ar <- Xr[,4]
# non-interaction version for comparison/conditioning
mod2r <- lm( y ~ hla1r + hla2r + hla3r + clec16ar + PCs )
summary( mod2r )
##              Estimate Std. Error t value Pr(>|t|)    
## hla1r       -0.060064   0.005993 -10.022  < 2e-16 ***
## hla2r        0.064710   0.011353   5.700 1.28e-08 ***
## hla3r       -0.117536   0.046956  -2.503   0.0123 *  
## clec16ar    -0.036785   0.007216  -5.098 3.57e-07 ***
r2_2r <- summary( mod2r )$r.squared
# 0.1173503
( r2_2r - r2_0 ) / ( 1 - r2_0 )
# 0.03647912
obj <- anova( mod0, mod2r )
obj$`Pr(>F)`
# [1]          NA 9.07548e-35

mod2ir <- lm( y ~ (hla1r + hla2r + hla3r + clec16ar)^2 + PCs )
summary( mod2ir )
##                  Estimate Std. Error t value Pr(>|t|)    
## hla1r          -0.0586577  0.0068397  -8.576  < 2e-16 ***
## hla2r           0.0809500  0.0140559   5.759 9.02e-09 ***
## hla3r          -0.1206309  0.0525138  -2.297  0.02166 *  
## clec16ar       -0.0401040  0.0088699  -4.521 6.30e-06 ***
## hla1r:hla2r    -0.0355470  0.0135118  -2.631  0.00855 ** 
## hla1r:hla3r     0.0551749  0.1039893   0.531  0.59573    
## hla1r:clec16ar  0.0054655  0.0074161   0.737  0.46117    
## hla2r:hla3r            NA         NA      NA       NA    
## hla2r:clec16ar -0.0008566  0.0139869  -0.061  0.95117    
## hla3r:clec16ar -0.0089067  0.0668296  -0.133  0.89398    
r2_2ir <- summary( mod2ir )$r.squared
# [1] 0.1188987
( r2_2ir - r2_0 ) / ( 1 - r2_0 )
# [1] 0.03816934
obj <- anova( mod0, mod2ir )
obj$`Pr(>F)`
# [1]           NA 1.233692e-32

# also calculate conditional version: inter|add
( r2_2ir - r2_2r ) / ( 1 - r2_2r )
# [1] 0.001754214
obj <- anova( mod2r, mod2ir )
obj$`Pr(>F)`
# [1]        NA 0.1665846

# dominant, better than recessive but worse than the original additive one
hla1d <- Xd[,1]
hla2d <- Xd[,2]
hla3d <- Xd[,3]
clec16ad <- Xd[,4]
# non-interaction version for comparison/conditioning
mod2d <- lm( y ~ hla1d + hla2d + hla3d + clec16ad + PCs )
summary( mod2d )
##              Estimate Std. Error t value Pr(>|t|)    
## hla1d       -0.085261   0.005895 -14.463  < 2e-16 ***
## hla2d        0.047498   0.005306   8.952  < 2e-16 ***
## hla3d       -0.057863   0.009471  -6.110 1.08e-09 ***
## clec16ad    -0.021526   0.005569  -3.865 0.000112 ***
r2_2d <- summary( mod2d )$r.squared
# 0.1521916
( r2_2d - r2_0 ) / ( 1 - r2_0 )
# 0.0745127
obj <- anova( mod0, mod2d )
obj$`Pr(>F)`
# [1]           NA 1.977647e-73

# dominant with interactions
mod2id <- lm( y ~ (hla1d + hla2d + hla3d + clec16ad)^2 + PCs )
summary( mod2id )
##                 Estimate Std. Error t value Pr(>|t|)    
## hla1d          -0.043608   0.010490  -4.157 3.29e-05 ***
## hla2d           0.110821   0.011668   9.498  < 2e-16 ***
## hla3d          -0.028341   0.020457  -1.385  0.16600    
## clec16ad        0.007472   0.011385   0.656  0.51165    
## hla1d:hla2d    -0.031809   0.005849  -5.439 5.66e-08 ***
## hla1d:hla3d     0.008072   0.009512   0.849  0.39611    
## hla1d:clec16ad -0.012193   0.005845  -2.086  0.03703 *  
## hla2d:hla3d    -0.027220   0.009580  -2.841  0.00451 ** 
## hla2d:clec16ad -0.010849   0.005313  -2.042  0.04121 *  
## hla3d:clec16ad -0.013500   0.009776  -1.381  0.16737    
r2_2id <- summary( mod2id )$r.squared
# [1] 0.1613523
( r2_2id - r2_0 ) / ( 1 - r2_0 )
# [1] 0.08451265
obj <- anova( mod0, mod2id )
obj$`Pr(>F)`
# [1]           NA 2.533478e-78

# also calculate conditional version: inter|add
( r2_2id - r2_2d ) / ( 1 - r2_2d )
# 0.01080506
obj <- anova( mod2d, mod2id )
obj$`Pr(>F)`
# [1]           NA 9.984967e-09

# conclusion: currently, interaction model is better under the additive model (for top loci, marginally additive was also better, perhaps circularly so).  Interactions between hla loci are strong, and with clec16 a bit less so but barely there!
