# uses top 4 loci from main/conditional analysis to test for interaction models (epistasis)

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

# try version without interactions
mod2 <- lm( y ~ X + PCs )
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

# try version with interactions
# better to name each part
hla1 <- X[,1]
hla2 <- X[,2]
hla3 <- X[,3]
clec16a <- X[,4]
# here R2 always increases with interactions, so that's not meaningful unless it were a huge difference.  The per-interaction p-values are more useful, sadly none are significant
mod2i <- lm( y ~ (hla1 + hla2 + hla3 + clec16a)^2 + PCs )
summary( mod2i )
##              Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   1.54076    0.07374  20.894   <2e-16 ***
## hla1          0.08574    0.04894   1.752   0.0804 .  
## hla2          0.11561    0.08072   1.432   0.1527    
## hla3         -0.44612    0.31567  -1.413   0.1582    
## clec16a      -0.16941    0.07549  -2.244   0.0253 *  
## hla1:hla2     0.01569    0.05146   0.305   0.7605    
## hla1:hla3     0.33470    0.17499   1.913   0.0564 .  
## hla1:clec16a  0.03773    0.04691   0.804   0.4216    
## hla2:hla3    -0.10611    0.17863  -0.594   0.5528    
## hla2:clec16a  0.01886    0.05512   0.342   0.7323    
## hla3:clec16a -0.02868    0.14483  -0.198   0.8431    
r2_2i <- summary( mod2i )$r.squared
# [1] 0.1312978
( r2_2i - r2_0 ) / ( 1 - r2_0 )
# [1] 0.09564643

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
mod2ir <- lm( y ~ (hla1r + hla2r + hla3r + clec16ar)^2 + PCs )
summary( mod2ir )
##                 Estimate Std. Error t value Pr(>|t|)    
## hla1r          -0.082032   0.039598  -2.072   0.0388 *  
## hla2r           0.147540   0.061815   2.387   0.0174 *  
## hla3r                 NA         NA      NA       NA    
## clec16ar       -0.050277   0.039198  -1.283   0.2002    
## hla1r:hla2r    -0.009769   0.063430  -0.154   0.8777    
## hla1r:hla3r           NA         NA      NA       NA    
## hla1r:clec16ar -0.033858   0.058760  -0.576   0.5647    
## hla2r:hla3r           NA         NA      NA       NA    
## hla2r:clec16ar -0.230807   0.123695  -1.866   0.0626 .  
## hla3r:clec16ar        NA         NA      NA       NA    
r2_2ir <- summary( mod2ir )$r.squared
# [1] 0.0727915
( r2_2ir - r2_0 ) / ( 1 - r2_0 )
# [1] 0.03473907

# dominant, better than recessive but worse than the original additive one
hla1d <- Xd[,1]
hla2d <- Xd[,2]
hla3d <- Xd[,3]
clec16ad <- Xd[,4]
mod2id <- lm( y ~ (hla1d + hla2d + hla3d + clec16ad)^2 + PCs )
summary( mod2id )
##                Estimate Std. Error t value Pr(>|t|)    
## hla1d          -0.04640    0.03419  -1.357   0.1754    
## hla2d           0.06727    0.03424   1.965   0.0500 .  
## hla3d           0.11952    0.07144   1.673   0.0950 .  
## clec16ad       -0.07036    0.03566  -1.973   0.0491 *  
## hla1d:hla2d    -0.01206    0.02074  -0.582   0.5610    
## hla1d:hla3d    -0.08216    0.04446  -1.848   0.0652 .  
## hla1d:clec16ad -0.01113    0.02031  -0.548   0.5841    
## hla2d:hla3d    -0.02418    0.04522  -0.535   0.5930    
## hla2d:clec16ad  0.02039    0.02050   0.995   0.3204    
## hla3d:clec16ad -0.01130    0.04033  -0.280   0.7794    
r2_2id <- summary( mod2id )$r.squared
# [1] 0.121662
( r2_2id - r2_0 ) / ( 1 - r2_0 )
# [1] 0.0856152

# conclusion: there's only one sort-of marginally significant interaction, and it's between two hla loci (1 and 3), and the interaction is strongest in the additive rathern than dominant mode (recessive couldn't be tested here because hla3 is fixed in recessive)
