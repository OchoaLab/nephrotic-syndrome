# uses 8 top C+T loci to demonstrate that recessive model is better, and by how much?

# need to calculate R2 as in the paper!  Need PCs, etc
# based on ldpred-02-score.R

# PACE is so bad right now, don't have readr/genio, so have to rely on other packages
library(BEDMatrix)

# constants
name_data <- 'mac20'
base <- 'base'
train <- 'train'
test <- 'test'

# all processing happens in subdirectory
setwd( test )
# combine base and train in new setup
base_train <- paste0( base, '-', train )

# add a dash to separate parts of path as needed

message( 'Loading testing dataset' )

# need to map SNPs from 'train' into 'test' here!
# load filtered sumstats `df_beta`!
file_df_beta <- paste0( 'betas-', base_train, '-clean-matched.txt.gz' )
df_beta <- read.table( file_df_beta, header = TRUE )

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

# focus on C+T method only
name <- '-ldpred2-ct-best'
file_in <- paste0( '../', train, '/betas-', base, name, '.txt.gz' )
# load input
betas <- as.numeric( readLines( file_in ) )
# subset betas using precalculated map of SNPs from 'train' into 'test'!
# NOTE "X" added by read.table (read_tsv doesn't do this!)
betas <- betas[ df_beta[["X_NUM_ID_.ss"]] ]

# what we actually want to do here is extract the sparse subset found before
indexes <- which( betas != 0 )
# had not appreciated that the betas for this model are identical to the initial GWAS betas, which is good for huge models, but in this very sparse one we can definitely refit coeffs and do better
df_beta <- df_beta[ indexes, ]
# get their genotypes from the X matrix
X <- X[ , match( df_beta$rsid, colnames(X) ) ]
# calculate score that assumes coefficients are good
preds <- drop( X %*% df_beta$beta )

# calculate fit in linear model
mod0 <- lm( y ~ PCs )
r2_0 <- summary( mod0 )$r.squared
# [1] 0.03942191
mod1 <- lm( y ~ preds + PCs )
summary( mod1 )
##             Estimate Std. Error t value Pr(>|t|)    
## (Intercept)  1.42523    0.05271  27.039  < 2e-16 ***
## preds        0.05252    0.01042   5.043  6.4e-07 ***
## PCs1        -0.72980    0.45478  -1.605   0.1092    
## PCs2         0.83370    0.45500   1.832   0.0675 .  
## PCs3         0.99156    0.45460   2.181   0.0296 *  
## PCs4        -0.60194    0.45569  -1.321   0.1871    
## PCs5        -0.27618    0.45494  -0.607   0.5441    
## PCs6         1.01473    0.45457   2.232   0.0260 *  
## PCs7         0.59290    0.46047   1.288   0.1985    
## PCs8        -0.58080    0.45593  -1.274   0.2033    
## PCs9        -0.66919    0.45467  -1.472   0.1417    
## PCs10        0.54916    0.45461   1.208   0.2276    
r2_1 <- summary( mod1 )$r.squared
# [1] 0.08547725
( r2_1 - r2_0 ) / ( 1 - r2_0 )
# [1] 0.04794543

# try version that refits coefficients.  Small improvement in R2
# coefficients are very different! Only one SNP is any good apparently, haha
mod2 <- lm( y ~ X + PCs )
summary( mod2 )
##                     Estimate Std. Error t value Pr(>|t|)    
## (Intercept)         1.356321   0.084785  15.997  < 2e-16 ***
## Xchr6:30770823:A:G -0.007187   0.095357  -0.075  0.93995    
## Xchr6:31151645:C:T -0.023645   0.040953  -0.577  0.56394    
## Xchr6:31633043:A:G  0.009274   0.049598   0.187  0.85176    
## Xchr6:32193589:T:C  0.058892   0.045682   1.289  0.19794    
## Xchr6:32593082:A:C  0.097738   0.037100   2.634  0.00869 ** 
## Xchr6:32658925:A:G  0.015877   0.035214   0.451  0.65229    
## Xchr6:32666541:C:A  0.045798   0.041376   1.107  0.26889    
## Xchr6:32669903:T:C  0.087063   0.051264   1.698  0.09007 .  
## PCs1               -0.839431   0.496466  -1.691  0.09150 .  
## PCs2                0.942120   0.486760   1.935  0.05350 .  
## PCs3                0.933565   0.485956   1.921  0.05529 .  
## PCs4               -0.710612   0.458260  -1.551  0.12161    
## PCs5               -0.267718   0.458839  -0.583  0.55984    
## PCs6                1.117551   0.459516   2.432  0.01537 *  
## PCs7                0.682950   0.469243   1.455  0.14618    
## PCs8               -0.521360   0.458937  -1.136  0.25650    
## PCs9               -0.758475   0.457196  -1.659  0.09775 .  
## PCs10               0.606045   0.460283   1.317  0.18855    
r2_2 <- summary( mod2 )$r.squared
# [1] 0.09766402
( r2_2 - r2_0 ) / ( 1 - r2_0 )
# [1] 0.06063235

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

# recessive first
mod2r <- lm( y ~ Xr + PCs )
# yey, now two SNPs are significant!!!  However, R2 got worse :(
summary( mod2r )
##                     Estimate Std. Error t value Pr(>|t|)    
## (Intercept)          1.60331    0.02598  61.713  < 2e-16 ***
## Xrchr6:30770823:A:G  0.10556    0.24124   0.438 0.661893    
## Xrchr6:31151645:C:T  0.01888    0.06408   0.295 0.768442    
## Xrchr6:31633043:A:G -0.01972    0.04845  -0.407 0.684195    
## Xrchr6:32193589:T:C  0.01363    0.04131   0.330 0.741642    
## Xrchr6:32593082:A:C  0.08649    0.03072   2.815 0.005064 ** 
## Xrchr6:32658925:A:G  0.03790    0.03136   1.209 0.227364    
## Xrchr6:32666541:C:A  0.10764    0.03195   3.369 0.000813 ***
## Xrchr6:32669903:T:C -0.04658    0.07494  -0.622 0.534520    
## PCs1                -0.86312    0.47087  -1.833 0.067397 .  
## PCs2                 0.70185    0.46501   1.509 0.131846    
## PCs3                 1.02891    0.46868   2.195 0.028601 *  
## PCs4                -0.68299    0.45846  -1.490 0.136922    
## PCs5                -0.10098    0.45853  -0.220 0.825789    
## PCs6                 1.06717    0.46069   2.316 0.020939 *  
## PCs7                 0.44321    0.46754   0.948 0.343612    
## PCs8                -0.51924    0.45934  -1.130 0.258843    
## PCs9                -0.71359    0.46035  -1.550 0.121750    
## PCs10                0.48662    0.45874   1.061 0.289308    
r2_2r <- summary( mod2r )$r.squared
# [1] 0.08435824
( r2_2r - r2_0 ) / ( 1 - r2_0 )
# [1] 0.0467805

# dominant, better than recessive but worse than the original additive one
mod2d <- lm( y ~ Xd + PCs )
summary( mod2d )
##                      Estimate Std. Error t value Pr(>|t|)    
## (Intercept)          1.566960   0.054958  28.512  < 2e-16 ***
## Xdchr6:30770823:A:G -0.001771   0.049404  -0.036  0.97142    
## Xdchr6:31151645:C:T -0.011659   0.023324  -0.500  0.61738    
## Xdchr6:31633043:A:G  0.012317   0.027409   0.449  0.65336    
## Xdchr6:32193589:T:C  0.042736   0.028605   1.494  0.13581    
## Xdchr6:32593082:A:C  0.047363   0.025843   1.833  0.06744 .  
## Xdchr6:32658925:A:G  0.009521   0.024245   0.393  0.69469    
## Xdchr6:32666541:C:A  0.007441   0.027003   0.276  0.78299    
## Xdchr6:32669903:T:C -0.073805   0.028498  -2.590  0.00988 ** 
## PCs1                -0.686130   0.491593  -1.396  0.16342    
## PCs2                 0.899454   0.485141   1.854  0.06433 .  
## PCs3                 0.928014   0.482694   1.923  0.05510 .  
## PCs4                -0.732589   0.462599  -1.584  0.11391    
## PCs5                -0.333961   0.463980  -0.720  0.47200    
## PCs6                 1.061658   0.460579   2.305  0.02157 *  
## PCs7                 0.665440   0.468027   1.422  0.15571    
## PCs8                -0.507950   0.460334  -1.103  0.27037    
## PCs9                -0.701604   0.459595  -1.527  0.12750    
## PCs10                0.697420   0.463220   1.506  0.13281    
r2_2d <- summary( mod2d )$r.squared
# [1] 0.08974151
( r2_2d - r2_0 ) / ( 1 - r2_0 )
# [1] 0.0523847

# ok, so clearly this model was originally overfit in favor of the additive encoding, sadly we can't include these random bunches of SNPs and expect good behavior
