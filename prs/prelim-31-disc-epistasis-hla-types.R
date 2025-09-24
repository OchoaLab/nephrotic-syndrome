# uses top HLA types (individually, not haplotype) to test for interaction models (epistasis)
# (haplotype has no epistasis because it's a single variant; would have to match against clec16a separately somehow).

# need to calculate R2 as in the paper!  Need PCs, etc
# based on ldpred-02-score.R

library(genio)
library(readr)
library(dplyr)
library(tidyr)

# constants
name_data <- 'mac20'

# we need the data at the base of this dir
setwd( '/datacommons/ochoalab/ssns_gwas/imputed' )

# load haplotypes
# unfortunately this is the only table I have, but we can extract data from here
haps <- read_tsv( 'haplotype_HapsbySubject.txt', show_col_types = FALSE )
haps <- haps %>% rename( id = SAMPLE.ID, hap1 = "DQA1~DQB1~DRB1.Hap1", hap2 = "DQA1~DQB1~DRB1.Hap2" )
# split each haplotype
haps <- haps %>% separate_wider_delim( c( 'hap1', 'hap2' ), delim = '~', names = c('dqa1', 'dqb1', 'drb1'), names_sep = '_' )
# process to get dosage of top types; this was a very manual process because of the way it's organized
# use the best independent set only
haps$dqb121 <- (haps$hap1_dqb1 == '02:01') + (haps$hap2_dqb1 == '02:01')
haps$drb171 <- (haps$hap1_drb1 == '07:01') + (haps$hap2_drb1 == '07:01')
haps$dqb162 <- (haps$hap1_dqb1 == '06:02') + (haps$hap2_dqb1 == '06:02')
haps$drb1151 <- (haps$hap1_drb1 == '15:01') + (haps$hap2_drb1 == '15:01')
haps$drb1131 <- (haps$hap1_drb1 == '13:01') + (haps$hap2_drb1 == '13:01')
haps$dqb151 <- (haps$hap1_dqb1 == '05:01') + (haps$hap2_dqb1 == '05:01')

# load PCs
obj <- read_eigenvec( name_data )
fam <- obj$fam
PCs <- obj$eigenvec
# and phenotypes
data <- read_tsv( 'patient-data.txt.gz', show_col_types = FALSE )
# haps isn't aligned to anything else, fix by merging into data!
data <- full_join( data, haps )
# data and fam are not aligned, fix!
indexes <- match( fam$id, data$id )
data <- data[ indexes, ]
stopifnot( all( data$id == fam$id ) ) # confirm alignment
# trait here is SSNS or not (combines controls and SRNS)
# treat 14 unclassifieds as NA though
y <- data$diagnosis == 'SSNS'
y[ data$diagnosis == 'NS UNCLASSIFIED' ] <- NA

# extract aligned types
dqb121 <- data$dqb121
drb171 <- data$drb171
dqb162 <- data$dqb162
drb1151 <- haps$drb1151
drb1131 <- haps$drb1131
dqb151 <- haps$dqb151

# calculate fit in linear model
mod0 <- lm( y ~ PCs )
r2_0 <- summary( mod0 )$r.squared
# [1] 0.08393301

# try version without interactions
# reorder variables by their paper p-value
mod2 <- lm( y ~ dqb121 + drb171 + dqb162 + drb1151 + drb1131 + dqb151 + PCs )
summary( mod2 )
##              Estimate Std. Error t value Pr(>|t|)    
## dqb121       0.065326   0.013577   4.811 1.55e-06 ***
## drb171       0.141822   0.010609  13.368  < 2e-16 ***
## dqb162      -0.058105   0.012497  -4.650 3.42e-06 ***
## drb1151     -0.066573   0.013115  -5.076 4.01e-07 ***
## drb1131     -0.088667   0.018683  -4.746 2.14e-06 ***
## dqb151      -0.060720   0.011399  -5.327 1.05e-07 ***
r2_2 <- summary( mod2 )$r.squared
# [1] 0.1523697
( r2_2 - r2_0 ) / ( 1 - r2_0 )
# [1] 0.07470703
# get anova p-value
obj <- anova( mod0, mod2 )
obj$`Pr(>F)`
# [1]           NA 1.123226e-71

# try version with interactions
# here R2 always increases with interactions, so that's not meaningful unless it were a huge difference.  The per-interaction p-values are more useful, some are significant!
mod2i <- lm( y ~ (dqb121 + drb171 + dqb162 + drb1151 + drb1131 + dqb151)^2 + PCs )
summary( mod2i )
##                  Estimate Std. Error t value Pr(>|t|)    
## dqb121           0.042560   0.017935   2.373  0.01769 *  
## drb171           0.173703   0.013203  13.156  < 2e-16 ***
## dqb162          -0.071355   0.016992  -4.199 2.73e-05 ***
## drb1151         -0.021253   0.018634  -1.141  0.25411    
## drb1131         -0.072556   0.025236  -2.875  0.00406 ** 
## dqb151          -0.019257   0.015558  -1.238  0.21587    
## dqb121:drb171    0.225801   0.031446   7.181 8.10e-13 ***
## dqb121:dqb162   -0.011681   0.044409  -0.263  0.79254    
## dqb121:drb1151  -0.102843   0.039967  -2.573  0.01011 *  
## dqb121:drb1131  -0.112730   0.055134  -2.045  0.04095 *  
## dqb121:dqb151   -0.072692   0.034274  -2.121  0.03399 *  
## drb171:dqb162   -0.071683   0.040381  -1.775  0.07594 .  
## drb171:drb1151  -0.184983   0.027843  -6.644 3.43e-11 ***
## drb171:drb1131  -0.129506   0.045985  -2.816  0.00488 ** 
## drb171:dqb151   -0.134290   0.025617  -5.242 1.66e-07 ***
## dqb162:drb1151   0.039136   0.024177   1.619  0.10557    
## dqb162:drb1131   0.091636   0.044412   2.063  0.03914 *  
## dqb162:dqb151   -0.006867   0.024464  -0.281  0.77896    
## drb1151:drb1131  0.032546   0.074370   0.438  0.66168    
## drb1151:dqb151   0.007281   0.044092   0.165  0.86886    
## drb1131:dqb151   0.055552   0.047896   1.160  0.24617    
r2_2i <- summary( mod2i )$r.squared
# [1] 0.1854231
( r2_2i - r2_0 ) / ( 1 - r2_0 )
# [1] 0.1107889
obj <- anova( mod0, mod2i )
obj$`Pr(>F)`
# [1]           NA 3.154975e-97

# also calculate conditional version: inter|add
( r2_2i - r2_2 ) / ( 1 - r2_2 )
# [1] 0.03899507
obj <- anova( mod2, mod2i )
obj$`Pr(>F)`
# [1]           NA 1.038682e-29

# now encode dominant and recessive models
dom <- function( x ) { x[ x == 1 ] <- 2; x }
rec <- function( x ) { x[ x == 1 ] <- 0; x }
dqb121r <- rec( dqb121 )
dqb121d <- dom( dqb121 )
drb171r <- rec( drb171 )
drb171d <- dom( drb171 )
dqb162r <- rec( dqb162 )
dqb162d <- dom( dqb162 )
drb1151r <- rec( drb1151 )
drb1151d <- dom( drb1151 )
drb1131r <- rec( drb1131 )
drb1131d <- dom( drb1131 )
dqb151r <- rec( dqb151 )
dqb151d <- dom( dqb151 )

# recessive first without interactions
mod2r <- lm( y ~ dqb121r + drb171r + dqb162r + drb1151r + drb1131r + dqb151r + PCs )
summary( mod2r )
##               Estimate Std. Error t value Pr(>|t|)    
## dqb121r      2.338e-05  3.083e-02   0.001  0.99940    
## drb171r      6.714e-02  1.689e-02   3.974 7.17e-05 ***
## dqb162r     -6.105e-02  1.905e-02  -3.205  0.00136 ** 
## drb1151r    -1.735e-02  2.398e-02  -0.723  0.46952    
## drb1131r    -7.944e-02  5.581e-02  -1.424  0.15464    
## dqb151r     -3.485e-02  2.112e-02  -1.650  0.09901 .  
r2_2r <- summary( mod2r )$r.squared
# [1] 0.09054688
( r2_2r - r2_0 ) / ( 1 - r2_0 )
# [1] 0.00721985
obj <- anova( mod0, mod2r )
obj$`Pr(>F)`
## [1]          NA 1.43666e-05

# then with interactions; absolutely sucks (too many zeroes make fixed products)
# overall, worse than additive with interactions
mod2ir <- lm( y ~ (dqb121r + drb171r + dqb162r + drb1151r + drb1131r + dqb151r)^2 + PCs )
summary( mod2ir )
##                     Estimate Std. Error t value Pr(>|t|)    
## dqb121r            6.409e-05  3.084e-02   0.002  0.99834    
## drb171r            6.879e-02  1.697e-02   4.054 5.12e-05 ***
## dqb162r           -6.126e-02  1.990e-02  -3.079  0.00209 ** 
## drb1151r          -1.302e-02  2.560e-02  -0.509  0.61107    
## drb1131r          -7.941e-02  5.582e-02  -1.423  0.15489    
## dqb151r           -3.544e-02  2.143e-02  -1.654  0.09824 .  
## dqb121r:drb171r           NA         NA      NA       NA    
## dqb121r:dqb162r           NA         NA      NA       NA    
## dqb121r:drb1151r          NA         NA      NA       NA    
## dqb121r:drb1131r          NA         NA      NA       NA    
## dqb121r:dqb151r           NA         NA      NA       NA    
## drb171r:dqb162r           NA         NA      NA       NA    
## drb171r:drb1151r  -9.501e-02  8.938e-02  -1.063  0.28784    
## drb171r:drb1131r          NA         NA      NA       NA    
## drb171r:dqb151r           NA         NA      NA       NA    
## dqb162r:drb1151r  -3.847e-03  3.943e-02  -0.098  0.92228    
## dqb162r:drb1131r          NA         NA      NA       NA    
## dqb162r:dqb151r    1.129e-02  6.394e-02   0.177  0.85989    
## drb1151r:drb1131r         NA         NA      NA       NA    
## drb1151r:dqb151r          NA         NA      NA       NA    
## drb1131r:dqb151r          NA         NA      NA       NA    
r2_2ir <- summary( mod2ir )$r.squared
# [1] 0.09078471
( r2_2ir - r2_0 ) / ( 1 - r2_0 )
# 0.007479472
obj <- anova( mod0, mod2ir )
obj$`Pr(>F)`
# [1]           NA 0.0001118132

# also calculate conditional version: inter|add
( r2_2ir - r2_2r ) / ( 1 - r2_2r )
# [1] 0.0002615103
obj <- anova( mod2r, mod2ir )
obj$`Pr(>F)`
# [1]        NA 0.7615871



# dominant, first without interactions
mod2d <- lm( y ~ dqb121d + drb171d + dqb162d + drb1151d + drb1131d + dqb151d + PCs )
summary( mod2d )
##              Estimate Std. Error t value Pr(>|t|)    
## dqb121d      0.035947   0.007238   4.967 7.07e-07 ***
## drb171d      0.083746   0.006013  13.927  < 2e-16 ***
## dqb162d     -0.032977   0.007262  -4.541 5.75e-06 ***
## drb1151d    -0.040314   0.007355  -5.481 4.46e-08 ***
## drb1131d    -0.046678   0.009725  -4.800 1.64e-06 ***
## dqb151d     -0.035265   0.006332  -5.569 2.71e-08 ***
r2_2d <- summary( mod2d )$r.squared
# [1] 0.1570518
( r2_2d - r2_0 ) / ( 1 - r2_0 )
# [1] 0.07981821
obj <- anova( mod0, mod2d )
obj$`Pr(>F)`
## [1]           NA 5.627074e-77


# and dominant with interactions
mod2id <- lm( y ~ (dqb121d + drb171d + dqb162d + drb1151d + drb1131d + dqb151d)^2 + PCs )
summary( mod2id )
##                    Estimate Std. Error t value Pr(>|t|)    
## dqb121d            0.027765   0.010329   2.688 0.007212 ** 
## drb171d            0.104779   0.007956  13.169  < 2e-16 ***
## dqb162d           -0.039909   0.010189  -3.917 9.11e-05 ***
## drb1151d          -0.019113   0.010790  -1.771 0.076577 .  
## drb1131d          -0.037831   0.013712  -2.759 0.005821 ** 
## dqb151d           -0.013114   0.009041  -1.450 0.146992    
## dqb121d:drb171d    0.045555   0.008370   5.443 5.53e-08 ***
## dqb121d:dqb162d   -0.003120   0.011447  -0.273 0.785227    
## dqb121d:drb1151d  -0.028761   0.011327  -2.539 0.011148 *  
## dqb121d:drb1131d  -0.028803   0.013954  -2.064 0.039066 *  
## dqb121d:dqb151d   -0.019583   0.009609  -2.038 0.041618 *  
## drb171d:dqb162d   -0.021215   0.010466  -2.027 0.042720 *  
## drb171d:drb1151d  -0.055084   0.008844  -6.229 5.14e-10 ***
## drb171d:drb1131d  -0.036290   0.012280  -2.955 0.003140 ** 
## drb171d:dqb151d   -0.038318   0.007587  -5.051 4.58e-07 ***
## dqb162d:drb1151d   0.019120   0.008459   2.260 0.023844 *  
## dqb162d:drb1131d   0.023258   0.013245   1.756 0.079163 .  
## dqb162d:dqb151d   -0.001351   0.008000  -0.169 0.865869    
## drb1151d:drb1131d  0.013879   0.018837   0.737 0.461312    
## drb1151d:dqb151d   0.004248   0.011587   0.367 0.713930    
## drb1131d:dqb151d   0.018503   0.012494   1.481 0.138686    
r2_2id <- summary( mod2id )$r.squared
# 0.1862805
( r2_2id - r2_0 ) / ( 1 - r2_0 )
# 0.1117249
obj <- anova( mod0, mod2id )
obj$`Pr(>F)`
## [1]           NA 3.298709e-98

# also calculate conditional version: inter|add
( r2_2id - r2_2d ) / ( 1 - r2_2d )
# [1] 0.03467438
obj <- anova( mod2d, mod2id )
obj$`Pr(>F)`
# [1]           NA 1.031166e-25


# conclusions: dominant with interactions is better than additive with interactions.  In this case interactions are also individually significant
# I had neglected to notice before that dqb162 is protective (both here and in paper), which is why interactions with it are also negative, whereas the two risk types have a positive interaction estimate

# bonferroni threshold for 21 tests:
0.01 / 21
# [1] 0.0004761905
