# uses top HLA types (individually, not as haplotype) to demonstrate that recessive model is better, and by how much?
# this version uses the discovery data for greater power

# need to calculate R2 as in the paper!  Need PCs, etc
# based on ldpred-02-score.R

library(genio)
library(readr)
library(dplyr)
library(tidyr)

# constants
name_data <- 'mac20'
# types of interest
## HLA_DRB1*07:01
## HLA_DQA1*02:01
## HLA_DQB1*02:01
## HLA_DQB1*06:02
## HLA_DQA1*01:02

# we need the data at the base of this dir
setwd( '/datacommons/ochoalab/ssns_gwas/imputed' )

# load haplotypes
# unfortunately this is the only table I have, but we can extract data from here
haps <- read_tsv( 'haplotype_HapsbySubject.txt', show_col_types = FALSE )
haps <- haps %>% rename( id = SAMPLE.ID, hap1 = "DQA1~DQB1~DRB1.Hap1", hap2 = "DQA1~DQB1~DRB1.Hap2" )
# split each haplotype
haps <- haps %>% separate_wider_delim( c( 'hap1', 'hap2' ), delim = '~', names = c('dqa1', 'dqb1', 'drb1'), names_sep = '_' )
# process to get dosage of top types; this was a very manual process because of the way it's organized
haps$drb171 <- (haps$hap1_drb1 == '07:01') + (haps$hap2_drb1 == '07:01')
haps$dqa121 <- (haps$hap1_dqa1 == '02:01') + (haps$hap2_dqa1 == '02:01')
haps$dqa112 <- (haps$hap1_dqa1 == '01:02') + (haps$hap2_dqa1 == '01:02')
haps$dqb121 <- (haps$hap1_dqb1 == '02:01') + (haps$hap2_dqb1 == '02:01')
haps$dqb162 <- (haps$hap1_dqb1 == '06:02') + (haps$hap2_dqb1 == '06:02')

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
drb171 <- data$drb171
dqa121 <- data$dqa121
dqa112 <- data$dqa112
dqb121 <- data$dqb121
dqb162 <- data$dqb162

# calculate fit in linear model
mod0 <- lm( y ~ PCs )
r2_0 <- summary( mod0 )$r.squared
# [1] 0.08393301

# try additive dosage
# reorder variables by their paper p-value
mod2 <- lm( y ~ drb171 + dqa121 + dqb121 + dqb162 + dqa112 + PCs )
summary( mod2 )
##               Estimate Std. Error t value Pr(>|t|)    
## drb171       0.1547861  0.0319466   4.845 1.31e-06 ***
## dqa121      -0.0027525  0.0317783  -0.087  0.93098    
## dqb121       0.0751663  0.0137510   5.466 4.85e-08 ***
## dqb162      -0.0651848  0.0161670  -4.032 5.62e-05 ***
## dqa112       0.0006269  0.0124788   0.050  0.95994    
r2_2 <- summary( mod2 )$r.squared
# [1] 0.1398235
( r2_2 - r2_0 ) / ( 1 - r2_0 )
# [1] 0.06101132

# now try the same with recessive encodings!
# all counted alleles are minor (though here it's all multiallelic, all are minor?); here the meaning is not symmetric, so we want to use counted allele regardless of freq
mean( drb171 )/2 # [1] 0.1526198
mean( dqa121 )/2 # [1] 0.1571906
mean( dqa112 )/2 # [1] 0.1921962
mean( dqb121 )/2 # [1] 0.07993311
mean( dqb162 )/2 # [1] 0.1

# now encode recessive model
# for completeness do dominant one as well, confirm that it is worse
drb171r <- drb171
drb171r[ drb171r == 1 ] <- 0
drb171d <- drb171
drb171d[ drb171d == 1 ] <- 2

dqa121r <- dqa121
dqa121r[ dqa121r == 1 ] <- 0
dqa121d <- dqa121
dqa121d[ dqa121d == 1 ] <- 2

dqa112r <- dqa112
dqa112r[ dqa112r == 1 ] <- 0
dqa112d <- dqa112
dqa112d[ dqa112d == 1 ] <- 2

dqb121r <- dqb121
dqb121r[ dqb121r == 1 ] <- 0
dqb121d <- dqb121
dqb121d[ dqb121d == 1 ] <- 2

dqb162r <- dqb162
dqb162r[ dqb162r == 1 ] <- 0
dqb162d <- dqb162
dqb162d[ dqb162d == 1 ] <- 2

# recessive first.  It's a worse model than additive :(
mod2r <- lm( y ~ drb171r + dqa121r + dqb121r + dqb162r + dqa112r + PCs )
summary( mod2r )
##               Estimate Std. Error t value Pr(>|t|)    
## drb171r      0.0286744  0.0729946   0.393   0.6945    
## dqa121r      0.0398174  0.0723420   0.550   0.5821    
## dqb121r      0.0008525  0.0308444   0.028   0.9780    
## dqb162r     -0.0544798  0.0233169  -2.336   0.0195 *  
## dqa112r     -0.0083182  0.0148181  -0.561   0.5746    
r2_2r <- summary( mod2r )$r.squared
# [1] 0.08961354
( r2_2r - r2_0 ) / ( 1 - r2_0 )
# [1] 0.006200998

# dominant, finally something better than both additive and recessive! though difference is small
mod2d <- lm( y ~ drb171d + dqa121d + dqb121d + dqb162d + dqa112d + PCs )
summary( mod2d )
##              Estimate Std. Error t value Pr(>|t|)    
## drb171d      0.086796   0.016653   5.212 1.95e-07 ***
## dqa121d      0.001482   0.016545   0.090 0.928648    
## dqb121d      0.039735   0.007328   5.423 6.19e-08 ***
## dqb162d     -0.031851   0.009131  -3.488 0.000491 ***
## dqa112d     -0.007464   0.007424  -1.005 0.314736    
r2_2d <- summary( mod2d )$r.squared
# [1] 0.1433968
( r2_2d - r2_0 ) / ( 1 - r2_0 )
# [1] 0.06491207

# here we may as well have a conditional test for dominant over additive, to get its p-value
# maybe this is too muddled because everything is together, only one dominant case is actually significant
mod2ad <- lm( y ~ drb171 + dqa121 + dqb121 + dqb162 + dqa112 + drb171d + dqa121d + dqb121d + dqb162d + dqa112d + PCs )
summary( mod2ad )
##              Estimate Std. Error t value Pr(>|t|)    
## drb171       0.005969   0.141931   0.042  0.96646    
## dqa121       0.021267   0.140743   0.151  0.87990    
## dqb121      -0.024803   0.061172  -0.405  0.68515    
## dqb162      -0.078615   0.045876  -1.714  0.08666 .  
## dqa112       0.071108   0.029745   2.391  0.01686 *  
## drb171d      0.082895   0.074009   1.120  0.26275    
## dqa121d     -0.008691   0.073404  -0.118  0.90575    
## dqb121d      0.053834   0.032805   1.641  0.10087    
## dqb162d      0.007846   0.026174   0.300  0.76438    
## dqa112d     -0.045798   0.017768  -2.577  0.00998 ** 
r2_2ad <- summary( mod2ad )$r.squared
# [1] 0.1447027
( r2_2ad - r2_0 ) / ( 1 - r2_0 )
# [1] 0.06633757

# pair off each type additive vs dominant
# yey, in every single case dominant is better than significant, and in 3 cases it's significant for dominance but not for additive!

mod2ad1 <- lm( y ~ drb171 + drb171d + PCs )
summary( mod2ad1 )
##             Estimate Std. Error t value Pr(>|t|)    
## drb171       0.01871    0.03394   0.551  0.58143    
## drb171d      0.08238    0.01935   4.257 2.12e-05 ***
r2_2ad1 <- summary( mod2ad1 )$r.squared
# 0.1306884
( r2_2ad1 - r2_0 ) / ( 1 - r2_0 )
# 0.05103927

mod2ad2 <- lm( y ~ dqa121 + dqa121d  + PCs )
summary( mod2ad2 )
##              Estimate Std. Error t value Pr(>|t|)    
## dqa121       0.027059   0.033732   0.802 0.422496    
## dqa121d      0.071655   0.019234   3.725 0.000197 ***
r2_2ad2 <- summary( mod2ad2 )$r.squared
# 0.1254311
( r2_2ad2 - r2_0 ) / ( 1 - r2_0 )
# 0.04530022

mod2ad3 <- lm( y ~ dqb121 + dqb121d + PCs )
summary( mod2ad3 )
##              Estimate Std. Error t value Pr(>|t|)    
## dqb121      -0.066706   0.062870  -1.061  0.28874    
## dqb121d      0.074251   0.033756   2.200  0.02788 *  
r2_2ad3 <- summary( mod2ad3 )$r.squared
# 0.08987281
( r2_2ad3 - r2_0 ) / ( 1 - r2_0 )
# 0.006484028

mod2ad4 <- lm( y ~ dqb162 + dqb162d + PCs )
summary( mod2ad4 )
##             Estimate Std. Error t value Pr(>|t|)    
## dqb162      -0.04686    0.03931  -1.192   0.2333    
## dqb162d     -0.03268    0.02292  -1.426   0.1539    
r2_2ad4 <- summary( mod2ad4 )$r.squared
# 0.09723863
( r2_2ad4 - r2_0 ) / ( 1 - r2_0 )
# 0.01452472

mod2ad5 <- lm( y ~ dqa112 + dqa112d + PCs )
summary( mod2ad5 )
##              Estimate Std. Error t value Pr(>|t|)    
## dqa112      -0.001282   0.025175  -0.051  0.95938    
## dqa112d     -0.046231   0.015531  -2.977  0.00293 ** 
r2_2ad5 <- summary( mod2ad5 )$r.squared
# 0.09716229
( r2_2ad5 - r2_0 ) / ( 1 - r2_0 )
# 0.01444139
