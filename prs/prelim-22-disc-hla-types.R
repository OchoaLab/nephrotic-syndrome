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
## HLA_DQB1*02:01	2.622	1.67E-45
## HLA_DRB1*07:01	2.965	4.44E-43
## HLA_DQA1*02:01	2.498	2.75E-34
## HLA_DQB1*06:02	0.2762	2.60E-15
## HLA_DQA1*01:02	0.5133	3.71E-13
## HLA_DQA1*01:03	0.4446	3.31E-11
## HLA_DQA1*01:01	0.4942	1.81E-10
## HLA_DRB1*15:01	0.3908	6.46E-10
## HLA_DRB1*13:01	0.2522	1.03E-09
## HLA_DQB1*06:03	0.1985	1.99E-09
## HLA_DQB1*05:01	0.5157	1.15E-08

# we need the data at the base of this dir
setwd( '/datacommons/ochoalab/ssns_gwas/imputed' )

# load haplotypes
# unfortunately this is the only table I have, but we can extract data from here
haps <- read_tsv( 'haplotype_HapsbySubject.txt', show_col_types = FALSE )
haps <- haps %>% rename( id = SAMPLE.ID, hap1 = "DQA1~DQB1~DRB1.Hap1", hap2 = "DQA1~DQB1~DRB1.Hap2" )
# split each haplotype
haps <- haps %>% separate_wider_delim( c( 'hap1', 'hap2' ), delim = '~', names = c('dqa1', 'dqb1', 'drb1'), names_sep = '_' )
# process to get dosage of top types; this was a very manual process because of the way it's organized
haps$dqb121 <- (haps$hap1_dqb1 == '02:01') + (haps$hap2_dqb1 == '02:01')
haps$drb171 <- (haps$hap1_drb1 == '07:01') + (haps$hap2_drb1 == '07:01')
haps$dqa121 <- (haps$hap1_dqa1 == '02:01') + (haps$hap2_dqa1 == '02:01')
haps$dqb162 <- (haps$hap1_dqb1 == '06:02') + (haps$hap2_dqb1 == '06:02')
haps$dqa112 <- (haps$hap1_dqa1 == '01:02') + (haps$hap2_dqa1 == '01:02')
# new extras
haps$dqa113 <- (haps$hap1_dqa1 == '01:03') + (haps$hap2_dqa1 == '01:03')
haps$dqa111 <- (haps$hap1_dqa1 == '01:01') + (haps$hap2_dqa1 == '01:01')
haps$drb1151 <- (haps$hap1_drb1 == '15:01') + (haps$hap2_drb1 == '15:01')
haps$drb1131 <- (haps$hap1_drb1 == '13:01') + (haps$hap2_drb1 == '13:01')
haps$dqb163 <- (haps$hap1_dqb1 == '06:03') + (haps$hap2_dqb1 == '06:03')
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
dqa121 <- data$dqa121
dqb162 <- data$dqb162
dqa112 <- data$dqa112
dqa113 <- haps$dqa113
dqa111 <- haps$dqa111
drb1151 <- haps$drb1151
drb1131 <- haps$drb1131
dqb163 <- haps$dqb163
dqb151 <- haps$dqb151

# calculate fit in linear model
mod0 <- lm( y ~ PCs )
r2_0 <- summary( mod0 )$r.squared
# [1] 0.08393301

# try additive dosage
# order variables by their paper p-value
# this is the set that survives significance adding one at the time, dropping insignificant cases (either immediately or in a later round, removing most insignificant when two or more arise)
mod2 <- lm( y ~ dqb121 + drb171 + dqb162 + drb1151 + drb1131 + dqb151 + PCs )
# dropped: dqa121, dqa112, dqa113, dqb163, dqa111
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

# now try the same with recessive encodings!
# all counted alleles are minor (though here it's all multiallelic, all are minor?); here the meaning is not symmetric, so we want to use counted allele regardless of freq
mean( dqb121 )/2 # [1] 0.07993311
mean( drb171 )/2 # [1] 0.1526198
mean( dqb162 )/2 # [1] 0.1
mean( drb1151 )/2 # 0.07915273
mean( drb1131 )/2 # 0.03812709
mean( dqb151 )/2 # 0.1124861

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

# recessive first.  It's a worse model than additive :(
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
# [1]          NA 1.43666e-05

# dominant, finally something better than both additive and recessive! though difference is small
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

# need combined additive + dominant model, only to calculate conditionals with it
mod2ad <- lm( y ~ dqb121 + drb171 + dqb162 + drb1151 + drb1131 + dqb151 +
                  dqb121d + drb171d + dqb162d + drb1151d + drb1131d + dqb151d + PCs )
#summary( mod2ad )
r2_2ad <- summary( mod2ad )$r.squared
# [1] 0.1574513
# dom|add
( r2_2ad - r2_2 ) / ( 1 - r2_2 )
# 0.005995126
obj <- anova( mod2, mod2ad )
obj$`Pr(>F)`
# [1]           NA 0.0001606887
# add|dom
( r2_2ad - r2_2d ) / ( 1 - r2_2d )
# 0.0004738962
obj <- anova( mod2d, mod2ad )
obj$`Pr(>F)`
# [1]       NA 0.909373

# pair off each type additive vs dominant
# yey, in all but one case dominant is better than additive, but only in 1 case it's significant for dominance but not for additive!

mod2ad1 <- lm( y ~ dqb121 + dqb121d + PCs )
summary( mod2ad1 )
##              Estimate Std. Error t value Pr(>|t|)    
## dqb121      -0.066706   0.062870  -1.061  0.28874    
## dqb121d      0.074251   0.033756   2.200  0.02788 *  
r2_2ad1 <- summary( mod2ad1 )$r.squared
# 0.08987281
( r2_2ad1 - r2_0 ) / ( 1 - r2_0 )
# 0.006484028

mod2ad2 <- lm( y ~ drb171 + drb171d + PCs )
summary( mod2ad2 )
##             Estimate Std. Error t value Pr(>|t|)    
## drb171       0.01871    0.03394   0.551  0.58143    
## drb171d      0.08238    0.01935   4.257 2.12e-05 ***
r2_2ad2 <- summary( mod2ad2 )$r.squared
# 0.1306884
( r2_2ad2 - r2_0 ) / ( 1 - r2_0 )
# 0.05103927

mod2ad3 <- lm( y ~ dqb162 + dqb162d + PCs )
summary( mod2ad3 )
##             Estimate Std. Error t value Pr(>|t|)    
## dqb162      -0.04686    0.03931  -1.192   0.2333    
## dqb162d     -0.03268    0.02292  -1.426   0.1539    
r2_2ad3 <- summary( mod2ad3 )$r.squared
# 0.09723863
( r2_2ad3 - r2_0 ) / ( 1 - r2_0 )
# 0.01452472

mod2ad4 <- lm( y ~ drb1151 + drb1151d + PCs )
summary( mod2ad4 )
##              Estimate Std. Error t value Pr(>|t|)    
## drb1151      0.039476   0.049596   0.796  0.42610    
## drb1151d    -0.069336   0.027912  -2.484  0.01302 *  
r2_2ad4 <- summary( mod2ad4 )$r.squared
# 0.09245018
( r2_2ad4 - r2_0 ) / ( 1 - r2_0 )
# 0.009297544

mod2ad5 <- lm( y ~ drb1131 + drb1131d + PCs )
summary( mod2ad5 )
##              Estimate Std. Error t value Pr(>|t|)    
## drb1131     -0.071653   0.113246  -0.633   0.5270    
## drb1131d    -0.011913   0.059190  -0.201   0.8405    
r2_2ad5 <- summary( mod2ad5 )$r.squared
# 0.08882465
( r2_2ad5 - r2_0 ) / ( 1 - r2_0 )
# 0.005339825

mod2ad6 <- lm( y ~ dqb151 + dqb151d + PCs )
summary( mod2ad6 )
##              Estimate Std. Error t value Pr(>|t|)    
## dqb151      -0.011053   0.043526  -0.254  0.79955    
## dqb151d     -0.033100   0.024294  -1.362  0.17312    
r2_2ad6 <- summary( mod2ad6 )$r.squared
# 0.09131131
( r2_2ad6 - r2_0 ) / ( 1 - r2_0 )
# 0.008054321
