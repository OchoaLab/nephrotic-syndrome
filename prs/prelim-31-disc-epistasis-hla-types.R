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

# try version without interactions
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

# try version with interactions
# here R2 always increases with interactions, so that's not meaningful unless it were a huge difference.  The per-interaction p-values are more useful, some are significant!
mod2i <- lm( y ~ (drb171 + dqa121 + dqb121 + dqb162 + dqa112)^2 + PCs )
summary( mod2i )
##                Estimate Std. Error t value Pr(>|t|)    
## drb171         0.176271   0.041270   4.271 1.98e-05 ***
## dqa121         0.053046   0.039757   1.334  0.18220    
## dqb121         0.027085   0.017312   1.565  0.11776    
## dqb162        -0.059795   0.035775  -1.671  0.09470 .  
## dqa112         0.022272   0.013979   1.593  0.11116    
## drb171:dqa121 -0.061309   0.019569  -3.133  0.00174 ** 
## drb171:dqb121  0.029067   0.145006   0.200  0.84113    
## drb171:dqb162 -0.181761   0.136062  -1.336  0.18166    
## drb171:dqa112  0.145683   0.120524   1.209  0.22683    
## dqa121:dqb121  0.192962   0.144125   1.339  0.18069    
## dqa121:dqb162  0.168704   0.133503   1.264  0.20641    
## dqa121:dqa112 -0.295140   0.119970  -2.460  0.01393 *  
## dqb121:dqb162 -0.001931   0.066701  -0.029  0.97691    
## dqb121:dqa112 -0.033537   0.053428  -0.628  0.53023    
## dqb162:dqa112 -0.005412   0.022441  -0.241  0.80945    
r2_2i <- summary( mod2i )$r.squared
# [1] 0.1607083
( r2_2i - r2_0 ) / ( 1 - r2_0 )
# [1] 0.08380966

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

# recessive first without interactions
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

# then with interactions; absolutely sucks (too many zeroes make fixed products)
# overall, worse than additive with interactions
mod2ir <- lm( y ~ (drb171r + dqa121r + dqb121r + dqb162r + dqa112r)^2 + PCs )
summary( mod2ir )
##                   Estimate Std. Error t value Pr(>|t|)    
## drb171r         -0.0290773  0.1248467  -0.233   0.8158    
## dqa121r          0.0109886  0.0883640   0.124   0.9010    
## dqb121r          0.0009001  0.0308411   0.029   0.9767    
## dqb162r          0.1058150  0.1019291   1.038   0.2993    
## dqa112r         -0.0048986  0.0149758  -0.327   0.7436    
## drb171r:dqa121r  0.0438427  0.0769533   0.570   0.5689    
## drb171r:dqb121r         NA         NA      NA       NA    
## drb171r:dqb162r         NA         NA      NA       NA    
## drb171r:dqa112r         NA         NA      NA       NA    
## dqa121r:dqb121r         NA         NA      NA       NA    
## dqa121r:dqb162r         NA         NA      NA       NA    
## dqa121r:dqa112r         NA         NA      NA       NA    
## dqb121r:dqb162r         NA         NA      NA       NA    
## dqb121r:dqa112r         NA         NA      NA       NA    
## dqb162r:dqa112r -0.0845877  0.0523572  -1.616   0.1063    
r2_2ir <- summary( mod2ir )$r.squared
# [1] 0.09021391
( r2_2ir - r2_0 ) / ( 1 - r2_0 )
# 0.006856373

# dominant, first without interactions
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

# and dominant with interactions
# the final r2 is a hair worse than additive with interactions
mod2id <- lm( y ~ (drb171d + dqa121d + dqb121d + dqb162d + dqa112d)^2 + PCs )
summary( mod2id )
##                   Estimate Std. Error t value Pr(>|t|)    
## drb171d          0.0677028  0.0294024   2.303  0.02135 *  
## dqa121d          0.0081392  0.0245896   0.331  0.74066    
## dqb121d          0.0133142  0.0098056   1.358  0.17459    
## dqb162d         -0.0455138  0.0522421  -0.871  0.38369    
## dqa112d          0.0088380  0.0086790   1.018  0.30858    
## drb171d:dqa121d  0.0056725  0.0177924   0.319  0.74988    
## drb171d:dqb121d  0.0005390  0.0366064   0.015  0.98825    
## drb171d:dqb162d -0.0484268  0.0341712  -1.417  0.15650    
## drb171d:dqa112d  0.0354606  0.0303546   1.168  0.24278    
## dqa121d:dqb121d  0.0529692  0.0363731   1.456  0.14539    
## dqa121d:dqb162d  0.0473680  0.0335377   1.412  0.15791    
## dqa121d:dqa112d -0.0736655  0.0302447  -2.436  0.01490 *  
## dqb121d:dqb162d  0.0006349  0.0167374   0.038  0.96974    
## dqb121d:dqa112d -0.0079652  0.0137712  -0.578  0.56303    
## dqb162d:dqa112d  0.0052617  0.0262111   0.201  0.84091    
r2_2id <- summary( mod2id )$r.squared
# 0.1602623
( r2_2id - r2_0 ) / ( 1 - r2_0 )
# 0.08332285


########################### CLEANER
# restarted with fewer types, since some are not even conditionally significant on the others (without interactions)

# additive without interactions, conditionally significant cases only
mod2 <- lm( y ~ drb171 + dqb121 + dqb162 + PCs )
summary( mod2 )
##              Estimate Std. Error t value Pr(>|t|)    
## drb171       0.152069   0.010608  14.336  < 2e-16 ***
## dqb121       0.075102   0.013614   5.516 3.66e-08 ***
## dqb162      -0.064607   0.012485  -5.175 2.38e-07 ***
r2_2 <- summary( mod2 )$r.squared
# 0.1398214
( r2_2 - r2_0 ) / ( 1 - r2_0 )
# 0.06100909

# additive with interactions; wow things are very significant and huge effect sizes just with these
mod2i <- lm( y ~ (drb171 + dqb121 + dqb162)^2 + PCs )
summary( mod2i )
##                Estimate Std. Error t value Pr(>|t|)    
## drb171         0.129845   0.011484  11.306  < 2e-16 ***
## dqb121         0.019359   0.016186   1.196  0.23175    
## dqb162        -0.055169   0.013265  -4.159 3.26e-05 ***
## drb171:dqb121  0.256939   0.031702   8.105 6.76e-16 ***
## drb171:dqb162 -0.107657   0.040744  -2.642  0.00826 ** 
## dqb121:dqb162 -0.027976   0.044878  -0.623  0.53307    
r2_2i <- summary( mod2i )$r.squared
# 0.1550253
( r2_2i - r2_0 ) / ( 1 - r2_0 )
# 0.07760601

# skip recessive version (nothing was significant)

# dominant without interactions; this one is better than additive without interactions
# interesting that the two same dqa1 types drop in both additive and dominant!
mod2d <- lm( y ~ drb171d + dqb121d + dqb162d + PCs )
summary( mod2d )
##              Estimate Std. Error t value Pr(>|t|)    
## drb171d      0.089253   0.006026  14.812  < 2e-16 ***
## dqb121d      0.040612   0.007271   5.585 2.47e-08 ***
## dqb162d     -0.037461   0.007257  -5.162 2.55e-07 ***
r2_2d <- summary( mod2d )$r.squared
# [1] 0.1431986
( r2_2d - r2_0 ) / ( 1 - r2_0 )
# [1] 0.06469565

# dominant with interactions
mod2id <- lm( y ~ (drb171d + dqb121d + dqb162d)^2 + PCs )
summary( mod2id )
##                  Estimate Std. Error t value Pr(>|t|)    
## drb171d          0.077588   0.006803  11.406  < 2e-16 ***
## dqb121d          0.011648   0.009057   1.286  0.19851    
## dqb162d         -0.029908   0.008034  -3.723  0.00020 ***
## drb171d:dqb121d  0.057721   0.008346   6.916 5.30e-12 ***
## drb171d:dqb162d -0.031064   0.010535  -2.949  0.00321 ** 
## dqb121d:dqb162d -0.005785   0.011579  -0.500  0.61738    
r2_2id <- summary( mod2id )$r.squared
# 0.155557
( r2_2id - r2_0 ) / ( 1 - r2_0 )
# 0.07818644

# conclusions: dominant with interactions is better than additive with interactions after we cut out some noisy types that probably inflate the measure.  In this case interactions are also individually significant
# I had neglected to notice before that dqb162 is protective (both here and in paper), which is why interactions with it are also negative, whereas the two risk types have a positive interaction estimate
