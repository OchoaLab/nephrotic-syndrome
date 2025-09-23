# uses top haplotype to demonstrate that recessive model is better, and by how much?
# this version uses the discovery data for greater power

# need to calculate R2 as in the paper!  Need PCs, etc
# based on ldpred-02-score.R

library(genio)
library(readr)
library(dplyr)

# constants
name_data <- 'mac20'
# ordered as in table, not as in paper (DRB1 listed first in paper, last here)
haplo <- '02:01~02:02~07:01'

# we need the data at the base of this dir
setwd( '/datacommons/ochoalab/ssns_gwas/imputed' )

# load haplotypes
haps <- read_tsv( 'haplotype_HapsbySubject.txt', show_col_types = FALSE )
haps <- haps %>% rename( id = SAMPLE.ID, hap1 = "DQA1~DQB1~DRB1.Hap1", hap2 = "DQA1~DQB1~DRB1.Hap2" )
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
# process to get dosage of top haplotype
data$dosage <- (data$hap1 == haplo) + (data$hap2 == haplo)
dosage <- data$dosage

# calculate fit in linear model
mod0 <- lm( y ~ PCs )
r2_0 <- summary( mod0 )$r.squared
# [1] 0.08393301

# try additive dosage
mod2 <- lm( y ~ dosage + PCs )
summary( mod2 )
##              Estimate Std. Error t value Pr(>|t|)    
## dosage       0.199002   0.011781  16.892  < 2e-16 ***
r2_2 <- summary( mod2 )$r.squared
# [1] 0.1390261
( r2_2 - r2_0 ) / ( 1 - r2_0 )
# [1] 0.06014085
coef( summary( mod2 ) )[ 'dosage', 'Pr(>|t|)' ]
# [1] 4.262415e-62

# now try the same with recessive encodings!
# note counted allele is minor (though here it's all multiallelic, all are minor?)
mean( dosage ) / 2
# [1] 0.1115942
# now encode recessive model
Xr <- dosage
Xr[ Xr == 1 ] <- 0
# for completeness do dominant one as well, confirm that it is worse
Xd <- dosage
Xd[ Xd == 1 ] <- 2

# recessive first.  It's a worse model than additive :(
mod2r <- lm( y ~ Xr + PCs )
summary( mod2r )
##              Estimate Std. Error t value Pr(>|t|)    
## Xr           0.096414   0.024209   3.983 6.92e-05 ***
r2_2r <- summary( mod2r )$r.squared
# [1] 0.08718001
( r2_2r - r2_0 ) / ( 1 - r2_0 )
# [1] 0.003544505
coef( summary( mod2r ) )[ 'Xr', 'Pr(>|t|)' ]
# [1] 6.924936e-05

# dominant, finally something better than both additive and recessive! though difference is small
mod2d <- lm( y ~ Xd + PCs )
summary( mod2d )
##              Estimate Std. Error t value Pr(>|t|)    
## Xd           0.110932   0.006408  17.311  < 2e-16 ***
r2_2d <- summary( mod2d )$r.squared
# [1] 0.1416196
( r2_2d - r2_0 ) / ( 1 - r2_0 )
# [1] 0.06297198
coef( summary( mod2d ) )[ 'Xd', 'Pr(>|t|)' ]
# [1] 4.996445e-65

# here we may as well have a conditional test for dominant over additive, to get its p-value
# wow, dominant is clearly the stronger signal here!  There's practically no gain to adding both (i.e. making the dominance deviation more flexible)
mod2ad <- lm( y ~ dosage + Xd + PCs )
summary( mod2ad )
##              Estimate Std. Error t value Pr(>|t|)    
## dosage       0.026357   0.047986   0.549 0.582854    
## Xd           0.097012   0.026141   3.711 0.000209 ***
r2_2ad <- summary( mod2ad )$r.squared
# [1] 0.1416776
( r2_2ad - r2_0 ) / ( 1 - r2_0 )
# [1] 0.06303539
# dom|add
( r2_2ad - r2_2 ) / ( 1 - r2_2 )
# [1] 0.003079752
# add|dom
( r2_2ad - r2_2d ) / ( 1 - r2_2d )
# [1] 6.766886e-05
