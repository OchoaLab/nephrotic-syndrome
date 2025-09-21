# uses top HLA haplotype to demonstrate that recessive model is better, and by how much?

# need to calculate R2 as in the paper!  Need PCs, etc
# based on ldpred-02-score.R

# constants
name_data <- 'mac20'
test <- 'test'
# ordered as in table, not as in paper (DRB1 listed first in paper, last here)
haplo <- '02:01~02:02~07:01'

# all processing happens in subdirectory
setwd( test )

# load haplotypes
data <- read.table('Bristol_hla_haplotypes_by_Subject.txt', header = TRUE)
colnames( data ) <- c('id', 'hap1', 'hap2')
# process to get dosage of top haplotype
data$dosage <- (data$hap1 == haplo) + (data$hap2 == haplo)

# load phenotype
fam <- read.table( paste0( name_data, '.fam' ) )
colnames( fam ) <- c('fam', 'id', 'pat', 'mat', 'sex', 'pheno')
y <- fam$pheno
y[ y == 0 ] <- NA # in plink format, zeros are missing, translate appropriately here!
# and PCs
PCs <- read.table( paste0( name_data, '.eigenvec' ) )
PCs <- PCs[ , -(1:2) ]
PCs <- as.matrix( PCs )
colnames( PCs ) <- NULL

# get dosage aligned with fam/y/PCs
indexes <- match( fam$id, data$id )
data <- data[ indexes, ]
stopifnot( all( data$id == fam$id ) )
dosage <- data$dosage

# calculate fit in linear model
mod0 <- lm( y ~ PCs )
r2_0 <- summary( mod0 )$r.squared
# [1] 0.03942191

# try additive dosage
# worse than R2 of top loci, but meh
mod2 <- lm( y ~ dosage + PCs )
summary( mod2 )
##             Estimate Std. Error t value Pr(>|t|)    
## dosage       0.15155    0.04106   3.691 0.000248 ***
r2_2 <- summary( mod2 )$r.squared
# [1] 0.06465434
( r2_2 - r2_0 ) / ( 1 - r2_0 )
# [1] 0.02626797

# now try the same with recessive encodings!
# note counted allele is minor (though here it's all multiallelic, all are minor?)
mean( dosage ) / 2
# [1] 0.1779497
# now encode recessive model
Xr <- dosage
Xr[ Xr == 1 ] <- 0
# for completeness do dominant one as well, confirm that it is worse
Xd <- dosage
Xd[ Xd == 1 ] <- 2

# recessive first.  It's a worse model than additive :(
mod2r <- lm( y ~ Xr + PCs )
summary( mod2r )
##             Estimate Std. Error t value Pr(>|t|)    
## Xr           0.14667    0.10547   1.391   0.1650    
r2_2r <- summary( mod2r )$r.squared
# [1] 0.04308611
( r2_2r - r2_0 ) / ( 1 - r2_0 )
# [1] 0.003814571

# dominant, better than recessive but worse than the original additive one
mod2d <- lm( y ~ Xd + PCs )
summary( mod2d )
##             Estimate Std. Error t value Pr(>|t|)    
## Xd           0.07687    0.02152   3.572 0.000388 ***
r2_2d <- summary( mod2d )$r.squared
# [1] 0.06309748
( r2_2d - r2_0 ) / ( 1 - r2_0 )
# [1] 0.02464721
