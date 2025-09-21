# uses top HLA type (individually, not as haplotype) to demonstrate that recessive model is better, and by how much?

# need to calculate R2 as in the paper!  Need PCs, etc
# based on ldpred-02-score.R

# constants
name_data <- 'mac20'
test <- 'test'
# types of interest
## HLA_DRB1*07:01
## HLA_DQA1*02:01
## HLA_DQB1*02:01
## HLA_DQB1*06:02
## HLA_DQA1*01:02

# all processing happens in subdirectory
setwd( test )

# load haplotypes
data <- read.table('bristol.hla.7loci.calls.4haplo.txt', header = TRUE)
# process to get dosage of top types; this was a very manual process because of the way it's organized
data$drb171 <- (data$DRB1 == '07:01') + (data$DRB1.1 == '07:01')
data$dqa121 <- (data$DQA1 == '02:01') + (data$DQA1.1 == '02:01')
data$dqa112 <- (data$DQA1 == '01:02') + (data$DQA1.1 == '01:02')
data$dqb121 <- (data$DQB1 == '02:01') + (data$DQB1.1 == '02:01')
data$dqb162 <- (data$DQB1 == '06:02') + (data$DQB1.1 == '06:02')

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

# get dosages aligned with fam/y/PCs
names( data )[ names( data ) == 'SubjectID' ] <- 'id'
indexes <- match( fam$id, data$id )
data <- data[ indexes, ]
stopifnot( all( data$id == fam$id ) )
drb171 <- data$drb171
dqa121 <- data$dqa121
dqa112 <- data$dqa112
dqb121 <- data$dqb121
dqb162 <- data$dqb162

# calculate fit in linear model
mod0 <- lm( y ~ PCs )
r2_0 <- summary( mod0 )$r.squared
# [1] 0.03942191

# try additive dosage
# reorder variables by their paper p-value
mod2 <- lm( y ~ drb171 + dqa121 + dqb121 + dqb162 + dqa112 + PCs )
summary( mod2 )
##             Estimate Std. Error t value Pr(>|t|)    
## drb171      -0.68649    0.27063  -2.537  0.01149 *  
## dqa121       0.81313    0.27148   2.995  0.00288 ** 
## dqb121       0.06082    0.03834   1.586  0.11329    
## dqb162       0.06068    0.09086   0.668  0.50455    
## dqa112      -0.13861    0.07147  -1.939  0.05301 .  
r2_2 <- summary( mod2 )$r.squared
# [1] 0.09283914
( r2_2 - r2_0 ) / ( 1 - r2_0 )
# [1] 0.05560946

# now try the same with recessive encodings!
# all counted alleles are minor (though here it's all multiallelic, all are minor?); here the meaning is not symmetric, so we want to use counted allele regardless of freq
mean( drb171 )/2 # [1] 0.212766
mean( dqa121 )/2 # [1] 0.2098646
mean( dqa112 )/2 # [1] 0.09574468
mean( dqb121 )/2 # [1] 0.2040619
mean( dqb162 )/2 # [1] 0.05222437

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
##             Estimate Std. Error t value Pr(>|t|)    
## drb171r      0.08310    0.08413   0.988   0.3238    
## dqa121r           NA         NA      NA       NA    
## dqb121r      0.03832    0.06331   0.605   0.5453    
## dqb162r      0.31846    0.15151   2.102   0.0361 *  
## dqa112r     -0.19356    0.10009  -1.934   0.0537 .  
r2_2r <- summary( mod2r )$r.squared
# [1] 0.05150205
( r2_2r - r2_0 ) / ( 1 - r2_0 )
# [1] 0.01257591

# dominant, better than recessive but worse than the original additive one
mod2d <- lm( y ~ drb171d + dqa121d + dqb121d + dqb162d + dqa112d + PCs )
summary( mod2d )
##              Estimate Std. Error t value Pr(>|t|)    
## drb171d     -0.346380   0.135465  -2.557  0.01085 *  
## dqa121d      0.409097   0.135989   3.008  0.00276 ** 
## dqb121d      0.031139   0.021283   1.463  0.14407    
## dqb162d     -0.001223   0.049707  -0.025  0.98038    
## dqa112d     -0.065600   0.040054  -1.638  0.10210    
r2_2d <- summary( mod2d )$r.squared
# [1] 0.0922656
( r2_2d - r2_0 ) / ( 1 - r2_0 )
# [1] 0.05501238
