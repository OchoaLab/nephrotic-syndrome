library(genio)
library(readr)

# decide how many individuals are in each subset
# only used for LDpred2-grid and lassosum2

set.seed(1)

# get number of individuals
n_ind <- count_lines( 'data.fam' )

# randomly select 70% to be training
ind_train <- sample( n_ind, n_ind * 0.7 )

# the rest will be testint subset
ind_test <- setdiff( 1L : n_ind, ind_train )

# NEWEST:
length( ind_train ) # [1] 408
length( ind_test )  # [1] 176

# save values for later loading
write_lines( ind_train, 'ind-train.txt.gz' )
write_lines( ind_test, 'ind-test.txt.gz' )
