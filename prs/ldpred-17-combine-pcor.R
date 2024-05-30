# combines odds ratios using a crude form of meta analysis

library(genio)
library(readr)
library(testthat)

#library(metafor)

# constants
# needed to replicate CI calculations
# includes intercept (1) plus 10 PCs
df_null <- 11
# process all of these
models <- c( 'grid-h0.1-best', 'lassosum-best', 'ct-best', 'ct-stacked', 'inf-best', 'auto-h0.1' )
# combine the first two, into the last one
tests <- c( 'test', 'test-curegn', 'test-bristol-curegn' )
# consider only this combination of base and training datasets
type_in <- 'base-train'

# sample sizes are essentially fixed in each of the two datasets, just look at the fam tables!
df1 <- count_lines( paste0( tests[1], '/mac20.fam' ) ) - df_null
df2 <- count_lines( paste0( tests[2], '/mac20.fam' ) ) - df_null
# the output combines degrees of freedom (yes, df_null needs to be subtracted twice from the combined sample sizes, because they were fit separately, which is inefficient but oh well, that's what we have)
df3 <- df1 + df2

# separate function for CIs only, given r and df, which is more amenable to combining the way I want to
# (code pulled from bigsnpr::pcor)
# also agrees with these formulas:
# https://en.wikipedia.org/wiki/Fisher_transformation
# https://en.wikipedia.org/wiki/Partial_correlation
pcor_CI <- function( r, df, alpha = 0.05 ) {
    # some transformation of this mean to sort of normalize?
    z <- (log(1 + r) - log(1 - r))/2
    # a kind of CI in a normal scale, which maybe hides a sense of variance?
    rad <- stats::qnorm(alpha/2, lower.tail = FALSE)/sqrt(df - 2)
    # CIs in the final, desired scale
    return( tanh( z + c( -rad, rad ) ) )
}

# read the two datasets
for ( model in models ) {
    message( model )
    # the part shared across datasets
    file <- paste0( '/cor-', type_in, '-ldpred2-', model, '.txt.gz' )
    cors1 <- as.numeric( read_lines( paste0( tests[1], file ) ) )
    cors2 <- as.numeric( read_lines( paste0( tests[2], file ) ) )
    # this just validates that we perfectly reproduce CIs here
    expect_equal( pcor_CI( cors1[1], df1 ), cors1[2:3] )
    expect_equal( pcor_CI( cors2[1], df2 ), cors2[2:3] )
    # now do the combining, using a weighted average (agrees with variance weighting)
    # this doesn't quite agree with metafor's estimates, but it's close enough (there is an additional "tau" correction)
    r3 <- ( cors1[1] * df1 + cors2[1] * df2 ) / df3
    # now get CIs using combined sample sizes
    ci3 <- pcor_CI( r3, df3 )
    # write that to output!
    cors3 <- c( r3, ci3 )
    write_lines( cors3, paste0( tests[3], file ) )
}


## # the input data is just these two things:
## # mean estimates of correlation
## ri <- c(0.233847540095574, 0.264317756196715)
## # cor-base-train-ldpred2-grid-h0.1-best.txt.gz
## ## 0.233847540095574
## ## 0.149815795826082
## ## 0.314525134818934
## #
## ## 0.264317756196715
## ## 0.171754021136691
## ## 0.352255285413333
## # and sample sizes
## # 'prs-base-train-ldpred2-grid-h0.1-best.txt.gz'
## ni <- c(517, 419)

## # this calculates variances given a correlation model
## # (there's also UCOR and ZCOR, each with non-default option `vtype = "UB"`)
## x <- escalc( 'COR', ri = ri, ni = ni - df_null  )
## ##       yi     vi 
## ## 1 0.2338 0.0018 
## ## 2 0.2643 0.0021 

## # this actually performs meta-analysis
## # default was random effects, this is fixed version (in example, the two versions gave the same mean and CIs)
## obj <- rma.uni( yi = x$yi, vi = x$vi, method = 'FE' )
## ## estimate      se    zval    pval   ci.lb   ci.ub      
## ##   0.2477  0.0311  7.9707  <.0001  0.1868  0.3086  *** 

## obj$beta
## obj$ci.lb
## obj$ci.ub

## # attempts to manually construct combined mean estimate
## # various weighted and unweighted versions, direct and transformed
## df <- ni - df_null
## df1 <- df[1]
## df2 <- df[2]
## df3 <- df1+df2

## z_xf <- function(r) (log(1 + r) - log(1 - r))/2

## #( ri[1] + ri[2] ) /2                                                                          # 0.2490826
## ( ri[1] * df1 + ri[2] * df2 ) / df3                                                           # 0.2474491 *
## ( ri[1] * (df1 - 2) + ri[2] * ( df2 - 2 ) ) / ( df3 - 4 )                                     # 0.2474419
## #( ri[1] * sqrt(df1) + ri[2] * sqrt(df2) ) / (sqrt(df1) + sqrt(df2) )                          # 0.2482635
## #tanh( ( z_xf( ri[1] ) + z_xf( ri[2] ) ) / 2 )                                                 # 0.2491443
## tanh( ( z_xf( ri[1] ) * df1 + z_xf( ri[2] ) * df2 ) / df3 )                                   # 0.24751 *
## tanh( ( z_xf( ri[1] ) * (df1-2) + z_xf( ri[2] ) * (df2-2) ) / (df3-4) )                       # 0.2475028
## #tanh( ( z_xf( ri[1] ) * sqrt(df1) + z_xf( ri[2] ) * sqrt(df2) ) / ( sqrt(df1) + sqrt(df2) ) ) # 0.2483249


## # comparison for a single value
## # the CIs here are narrower than what we actually have
## # it's because pcor corrects for PCs, but here that's not used or taken into account
## # PCs reduce the effective sample size

## # COR
## ##   0.2338  0.0416  5.6193  <.0001  0.1523  0.3154  ***
## # UCOR
## ##   0.2341  0.0416  5.6251  <.0001  0.1525  0.3156  *** 
## # ZCOR
## ##   0.2383  0.0441  5.4016  <.0001  0.1518  0.3247  *** 

## x <- escalc( 'COR', ri = ri[1], ni = ni[1]-30  )
## obj <- rma.uni( yi = x$yi, vi = x$vi, method = 'FE' )
## obj
## ##   0.2338  0.0420  5.5646  <.0001  0.1515  0.3162  *** # -10
## ##   0.2338  0.0424  5.5093  <.0001  0.1507  0.3170  *** # -20
## ##   0.2338  0.0429  5.4535  <.0001  0.1498  0.3179  *** # -30

## x <- escalc( 'COR', ri = ri[2], ni = ni[2]-30  )
## obj <- rma.uni( yi = x$yi, vi = x$vi, method = 'FE' )
## obj
## ##   0.2643  0.0472  5.5975  <.0001  0.1718  0.3569  *** # -30


## # NOTES
## # this is what bigsnpr::pcor actually does, simplified to exclude checks and special cases
## # makes me wonder if the CIs are right, have to think about it...
## pcor <- function (x, y, z, alpha = 0.05) {
##     m <- stats::model.matrix(~., data = cbind.data.frame(x, y, z))
##     mod1 <- stats::lm.fit(x = m[, -(2:3), drop = FALSE], y = m[, 2])
##     mod2 <- stats::lm.fit(x = m[, -(2:3), drop = FALSE], y = m[, 3])
##     # so the mean correlation is just a partial correlation (after removing the effect of PCs (z) on both x and y)
##     r <- stats::cor(mod1$residuals, mod2$residuals)
##     # some transformation of this mean to sort of normalize?
##     z <- (log(1 + r) - log(1 - r))/2
##     # a kind of CI in a normal scale, which maybe hides a sense of variance?
##     rad <- stats::qnorm(alpha/2, lower.tail = FALSE)/sqrt(mod2$df.residual - 2)
##     # mean and CIs in the final, desired scale
##     c(r, tanh(z - rad), tanh(z + rad))
## }

# Idea: recover the internal variance of this model, provide that to rma.uni instead?

