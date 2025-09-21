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

pcor_meta <- function( r1, df1, r2, df2 ) {
    # the output combines degrees of freedom (yes, df_null needs to be subtracted twice from the combined sample sizes, because they were fit separately, which is inefficient but oh well, that's what we have)
    df3 <- df1 + df2
    # now do the combining, using a weighted average (agrees with variance weighting)
    # this doesn't quite agree with metafor's estimates, but it's close enough (there is an additional "tau" correction)
    r3 <- ( r1 * df1 + r2 * df2 ) / df3
    # now get CIs using combined sample sizes
    ci3 <- pcor_CI( r3, df3 )
    # this is the return value we want!
    cors3 <- c( r3, ci3 )
    return( cors3 )
}

