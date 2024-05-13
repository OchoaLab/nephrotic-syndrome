# the hardest plot is left, calculating and plotting ORs for each quantile
# the data just has to be made the hard way
# try quartiles, which is what Debo suggested given our sample sizes
# NOTES:
# - has a lot of hardcoded choices we might want to change in the future, like
#   - predictor is "PRS"
#   - response is "Type", and the non-default value is "SSNS"
prs_quant <- function( data, n_groups = 4 ) {
    # first get quantiles
    breaks <- quantile( data$PRS, probs = ( 0 : n_groups ) / n_groups )
    # now "cut" data
    data$Quantile <- cut( data$PRS, breaks, include.lowest = TRUE )

    # functions that reorganize data as needed
    # get unique quantiles
    quants <- levels( data$Quantile )
    # make parallel data
    ORs <- vector('numeric', length( quants ) )
    Ls <- vector('numeric', length( quants ) )
    Us <- vector('numeric', length( quants ) )
    Ps <- vector('numeric', length( quants ) )
    for ( quant in quants ) {
        i <- which( quant == quants )
        # this produces the input table we desire!
        x <- table( data$Type, data$Quantile == quant )
        # this is the report
        obj <- oddsratio.wald( x )
        # store the data of interest
        ORs[i] <- obj$measure['SSNS', 'estimate']
        Ls[i] <- obj$measure['SSNS', 'lower']
        Us[i] <- obj$measure['SSNS', 'upper']
        Ps[i] <- obj$p.value['SSNS', 'fisher.exact']
    }

    # gather into tibble
    data_quant <- tibble(
        Quantile = factor( quants, levels = quants ),
        OR = ORs,
        lower = Ls,
        upper = Us,
        pvalue = Ps
    )
    # return the tibble!
    return( data_quant )
}
