# the hardest plot is left, calculating and plotting ORs for each quantile
# try quartiles, which is what Debo suggested given our sample sizes
# NOTES:
# - has a lot of hardcoded choices we might want to change in the future, like
#   - predictor is "PRS" (binned)
#   - response is "Type", and the non-default value is "SSNS" (because it's alphabetically after SRNS)
prs_quant <- function( data, n_groups = 4 ) {
    # first get quantiles
    breaks <- quantile( data$PRS, probs = ( 0 : n_groups ) / n_groups )
    # now "cut" data into bins
    data$bin <- cut( data$PRS, breaks, include.lowest = TRUE )

    # functions that reorganize data as needed
    quantiles <- paste0( 'Q', 1 : n_groups )
    # get unique quantile bins
    bins <- levels( data$bin )

    # this is prefered input/orientation
    x <- table( data$bin, data$Type )
    # this is the report
    obj <- oddsratio( x )
    
    # this gathers the data exactly as desired
    data_quant <- as_tibble( obj$measure, rownames = 'Bin' )
    # add a few more things
    data_quant$pvalue <- obj$p.value[, 'midp.exact']
    data_quant$Quantile <- quantiles
    # now just rename one column
    data_quant <- data_quant %>% rename( OR = estimate )
    # first row (ref) is trivial exclude from plot
    data_quant <- data_quant[ -1, ]
    
    # return the tibble!
    return( data_quant )
}
