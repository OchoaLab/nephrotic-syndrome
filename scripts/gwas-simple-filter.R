library(readr)
library(qqman)
library(ochoalabtools)
library(tibble)
library(dplyr)

# reads some noisy gwas data (from plink pca logistic), apply a simple cleanup filter

# for now let's just stick with the standard threshold
p_cut_gwas <- 5e-8 # significant
p_cut_sugg <- 1e-5 # suggestive

# switch to data directory
setwd( '../data/' )

# file to read
file_in <- 'pca_combine3_maf10.PHENO1.glm.logistic'
file_out_pvals <- paste0( file_in, '-pval-hist' )
file_out_np <- paste0( file_in, '-pval-neighbors' )
file_out_nd <- paste0( file_in, '-pval-neighbors-dist' )
file_out_manhattan_all <- paste0( file_in, '-manhattan-all' )
file_out_manhattan_clean <- paste0( file_in, '-manhattan-clean' )

data <- read_tsv( file_in )
# fix first name, starts with a comment symbol
# other renames are for coherence with qqman::manhattan
names(data)[ names(data) == '#CHROM' ] <- 'CHR'
names(data)[ names(data) == 'POS' ] <- 'BP'
names(data)[ names(data) == 'ID' ] <- 'SNP'
# "P" is the only other column that matters, and it's good

# remove some useless columns to simplify inspection (doesn't matter otherwise)
data$A1 <- NULL
data$TEST <- NULL

# there are only a few NAs, but for the logic of our script to make sense, it's best to remove them immediately
## sum( is.na( data$P ))
## [1] 14
data <- data[ !is.na( data$P ), ]

# visualize p-values
fig_start( file_out_pvals )
hist(
    data$P,
    freq = FALSE,
    xlab = 'p-value',
    main = ''
)
abline( h = 1, lty = 2, col = 'red' )
fig_end()
# I see an overall skew, so it's not just very small p-values that are a problem

# list of all significant cases
indexes_signif <- which( data$P < p_cut_gwas )
# there are fewer than I remembered (I guess we've filtered a lot to get here)
m_sig <- length( indexes_signif )
# [1] 567

# overall patterns:
# chromosomes:
## table( data$CHR[ indexes_signif ] )
##  1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20 
## 26  23  49  19  22 290  13  12  12  18   9  13   9   9   5   9  10   6   5   6 
## 22 
##  2 

# let's gather some statistics for these cases
neigh_dist <- vector( 'numeric', m_sig )
neigh_p <- vector( 'numeric', m_sig )
# j is index in table of significant loci
for ( j in 1 : m_sig ) {
    # row index in big table
    i <- indexes_signif[ j ]
    # extract row
    # also get previous and next rows
    data_i <- data[ i + (-1):1, ]
    
    # it's extremely unlikely that snp is at the edges of a chromosome, but let's check (code doesn't handle that case)
    if ( length( unique( data_i$CHR ) ) != 1 )
        stop( 'Got a SNP on an edge!' )

    # pick the most significant neighbor
    # first is TRUE if its p-value was smaller than last, FALSE otherwise
    first <- data_i$P[1] < data_i$P[3]

    # extract that row ("best neighbor") and save its stats
    data_neigh <- data_i[ if (first) 1 else 3, ]
    neigh_dist[ j ] <- abs( data_neigh$BP - data_i$BP[2] )
    neigh_p[ j ] <- data_neigh$P
}

# gather vectors into data frame
# copy subset first, to inherit some of the other columns
data_new <- data[ indexes_signif, ]
# add new columns
data_new$nd <- neigh_dist
data_new$np <- neigh_p
# remove other less useful columns for clarity
#data_new$`LOG(OR)_SE` <- NULL

# this plot shows the two lines, one of strong p-value correlation as expected for neighbors in LD, and the other with no correlation
# the second line is well separated by a p < 1e-5 requirement on both neighbors
fig_start( file_out_np )
plot(
    data_new$P,
    data_new$np,
    pch = '.',
    log = 'xy',
    xlab = 'p-value (GWAS significant)',
    ylab = 'p-value best neighbor'
)
abline( 0, 1, lty = 2, col = 'gray' )
abline( h = p_cut_sugg, lty = 2, col = 'gray', untf = TRUE )
abline( v = p_cut_sugg, lty = 2, col = 'gray', untf = TRUE )
fig_end()

# how does distance fit in?
fig_start( file_out_nd )
plot(
    data_new$nd,
    data_new$np,
    xlab = 'Distance best neighbor',
    ylab = 'p-value best neighbor',
    pch = '.',
    log = 'y'
)
abline( h = p_cut_sugg, lty = 2, col = 'gray', untf = TRUE )
fig_end()
# looks like the expected inverse correlation
# Most of the strong neighbors are close, < 25KB
# Most of the bad neighbors are about as close, though maybe less so on average, but a few are much farther, up to 150KB away, which may be false negatives after my filter but it's hard to tell.
# these are still closer than other data (a more typical 1 MB distance threshold)

## # out of curiosity, what if we did this for the entire dataset?
## # the log scale distorts things in favor of strong hits (good), so poor correlation between bad p-values won't stand out
## # let's gather some statistics for these cases
## m <- nrow( data )
## data$np <- NA
## # handle two edge cases (only one neighbor
## data$np[1] <- data$P[2]
## data$np[m] <- data$P[m-1]
## # navigate all loci (not just significant cases as before)
## for ( i in 2 : (m-1) ) {
##     # handle chromosome changes
##     # assumes no chromosome has a single row
##     chr_i <- data$CHR[ i ] # current chromosome
##     if ( data$CHR[ i - 1 ] != chr_i ) {
##         # only acceptable neighbor is next row
##         data$np[ i ] <- data$P[ i + 1 ]
##     } else if ( data$CHR[ i + 1 ] != chr_i ) {
##         # only acceptable neighbor is previous row
##         data$np[ i ] <- data$P[ i - 1 ]
##     } else {
##         # have to compare both rows, pick minimum
##         data$np[ i ] <- min( data$P[ c( i - 1, i + 1 ) ] )
##     }
## }

## # repeat neighbor p-value plot
## plot( data$P, data$np, log='xy', pch = '.' )
## abline( 0, 1, lty = 2, col = 'gray' )
## abline( h = p_cut_sugg, lty = 2, col = 'gray', untf = TRUE )
## abline( v = p_cut_sugg, lty = 2, col = 'gray', untf = TRUE )
## # that was surprisingly boring, I'm sad I wasted time making that

# now it's time to flag the SNPs to be removed on the basis of their direct neighbors (sure to be in LD) not having good p-values)
# start from scratch in case I don't want to regenerate all the previous work
data$keep <- TRUE # most things are kept
# navigate significant loci only
# j is index in table of significant loci
for ( j in 1 : m_sig ) {
    # row index in big table
    i <- indexes_signif[ j ]
    # extract row
    # also get previous and next rows
    data_i <- data[ i + (-1):1, ]
    
    # it's extremely unlikely that snp is at the edges of a chromosome, but let's check (code doesn't handle that case)
    if ( length( unique( data_i$CHR ) ) != 1 )
        stop( 'Got a SNP on an edge!' )

    # pick the most significant neighbor
    # first is TRUE if its p-value was smaller than last, FALSE otherwise
    p2 <- min( data_i$P[ c(1, 3) ] )
    if ( p2 > p_cut_sugg ) {
        # neighbors were too insignificant, flag to remove
        data$keep[ i ] <- FALSE
    }
}

# these many are getting removed now
# sum( !data$keep )
# [1] 288

# create filtered data for plotting
data_clean <- data[ data$keep, ]

# for plots, it's also helpful to visualize previous data with a cap so it's in the same scale as the cleaned data
p_cap <- min( data_clean$P )

# another copy of original data
data_capped <- data
# apply cap to all p-values (by construction, only bad cases are affected)
data_capped$P[ data_capped$P < p_cap ] <- p_cap

fig_start(
    file_out_manhattan_all,
    width = 12
)
manhattan(
    data_capped,
    cex.axis = 0.5
)
fig_end()

fig_start(
    file_out_manhattan_clean,
    width = 12
)
manhattan(
    data_clean,
    cex.axis = 0.5
)
fig_end()

# zoom into interesting chromosomes (clean only)
## table( data_clean$CHR[ data_clean$P < p_cut_gwas ] )
##  3   6 
## 24 255 
chrs_hits <- unique( data_clean$CHR[ data_clean$P < p_cut_gwas ] )

for ( chr in chrs_hits ) {
    fig_start(
        paste0( file_out_manhattan_clean, '_chr', chr ),
        width = 12
    )
    manhattan(
        data_clean[ data_clean$CHR == chr, ]
    )
    fig_end()
}

# a second pass reports position ranges per peak (easy cause it's one per chr)
range_all <- NULL
for ( chr in chrs_hits ) {
    # this chromosome only
    data_clean_chr <- data_clean[ data_clean$CHR == chr, ]
    # and significant loci only
    data_clean_chr_sig <- data_clean_chr[ data_clean_chr$P < p_cut_gwas, ]
    # report range now
    range_chr <- range( data_clean_chr_sig$BP )
    # turn into tibble
    range_chr <- tibble(
        chr = chr,
        pos_min = range_chr[1],
        pos_max = range_chr[2],
    )
    # add to complete data structure
    range_all <- bind_rows(
        range_all,
        range_chr
    )
}

## range_all
##     chr  pos_min  pos_max
## 1     3 86793259 86917711 # 86.8 86.9 Mb
## 2     6 29938571 32681568 # 29.9 32.7 Mb
