# test whether dominance is present in top hits
# look at most powerful data for this question, ssns-vs-control

library(tidyverse)
library(BEDMatrix)
library(ochoalabtools)

# go where the data is
setwd( '/datacommons/ochoalab/ssns_gwas/imputed' )

# load phenotype and other covariates
data <- read_tsv( 'patient-data.txt.gz', show_col_types = FALSE )

# load genotypes
X <- BEDMatrix( 'ssns_ctrl/mac20', simple_names = TRUE )

# align phenotype data
indexes <- match( rownames( X ), data$id )
data <- data[ indexes, ]
stopifnot( all( rownames( X ) == data$id ) )

# now go to PRS data
setwd( 'prs-new/train/' )

# load C+T report containing top independent SNPs to consider
# NOTES:
# - this only has HLA SNPs, we should also look at clec16a
# - SNPs appear in pos order, not ranked by p-value
top_snps <- read_tsv( 'report-nonzero-betas-base-ldpred2-ct-best.txt.gz', show_col_types = FALSE )
ids <- top_snps$rsid

# extract these variants now!
X2 <- X[ , ids ]

# save key data to revisit more quickly next time
save( X2, data, top_snps, file = 'dominance.RData' )
#load( 'dominance.RData' )

# the response is the same for all SNPs
y <- data$ssns_ctrl

plot_geno_case_ctrl <- function( x, y, name_out ) {
    # look at counts
    c <- table( x, y )
    n <- rowSums( c )
    # perform calculation of point estimates and CIs
    obj <- binconf( x = c[, '1'], n = n, method = 'exact' )
    
    # visualize!
    data2 <- tibble(
        x = as.numeric( names( n ) ),
        p = obj[,'PointEst'],
        lower = obj[,'Lower'],
        upper = obj[,'Upper']
    )

    # for info, fit a linear and logistic regressions
    obj_lm <- lm( y ~ x )
    coef_lm <- coef( obj_lm )
    obj_glm <- glm( y ~ x, family = binomial() )
    coef_glm <- coef( obj_glm )
    # to get nice curve, predict this way
    xp <- ( 0 : 100 ) / 100 * 3 - 0.5
    yp <- predict( obj_glm, tibble( x = xp ), type = 'response' )
    data3 <- tibble(
        x = xp,
        y = yp
    )

    fig_start( name_out )
    print(
        ggplot( data2, aes( x = x, y = p ) ) +
        geom_errorbar( aes( ymin = lower, ymax = upper ), width = .5 ) +
        geom_point( ) +
        #    expand_limits( y = c(0,1) ) + 
        theme_classic() +
        labs( x = 'Genotype', y = 'Prob(case)' ) +
        geom_abline( intercept = coef_lm[1], slope = coef_lm[2], col = 'red' ) +
        geom_line( aes( x = x, y = y ), data = data3, col = 'blue' )
    )
    fig_end()
}

dom_test <- function( x, y ) {
    # the het encoding version
    h <- ifelse( x == 1, 2, 0 )
    obj <- glm( y ~ x + h, family = binomial() )
    objsum <- summary( obj )
    coefs <- coef( objsum )
    # extract desired data and return that only
    return(
        tibble(
            bx = coefs[ 'x', 'Estimate' ],
            bh = coefs[ 'h', 'Estimate' ],
            px = coefs[ 'x', 'Pr(>|z|)' ],
            ph = coefs[ 'h', 'Pr(>|z|)' ]
        )
    )
}

# keep track of dominance tests
dt <- NULL

# repeat analysis for all SNPs
for ( i in 1 : nrow( top_snps ) ) {
    ## start with the top SNP
    #i <- which.min( top_snps$p )
    id <- top_snps$rsid[ i ]
    x <- X2[ , id ]

    # plot
    name_out <- paste0( 'dominance-', top_snps$chr[i], '-', top_snps$pos[i] )
    plot_geno_case_ctrl( x, y, name_out )

    # proper statistical test of dominance
    dti <- dom_test( x, y )
    # add info to row before concatenating
    dti$id <- id
    dt <- bind_rows( dt, dti )
}

# wow, one of these SNPs has a highly significance dominance result!
dt
##      bx      bh       px       ph id               
## 1 1.05   0.154  6.24e- 6 2.24e- 1 chr6:30770823:A:G
## 2 0.232  0.160  2.91e- 2 1.11e- 2 chr6:31151645:C:T
## 3 0.789  0.0524 7.60e-19 3.48e- 1 chr6:31633043:A:G
## 4 0.414  0.101  1.38e- 8 3.49e- 2 chr6:32193589:T:C
## 5 0.645  0.277  1.53e-21 2.81e-10 chr6:32593082:A:C
## 6 0.471  0.152  1.79e-14 3.18e- 4 chr6:32658925:A:G
## 7 1.06  -0.0191 9.56e-50 6.90e- 1 chr6:32666541:C:A
## 8 0.921 -0.0734 1.09e- 5 5.37e- 1 chr6:32669903:T:C

# a manual version of the SNPs from the main conditional analysis (table 1 from paper)
tab1 <- tibble(
    chr = c(6, 6, 6, 10, 16),
    pos = c(32652506, 32689478, 31361670, 28810849, 11077745),
    ref = c('C', 'C', 'A', 'C', 'A'),
    alt = c('T', 'T', 'G', 'T', 'G')
)
tab1 <- tab1 %>% mutate( id = paste( paste0('chr', chr), pos, ref, alt, sep = ':' ) )

# extract these variants now!
# order automatically matches tab1
X3 <- X[ , tab1$id ]

# save key data to revisit more quickly next time
save( X3, data, tab1, file = 'dominance-tab1.RData' )

# keep track of dominance tests
dt2 <- NULL

# repeat analysis for all SNPs
for ( i in 1 : nrow( tab1 ) ) {
    ## start with the top SNP
    #i <- which.min( top_snps$p )
    id <- tab1$id[ i ]
    x <- X3[ , id ]

    # plot
    name_out <- paste0( 'dominance-tab1-', tab1$chr[i], '-', tab1$pos[i] )
    plot_geno_case_ctrl( x, y, name_out )

    # proper statistical test of dominance
    dti <- dom_test( x, y )
    # add info to row before concatenating
    dti$id <- id
    dt2 <- bind_rows( dt2, dti )
}

# sad, this is way less significant, but sort of for one of the chr6 loci
dt2
##       bx      bh       px     ph id                
## 1  0.780 -0.0525 3.50e-35 0.240  chr6:32652506:C:T 
## 2  0.711  0.117  1.71e-21 0.0144 chr6:32689478:C:T 
## 3 -6.50   3.14   9.58e- 1 0.959  chr6:31361670:A:G 
## 4 -6.45   3.74   9.63e- 1 0.957  chr10:28810849:C:T
## 5 -0.245  0.0836 1.31e- 4 0.0581 chr16:11077745:A:G

# combine dominance tables to save for later
dt$dataset <- 'ct-best'
dt2$dataset <- 'conditional'
dt <- bind_rows( dt2, dt )
# add dominance factor
dt <- dt %>% mutate( d = 0.5 + bh/bx )

##        bx      bh       px       ph id                 dataset           d
##  1  0.780 -0.0525 3.50e-35 2.40e- 1 chr6:32652506:C:T  conditional  0.433 
##  2  0.711  0.117  1.71e-21 1.44e- 2 chr6:32689478:C:T  conditional  0.664 
##  3 -6.50   3.14   9.58e- 1 9.59e- 1 chr6:31361670:A:G  conditional  0.0169
##  4 -6.45   3.74   9.63e- 1 9.57e- 1 chr10:28810849:C:T conditional -0.0805
##  5 -0.245  0.0836 1.31e- 4 5.81e- 2 chr16:11077745:A:G conditional  0.159 
##  6  1.05   0.154  6.24e- 6 2.24e- 1 chr6:30770823:A:G  ct-best      0.646 
##  7  0.232  0.160  2.91e- 2 1.11e- 2 chr6:31151645:C:T  ct-best      1.19  
##  8  0.789  0.0524 7.60e-19 3.48e- 1 chr6:31633043:A:G  ct-best      0.566 
##  9  0.414  0.101  1.38e- 8 3.49e- 2 chr6:32193589:T:C  ct-best      0.744 
## 10  0.645  0.277  1.53e-21 2.81e-10 chr6:32593082:A:C  ct-best      0.929 
## 11  0.471  0.152  1.79e-14 3.18e- 4 chr6:32658925:A:G  ct-best      0.822 
## 12  1.06  -0.0191 9.56e-50 6.90e- 1 chr6:32666541:C:A  ct-best      0.482 
## 13  0.921 -0.0734 1.09e- 5 5.37e- 1 chr6:32669903:T:C  ct-best      0.420 

write_tsv( dt, 'dominance-pvals.txt.gz' )

## colMeans(X3)/2
##  chr6:32652506:C:T  chr6:32689478:C:T  chr6:31361670:A:G chr10:28810849:C:T 
##         0.50035063         0.24076671         0.04301075         0.03085554 
## chr16:11077745:A:G 
##         0.39703132 
