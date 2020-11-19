library(genio)
library(BEDMatrix)
library(tibble)
library(readr)
library(ochoalabtools)

# get name from command line
args <- args_cli()
name <- args[1]
if (is.na(name))
    stop('provide BED file name as first argument!')

#name <- 'combine2' # OLD
#name <- '/dmpi/analysis/Rasheed/ssns_gwas/merge_tgp/combine'   # NEW

fam <- read_fam(name)
bim <- read_bim(name)
X <- BEDMatrix(name, simple_names = TRUE)

## indexes_remove <- ( bim$ref == 'A' & bim$alt == 'T') | ( bim$ref == 'T' & bim$alt == 'A' ) | ( bim$ref == 'C' & bim$alt == 'G' ) | ( bim$ref == 'G' & bim$alt == 'C' )
## ids_remove <- bim$id[ indexes_remove ]

## mean( indexes_remove ) 
## # [1] 0.1157676 # OLD
## # [1] 0         # NEW

icase <- fam$pheno == 2
ictrl <- fam$pheno == 1

# all in mem version
pcase <- colMeans( X[ icase, ], na.rm = TRUE ) / 2
pctrl <- colMeans( X[ ictrl, ], na.rm = TRUE ) / 2

# put in a tibble
data <- tibble(
    chr = bim$chr,
    pos = bim$pos,
    id = bim$id,
    ref = bim$ref,
    alt = bim$alt,
#    rm = indexes_remove, # OLD
    pcase = pcase,
    pctrl = pctrl,
    pdiff = abs(pcase - pctrl)
)

# save analysis
#name <- 'OLD-combine2-transversions' # OLD
name <- paste0( name, '.af-case-ctrl' ) # data name first
write_tsv( data, paste0( name, '.txt.gz') )

# # read back locally:
# data <- read_tsv( paste0( name, '.txt.gz') )

# plot frequencies case vs ctrl
fig_start( name )
plot(
    data$pctrl,
    data$pcase,
    xlab = 'AF controls',
    ylab = 'AF cases',
    xlim = c(0, 1), # only this one is needed, most AFs are < 0.5, a few over 0.6 though
    ylim = c(0, 1),
    pch = '.'
)
abline( 0, 1, lty = 2, col = 'gray' )
fig_end()

## quantile( data$pdiff )
## # OLD
## ##          0%         25%         50%         75%        100% 
## ## 0.000000000 0.002134924 0.006561533 0.021533894 0.999800319 
## # NEW
## ##          0%         25%         50%         75%        100% 
## ## 0.000000000 0.002136766 0.006662031 0.021286048 0.999800319 


## ##############################################
## # OLD versions
## # require `rm` column to be non-trivial


## mean( data$pdiff[ data$rm ] )
## # [1] 0.03239746
## mean( data$pdiff[ !data$rm ] )
## # [1] 0.01855753

## median( data$pdiff[ data$rm ] )
## # [1] 0.005679507
## median( data$pdiff[ !data$rm ] )
## # [1] 0.006697584

## max( data$pdiff[ data$rm ] )
## # [1] 0.9998003
## max( data$pdiff[ !data$rm ] )
## # [1] 0.9998003

## # no huge differences by RM status, or by CHR, huge outliers across all cases!
## boxplot( data$pdiff ~ data$rm )
## boxplot( data$pdiff ~ data$chr )

## data[ data$pdiff > 0.9, ]
## ## # A tibble: 2,766 x 9
## ##      chr     pos id          ref   alt   rm    pcase   pctrl pdiff
## ##  1     1 1302327 rs112122411 C     A     FALSE 0.932 0.00759 0.925
## ##  2     1 1799851 rs556554846 C     G     TRUE  1     0.00479 0.995
## ##  3     1 2219479 rs149672889 G     C     TRUE  1     0.00180 0.998
## ##  4     1 3034557 rs373023536 T     C     FALSE 1     0.00879 0.991
## ##  5     1 3645289 rs181718533 A     G     FALSE 1     0.00499 0.995
## ##  6     1 5715710 rs114321438 C     T     FALSE 1     0.0198  0.980
## ##  7     1 5731935 rs77945749  A     G     FALSE 1     0.0483  0.952
## ##  8     1 5743801 rs3913320   C     G     TRUE  0.970 0.0455  0.924
## ##  9     1 6855514 rs144409026 A     G     FALSE 1     0.00220 0.998
## ## 10     1 7993901 rs9657983   T     A     TRUE  1     0.00399 0.996
## ## # … with 2,756 more rows

## setwd('~/dbs/ssns_gwas')
## bim <- read_bim('ssns.gwas.hg37.bim')

## bim[ bim$chr == 1 & bim$pos == 1302327, ]
## ##   chr   id             posg     pos ref   alt  
## ## 1 1     JHU_1.1302326     0 1302327 A     C    
## ## 1 1     1:1799851-C-G     0 1799851 0     C    
## ## 1 1     rs149672889   0.202 2219479 C     G    
## ## 1 1     JHU_1.3034556  2.47 3034557 0     T    
## ## 1 1     rs181718533    4.17 3645289 0     A    
## ## 1 1     JHU_1.5715709  9.32 5715710 0     C    
## ## 1 1     JHU_1.5731934  9.36 5731935 0     A    
## ## 1 1     rs3913320      9.38 5743801 G     C    
## ## 1 1     rs144409026    12.0 6855514 0     A    
## ## 1 1     rs9657983      13.8 7993901 A     T    

## mean( bim$ref == '0' )
## ## [1] 0.08419848
## mean( bim$alt == '0' )
## ## [1] 0.0009049049

## table( bim$ref )
## ##      0      A      C      D      G      I      T 
## ## 147200 469038 326756  12340 326008   4488 462420 
## table( bim$alt )
## ##      0      A      C      D      G      I      T 
## ##   1582 356620 507116   6580 500595  17869 357888 
