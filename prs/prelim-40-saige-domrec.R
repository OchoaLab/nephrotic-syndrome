# look at newest "domrec + saige" results, determine which model overall seems more significant, how the top loci shift

library(readr)
library(dplyr)

# constants
pcut <- 5e-8

# these are all "base" gwas runs for PRS devel
setwd( 'base' )

# load results from all models
add <- read_tsv( "saige_output.txt.gz", show_col_types = FALSE )
rec <- read_tsv( "mac20-rec_saige.txt.gz", show_col_types = FALSE )
dom <- read_tsv( "mac20-dom_saige.txt.gz", show_col_types = FALSE )

# rec has more fixed cases, which don't produce outputs; dom has same snps as add!
nrow( add ) # 20511795
nrow( rec ) # 15074085
nrow( dom ) # 20511795

add_i <- which.min( add$p.value )
rec_i <- which.min( rec$p.value )
dom_i <- which.min( dom$p.value )

# top locus is exactly the same for add and rec, but it is less significant for rec
# It's a different locus for dom, but sadly significance is also not greater overall (they are close though)
add[ add_i, ] %>% select( MarkerID, BETA, p.value )
##   MarkerID           BETA  p.value
## 1 chr6:32652506:C:T 0.927 1.32e-32
rec[ rec_i, ] %>% select( MarkerID, BETA, p.value )
##   MarkerID           BETA  p.value
## 1 chr6:32652506:C:T 0.644 1.30e-30
dom[ dom_i, ] %>% select( MarkerID, BETA, p.value )
##   MarkerID           BETA  p.value
## 1 chr6:32714062:A:G 0.644 6.18e-30

##################
### ADD vs REC ###
##################

# make a merged version to make better statements about all other significant loci
# omit things unless they're in both
ar <- inner_join(
    add %>% select( id = MarkerID, pa = p.value ),
    rec %>% select( id = MarkerID, pr = p.value )
)
# now look at significant subset (significant in at least one)
ars <- ar %>% filter( pa < pcut | pr < pcut )
# overall, how often are p-values smaller in additive?  The vast majority of the time!
mean( ars$pa < ars$pr )
# [1] 0.9721497

# what about the cases where rec is better?
ars_rb <- ars %>% filter( pa > pr )
ars_rb <- ars_rb %>% mutate( fac = pa/pr )
ars_rb <- ars_rb %>% arrange( desc( fac ) )
# one chr15 locus becomes crazy significant, and also 2 from chr9, after that it's chr6 and factors become more modest (relatively speaking)
##    id                       pa       pr     fac
##  1 chr15:78962606:G:A 0.199    9.22e-29 2.16e27
##  2 chr9:81736734:G:C  0.529    3.28e-10 1.61e 9
##  3 chr9:81738631:A:G  0.514    3.28e-10 1.56e 9
##  4 chr6:32603321:A:G  0.000118 2.68e- 9 4.41e 4
##  5 chr6:32601786:C:T  0.000117 2.84e- 9 4.11e 4
##  6 chr6:32601795:G:A  0.000117 2.84e- 9 4.11e 4
##  7 chr6:32601945:G:T  0.000117 2.84e- 9 4.11e 4
##  8 chr6:32602015:T:C  0.000117 2.84e- 9 4.11e 4
##  9 chr6:32602484:AG:A 0.000117 2.84e- 9 4.11e 4
## 10 chr6:32603182:C:T  0.000117 2.84e- 9 4.11e 4

# in fact, those three weird loci are the only ones outside chr6.  Since the first one is a singleton, I believe it's more likely a false positive; the second two maybe too but less clear.
ars_rb %>% filter( !grepl( 'chr6', id ) )
##   id                    pa       pr     fac
## 1 chr15:78962606:G:A 0.199 9.22e-29 2.16e27
## 2 chr9:81736734:G:C  0.529 3.28e-10 1.61e 9
## 3 chr9:81738631:A:G  0.514 3.28e-10 1.56e 9

# conclusions: additive is a better fit than recessive for chr6 and everywhere else, except 3 weirdo loci (one in chr15, two neighbors in chr9)
# their significance can be re-assessed later with conditional analyses or looking at genes/annotations

##################
### ADD vs DOM ###
##################

# make a merged version to make better statements about all other significant loci
# yey these are already perfectly aligned!
stopifnot( all( add$MarkerID == dom$MarkerID ) )
# lazy merge for this special case
ad <- bind_cols(
    add %>% select( id = MarkerID, pa = p.value ),
    dom %>% select( pd = p.value )
)
# now look at significant subset (significant in at least one)
ads <- ad %>% filter( pa < pcut | pd < pcut )
# overall, how often are p-values smaller in additive?  Unlike recessive, here p-values are comparable in bulk sign
mean( ads$pa < ads$pd )
# [1] 0.5118197

# what about the cases where dom is better?
ads_db <- ads %>% filter( pa > pd )
ads_db <- ads_db %>% mutate( fac = pa/pd )
ads_db <- ads_db %>% arrange( desc( fac ) )
# these seem largely comparable, there's only one chr6 locus that is way more significant under dominance, the rest are not big factors
##    id                          pa       pd      fac
##  1 chr6:32634840:C:T 0.0000000554 2.58e-14 2146772.
##  2 chr6:32632540:T:C 0.0000211    2.80e- 9    7556.
##  3 chr6:32632818:T:C 0.0000211    2.80e- 9    7525.
##  4 chr6:32631140:A:C 0.0000210    2.80e- 9    7488.
##  5 chr6:32631169:T:G 0.0000210    2.80e- 9    7488.
##  6 chr6:32631403:G:A 0.0000210    2.80e- 9    7488.
##  7 chr6:32631413:T:A 0.0000210    2.80e- 9    7488.
##  8 chr6:32631488:G:A 0.0000210    2.80e- 9    7488.
##  9 chr6:32631708:T:G 0.0000210    2.80e- 9    7488.
## 10 chr6:32632141:T:C 0.0000210    2.80e- 9    7488.

# only three loci outside chr6 are more significant under the dominant model.  all are on chr10 and appear to be neighbors, so that seems like a good signal!  significance is only marginally better though
ads_db %>% filter( !grepl( 'chr6', id ) )
##   id                       pa       pd   fac
## 1 chr10:28810849:C:T 8.70e-12 4.73e-13 18.4 
## 2 chr10:28809695:C:T 5.08e- 8 5.06e- 9 10.0 
## 3 chr10:28815914:C:T 1.59e- 7 2.30e- 8  6.90

# conclusions:
# - additive is a better fit than recessive for chr6 and everywhere else, except 3 weirdo loci (one in chr15, two neighbors in chr9)
#   - their significance can be re-assessed later with conditional analyses or looking at genes/annotations
# - dominant is mostly comparable to additive, at least on this broad view of the data
