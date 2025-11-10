library(tidyverse)
library(PRROC)
library(ochoalabtools)
#library(ggpubr)

### DATA: copy of previous net-benefit analysis ###

# constants
#prevalence_srns <- 0.2
# ordered as in table, not as in paper (DRB1 listed first in paper, last here)
haplo <- '02:01~02:02~07:01'

setwd( 'test' )

# load clean data I had gathered for Debo previously
data <- read_tsv( 'bristol-curegn-info-pcs.txt.gz', show_col_types = FALSE )

# load HLA haplotype, which is in different locations
haps <- read_tsv( 'hla-haplotype.txt', col_types = 'ccc' )
haps2 <- read_tsv( '../test-curegn/hla-haplotype.txt', col_types = 'ccc' )
haps <- bind_rows( haps, haps2 )
colnames( haps ) <- c('id', 'hap1', 'hap2')
# process to get dosage of top haplotype
haps$dosage <- (haps$hap1 == haplo) + (haps$hap2 == haplo)

# numbers of people are wildly different
nrow( haps ) #[1] 2434
nrow( data ) #[1] 936
# everybody in data is also in haps, makes sense
stopifnot( all( data$id %in% haps$id ) )
# visually, it looks like the majority are extra cureGN samples (probably adults or non-NS)
# use inner join!
data <- inner_join( data, haps )

# make sure SSNS is treated as disease (because steroids treats it)
data <- data %>% mutate( y = as.numeric( pheno == 'SRNS' ) )

# toss individuals without age, much more complicated keeping them conditionally (tosses 15 individuals only)
sum( is.na( data$age ) )
# [1] 15
data <- data %>% filter( !is.na( age ) )

# report final sample size
nrow( data )
# [1] 921


### MODEL FITS ###

## # calculate weights to adjust for prevalence in fitting GLMs!
## x <- table( data$y )
## p <- x/sum(x)
## ##         0         1 
## ## 0.6471227 0.3528773 
## # we want p = prevalence, these are the weights we need.  Note SSNS must be upweighted and SRNS downweighted
## w <- c( 1 - prevalence_srns, prevalence_srns ) / p
## ##         0         1 
## ## 1.2362416 0.5667692 
## # spread over data frame
## data <- data %>% mutate( weight = ifelse( y == 1, w[2], w[1] ) )

# calibration isn't important here, but this brings uniformity to various single predictors (so we don't have to reverse signs) and especially for combined predictors!
# unlike net benefit analysis, don't use weights here to match what Debo calculated before
mod_ldpred <- glm( y ~ prs_ldpred, data = data, family = binomial() )
mod_hla <- glm( y ~ dosage, data = data, family = binomial() )
mod_age <- glm( y ~ age, data = data, family = binomial() )
mod_comb1 <- glm( y ~ prs_ldpred + dosage, data = data, family = binomial() )
mod_comb2 <- glm( y ~ prs_ldpred + dosage + age, data = data, family = binomial() )

# apply calibration here.  this overfits, but right now we don't have a better way!
data <- data %>% mutate(
                     prob_ldpred = broom::augment( mod_ldpred, type.predict = "response" ) %>% pull( ".fitted" ),
                     prob_hla = broom::augment( mod_hla, type.predict = "response" ) %>% pull( ".fitted" ),
                     prob_age = broom::augment( mod_age, type.predict = "response" ) %>% pull( ".fitted" ),
                     prob_comb1 = broom::augment( mod_comb1, type.predict = "response" ) %>% pull( ".fitted" ),
                     prob_comb2 = broom::augment( mod_comb2, type.predict = "response" ) %>% pull( ".fitted" )
                 )

### ROC/AUC ###

# let's do unweighted version
# same for all cases, true classifications
y <- data$y

labels <- list(
    prob_age = 'Age',
    prob_hla = 'Haplotype',
    prob_ldpred = 'PRS',
    prob_comb1 = 'Haplotype + PRS',
    prob_comb2 = 'Haplotype + PRS + Age'
)

# keep list of AUCs
aucs <- list()
# plotting data
datap <- NULL

for ( name in names( labels ) ) {
    obj <- roc.curve( data[[ name ]], weights.class0 = y, curve = TRUE )
    # store AUCs for later calculations
    aucs[[ name ]] <- obj$auc
    # form AUC plotting data
    # columns: fpr, recall/sensitivity, score
    x <- obj$curve
    colnames( x ) <- c( 'FPR', 'Sensitivity', 'score' )
    datapi <- as_tibble( x )
    datapi$Model <- labels[[ name ]]
    datap <- bind_rows( datap, datapi )
}

# save AUCs
aucs <- pivot_longer( as_tibble( aucs ), cols = everything(), names_to = 'model', values_to = 'auc' )
aucs <- aucs %>% mutate( model = unlist( labels[ model ] ) )
write_tsv( aucs, 'delta-auc.txt' )

### PLOTS ###

wh <- fig_scale( 4/3 )
fig_start( 'delta-auc', width = wh[1], height = wh[2] )
ggplot( datap, aes( x = FPR, y = Sensitivity, col = Model ) ) +
    geom_line() +
    theme_classic() +
    geom_abline( slope = 1, intercept = 0, linetype = 'dashed', color = 'gray' ) +
    labs( x = '1 - Specificity' )
fig_end()

