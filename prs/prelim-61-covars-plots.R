# a basic covariate analysis motivating aim3 (probabilistic prediction of SSNS vs SRNS using age, sex, race, etc)

library(tidyverse)
library(ochoalabtools)

# go where the main data is
setwd( '/datacommons/ochoalab/ssns_gwas/' )
data <- read_tsv( 'covars-ssns-srns-discov-bristol-curegn.txt.gz' )
# make response numeric (binary)
data <- data %>% mutate( y = diagnosis == 'SRNS' ) %>% select( -diagnosis )
data <- data %>% mutate( y2 = diagnosis2 == 'SRNS' ) %>% select( -diagnosis2 )

# simple models: single variables (always toss intercept and sexunknown)
summary( lm( y ~ sex, data = data ) )
##              Estimate Std. Error t value Pr(>|t|)    
## sexmale     -0.064483   0.021052  -3.063  0.00222 ** 
summary( lm( y ~ age, data = data ) )
##              Estimate Std. Error t value Pr(>|t|)    
## age         0.0021828  0.0006762   3.228  0.00127 ** 
summary( lm( y ~ race, data = data ) )
##              Estimate Std. Error t value Pr(>|t|)    
## raceBlack     0.27313    0.02989   9.139  < 2e-16 ***
## raceHispanic  0.33578    0.06354   5.285 1.40e-07 ***
## raceOther     0.23702    0.04078   5.812 7.22e-09 ***
## raceWhite     0.26619    0.02447  10.878  < 2e-16 ***
summary( lm( y ~ dataset, data = data ) )
##                  Estimate Std. Error t value Pr(>|t|)    
## datasetCureGN     0.07192    0.02916   2.466   0.0137 *  
## datasetDiscovery -0.12073    0.02404  -5.022 5.58e-07 ***

# since they're all strong separately, let's just consider them all together
summary( lm( y ~ sex + age + race + dataset, data = data ) )
##                    Estimate Std. Error t value Pr(>|t|)    
## sexmale          -0.0563411  0.0205036  -2.748 0.006057 ** 
## age               0.0003176  0.0007691   0.413 0.679718    
## raceBlack         0.2535681  0.0307873   8.236 3.33e-16 ***
## raceHispanic      0.3681251  0.0663288   5.550 3.27e-08 ***
## raceOther         0.1544986  0.0436950   3.536 0.000417 ***
## raceWhite         0.2006926  0.0278262   7.212 7.99e-13 ***
## datasetCureGN     0.0603965  0.0312378   1.933 0.053334 .  
## datasetDiscovery -0.0739610  0.0302724  -2.443 0.014652 *  

# lose age signal, try making it non-linear: yey that works excellently!
summary( lm( y ~ sex + age + I(age^2) + race + dataset, data = data ) )
##                    Estimate Std. Error t value Pr(>|t|)    
## sexmale          -5.054e-02  2.017e-02  -2.506 0.012287 *  
## age               1.684e-02  2.173e-03   7.751 1.49e-14 ***
## I(age^2)         -2.675e-04  3.298e-05  -8.112 8.99e-16 ***
## raceBlack         2.229e-01  3.050e-02   7.309 3.99e-13 ***
## raceHispanic      3.547e-01  6.521e-02   5.439 6.06e-08 ***
## raceOther         1.450e-01  4.296e-02   3.376 0.000751 ***
## raceWhite         1.903e-01  2.738e-02   6.950 5.04e-12 ***
## datasetCureGN     4.318e-02  3.078e-02   1.403 0.160792    
## datasetDiscovery -4.878e-02  2.992e-02  -1.631 0.103163    

# conditional on age and race, sex is less important, and the dataset differences vanish!

# try expanded curegn (y2) that uses adults too.  Is this very different?
# YES! curegn as a dataset becomes a significant variable, probably means bad ascertainment effects!!!
summary( lm( y2 ~ sex + age + I(age^2) + race + dataset, data = data ) )
##                    Estimate Std. Error t value Pr(>|t|)    
## sexmale          -0.0165760  0.0186032  -0.891 0.373005    
## age               0.0163846  0.0017149   9.554  < 2e-16 ***
## I(age^2)         -0.0002006  0.0000242  -8.287  < 2e-16 ***
## raceBlack         0.2208646  0.0291713   7.571 5.31e-14 ***
## raceHispanic      0.3483156  0.0673421   5.172 2.51e-07 ***
## raceOther         0.1485807  0.0408408   3.638 0.000281 ***
## raceWhite         0.1605643  0.0261565   6.139 9.77e-10 ***
## datasetCureGN     0.1603677  0.0250265   6.408 1.78e-10 ***
## datasetDiscovery  0.0008025  0.0288679   0.028 0.977825    

# ok, let's stick with pediatric curegn (y)

###############
### DATASET ###
###############

# I want a simple table showing me first the counts of ssns vs srns by dataset
x <- with( data, table( y, dataset ) )
# marginal counts per dataset
n <- colSums( x )
# normalized counts, want for SRNS only (y=TRUE)
p <- x[ 'TRUE', ] / n
# and error bars:
# binomial variance of count is n*p*(1-p), after normalizing it becomes p*(1-p)/n.  square root of this to have SD?  Yes, the internet agrees
se <- sqrt( p * ( 1 - p ) / n )

# make a plot of this
data_dataset <- tibble(
    Dataset = names( p ),
    SRNS = p,
    se = se
)
data_dataset$Dataset <- factor( data_dataset$Dataset, levels = c('Discovery', 'Bristol', 'CureGN') )

fig_start( 'dataset' )
ggplot( data_dataset, aes( x = Dataset, y = SRNS ) ) +
    geom_col() +
    geom_errorbar( aes( ymin = SRNS - se, ymax = SRNS + se ), width = 0.5 ) +
    theme_classic()
fig_end()

############
### RACE ###
############

# let's to the same but for race
x <- with( data, table( y, race ) )
n <- colSums( x )
p <- x[ 'TRUE', ] / n
se <- sqrt( p * ( 1 - p ) / n )
data_race <- tibble(
    Race = names( p ),
    SRNS = p,
    se = se
)
data_race$Race <- factor( data_race$Race, levels = c('Black', 'White', 'Asian', 'Hispanic', 'Other') )

fig_start( 'race' )
ggplot( data_race, aes( x = Race, y = SRNS ) ) +
    geom_col() +
    geom_errorbar( aes( ymin = SRNS - se, ymax = SRNS + se ), width = 0.5 ) +
    theme_classic()
fig_end()

# do it by race now, but stratifying by dataset still
# this is a more complicated table
# glad as_tibble does a very good job giving me what I want!!!
c <- with( data %>% filter( y ), table( race, dataset ) ) %>% as_tibble %>% rename( c = n )
n <- with( data, table( race, dataset ) ) %>% as_tibble
x <- full_join( c, n )
x <- x %>% mutate( p = c/n, se = sqrt( p * ( 1 - p ) / n ) )
x <- x %>% mutate(
               race = factor( race, levels = c('Black', 'White', 'Asian', 'Hispanic', 'Other') ),
               dataset = factor( dataset, levels = c('Discovery', 'Bristol', 'CureGN') )
           )
x <- x %>% rename( SRNS = p, Race = race, Dataset = dataset )

pd <- position_dodge(width = 0.9) # Define a common dodge width

fig_start( 'race-dataset', width = 6 )
ggplot( x, aes( x = Race, y = SRNS, fill = Dataset ) ) +
    geom_col( position = pd ) +
    geom_errorbar( aes( ymin = SRNS - se, ymax = SRNS + se ), width = 0.5, position = pd ) +
    theme_classic()
fig_end()

###########
### AGE ###
###########

# do something nice as in the previous plot, or maybe something more like the current analyses where we estimate proportions with some smooth model

# TODO: fold large ages or exclude them???
