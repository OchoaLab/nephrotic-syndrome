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
## Multiple R-squared:  0.004856,	Adjusted R-squared:  0.003824 
## F-statistic: 4.708 on 2 and 1930 DF,  p-value: 0.009123

summary( lm( y ~ age, data = data ) )
##              Estimate Std. Error t value Pr(>|t|)    
## age         0.0021828  0.0006762   3.228  0.00127 **
## Multiple R-squared:  0.005595,	Adjusted R-squared:  0.005058 
## F-statistic: 10.42 on 1 and 1852 DF,  p-value: 0.001269

summary( lm( y ~ age + I(age^2), data = data ) )
##               Estimate Std. Error t value Pr(>|t|)    
## age          2.279e-02  2.099e-03  10.858  < 2e-16 ***
## I(age^2)    -3.396e-04  3.285e-05 -10.338  < 2e-16 ***
## Multiple R-squared:  0.05988,	Adjusted R-squared:  0.05886 
## F-statistic: 58.95 on 2 and 1851 DF,  p-value: < 2.2e-16
anova( lm( y ~ 1, data = data %>% filter( !is.na( age ) ) ), lm( y ~ age + I(age^2), data = data ) )$`Pr(>F)`
# [1]           NA 1.518723e-25

summary( lm( y ~ race, data = data ) )
##              Estimate Std. Error t value Pr(>|t|)    
## raceBlack     0.27313    0.02989   9.139  < 2e-16 ***
## raceHispanic  0.33578    0.06354   5.285 1.40e-07 ***
## raceOther     0.23702    0.04078   5.812 7.22e-09 ***
## raceWhite     0.26619    0.02447  10.878  < 2e-16 ***
## Multiple R-squared:  0.06867,	Adjusted R-squared:  0.06674 
## F-statistic: 35.54 on 4 and 1928 DF,  p-value: < 2.2e-16
anova( lm( y ~ 1, data = data ), lm( y ~ race, data = data ) )$`Pr(>F)`
# [1]           NA 1.102155e-28

summary( lm( y ~ dataset, data = data ) )
##                  Estimate Std. Error t value Pr(>|t|)    
## datasetCureGN     0.07192    0.02916   2.466   0.0137 *  
## datasetDiscovery -0.12073    0.02404  -5.022 5.58e-07 ***
## Multiple R-squared:  0.03166,	Adjusted R-squared:  0.03065 
## F-statistic: 31.55 on 2 and 1930 DF,  p-value: 3.296e-14

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
## Multiple R-squared:  0.08714,	Adjusted R-squared:  0.08268 
## F-statistic: 19.56 on 9 and 1844 DF,  p-value: < 2.2e-16

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
## Multiple R-squared:  0.1186,	Adjusted R-squared:  0.1138 
## F-statistic:  24.8 on 10 and 1843 DF,  p-value: < 2.2e-16

anova( lm( y ~ sex + age + I(age^2) + race + dataset, data = data ) )
##             Df Sum Sq Mean Sq  F value    Pr(>F)    
## sex          2   2.04  1.0221   5.7462 0.0032525 ** 
## age          1   2.02  2.0183  11.3464 0.0007715 ***
## I(age^2)     1  19.72 19.7180 110.8512 < 2.2e-16 ***
## race         4  18.41  4.6015  25.8688 < 2.2e-16 ***
## dataset      2   1.93  0.9647   5.4233 0.0044834 ** 

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
## Multiple R-squared:  0.1691,	Adjusted R-squared:  0.1655 
## F-statistic: 46.92 on 10 and 2305 DF,  p-value: < 2.2e-16

# ok, let's stick with pediatric curegn (y)

###################
### CONDITIONAL ###
###################

# need to do things a different way to get the conditional values I'd like to have on the table
mod_full <- lm( y ~ sex + age + I(age^2) + race + dataset, data = data )
mod_minus_sex <- lm( y ~ age + I(age^2) + race + dataset, data = data )
mod_minus_age <- lm( y ~ sex + race + dataset, data = data %>% filter( !is.na( age ) ) ) # ages are missing, it's important to match them up this way or ANOVA isn't happy
mod_minus_race <- lm( y ~ sex + age + I(age^2) + dataset, data = data )
mod_minus_dataset <- lm( y ~ sex + age + I(age^2) + race, data = data )

# conditional R2 values for each set
( summary( mod_full )$r.squared - summary( mod_minus_sex )$r.squared ) / ( 1 - summary( mod_minus_sex )$r.squared )
# 0.0041969
( summary( mod_full )$r.squared - summary( mod_minus_age )$r.squared ) / ( 1 - summary( mod_minus_age )$r.squared )
# 0.03456379
( summary( mod_full )$r.squared - summary( mod_minus_race )$r.squared ) / ( 1 - summary( mod_minus_race )$r.squared )
# 0.04287696
( summary( mod_full )$r.squared - summary( mod_minus_dataset )$r.squared ) / ( 1 - summary( mod_minus_dataset )$r.squared )
# 0.005850832

# p-values are easier
anova( mod_minus_sex, mod_full )$`Pr(>F)`
# [1]         NA 0.02074227
anova( mod_minus_age, mod_full )$`Pr(>F)`
# [1]           NA 8.371106e-15
anova( mod_minus_race, mod_full )$`Pr(>F)`
# [1]           NA 1.173209e-16
anova( mod_minus_dataset, mod_full )$`Pr(>F)`
# [1]         NA 0.00448339


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

# change language to self-reported ancestry
data$race[ data$race == 'Black' ] <- 'African'
data$race[ data$race == 'White' ] <- 'European'

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
data_race$Race <- factor( data_race$Race, levels = c('African', 'European', 'Asian', 'Hispanic', 'Other') )

fig_start( 'race', width = 4 )
ggplot( data_race, aes( x = Race, y = SRNS ) ) +
    geom_col() +
    geom_errorbar( aes( ymin = SRNS - se, ymax = SRNS + se ), width = 0.5 ) +
    theme_classic() +
    labs( x = 'Self-reported Ancestry' )
fig_end()

# do it by race now, but stratifying by dataset still
# this is a more complicated table
# glad as_tibble does a very good job giving me what I want!!!
c <- with( data %>% filter( y ), table( race, dataset ) ) %>% as_tibble %>% rename( c = n )
n <- with( data, table( race, dataset ) ) %>% as_tibble
x <- full_join( c, n )
x <- x %>% mutate( p = c/n, se = sqrt( p * ( 1 - p ) / n ) )
x <- x %>% mutate(
               race = factor( race, levels = c('African', 'European', 'Asian', 'Hispanic', 'Other') ),
               dataset = factor( dataset, levels = c('Discovery', 'Bristol', 'CureGN') )
           )
x <- x %>% rename( SRNS = p, Race = race, Dataset = dataset )

pd <- position_dodge(width = 0.9) # Define a common dodge width

fig_start( 'race-dataset', width = 6 )
ggplot( x, aes( x = Race, y = SRNS, fill = Dataset ) ) +
    geom_col( position = pd ) +
    geom_errorbar( aes( ymin = SRNS - se, ymax = SRNS + se ), width = 0.5, position = pd ) +
    theme_classic() +
    labs( x = 'Self-reported ancestry', y = 'SRNS prevalence', fill = 'Cohort' )
fig_end()

###########
### AGE ###
###########

# do something nice as in the previous plot, or maybe something more like the current analyses where we estimate proportions with some smooth model

# Bristol is the only case with larger values, though they're all supposed to be pediatric.  Let's toss the likely bad/corrupted age data
data2 <- data %>% filter( age <= 21 )
## # copy data to have a combined dataset with a more stable curve
## data3 <- data2
## data3$dataset <- 'Combined'
## data2 <- bind_rows( data2, data3 )
## data2$dataset <- factor( data2$dataset, levels = c('Combined', 'Discovery', 'Bristol', 'CureGN') )
data2$dataset <- factor( data2$dataset, levels = c('Discovery', 'Bristol', 'CureGN') )

# plot smooth fits, default ggplot stuff!
fig_start( 'age', width = 6 )
ggplot( data2, aes( x = age, y = as.numeric( y ), color = dataset ) ) +
    geom_smooth() +
    theme_classic() +
    labs( x = 'Age at onset of NS', y = 'SRNS prevalence', color = 'Cohort' )
fig_end()
