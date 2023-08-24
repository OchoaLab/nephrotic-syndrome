# Run simple trait ~ sex + race model
# conclusions:
# - sex is significant explanatory variable for SSNS but not SRNS (may be due to small sample size of latter)
# - race is always significant

library(tidyverse)

# this is where our data is
setwd( '/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/' )

# load covariates file
data <- read_tsv( 'patient-data-merged-tgp.txt.gz', col_types = 'ccccdiiii' )

# regress ns-vs-control to covariates
summary( glm( ns_ctrl ~ sex + race, data = data, family = binomial ) )
##               Estimate Std. Error z value Pr(>|z|)    
## (Intercept)   -1.53963    0.07634 -20.168  < 2e-16 ***
## sexmale        0.42073    0.07639   5.508 3.63e-08 ***
## sexunknown    14.10569  187.49088   0.075  0.94003    
## raceBlack     -0.20855    0.09255  -2.253  0.02424 *  
## raceHispanic  -0.50070    0.15522  -3.226  0.00126 ** 
## raceOther      2.85222    0.49486   5.764 8.23e-09 ***
## raceWhite      0.13890    0.09519   1.459  0.14450    

## repeat for subtypes!  Missing responses are automatically ignored, so no need to explicitly subset rows!

summary( glm( ssns_ctrl ~ sex + race, data = data, family = binomial ) )
##               Estimate Std. Error z value Pr(>|z|)    
## (Intercept)   -1.62249    0.08144 -19.922  < 2e-16 ***
## sexmale        0.48419    0.08560   5.656 1.55e-08 ***
## sexunknown    14.18855  229.62848   0.062   0.9507    
## raceBlack     -0.52958    0.10242  -5.171 2.33e-07 ***
## raceHispanic  -1.13459    0.20298  -5.590 2.27e-08 ***
## raceOther      2.62862    0.50850   5.169 2.35e-07 ***
## raceWhite     -0.23218    0.10652  -2.180   0.0293 *  

summary( glm( srns_ctrl ~ sex + race, data = data, family = binomial ) )
##              Estimate Std. Error z value Pr(>|z|)    
## (Intercept)   -4.5823     0.2747 -16.683  < 2e-16 ***
## sexmale        0.2012     0.1514   1.329    0.184    
## sexunknown    17.1484   324.7438   0.053    0.958    
## raceBlack      1.7097     0.2875   5.947 2.73e-09 ***
## raceHispanic   1.7162     0.3404   5.042 4.62e-07 ***
## raceOther      4.2244     0.7201   5.866 4.46e-09 ***
## raceWhite      2.2119     0.2849   7.765 8.15e-15 ***

summary( glm( ssns_srns ~ sex + race, data = data, family = binomial ) )
##              Estimate Std. Error z value Pr(>|z|)    
## (Intercept)   -2.9446     0.2824 -10.427  < 2e-16 ***
## sexmale       -0.3157     0.1766  -1.787   0.0739 .  
## sexunknown     2.2514     1.2569   1.791   0.0732 .  
## raceBlack      2.2537     0.2999   7.515 5.68e-14 ***
## raceHispanic   2.8587     0.3873   7.380 1.58e-13 ***
## raceOther      1.6134     0.6113   2.639   0.0083 ** 
## raceWhite      2.4424     0.2973   8.215  < 2e-16 ***
