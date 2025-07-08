library(tidyverse)
library(genio)

# gathers individual data for Bristol and CureGN to share with Debo

# constants
name_data <- 'mac20'
names <- paste0( 'base-train-ldpred2-', c( 'grid-h0.1', 'ct' ), '-best' )
files_prs <- paste0( 'prs-', names, '.txt.gz' )

# all processing happens in subdirectory
setwd( 'test-curegn' )

# load data the same way in both cases
load_data <- function( dataset, info_path ) {
    # load main fam table, with IDs, sex, true phenotype
    # immediately throw away some useless columns (there's no family structure)
    data <- read_fam( name_data ) %>% select( id, sex, pheno )

    # load patient info
    info <- read_tsv( paste0( '../../../', info_path, '/patient-data.txt.gz' ), show_col_types = FALSE )
    # for curegn we can get rid of diagnoses (redundant with pheno, checked manually; sex also matches but for discovery format is different, easier to exclude while joining; match will be by ID, and race and age are added!)
    info <- info %>% select( id, race, age )

    # do inner join
    data <- inner_join( data, info )

    # make sex readable
    data <- data %>% mutate( sex = sex_to_char( sex ) )
    # mark dataset, in case we need that info later
    data$dataset <- dataset
    # load PRS per individual, attach to main table
    data$prs_ldpred <- as.numeric( read_lines( files_prs[1] ) )
    data$prs_ct <- as.numeric( read_lines( files_prs[2] ) )
    
    # load PCs as columns
    # in this case it's easier to use regular tidyverse instead of genio (because we want a data frame, not a matrix)
    pcs <- read_tsv( paste0( name_data, '.eigenvec' ), show_col_types = FALSE )
    # rename two columns before merging, but actually we don't need fam (it's just awkward with the hash)
    pcs <- pcs %>% rename( fam = '#FID', id = IID ) %>% select( -fam )
    # merge PCs into main table
    data <- full_join( data, pcs )
    
    return( data )
}

data <- load_data( 'CureGN', '../curegn' )

# now go back to test, load and merge, and outputs will be there too
setwd( '../test' )
# load as before
data2 <- load_data( 'Bristol', 'array' )

# concatenate, putting Bristol first
data <- bind_rows( data2, data )
# make trait a string, for clarity
data <- data %>% mutate( pheno = ifelse( pheno == 2, 'SSNS', 'SRNS' ) )

# minor harmonization
data$race[ data$race == 'Multiracial' ] <- 'Mixed'

# save this table!
write_tsv( data, 'bristol-curegn-info-pcs.txt.gz' )

## # do some test regressions
## data$pheno <- as.numeric( as.factor( data$pheno ) )

## lm( pheno ~ sex, data = data )
## ##             Estimate Std. Error t value Pr(>|t|)    
## ## (Intercept)  1.60105    0.02460  65.083   <2e-16 ***
## ## sexM         0.06201    0.03195   1.941   0.0525 .  
## ## Residual standard error: 0.4802 on 934 degrees of freedom
## ## Multiple R-squared:  0.004018,	Adjusted R-squared:  0.002952 
## ## F-statistic: 3.768 on 1 and 934 DF,  p-value: 0.05254
## summary( lm( pheno ~ race, data = data ) )
## ##             Estimate Std. Error t value Pr(>|t|)    
## ## (Intercept)  1.75000    0.04846  36.111  < 2e-16 ***
## ## raceBlack   -0.28097    0.06591  -4.263 2.22e-05 ***
## ## raceMixed    0.05000    0.09376   0.533  0.59395    
## ## raceNatAmr  -0.58333    0.19981  -2.919  0.00359 ** 
## ## racePacific -0.25000    0.33923  -0.737  0.46133    
## ## raceUnknown -0.11250    0.07188  -1.565  0.11790    
## ## raceWhite   -0.10265    0.05217  -1.968  0.04942 *  
## ## Residual standard error: 0.4748 on 929 degrees of freedom
## ## Multiple R-squared:  0.03132,	Adjusted R-squared:  0.02507 
## ## F-statistic: 5.007 on 6 and 929 DF,  p-value: 4.623e-05
## summary( lm( pheno ~ prs_ldpred, data = data ) )
## ##             Estimate Std. Error t value Pr(>|t|)    
## ## (Intercept)  1.40792    0.03486   40.39  < 2e-16 ***
## ## prs_ldpred   0.27630    0.03765    7.34 4.66e-13 ***
## ## Residual standard error: 0.4678 on 934 degrees of freedom
## ## Multiple R-squared:  0.05453,	Adjusted R-squared:  0.05352 
## ## F-statistic: 53.87 on 1 and 934 DF,  p-value: 4.657e-13
## summary( lm( pheno ~ prs_ct, data = data ) )
## ##             Estimate Std. Error t value Pr(>|t|)    
## ## (Intercept)  1.38593    0.03987  34.760  < 2e-16 ***
## ## prs_ct       0.05435    0.00794   6.845 1.38e-11 ***
## ## Residual standard error: 0.4695 on 934 degrees of freedom
## ## Multiple R-squared:  0.04777,	Adjusted R-squared:  0.04675 
## ## F-statistic: 46.85 on 1 and 934 DF,  p-value: 1.383e-11
## summary( lm( pheno ~ PC1, data = data ) )
## ##             Estimate Std. Error t value Pr(>|t|)    
## ## (Intercept)  1.63782    0.01568 104.465   <2e-16 ***
## ## PC1         -0.81618    0.33917  -2.406   0.0163 *  
## ## Residual standard error: 0.4797 on 934 degrees of freedom
## ## Multiple R-squared:  0.006162,	Adjusted R-squared:  0.005098 
## ## F-statistic: 5.791 on 1 and 934 DF,  p-value: 0.0163
## summary( lm( pheno ~ PC1 * dataset, data = data ) )
## ##                   Estimate Std. Error t value Pr(>|t|)    
## ## (Intercept)        1.67118    0.02105  79.382   <2e-16 ***
## ## PC1               -0.65844    0.47868  -1.376   0.1693    
## ## datasetCureGN     -0.07452    0.03147  -2.368   0.0181 *  
## ## PC1:datasetCureGN -0.31548    0.67696  -0.466   0.6413    
## ## Residual standard error: 0.4787 on 932 degrees of freedom
## ## Multiple R-squared:  0.01234,	Adjusted R-squared:  0.009157 
## ## F-statistic:  3.88 on 3 and 932 DF,  p-value: 0.008991
## summary( lm( pheno ~ PC1, data = data %>% filter( dataset == 'Bristol' ) ) )
## ##             Estimate Std. Error t value Pr(>|t|)    
## ## (Intercept)  1.67118    0.02066  80.883   <2e-16 ***
## ## PC1         -0.65844    0.46980  -1.402    0.162    
## ## Residual standard error: 0.4698 on 515 degrees of freedom
## ## Multiple R-squared:  0.0038,	Adjusted R-squared:  0.001865 
## ## F-statistic: 1.964 on 1 and 515 DF,  p-value: 0.1617
## summary( lm( pheno ~ PC1, data = data %>% filter( dataset != 'Bristol' ) ) )
## ##             Estimate Std. Error t value Pr(>|t|)    
## ## (Intercept)  1.59666    0.02391   66.78   <2e-16 ***
## ## PC1         -0.97392    0.48942   -1.99   0.0473 *  
## ## Residual standard error: 0.4894 on 417 degrees of freedom
## ## Multiple R-squared:  0.009407,	Adjusted R-squared:  0.007031 
## ## F-statistic:  3.96 on 1 and 417 DF,  p-value: 0.04725
## summary( lm( pheno ~ PC1 : (dataset == 'Bristol'), data = data ) )
## ##                               Estimate Std. Error t value Pr(>|t|)    
## ## (Intercept)                    1.63782    0.01568 104.421   <2e-16 ***
## ## PC1:dataset == "Bristol"FALSE -0.97392    0.47986  -2.030   0.0427 *  
## ## PC1:dataset == "Bristol"TRUE  -0.65844    0.47986  -1.372   0.1704    
## ## Residual standard error: 0.4799 on 933 degrees of freedom
## ## Multiple R-squared:  0.006392,	Adjusted R-squared:  0.004262 
## ## F-statistic: 3.001 on 2 and 933 DF,  p-value: 0.05022
## summary( lm( pheno ~ age, data = data ) )
## ##              Estimate Std. Error t value Pr(>|t|)    
## ## (Intercept) 1.6102029  0.0213366  75.467   <2e-16 ***
## ## age         0.0021137  0.0008267   2.557   0.0107 *  
## ## Residual standard error: 0.4767 on 919 degrees of freedom
## ##   (15 observations deleted due to missingness)
## ## Multiple R-squared:  0.007063,	Adjusted R-squared:  0.005982 
## ## F-statistic: 6.537 on 1 and 919 DF,  p-value: 0.01073
