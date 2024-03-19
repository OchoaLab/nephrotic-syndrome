library(tidyverse)
library(genio)
# this is where our data is
dir = "/datacommons/ochoalab/ssns_gwas/imputed/prs-new/base/"
name = "mac20"

setwd(dir)

# load phenotype and fixed covariates file
data <- read_tsv( '/datacommons/ochoalab/ssns_gwas/imputed/patient-data.txt.gz' )
grm <- read_grm( name )
eigenvec <- read_eigenvec( name )

data <- data[ match( grm$fam$id, data$id ), ]
# now that data is aligned, include all PCs as a single, convenient covariate
data$PCs <- eigenvec$eigenvec

data_covar = data %>% select(ssns_ctrl, id, sex, race, PCs) %>% 
  mutate(sex = ifelse(sex == "female", 1, 0))
table(data_covar$ssns_ctrl)
table(data_covar$race)
table(data_covar$sex)
write.table(data_covar, "covar_ssns_ctrl.txt", 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = " ")

