library(tidyverse)
library(writexl)

# start by splitting discovery data
setwd( '/datacommons/ochoalab/ssns_gwas/imputed/prs-new/' )

# load data
# this is discovery and bristol together
data <- read_tsv( '../../array/patient-data.txt.gz', show_col_types = FALSE )

# exclude controls
data <- data %>% filter( diagnosis != 'Control' )

# subset to few columns of interest
data <- data %>% select( !bristol )

# many of these samples don't have ages, even though they should
sum( is.na( data$age ) ) # [1] 116

# reorder by age, to identify high cases easily
data <- data %>% arrange( age )

# write to file to share with Rasheed
write_xlsx( data, 'patient-data-ages.xlsx' )
