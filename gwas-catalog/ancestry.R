# load big table for discovery data only, to extract ancestry breakdown for each subanalysis

library(tidyverse)

# go where the data is
setwd( '/datacommons/ochoalab/ssns_gwas/imputed' )

# load the big table!
data <- read_tsv( 'patient-data.txt.gz' )

### NS ###

# report separately TGP and our internal cases and controls
x <- table( data$ancestry, data$ns_ctrl, data$dataset, useNA='i' )
## , ,  = array
##          0   1
##   afr  301 194 #
##   amr    0   0
##   eas    0   0
##   eur  274 212 #
##   sas  248 343 #
##   <NA> 226 183 #

## , ,  = tgp
##          0   1
##   afr  661   0 #
##   amr  347   0
##   eas  504   0
##   eur  503   0 #
##   sas  489   0 #
##   <NA>   0   0

# counts for array data only
rowSums( x[,, 1] )
## afr  amr  eas  eur  sas <NA> 
## 495    0    0  486  591  409 


### SSNS ###

# exclude all but SSNS or controls
x <- with( data[ !is.na(data$ssns_ctrl), ], table( ancestry, ssns_ctrl, dataset, useNA='i' ) )
## , , dataset = array
##         ssns_ctrl
## ancestry   0   1
##     afr  301 148 #
##     amr    0   0
##     eas    0   0
##     eur  274 136 #
##     sas  248 329 #
##     <NA> 226 112 #

## , , dataset = tgp
##         ssns_ctrl
## ancestry   0   1
##     afr  661   0 #
##     amr  347   0
##     eas  504   0
##     eur  503   0 #
##     sas  489   0 #
##     <NA>   0   0

# array only
rowSums( x[,, 1] )
## afr  amr  eas  eur  sas <NA> 
## 449    0    0  410  577  338 


### SRNS ###

# exclude all but SRNS or controls
x <- with( data[ !is.na(data$srns_ctrl), ], table( ancestry, srns_ctrl, dataset, useNA='i' ) )
## , , dataset = array
##         srns_ctrl
## ancestry   0   1
##     afr  301  41 #
##     amr    0   0
##     eas    0   0
##     eur  274  75 #
##     sas  248  14 #
##     <NA> 226  63 #

## , , dataset = tgp
##         srns_ctrl
## ancestry   0   1
##     afr  661   0 #
##     amr  347   0
##     eas  504   0
##     eur  503   0 #
##     sas  489   0 #
##     <NA>   0   0 #

# array only
rowSums( x[,, 1] )
## afr  amr  eas  eur  sas <NA> 
## 342    0    0  349  262  289 


### SSNS-SRNS ###

# exclude all but SSNS or SRNS
# only case without TGP data, which simplifies whole setup
x <- with( data[ !is.na(data$ssns_srns), ], table( ancestry, ssns_srns, useNA='i' ) )
##         ssns_srns
## ancestry   0   1
##     afr  148  41
##     eur  136  75
##     sas  329  14
##     <NA> 112  63

rowSums(x)
 ## afr  eur  sas <NA> 
 ## 189  211  343  175 
