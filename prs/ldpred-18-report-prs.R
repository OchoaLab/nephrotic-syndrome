library(tidyverse)
library(genio)

# right now only this case is of interest
setwd( 'discovery' )

# read individual info
fam <- read_fam( 'mac20' )
# only IDs are of interest, in fact the rest of this info is trivial, so remove
fam <- fam %>% select( id )

# read scores of interest
# these are already aligned to fam, so just merge
fam$prs_ldpred <- read_lines( 'prs-base-train-ldpred2-grid-h0.1-best.txt.gz' )
fam$prs_ct <- read_lines( 'prs-base-train-ldpred2-ct-best.txt.gz' )

# just save this now!
write_tsv( fam, 'prs-best.txt.gz' )
