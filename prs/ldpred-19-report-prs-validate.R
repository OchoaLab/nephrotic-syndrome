library(tidyverse)
library(ochoalabtools)

# want to combine these two
data <- read_tsv( 'patient-data.txt.gz' )
data2 <- read_tsv( 'prs-best.txt.gz' )

# merge tables now, overlap is complete so this is the same as other types of join
data <- inner_join( data, data2, by = 'id' )

fig_start( 'prs-best', width = 4 )
ggplot( data, aes( x = diagnosis, y = prs_ldpred ) ) +
    geom_boxplot() +
    theme_classic()
fig_end()
