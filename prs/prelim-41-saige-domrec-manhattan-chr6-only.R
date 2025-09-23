# look at newest "domrec + saige" results, make plots of significant chrs only
# this version for chr6 only is geared to use way less memory up-front, for a plot that needs a lot of tweaking

library(tidyverse)
library(ochoalabtools)

# constants
pcut <- 5e-8

# these are all "base" gwas runs for PRS devel
setwd( 'base' )

# load results from all models
add <- read_tsv( "saige_output.txt.gz", show_col_types = FALSE, col_select = c(MarkerID, CHR, POS, p.value) )
rec <- read_tsv( "mac20-rec_saige.txt.gz", show_col_types = FALSE, col_select = c(MarkerID, p.value) )
dom <- read_tsv( "mac20-dom_saige.txt.gz", show_col_types = FALSE, col_select = c(MarkerID, p.value) )

# keep chr6 only
add <- add %>% filter( CHR == 6 )
rec <- rec %>% filter( grepl( '^chr6:', MarkerID, perl = TRUE ) )
dom <- dom %>% filter( grepl( '^chr6:', MarkerID, perl = TRUE ) )

# merge all of these 
# yey these are already perfectly aligned!
stopifnot( all( add$MarkerID == dom$MarkerID ) )
# lazy merge for this special case
data <- bind_cols(
    add %>% rename( chr = CHR, pos = POS, id = MarkerID, pa = p.value ),
    dom %>% select( pd = p.value )
)

# now merge in recessive data, but keep all SNPs (don't care if they have missing rec)
data <- full_join(
    data,
    rec %>% rename( id = MarkerID, pr = p.value )
)

# need to pivot-longer now, so each p-value is a different point to plot
data <- data %>% pivot_longer( cols = pa:pr, names_to = 'model', names_prefix = 'p', values_to = 'p', values_drop_na = TRUE )
# assign nicer names to models
data$model[ data$model == 'a' ] <- 'Additive'
data$model[ data$model == 'd' ] <- 'Dominant'
data$model[ data$model == 'r' ] <- 'Recessive'
data$model <- factor( data$model, levels = c('Additive', 'Dominant', 'Recessive') )

# keep current order, plots a first, then d, then r (the last is the "smallest" in -log10p, the first two are mixed)

# zoom into chr6 HLA region only
data2 <- data %>% filter( pos >= 2.8e7 & pos <= 3.4e7 )
#6:28,510,120 to 33,480,577

# actual plot!
wh <- fig_scale( 2 )
fig_start( 'domrec-saige-hla', width = wh[1], height = wh[2]  )
ggplot( data2, aes( x = pos, y = -log10( p ), color = model ) ) +
    geom_hline(
        yintercept = -log10( 5e-8 ), color = "grey40",
        linetype = "dashed"
    ) +
    geom_point() +
    theme_classic() +
    labs(
        x = 'Chr 6 Position',
        y = expression(-log[10](p)),
        color = 'Model'
    )
fig_end()
