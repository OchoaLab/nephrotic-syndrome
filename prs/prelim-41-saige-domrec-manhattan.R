# look at newest "domrec + saige" results, make plots of significant chrs only

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

# most chromosomes don't have much going on, let's see
chrs_sig <- c()
for ( chr_i in 1 : 22 ) {
    # look at this chr only
    data_chr <- data %>% filter( chr == chr_i )
    # keep significant loci only, in at least one of the three tests
    data_chr_filt <- data_chr %>% filter( pa < pcut | pd < pcut | ( !is.na( pr ) & pr < pcut ) )
    # if we got any hits
    if ( nrow( data_chr_filt ) > 0 ) {
        # add to list
        chrs_sig <- c( chrs_sig, chr_i )
        # report number
        message( 'chr ', chr_i, ': ', nrow( data_chr_filt ) )
    }
}
## chr 6: 3652
## chr 9: 2
## chr 10: 3
## chr 15: 1
## chr 21: 22

# filter data to keep those chrs only
data2 <- data %>% filter( chr %in% chrs_sig )

# need to pivot-longer now, so each p-value is a different point to plot
data2 <- data2 %>% pivot_longer( cols = pa:pr, names_to = 'model', names_prefix = 'p', values_to = 'p', values_drop_na = TRUE )
# assign nicer names to models
data2$model[ data2$model == 'a' ] <- 'Additive'
data2$model[ data2$model == 'd' ] <- 'Dominant'
data2$model[ data2$model == 'r' ] <- 'Recessive'
data2$model <- factor( data2$model, levels = c('Additive', 'Dominant', 'Recessive') )

# keep current order, plots a first, then d, then r (the last is the "smallest" in -log10p, the first two are mixed)

# start plot-specific processing

# make chromosomes ordered factors
data2$chr <- factor( data2$chr, levels = chrs_sig )

# code adapted from here to make nicer plots
# https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/
# finds cumulative length of chromosomes
data_cum <- data2 |>
    group_by( chr ) |>
    summarise( max_pos = max( pos ) ) |>
    mutate( pos_add = lag( cumsum( as.numeric( max_pos ) ), default = 0 ) ) |>
    select( chr, pos_add )

# calculates shifted position `x` nicer for plot
data2 <- data2 |>
    inner_join( data_cum, by = "chr" ) |>
    mutate( x = pos + pos_add )

# center of each chromosome, for x-axis labels
axis_set <- data2 |>
    group_by( chr ) |>
    summarize( center = mean( x ) )

# actual plot!
wh <- fig_scale( 3/(2*1.5) )
fig_start( 'domrec-saige', width = wh[1], height = wh[2]  )
ggplot( data2, aes( x = x, y = -log10( p ), color = model ) ) +
    geom_hline(
        yintercept = -log10( 5e-8 ), color = "grey40",
        linetype = "dashed"
    ) +
    geom_point() +
    scale_x_continuous(
        label = axis_set$chr,
        breaks = axis_set$center
    ) +
    theme_classic() +
    labs(
        x = 'Chromosome, Position',
        y = expression(-log[10](p)),
        color = 'Model'
    )
#    facet_wrap( vars( type ), nrow = 4 )
fig_end()
