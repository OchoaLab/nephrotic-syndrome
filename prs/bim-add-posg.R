library(genio)
library(ochoalabtools)
library(simfam)

# script adds genetic distance, useful for informing relevant/dynamic window sizes for LD calculations

# support old data for now, expect ssns_ctrl or ssns_srns
args <- args_cli()
name <- args[1]
hgv <- args[2]

# process version first because it's faster and more limited
if ( is.na( hgv ) )
    stop( 'Usage: <name> <hg version>' )
if ( ! hgv %in% c(37, 38) )
    stop( 'Human genome version (second argument) must be 37 or 38, passed: ', hgv )
map <- if ( hgv == 38 ) recomb_map_hg38 else recomb_map_hg37

# load variant table to edit
bim <- read_bim( name )

# calculate and add posg column:
bim <- bim_add_posg( bim, map )

# instead of overwriting, make copy of original to allow for undoing
# assumes name is lacking extension, won't work otherwise!!!
file.rename( paste0( name, '.bim' ), paste0( name, '.bim~' ) )

# now write new output
write_bim( name, bim )
