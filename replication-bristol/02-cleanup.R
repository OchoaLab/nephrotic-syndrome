library(readr)
library(dplyr)
library(ochoalabtools)

# get arguments from command line
args <- args_cli()
# must have input and output
file_in <- args[1]
file_out <- args[2]
# complain if they are missing!
if ( is.na( file_in ) || is.na( file_out ) )
    stop( 'Usage: <file_in> <file_out>' )

# reload raw data
data <- read_tsv( file_in, show_col_types = FALSE )

# process again as in script 01-*.R, except keeping more data this time
# subset for simplicity to non-assoc stats, essentially keep BIM-like data
# purposefully list A2 before A2, aligns better visually with gnomad later
# (rsid is now optional, the random SNPs table doesn't have it!)
data <- select( data, CHR, POS, any_of( 'rsid' ), A2, A1 )
# there are many duplicates (because some were significant in several substudies), subset to unique
data <- distinct( data )
# resort by CHROM/POS, not strictly needed but it's just cleaner
data <- data[ order( data$CHR, data$POS ), ]
# add chr prefix for chromosomes, otherwise they're not matched to gnomad data
data$CHR <- paste0( 'chr', data$CHR )

# and gnomad data
gnomad <- read_tsv( file_out, show_col_types = FALSE )
# this one is already sorted correctly, and chr/pos should match `data`

# ref/alt appear to be aligned this way, pre-emptively align them by rename
data <- rename( data, REF = A2, ALT = A1 )

# easiest matching concatenates chr:pos:ref:alt
data$chrposrefalt <- paste0( data$CHR, ':', data$POS, ':', data$REF, ':', data$ALT )
gnomad$chrposrefalt <- paste0( gnomad$CHROM, ':', gnomad$POS, ':', gnomad$REF, ':', gnomad$ALT )

# there's one SNP that wasn't found directly on gnomad (i.e. matching chrposrefalt perfectly):
index_unmatched <- which( !( data$chrposrefalt %in% gnomad$chrposrefalt ) )
# inspection confirms that it is simply missing
## data[ index_unmatched, ]
## ##   CHR       POS rsid         REF   ALT   chrpos      
## ## 1 chr4  9172304 rs1319121362 T     C     chr4:9172304
## gnomad[ gnomad$CHROM == 'chr4', ] 
## ##   CHROM       POS ID       REF   ALT   AC_nfe AN_nfe AC_afr AN_afr AC_sas AN_sas
## ## 1 chr4  132276787 rs18749… G     A        648  68036     48  41438      2   4836
## ## 2 chr4  136790863 rs62313… T     C       2009  67984   1782  41430    178   4830
## ## 3 chr4  167203980 rs18493… T     C        686  68028     41  41450     13   4834
## ## 4 chr4  175927544 rs78857… G     T       1562  68040    178  41460     55   4836
## ## 5 chr4  176069183 rs14637… C     A       1873  67972    213  41382     72   4830
# this unmatched SNP truly is missing in gnomad (also searched by rsID)
# it was not actually significant in ssns_control_original (1.34E-05), but it rose to near-gwas significant in ssns_control_conditional_4 (4.47E-08)
# if it has LD buddies, maybe we can replace it with one of them?

# anyway, let's report deletion, then remove and move on
if ( length( index_unmatched ) > 0 ) {
    message( 'Removing loci unmatched in gnomad:' )
    print( data[ index_unmatched, ] )
    data <- data[ -index_unmatched, ]
}

# out of curiosity, let's look at the extras in gnomad
## index_extras_gnomad <- which( !( gnomad$chrposrefalt %in% data$chrposrefalt ) )
## View( gnomad[ index_extras_gnomad, ] )
## View( gnomad[ -index_extras_gnomad, ] )
# won't describe thoroughly, but inspection shows that:
# - excluded cases are generally rarer than included
# - SNPs can be multiallelic
# - other variants just overlap, but are unrelated (some are common, but who cares)
# - only other common excluded cases are weird (like kept is TA-T but excluded is T-TA, which appear the same except thinking about it I see they are not, one is a deletion of the second A, the second is an insertion of another A), but again in all cases we had a match in one directly already so meh, just have to assume it's the right one

# overall, I'm convinced that the 1-1 matching is the right one, so let's just do that
indexes <- match( data$chrposrefalt, gnomad$chrposrefalt )
# subset and reorder gnomad!
gnomad <- gnomad[ indexes, ]
# remove extra columns from gnomad
gnomad <- select( gnomad, -chrposrefalt )
# write gnomad data back out!
# just overwrite table, I'm pretty sure we exclusively removed junk
write_tsv( gnomad, file_out )
