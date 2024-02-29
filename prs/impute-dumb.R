library(genio)
library(ochoalabtools)

# get name from terminal
name <- args_cli()[1]

# assume file will be huge, so save some time by loading X only

# get dimensions from bim, fam
m <- count_lines( name, 'bim' )
n <- count_lines( name, 'fam' )

# load entire genotype matrix
X <- read_bed( name, m_loci = m, n_ind = n )

# move old file to replace it with new one
file.rename( paste0( name, '.bed' ), paste0( name, '_ORIG.bed' ) )

# begin dumb imputation
for ( i in 1 : m ) {
    # skip rows without missingness (most of them)
    if ( !anyNA( X[i,] ) )
        next
    # copy down to simplify the rest
    xi <- X[i,]
    # identify individuals with missingness at this SNP
    indexes <- is.na( xi )
    # get AF
    af <- mean( xi, na.rm = TRUE ) / 2
    # draw random replacements, store immediately
    X[ i, indexes ] <- rbinom( sum( indexes ), 2, af )
}

# write new version
write_bed( name, X )
