## prep array data for merging with TGP

# In between the previous block and here, you need to create `TGP-subset.{bed,bim,fam}`!  That in turn cannot be created until `ssns_gwas_maf_dedup` above is created.  See README.txt.

library(genio)
library(readr)

# since all are SNPs, complements are reverse complements (i.e. no need to reverse)
comp <- function(x) chartr("ATGC","TACG",x)

# on DCC we should be here
setwd( '/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/' )

# just need to compare annotations right now
arr <- read_bim( 'ssns_gwas_maf_dedup' )
tgp <- read_bim( 'TGP-subset' )

# keep an unedited copy of the original, to edit in a second pass
arr_orig <- arr

# confirm assumption that IDs are unique in each input
stopifnot( nrow( arr ) == length( unique( arr$id ) ) )
stopifnot( nrow( tgp ) == length( unique( tgp$id ) ) )
# TGP is a subset of array data by construction, confirm
stopifnot( all( tgp$id %in% arr$id ) )

# let's throw away SNPs only in array (highly unlikely to be real, or to be useful for imputation anyway)
indexes <- arr$id %in% tgp$id
ids_rm <- arr$id[ !indexes ]
length( ids_rm ) # 70337, fewer than previous report, but we're not done
# actually remove from working copy
arr <- arr[ indexes, ]
# because SNPs are ordered by chr:pos and there's no duplicates, tibbles should now be aligned, confirm!
stopifnot( all( arr$id == tgp$id ) )

nrow( arr ) # 762710

# by construction all columns match now except possible ref/alt, let's check those now
# many of these cases are mutually exclusive, but when ref/alt are reverse complements of each other then they can fit in multiple cases, the way it's coded those will tend to stay unchanged (we need to compare allele frequencies later)

# the transformations form a group!
# let's enumerate the equivalent cases for SNPs according to the two symmetries
# pairs from different cases are incompatible (error, throw away SNP)
# there are 4*3 = 12 cases total:

# ref alt
# A   C   # seed
# C   A   # flip
# T   G   # revcomp
# G   T   # revcomp + flip

# ref alt
# A   G   # seed
# G   A   # flip
# T   C   # revcomp
# C   T   # revcomp + flip

# ref alt
# A   T   # seed
# T   A   # flip, revcomp, revcomp + flip

# ref alt
# C   G   # seed
# G   C   # flip, revcomp, revcomp + flip

# there's no need to change rows where both ref and alt match, identify those and remove from working copies
# this is most of the data!
indexes <- arr$ref == tgp$ref & arr$alt == tgp$alt
mean( indexes ) # 0.8097022
arr <- arr[ !indexes, ]
tgp <- tgp[ !indexes, ]

nrow( arr ) # 145142

# now try flipping ref/alt
# this is the most common case left
indexes <- arr$alt == tgp$ref & arr$ref == tgp$alt
mean( indexes ) # 0.6521682
# remember how to edit this case
ids_flip <- arr$id[ indexes ]
arr <- arr[ !indexes, ]
tgp <- tgp[ !indexes, ]

nrow( arr ) # 50485

# two more possibilities, revcomp and that plus flip
arr$altcomp <- comp( arr$alt )
arr$refcomp <- comp( arr$ref )

# this is next most common case
indexes <- arr$refcomp == tgp$ref & arr$altcomp == tgp$alt
mean( indexes ) # 0.6330197
# remember how to edit this case
ids_revcomp <- arr$id[ indexes ]
arr <- arr[ !indexes, ]
tgp <- tgp[ !indexes, ]

nrow( arr ) # 18527

# last case, if these don't match then they are incompatible
indexes <- arr$altcomp == tgp$ref & arr$refcomp == tgp$alt
mean( indexes ) # 0.995628
# remember how to edit this case
ids_revcomp_flip <- arr$id[ indexes ]
arr <- arr[ !indexes, ]
tgp <- tgp[ !indexes, ]

nrow( arr ) # 81

# inspected what was left, agreed they are true incompatibilities, remove all of them!
ids_rm <- c( ids_rm, arr$id )

# final counts
length( ids_flip )         # 94657
length( ids_revcomp )      # 31958
length( ids_revcomp_flip ) # 18446
length( ids_rm )           # 70418

# write lists to pass to plink
# names match plink1 commands
# remove list is easy
write_lines( ids_rm, 'premerge-array-exclude.txt' )
# plink's "flip" is actually my revcomp, it doesn't changes X, just changes bim table!
write_lines( c( ids_revcomp, ids_revcomp_flip ), 'premerge-array-flip.txt' )

# my flips are actually changes of reference, in this case it's best to let plink handle it automatically, no new outputs needed (the TGP bim table will do).  This changes both bim and X, too messy for R I think when plink is way more reliable and fast!

#########################

# in this test file manually renamed B alleles to C so they're valid snps
## plink1 --keep-allele-order --bfile dummy-33-101-0.1 --flip flip.txt --make-bed --out test
# bim table is indeed revcomp'd:
## diff *.bim
## 6c6
## < 1	snp5	0	5	C	A
## ---
## > 1	snp5	0	5	G	T
## 11c11
## < 1	snp10	0	10	C	A
## ---
## > 1	snp10	0	10	G	T
## 16c16
## < 1	snp15	0	15	C	A
## ---
## > 1	snp15	0	15	G	T
# bed files are unchanged!
## diff -q *.bed

## # plink2 version
## --ref-allele TGP-subset.bim 6 2

## # plink2 test with ref change, no revcomp/flip
## plink2 --bfile dummy-33-101-0.1 --ref-allele newref.bim 6 2 --make-bed --out test2
## ## # bim files exactly as expected
## ## # one line diff
## ## diff dummy-33-101-0.1.bim test2.bim
## ## ## 2c2
## ## ## < 1	snp1	0	1	C	A
## ## ## ---
## ## ## > 1	snp1	0	1	A	C
## ## # identical
## ## diff -q newref.bim test2.bim
## ## diff -q dummy-33-101-0.1.bed test2.bed
## ## # Files dummy-33-101-0.1.bed and test2.bed differ
## # confirmed in R using BEDMatrix that only snp1 was changed as x <- 2-x, as desired!

## # included one revcomp first line, let's see if this works
## plink2 --bfile dummy-33-101-0.1 --ref-allele newref3.bim 6 2 --make-bed --out test3
## # ok that didn't work
## ## Warning: --ref-allele mismatch for biallelic variant 'snp0'.
