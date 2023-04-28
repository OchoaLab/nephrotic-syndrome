
library(tidyverse)
library(dplyr)
library(ligera)
library(popkin)
library(BEDMatrix)
library(genio)
library(qqman)
library(vroom)
library(MASS)

name <- '/datacommons/ochoalab/ssns_gwas/imputation/merge/ssns_tgp'
X <- BEDMatrix( name )
bim <- read_bim( name )
fam <- read_fam( name )
phen <- read_phen('/hpc/home/tt207/NS_data/pheno_filtered.txt')

name_phen <- '/hpc/home/tt207/NS_data/pheno_filtered.txt'
if ( !is.na( name_phen ) ) {
  # if the phenotype is in a separate PHEN file, load that
  phen <- read_phen( name_phen )
  # reorder phen to match fam, if needed
  phen <- phen[ match( fam$id, phen$id ), ]
  # sanity check
  # when phen has fewer individuals than fam, some are NAs after match above, so this takes care of those cases
  stopifnot( all(phen$id == fam$id, na.rm = TRUE) )
  # save trait now
  trait <- phen$pheno
} else {
  # if the trait of interest is in the fam file, extract it now
  trait <- fam$pheno
}


ligera_f <- function(
  X,
  trait,
  kinship,
  kinship_inv = NULL,
  covar = NULL,
  loci_on_cols = FALSE,
  mem_factor = 0.7,
  mem_lim = NA,
  m_chunk_max = 1000,
  # cgsolve options
  tol = 1e-15, # default 1e-6
  maxIter = 1e6 # default 1e3
) {
  # - supports missingness in trait (exact kinship matrix inverse in those cases)
  # TODO
  # - support true missingness in genotypes (inversion of matrix subsets, etc; right now only approximate)
  
  # informative errors when things are missing
  if ( missing( X ) )
    stop( 'Genotype matrix `X` is required!' )
  if ( missing( trait ) )
    stop( '`trait` is required!' )
  if ( missing( kinship ) )
    stop( '`kinship` is required!' )
  # function from popkin validates further (includes square matrix test)
  popkin::validate_kinship( kinship )
  
  # and even further, as unexpected NAs are a pain
  if ( anyNA( kinship ) )
    stop( '`kinship` must not have any missing values!' )
  if ( !is.null( kinship_inv ) && anyNA( kinship_inv ) )
    stop( '`kinship_inv` must not have any missing values!' )
  # NOTE: trait and X may have missing values
  
  # override this for BEDMatrix
  if ( 'BEDMatrix' %in% class(X) ) {
    loci_on_cols <- TRUE
  } else if (!is.matrix(X))
    stop('X has unsupported class: ', toString( class( X ) ) )
  
  # need these dimensions
  if (loci_on_cols) {
    m_loci <- ncol(X)
    n_ind <- nrow(X)
  } else {
    m_loci <- nrow(X)
    n_ind <- ncol(X)
  }
  
  # check dimensions of other items
  if ( length( trait ) != n_ind )
    stop('Number of individuals in `trait` (', length( trait ), ') does not match genotype matrix (', n_ind , ')')
  if ( nrow( kinship ) != n_ind )
    stop('Number of individuals in `kinship` (', nrow( kinship ), ') does not match genotype matrix (', n_ind , ')')
  if ( !is.null( kinship_inv ) && nrow( kinship_inv ) != n_ind )
    stop('Number of individuals in `kinship_inv` (', nrow( kinship_inv ), ') does not match genotype matrix (', n_ind , ')')
  if ( !is.null( covar ) ) {
    if ( nrow( covar ) != n_ind )
      stop('Number of individuals in `covar` (', nrow( covar ), ') does not match genotype matrix (', n_ind , ')')
  }
  
  # update kinship, etc, if the trait has missing values
  # the good thing is that this is shared across loci, so comparably it's not so expensive
  # this NULL means there are no filters to apply
  indexes_ind <- NULL
  if ( anyNA( trait ) ) {
    # indexes to keep (need to subset genotypes at load time)
    indexes_ind <- !is.na( trait )
    # subset trait
    trait <- trait[ indexes_ind ]
    # subset kinship matrix
    kinship <- kinship[ indexes_ind, indexes_ind ]
    # force recomputing inverse of kinship matrix (see further below)
    kinship_inv <- NULL
    # subset covariates, if present
    if ( !is.null( covar ) )
      covar <- covar[ indexes_ind, ]
    # reduce number of individuals, used in some calculations
    n_ind <- length( trait )
    # NOTE: only genotypes are left to filter with indexes_ind
  }
  
  # gather matrix of trait, intercept, and optional covariates
  Y1 <- cbind( trait, 1 )
  # add covariates, if present
  if ( !is.null( covar ) ) {
    # handle NAs now, so final Y has no missingness whatsoever
    covar <- covar_fix_na( covar )
    Y1 <- cbind( Y1, covar )
  }
  # compute inverse if needed
  if ( is.null( kinship_inv ) ) {
    # initialize matrix to fill
    Z1 <- matrix( NA, nrow = nrow(Y1), ncol = ncol(Y1) )
    for ( k in 1:ncol(Y1) ) {
      Z1[ , k ] <- drop( cPCG::cgsolve( kinship, Y1[ , k ], tol = tol, maxIter = maxIter ) )
    }
  } else {
    # use kinship inverse if given
    Z1 <- kinship_inv %*% Y1
  }
  
  # need null model too
  # just remove first column of these two
  Y0 <- Y1[ , -1 ]
  Z0 <- Z1[ , -1 ]
  
  # calculate other intermediate parts
  # all have dimensions n x k
  H1 <- Z1 %*% solve( crossprod( Y1, Z1 ) )
  H0 <- Z0 %*% solve( crossprod( Y0, Z0 ) )
  proj1 <- H1[ , 1 ] # projection for trait coefficient only
  # another recurrent product for getting residuals/stats
  # dimensions n x n
  R1 <- tcrossprod( Y1, H1 ) - diag( n_ind )
  R0 <- tcrossprod( Y0, H0 ) - diag( n_ind )
  
  # last thing is matrix that returns sums of residuals quickly
  # dimensions n x n
  # NOTE: here conjugate gradient probably isn't very efficient, as products are as big as explicit inverse :(  Better luck computing residuals more directly from genotypes!
  if ( is.null( kinship_inv ) ) {
    # initialize matrix to fill
    SSR1 <- matrix( NA, nrow = n_ind, ncol = n_ind )
    SSR0 <- matrix( NA, nrow = n_ind, ncol = n_ind )
    # first do one of these products
    for ( k in 1 : n_ind ) {
      SSR1[ , k ] <- drop( cPCG::cgsolve( kinship, R1[ , k ], tol = tol, maxIter = maxIter ) )
      SSR0[ , k ] <- drop( cPCG::cgsolve( kinship, R0[ , k ], tol = tol, maxIter = maxIter ) )
    }
    # complete product
    SSR1 <- crossprod( R1, SSR1 )
    SSR0 <- crossprod( R0, SSR0 )
  } else {
    # use kinship inverse if given
    SSR1 <- crossprod( R1, kinship_inv %*% R1 )
    SSR0 <- crossprod( R0, kinship_inv %*% R0 )
  }
  
  
  ##############################
  ### COEFFICIENT ESTIMATION ###
  ##############################
  
  # initialize output vectors
  # Do before get_mem_lim_m so free memory is accounted for properly
  beta <- vector('numeric', m_loci)
  f_stat <- vector('numeric', m_loci)
  df <- vector('numeric', m_loci)
  
  # this overcounts since there are logical branches, not all overlap, but meh seriously
  # as usual, this should be conservative
  # recall that ints count as 0.5, doubles as 1
  #
  # vec_m, int
  # # indexes_loci_chunk
  # vec_m, double
  # # drop( Xi %*% proj )
  # # res1, res0
  #
  # vec_n # int
  # # n_ind_no_NA
  #
  # mat_m_n int
  # # Xi
  # # M # unnamed
  # # ( Xi == 1 )
  # mat_m_n double
  # # Res1, Res0
  # add one more double copy of Xi, this happens in matrix operations, are only temporary unnamed matrices
  
  # estimating total memory usage in bytes
  data <- popkin:::solve_m_mem_lim(
    n = n_ind,
    m = m_loci,
    mat_m_n = 5,
    vec_m = 3.5,
    vec_n = 0.5,
    mem = mem_lim,
    mem_factor = mem_factor
  )
  m_chunk <- data$m_chunk
  # cap value to a nice performing value (very good speed, minimal memory)
  if ( m_chunk > m_chunk_max )
    m_chunk <- m_chunk_max
  
  # navigate chunks
  i_chunk <- 1 # start of first chunk (needed for matrix inputs only; as opposed to function inputs)
  while (TRUE) { # start an infinite loop, break inside as needed
    # this means all SNPs have been covered!
    if (i_chunk > m_loci)
      break
    
    # range of SNPs to extract in this chunk
    indexes_loci_chunk <- i_chunk : min(i_chunk + m_chunk - 1, m_loci)
    
    if ( !is.null( indexes_ind ) ) {
      # individuals get filtered here (indexes_ind; required when there's missingness in trait)
      if (loci_on_cols) {
        Xi <- t( X[ indexes_ind, indexes_loci_chunk, drop = FALSE ] ) # transpose for our usual setup
      } else {
        Xi <- X[ indexes_loci_chunk, indexes_ind, drop = FALSE ]
      }
    } else {
      # keep all individuals (is this faster in that case?)
      if (loci_on_cols) {
        Xi <- t( X[ , indexes_loci_chunk, drop = FALSE ] ) # transpose for our usual setup
      } else {
        Xi <- X[ indexes_loci_chunk, , drop = FALSE ]
      }
    }
    
    # to have good averages, we need the number of non-NA individuals per row
    n_ind_no_NA <- rowSums( !is.na(Xi) )
    # now we can turn all NAs to zeroes (as ints, lower mem)
    Xi[ is.na(Xi) ] <- 0L
    
    # the coefficients are simply the genotypes projected!
    # these are matrices though, a vector for every locus (entry for every covariate)
    # adjust for the NAs? (not sure if this is reasonable or not yet)
    # store trait coefficients
    beta[ indexes_loci_chunk ] <- drop( Xi %*% proj1 ) * n_ind / n_ind_no_NA
    # rest are for getting residuals
    # Xi %*% t(SSRX) are dims m_chunk x n
    # NOTE: here missing values are just not part of sums, so setting them to zero is perfectly fine
    rss1 <- rowSums( tcrossprod( Xi, SSR1 ) * Xi )
    rss0 <- rowSums( tcrossprod( Xi, SSR0 ) * Xi )
    # calculate F statistic now that all the parts are in place!
    # here NAs are accounted for in formula (to be normalized in the end!
    f_stat[ indexes_loci_chunk ] <- (rss0 - rss1) / rss1
    df[ indexes_loci_chunk ] <- n_ind_no_NA - ncol( Y1 )
    
    # update starting point for next chunk! (overshoots at the end, that's ok)
    i_chunk <- i_chunk + m_chunk
  }
  
  ################
  ### P-VALUES ###
  ################
  
  # normalize stats
  f_stat <- f_stat * df
  
  # Get p-values for F test!!! (should be exact)
  pval <- stats::pf( f_stat, 1, df, lower.tail = FALSE )
  
  # done, return quantities of interest (nice table!)
  return(
    tibble::tibble(
      pval = pval,
      beta = beta,
      f_stat = f_stat,
      df = df
    )
  )
}


kinship <- popkin( X)
kinship_inv <- ginv( kinship )
print('kinship finished')
tib_f <- ligera_f( X, trait, kinship, kinship_inv)
tib_f <- cbind( bim, tib_f )

file_out <- paste0( 'ligera_f_NS_tgp', '.txt' )
write_tsv( tib_f, file_out )


# tib_manplot = tib %>% select(SNP = id, CHR = chr, BP = pos, P = pval) %>% mutate(CHR = as.numeric(CHR))
# ligera_manplot = manhattan(tib_manplot, annotatePval = 0.01)
# 
# qq(tib_manplot$P)
# hist(tib_manplot$P)
# 
# plink_freq <- read_delim("/hpc/home/tt207/txa_files/plink2.afreq", delim = "\t",
#                          col_names = c("chr", "id", "ref", "alt", "alt_freqs", "obs_count"), skip = 1)
# 
# combined_data = merge(plink_freq, tib, by = c("chr", "id"), all.y = T) 
# ggplot(combined_data, aes(x=-log(pval), y=alt_freqs)) + geom_point()