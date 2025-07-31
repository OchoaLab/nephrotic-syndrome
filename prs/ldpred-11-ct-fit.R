library(ochoalabtools)
library(tidyverse)

# determine which type to run
args <- args_cli()
base <- args[1]
train <- args[2]
if ( is.na( train ) )
    stop( 'Usage: <base> <train>' )

# suffix shared by several in/out files
name_in <- paste0( base, '-ldpred2-ct' )
# check if output already exists, avoid recalculations
name_out <- paste0( 'eval-', name_in )
file_params <- paste0( name_out, '.txt.gz' )

if ( file.exists( paste0( train, '/', file_params ) ) ) {
    # load existing results
    setwd( train )
    params <- read_tsv( file_params, col_types = 'idiididdd' )
} else {
    library(bigsnpr)
    library(genio)

    # constants
    name_data <- 'mac20'
    NCORES <- 10
    
    # start at base
    setwd( base )

    # start by loading several previous calculations and other needed data

    # for outputs of second method that gets run together (not actually different input)
    name_in_stacked <- paste0( name_in, '-stacked' )

    # load previously calculated results
    load( paste0( 'all-keep-', name_in, '.RData' ) )
    # loads: all_keep
    # retrieve params from here
    # (params don't change between base and train)
    params <- attr( all_keep, "grid" )
    # rename because it will be modified (remapped) for train later
    all_keep_base <- all_keep
    rm( all_keep )

    # load filtered sumstats `df_beta_base`!
    df_beta_base <- read_tsv( paste0( 'betas-', base, '-clean-matched.txt.gz' ), show_col_types = FALSE )

    # rest always happens in training set
    setwd( paste0( '../', train ) )

    # load filtered sumstats `df_beta_train`!
    df_beta_train <- read_tsv( paste0( 'betas-', base, '-clean-matched.txt.gz' ), show_col_types = FALSE )

    # load training dataset
    # Attach the "bigSNP" object in R session
    obj.bigSNP <- snp_attach( paste0( name_data, '.rds' ) )
    G <- obj.bigSNP$genotypes
    y <- obj.bigSNP$fam$affection
    y[ y == 0 ] <- NA # in plink format, zeros are missing, translate appropriately here!
    # load PCs, to condition on when scoring with R2
    PCs <- read_eigenvec( name_data )$eigenvec

    # the best choice is to re-define both betas and lpvals using df_beta_train (it has all we need!)
    betas <- rep( NA, ncol( G ) )
    betas[ df_beta_train$`_NUM_ID_` ] <- df_beta_train$beta
    lpvals <- rep( NA, ncol( G ) )
    lpvals[ df_beta_train$`_NUM_ID_` ] <- -log10(df_beta_train$p)

    # only remaining challenge is to re-map kept SNPs from base to train indexes, or remove if they are no longer present
    # this is probably a memory-hungry solution, but it might be fastest and we can at least use to validate other algorithms
    # initialize whole vector for all possible positions we want to map
    index_map <- rep( NA, max( df_beta_base[["_NUM_ID_"]] ) )
    index_map[ df_beta_base$`_NUM_ID_`[ df_beta_train$`_NUM_ID_.ss` ] ] <- df_beta_train$`_NUM_ID_`

    # now apply map to all_keep_base!
    # fixed dimensions of this data, for loops and checks
    n_chr <- length( all_keep_base )
    n_params <- nrow( params )
    all_keep_train <- vector( 'list', n_chr )
    # keep some stats, since there's a worry about losing too many SNPs between datasets
    m_in <- vector( 'numeric', n_params )
    m_out <- vector( 'numeric', n_params )
    for ( chr in 1 : n_chr ) {
        all_keep_base_chr <- all_keep_base[[ chr ]]
        stopifnot( length( all_keep_base_chr ) == n_params )
        # initialize output structure for this chromosome
        all_keep_train_chr <- vector( 'list', n_params )
        for ( p in 1 : n_params ) {
            all_keep_base_chr_p <- all_keep_base_chr[[ p ]]
            stopifnot( !anyNA( all_keep_base_chr_p ) )
            # apply map!
            all_keep_train_chr_p <- index_map[ all_keep_base_chr_p ]
            # remove missing cases
            all_keep_train_chr_p <- all_keep_train_chr_p[ !is.na( all_keep_train_chr_p ) ]
            # store in output structure
            all_keep_train_chr[[ p ]] <- all_keep_train_chr_p
            # record removal stats
            m_in[ p ] <- m_in[ p ] + length( all_keep_base_chr_p )
            m_out[ p ] <- m_out[ p ] + length( all_keep_train_chr_p )
        }
        # add to bigger output structure
        all_keep_train[[ chr ]] <- all_keep_train_chr
    }

    # missingness stats
    # TODO: save to a file properly, and/or make a figure!
    ## 1 - m_out/m_in
    ##  [1] 0.4128010 0.3995251 0.3848091 0.3673955 0.3692467 0.3688268 0.3693886
    ##  [8] 0.3700783 0.3355775 0.3333294 0.3328007 0.3326910 0.2994635 0.2970717
    ## [15] 0.2960130 0.2958030 0.2440761 0.2434372 0.2427161 0.2424080 0.2232495
    ## [22] 0.2230807 0.2228071 0.2226956 0.2180747 0.2180062 0.2179124 0.2178683

    # calculate PRS for training individuals only
    # this pred_grid isn't a regular matrix, instead it is an FBM (matrix on disk), so watch out!  Can't save it as usual without some additional work anyway
    pred_grid <- snp_grid_PRS(
        G,
        all_keep_train,
        betas,
        lpvals,
        #ind.row = ind.train,
        #backingfile = "tmp-data/public-data-scores", 
        n_thr_lpS = 50,
        ncores = NCORES
    )
    ## dim( pred_grid ) # [1]   386 30800 # individuals by params?

    # add more info to params
    # this grows original table of 28 rows to 1400 (50x larger)
    params <- params %>%
        mutate( thr.lp = list(attr(pred_grid, "grid.lpS.thr")), id = row_number() ) %>%
        tidyr::unnest(cols = "thr.lp")
    # NOTE: pred_grid above is itself for 30800 params, which is 22x 1400, i.e. each chromosome appears separately!
    # get size of new table
    s <- nrow(params)

    # use training individuals to score grid values using correlation, adjusting for PCs, determine which is best
    # (apparently a.FUN in big_apply doesn't see any global variables, have to pass everything explicitly)
    data <- big_apply(
        pred_grid,
        # this setup sums over all chromosomes, for the same C+T parameters
        a.FUN = function( x, ind, s, y, n_chr, PCs, pcor ) {
            xs <- rowSums( x[, ind + s * (0:(n_chr-1))] )
            pcor( xs, y, PCs )
        },
        ind = 1 : s,
        s = s,
        y = y,
        n_chr = n_chr,
        PCs = PCs,
        pcor = pcor,
        a.combine = 'rbind',
        block.size = 1,
        ncores = NCORES
    )
    # incorporate back into tibble
    colnames( data ) <- c('cor', 'cor_lower', 'cor_upper')
    params <- bind_cols( params, as_tibble( data ) )

    # sort by correlation
    params <- params %>%
        arrange( desc( cor ) )

    # save this table of results!
    write_tsv( params, file_params )

    # pick out the best set of parameters, to use and score out of sample later!
    # first obtain the R2 set of SNPs (first threshold)
    loci_keep <- unlist( purrr::map( all_keep_train, params$id[1] ) )
    # then apply log p-value threshold (second threshold)
    loci_keep2 <- loci_keep[ lpvals[ loci_keep ] > params$thr.lp[1] ]
    ## # these are the coefficients in isolation, not the vector the same length as df_beta (like the rest)
    ## betas[loci_keep2]
    # this is starting vector we want, only difference is thresholds effectively set most of these to zero
    betas_best <- df_beta_train$beta
    # this sets all to zero except the best ones!
    betas_best[ ! df_beta_train$`_NUM_ID_` %in% loci_keep2 ] <- 0

    # save betas!
    # this is a simple vector
    file_out <- paste0( 'betas-', name_in, '-best.txt.gz' )
    write_lines( betas_best, file_out )


    ###################
    ### STACKED C+T ###
    ###################

    # fit a stacked model, as an alternative
    final_mod <- snp_grid_stacking( pred_grid, y-1, ncores = NCORES ) # , K = 4
    # extract betas to save
    # this is all of them
    betas_stacked <- final_mod$beta.G
    # subset to betas in data frame
    betas_stacked <- betas_stacked[ df_beta_train$`_NUM_ID_` ]
    file_out <- paste0( 'betas-', name_in_stacked, '.txt.gz' ) # no need for "best" here, there aren't other versions in any outputs
    write_lines( betas_stacked, file_out )

    # there's no other figures associated with stacked version, nothing worth showing anyway
}

# set all negatives to zero
params$cor[ params$cor < 0 ] <- 0
params$cor_lower[ params$cor_lower < 0 ] <- 0
params$cor_upper[ params$cor_upper < 0 ] <- 0

# plot results
fig_start( name_out, width = 6 )
ggplot( params, aes( x = 10^(-thr.lp), y = cor^2, color = as.factor( thr.r2 ) ) ) +
    theme_classic() +
    geom_point() +
    geom_line() +
    scale_x_log10() +
    geom_errorbar( aes( ymin = cor_lower^2, ymax = cor_upper^2 ), width = .5 ) +
    expand_limits( y = 0 ) + 
    labs( x = 'P-value threshold', y = expression(R^2 * " to trait"), color = "LD threshold" )
fig_end()

