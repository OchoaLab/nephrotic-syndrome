library(bigsnpr)
library(readr)
library(genio)
library(ochoalabtools)

# this script stratifies by ancesty, to confirm that accuracy is comparable across ancestries

# constants
name_data <- 'mac20'

args <- args_cli()
base <- args[1]
train <- args[2]
test <- args[3]
if ( is.na( test ) )
    stop( 'Usage: <base> <train> <test>' )

# load ancestry info as needed
if ( test == 'test' ) {
    # this is Bristol

    # shortcut function for loading subanalysis IDs
    read_ids <- function( anc ) 
        read_table( paste0( '../../replication/saige_results/bristol/AF/ns_ctrl/bristol_ns_', anc, '_allage.txt' ), col_types = 'cc', col_names = c('fam', 'id') )$id

    # NOTE: this main table doesn't have admixture results!  But it's base table over which we add this info...
    data <- read_tsv( '../../array/patient-data.txt.gz', show_col_types = FALSE )
    # let's remove all non-Bristol data
    data <- data[ data$bristol, ]
    # load subanalysis IDs
    ids_afr <- read_ids( 'afr' )
    ids_eur <- read_ids( 'eur' )
    ids_sas <- read_ids( 'sas' )
    # confirm they are disjoint
    stopifnot( !any( ids_afr %in% ids_eur ) )
    stopifnot( !any( ids_afr %in% ids_sas ) )
    stopifnot( !any( ids_eur %in% ids_sas ) )
    # check that all are indeed in Bristol
    stopifnot( all( ids_afr %in% data$id ) )
    stopifnot( all( ids_eur %in% data$id ) )
    stopifnot( all( ids_sas %in% data$id ) )
    # now massage data to look like previous cases... (Discovery and CureGN)
    # this is our new column with desired groupings; the default (cases not classified) are labeled as Other
    data$ancestry <- 'Other'
    # add labels
    data$ancestry[ data$id %in% ids_afr ] <- 'AFR'
    data$ancestry[ data$id %in% ids_eur ] <- 'EUR' 
    data$ancestry[ data$id %in% ids_sas ] <- 'SAS'
    # check correspondence with race; before (Discovery) they had to be of main race to also be of ancestry (i.e. to be EUR they had to be White), but here that's not true, admixture results were used regardless of race to group; here that seems like a distinct advantage because there are so many Unknowns that ended up being mostly unadmixed!  Results are as expected mostly I think:
    # - most AFR individuals are Black, but there are also Mixed/Unknown, and surprisingly, one Asian
    # - most EUR individuals are White, also Mixed/Unknown, but surprisingly one Asian (again)
    # - most SAS were Asian, also Mixed/Unknown, but surprisingly one Black and 7 White
    # table( data$race, data$ancestry )
    ##         Other AFR EUR SAS
    ## Asian       3   1   1  79
    ## Black       3  12   0   1
    ## Mixed       7   2   1   5
    ## Unknown    16   1  37  19
    ## White       6   0 389   7
    # a small proportion of individuals are Other, let's not label more finely
    # table( data$ancestry )
    ## Other     AFR     EUR     SAS 
    ##    35      16     428     111
    # since AFR is way too small to be useful, merge into Other...
    data$ancestry[ data$ancestry == 'AFR' ] <- 'Other'
} else {
    # all other cases so far are CureGN...
    data <- read_tsv( '../../../curegn/curegn_covar.txt', show_col_types = FALSE )
    # In theory I want to show as many results as possible, but sample size will limit 
    ## table( data$ancestry )
    ##   AFR   AFR_admix Asian_admix         EAS         EUR   EUR_admix 
    ##   168         108          11          88        1134         161 
    ## Other         SAS 
    ##   125          55 
    # merge smaller groups into Other, until they're more reasonable
    data$ancestry[ data$ancestry %in% c('Asian_admix', 'EAS', 'SAS', 'AFR_admix', 'EUR_admix') ] <- 'Other'
}

# all processing happens in subdirectory
setwd( test )
# combine base and train in new setup
base_train <- paste0( base, '-', train )

message( 'Loading testing dataset' )

# need to map SNPs from 'train' into 'test' here!
# load filtered sumstats `df_beta`!
file_df_beta <- paste0( 'betas-', base_train, '-clean-matched.txt.gz' )
df_beta <- read_tsv( file_df_beta, show_col_types = FALSE )

# load testing dataset
# Attach the "bigSNP" object in R session
obj.bigSNP <- snp_attach( paste0( name_data, '.rds' ) )
G <- obj.bigSNP$genotypes
y <- obj.bigSNP$fam$affection
y[ y == 0 ] <- NA # in plink format, zeros are missing, translate appropriately here!
# load PCs, to condition on when scoring with R2
PCs <- read_eigenvec( name_data )$eigenvec

# names of cases to score in testing data
names <- paste0( '-ldpred2-', c( 'inf-best', 'grid-h0.1-best', 'auto-h0.1', 'lassosum-best', 'ct-best', 'ct-stacked' ) )

# align ancestry annotations to genotypes; in all cases further subsetting has occurred, so re-evaluate sample sizes then
# in Bristol, only unclassified NS are removed
# in CureGN, most samples are removed because we only use MCD/FSGS, and in one case pediatric cases only
indexes <- match( obj.bigSNP$fam$sample.ID, data$id )
data <- data[ indexes, ]
# table( data$ancestry )
## ## test (Bristol)
## Other     AFR     EUR     SAS 
##    27      15     371     101
## ## test-curegn
## AFR   EUR Other
##  64   218   137

# process preexisting results
for ( name in names ) {
    # the following steps are shared by all ancestries
    
    # load PRS calculated previously from training
    # location of PRSs is training!
    file_in <- paste0( '../', train, '/betas-', base, name, '.txt.gz' )
    
    # skip silently if input is missing
    if ( !file.exists( file_in ) ) next
    # report what is being processed right now
    message( name )

    # load input
    betas <- as.numeric( read_lines( file_in ) )

    # subset betas using precalculated map of SNPs from 'train' into 'test'!
    betas <- betas[ df_beta[["_NUM_ID_.ss"]] ]

    # process each ancestry separately!
    for ( anc in unique( data$ancestry ) ) {
        # output now includes ancestry name
        file_out <- paste0( 'cor-', base_train, name, '-', anc, '.txt.gz' )
        # skip costly calculations if output already exists!
        if ( file.exists( file_out ) ) next

        # report what is being processed right now
        message( name, '-', anc )
        
        # determine indexes of subset of individuals to test
        indexes_inds <- which( data$ancestry == anc )
        
        # calculate PRS for test individuals now, requesting subset
        preds <- big_prodVec( G, betas, ind.col = df_beta[["_NUM_ID_"]], ind.row = indexes_inds )
        # calculate and save only correlation coefficient to truth, adjusting for PCs, subsetting appropriately
        cor <- pcor( preds, y[ indexes_inds ], PCs[ indexes_inds, ] )
        write_lines( cor, file_out )
    }
}
