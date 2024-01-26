library(GMMAT)
library(genio)
library(readr)
library(SeqArray)    # convert gds

# constants
name <- "mac20"
diagnosis_subtype <- 'ssns_ctrl'
covariates <- c("sex", "race", "PCs")
number_cores <- 20

# this version only runs on base data, so no options vs earlier cases

### ANALYSIS ###

# this is where our data is
setwd( '/datacommons/ochoalab/ssns_gwas/imputed/prs-new/base/' )

# load phenotype and fixed covariates file
data <- read_tsv( '../../patient-data.txt.gz', col_types = 'cccccdciiii' )

file_data <- paste0( name, '-model0.RData' )

if ( file.exists( file_data ) ) {
    # load precalculated "model0" from previous session
    load( file_data ) # loads model0
} else {
    # calculate model0 for the first time
    
    # load GRM and PCs
    grm <- read_grm( name )
    eigenvec <- read_eigenvec( name )
    
    # take for granted that grm and eigenvec are aligned to each other, both came from genetic data
    # phenotype/covars are in different order and may contain extras (in more general case of subanalyses), so subset and reorder now
    data <- data[ match( grm$fam$id, data$id ), ]
    # now that data is aligned, include all PCs as a single, convenient covariate
    data$PCs <- eigenvec$eigenvec
    
    # fit null model
    message( 'glmmkin' )

    model0 <- glmmkin(
        as.formula( paste( diagnosis_subtype, "~", paste( covariates, collapse = "+" ) ) ),
        data = data,
        kins = grm$kinship,
        id = "id",
        family = binomial( link = "logit" )
    )
    
    # to restart here if needed, save model
    save( model0, file = file_data )
}

# make GDS version if needed
file_gds <- paste0(name, ".gds")
if ( !file.exists( file_gds ) ) {
    message( 'seqBED2GDS' )
    seqBED2GDS(
        bed.fn = paste0(name, ".bed"),
        fam.fn = paste0(name, ".fam"),
        bim.fn = paste0(name, ".bim"),
        out.gdsfn = file_gds,
        verbose = TRUE
    )
}

# perform GWAS!
message( 'glmm.score' )
glmm.score(
    model0,
    # infile here is GDS
    infile = file_gds,
    outfile = paste0( name, "-glmm-score.txt" ),
    ncores = number_cores
)
