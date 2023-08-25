library(GMMAT)
library(genio)
library(readr)

# this is where our data is
setwd( '/datacommons/ochoalab/ssns_gwas/imputed/' )

# constants
name <- "mac20"

# load phenotype and fixed covariates file
data <- read_tsv( 'patient-data.txt.gz', col_types = 'ccccdiiii' )
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
    ns_ctrl ~ sex + race + PCs,
    data = data,
    kins = grm$kinship,
    id = "id",
    family = binomial( link = "logit" )
)

# to restart here if needed, save model
save( model0, file = paste0( name, '-model0.RData' ) )

# perform GWAS!
message( 'glmm.score' )
glmm.score(
    model0,
    infile = name,
    outfile = paste0( name, "-glmm-score.txt" ),
    ncores = 92
)
