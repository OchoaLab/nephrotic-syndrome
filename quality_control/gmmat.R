library(GMMAT)
library(genio)
library(readr)
library(optparse)    # for terminal options
library(SeqArray)    # convert gds
# constants
name <- "mac20"

############
### ARGV ###
############

# define options
option_list = list(
  make_option(c( "-d", "--diagnosis"), type = "character", default = 'ns_ctrl', 
              help = "Diagnosis subtype: ns_ctrl, ssns_ctrl, srns_ctrl, ssns_srns", metavar = "character"),
  make_option(c( "-a", "--ancestry"), type = "character", default = NA, 
              help = "Ancestry subanalysis: NA, afr, eur, sas", metavar = "character"),
  make_option(c( "-c", "--numberofcore"), type = "numeric", default = 30,
              help = "Number of cores: numeric value", metavar = "numeric")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# get values
diagnosis_subtype <- opt$diagnosis # 'ns_ctrl'
ancestry_subtype <- opt$ancestry 


#ancestry_subtype <- opt$ancestry # 'eur' # set to NA if not ancestry subanalysis
number_cores <- opt$numberofcore
message( 'diagnosis: ', diagnosis_subtype, '; ancestry: ', ancestry_subtype, "; cores: ", number_cores)

### ANALYSIS ###


# this is where our data is
setwd( '/datacommons/ochoalab/ssns_gwas/imputed/' )

# load phenotype and fixed covariates file
data <- read_tsv( 'patient-data.txt.gz', col_types = 'cccccdiiii' )

# go to desired subdirectories
if ( diagnosis_subtype != 'ns_ctrl' )
  setwd( diagnosis_subtype )
if ( ancestry_subtype != 'NA')
  setwd( ancestry_subtype )

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
  # TIFFANY:
  # - change ns_ctrl to the value of diagnosis_subtype
  # - exclude race from ancestry subanalyses
  
  if (ancestry_subtype != 'NA' ){
    factors <- c("sex", "PCs")
    message('factors: sex, PCs')
  } else {
    factors <- c("sex", "race", "PCs")
    message('factors: sex, race, PCs')
  }
  
  model0 <- glmmkin(
    as.formula(paste(diagnosis_subtype, "~", paste(factors, collapse = "+"))),
    #ns_ctrl ~ sex + race + PCs,
    data = data,
    kins = grm$kinship,
    id = "id",
    family = binomial( link = "logit" )
  )
  
  # to restart here if needed, save model
  save( model0, file = file_data )
}

# make GDS or whatever version
bed.fn <- paste0(name, ".bed")
fam.fn <- paste0(name, ".fam")
bim.fn <- paste0(name, ".bim")
library(SeqArray)
seqBED2GDS(bed.fn = bed.fn, fam.fn = fam.fn, bim.fn = bim.fn, 
           out.gdsfn = paste0(name, ".gds"), verbose = TRUE)

# perform GWAS!
message( 'glmm.score' )
glmm.score(
  model0,
  # infile here is GDS
  infile = paste0( name, ".gds" ),
  outfile = paste0( name, "test-glmm-score.txt" ),
  ncores = number_cores
)