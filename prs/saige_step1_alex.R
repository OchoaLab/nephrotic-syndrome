library(SAIGE)
library(optparse) 

# terminal inputs
option_list = list(
  make_option(c( "-d", "--diagnosis"), type = "character", default = 'ssns_ctrl', 
              help = "Diagnosis subtype: ns_ctrl, ssns_ctrl, srns_ctrl, ssns_srns", metavar = "character"),
  make_option("--bfile", type = "character", default = 'mac20', 
              help = "Base name of bed/bim/fam file", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
# get values
disease_subtype <- opt$diagnosis
name <- opt$bfile

fitNULLGLMM(
    plinkFile = name,
    phenoFile = "covar_ssns_ctrl.txt",
    phenoCol = disease_subtype,
    sampleIDColinphenoFile = 'id',
    traitType = 'binary',
    outputPrefix = name,
    covarColList = c('sex', "race",  "PCs.1","PCs.2","PCs.3","PCs.4","PCs.5","PCs.6","PCs.7","PCs.8","PCs.9","PCs.10"),
    qCovarCol = c('sex', 'race'),
    IsOverwriteVarianceRatioFile = TRUE,
    LOCO = FALSE,
    sexCol = "sex",
    minMAFforGRM = 0,
    maxMissingRateforGRM = 1
)
