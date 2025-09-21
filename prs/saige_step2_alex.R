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

SPAGMMATtest(
    bedFile = paste0( name, ".bed" ),
    bimFile = paste0( name, ".bim" ),
    famFile = paste0( name, ".fam" ),
    AlleleOrder = 'alt-first',
    is_imputed_data = TRUE,
    GMMATmodelFile = paste0( name, ".rda" ),
    varianceRatioFile = paste0( name, ".varianceRatio.txt" ),
    SAIGEOutputFile = paste0( name, "_saige.txt" ),
    is_output_moreDetails = TRUE,
    is_overwrite_output = TRUE,
    is_Firth_beta = TRUE, # for binary traits
    LOCO = FALSE,
    min_MAF = 0,
    min_MAC = 0.5,
    max_missing = 1,
    dosage_zerod_cutoff = 0,
    dosage_zerod_MAC_cutoff = 0
)
