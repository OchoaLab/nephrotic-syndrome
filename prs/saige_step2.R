library(SAIGE)
library(optparse) 

dir = "/datacommons/ochoalab/ssns_gwas/imputed/prs-new/base/"
name = "mac20"

# terminal inputs
option_list = list(
  make_option(c( "-d", "--diagnosis"), type = "character", default = 'ns_ctrl', 
              help = "Diagnosis subtype: ns_ctrl, ssns_ctrl, srns_ctrl, ssns_srns", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
# get values
disease_subtype <- opt$diagnosis # 'ns_ctrl'
print(disease_subtype)


# load saige inputs
#plinkFile= paste0(dir, 'imputed/', disease_subtype, "/", ancestry, '/mac20') 

plinkFile = paste0(dir, name)

GMMATmodelFile = "ssns_ctrl.rda"
varianceRatioFile = "ssns_ctrl.varianceRatio.txt"
SAIGEOutputFile ="saige_output.txt"


print( 'saige step 2')
print(plinkFile)
print(GMMATmodelFile)
print(varianceRatioFile)
print(SAIGEOutputFile)

SPAGMMATtest(bedFile=paste0(plinkFile, ".bed"),
             bimFile=paste0(plinkFile, ".bim"),
             famFile=paste0(plinkFile, ".fam"),
             AlleleOrder= 'alt-first',
             is_imputed_data=TRUE,
             #impute_method = opt$impute_method,
             GMMATmodelFile=GMMATmodelFile,
             varianceRatioFile=varianceRatioFile,
             SAIGEOutputFile=SAIGEOutputFile,
             is_output_moreDetails =TRUE,
             is_overwrite_output = TRUE,
             #SPAcutoff = opt$SPAcutoff, default 2
             is_Firth_beta = TRUE, # for binary traits
             #pCutoffforFirth = opt$pCutoffforFirth, # default 0.01
             LOCO = FALSE,
             min_MAF=0,
             min_MAC=0.5,
             max_missing = 1,
             dosage_zerod_cutoff = 0,
             dosage_zerod_MAC_cutoff = 0
)