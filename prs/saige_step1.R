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

plinkFile = paste0(dir, name)
phenoFile= "covar_ssns_ctrl.txt"
outputPrefix= "ssns_ctrl"
# covariates for main + conditional 
covars=c('sex', "race",  "PCs.1","PCs.2","PCs.3","PCs.4","PCs.5","PCs.6","PCs.7","PCs.8","PCs.9","PCs.10")
# categorical
qcovars=c('sex', 'race')


print('start saige step1')
print(plinkFile)
print(phenoFile)
print(outputPrefix)
print(qcovars)

phenoCol= disease_subtype
sampleIDColinphenoFile='id' 
traitType='binary'        
IsOverwriteVarianceRatioFile=TRUE

fitNULLGLMM(plinkFile = plinkFile,
            phenoFile = phenoFile,
            phenoCol = phenoCol,
            sampleIDColinphenoFile = sampleIDColinphenoFile,
            traitType = traitType,
            outputPrefix = outputPrefix,
            covarColList = covars,
            qCovarCol = qcovars,
            IsOverwriteVarianceRatioFile = TRUE,
            LOCO = FALSE,
            sexCol = "sex",
            minMAFforGRM = 0,
            maxMissingRateforGRM = 1
)
