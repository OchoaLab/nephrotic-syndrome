library(SAIGE)


plinkFile='/datacommons/ochoalab/ssns_gwas/imputed/srns_ctrl/mac20'
GMMATmodelFile='/datacommons/ochoalab/ssns_gwas/saige/srns_ctrl/srns_ctrl_binary.rda'
varianceRatioFile='/datacommons/ochoalab/ssns_gwas/saige/srns_ctrl/srns_ctrl_binary.varianceRatio.txt'
SAIGEOutputFile = '/datacommons/ochoalab/ssns_gwas/saige/srns_ctrl/saige_output.txt'
  
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
             LOCO = FALSE
)