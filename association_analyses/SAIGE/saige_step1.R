library(SAIGE)
plinkFile='/datacommons/ochoalab/ssns_gwas/imputed/srns_ctrl/mac20'
phenoFile='/datacommons/ochoalab/ssns_gwas/saige/srns_ctrl/covar_srns_ctrl.txt'
phenoCol='srns_ctrl'
covars=c('sex','race', "PCs.1","PCs.2","PCs.3","PCs.4","PCs.5","PCs.6","PCs.7","PCs.8","PCs.9","PCs.10")
# categorical
qcovars=c('sex', 'race')
sampleIDColinphenoFile='id' 
traitType='binary'        
outputPrefix='/datacommons/ochoalab/ssns_gwas/saige/srns_ctrl/srns_ctrl_binary' 
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
            sexCol = "sex"
)
