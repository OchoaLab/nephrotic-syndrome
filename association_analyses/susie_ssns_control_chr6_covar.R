library("susieR")
library(tidyverse)

library(BEDMatrix)
X <- BEDMatrix("/datacommons/ochoalab/tiffany_data/susie/GMMAT_NS/ssns_ctr_hwe_mac20_chr6_29000-34000")
X <- as.matrix(X)
storage.mode(X) <- 'double'
rownames(X) <- NULL
colnames(X) <- NULL

phenofile = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/ssns_ctr_pheno.txt", sep = " ", stringsAsFactors=FALSE, quote = "", header = FALSE) %>% mutate(V3 = ifelse(V3 == 1, 0, 1))
covarfile = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/ssns_ctr_covar.txt", sep = " ", stringsAsFactors=FALSE, quote = "", header = FALSE) 
eigenfile = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/split_filter/ssns_ctr.eigenvec", sep = " ", stringsAsFactors=FALSE, quote = "", header = FALSE) 
Z_mat = merge(covarfile, eigenfile, by = c("V1", "V2")) %>% 
  select(-V1, -V2) %>% mutate(V3.x = ifelse(V3.x == "male", 0, 1),
                              V4.x = as.factor(V4.x) %>% unclass()) %>% as.matrix() 

Z_mat_numeric = matrix(as.numeric(Z_mat), ncol = ncol(Z_mat))
library(Matrix)
remove.covariate.effects <- function (X, Z, y) {
  # include the intercept term
  if (any(Z[,1]!=1)) Z = cbind(1, Z)
  A   <- forceSymmetric(crossprod(Z))
  SZy <- as.vector(solve(A,c(y %*% Z)))
  SZX <- as.matrix(solve(A,t(Z) %*% X))
  y <- y - c(Z %*% SZy)
  X <- X - Z %*% SZX
  return(list(X = X,y = y,SZy = SZy,SZX = SZX))
}

out = remove.covariate.effects(X = X, Z = Z_mat_numeric, y = phenofile$V3)
fitted_adjusted_all = susie(out$X, out$y, L = 10, verbose = TRUE)
susie_plot(fitted_adjusted_all, y="PIP", b=NULL, add_legend = "topleft", main = "SSNS vs Control: all ancestry - covariates")


L1 = print(fitted_adjusted_all$sets)$cs$L1
L2 = print(fitted_adjusted_all$sets)$cs$L2
L3 = print(fitted_adjusted_all$sets)$cs$L3
L4 = print(fitted_adjusted_all$sets)$cs$L4
L5 = print(fitted_adjusted_all$sets)$cs$L6
L6 = print(fitted_adjusted_all$sets)$cs$L10

library(genio)
name <- '/datacommons/ochoalab/tiffany_data/susie/GMMAT_NS/ssns_ctr_hwe_mac20_chr6_29000-34000'
X_bim <- read_bim( name )

L1_bim = X_bim[L1,]
print('L1')
print(L1_bim)
L2_bim = X_bim[L2,]
print('L2')
print(L2_bim)
L3_bim = X_bim[L3,]
print('L3')
print(L3_bim)
L4_bim = X_bim[L4,]
print('L4')
print(L4_bim)
L5_bim = X_bim[L5,]
print('L5')
print(L5_bim)
L6_bim = X_bim[L6,]
print('L6')
print(L6_bim)


# PIP dataframe
L1_df = cbind(L1_bim, fitted_adjusted_all$pip[L1], "L1") %>% mutate(PIP = fitted_adjusted_all$pip[L1]) %>% rename(L = `"L1"`) 
L1_df = L1_df[,-7]
L2_df = cbind(L2_bim, fitted_adjusted_all$pip[L2], "L2") %>% mutate(PIP = fitted_adjusted_all$pip[L2]) %>% rename(L = `"L2"`) 
L2_df = L2_df[,-7]
L3_df = cbind(L3_bim, fitted_adjusted_all$pip[L3], "L3") %>% mutate(PIP = fitted_adjusted_all$pip[L3]) %>% rename(L = `"L3"`) 
L3_df = L3_df[,-7]
L4_df = cbind(L4_bim, fitted_adjusted_all$pip[L4], "L4") %>% mutate(PIP = fitted_adjusted_all$pip[L4]) %>% rename(L = `"L4"`) 
L4_df = L4_df[,-7]
L5_df = cbind(L5_bim, fitted_adjusted_all$pip[L5], "L5") %>% mutate(PIP = fitted_adjusted_all$pip[L5]) %>% rename(L = `"L5"`)
L5_df = L5_df[,-7]
L6_df = cbind(L6_bim, fitted_adjusted_all$pip[L6], "L6") %>% mutate(PIP = fitted_adjusted_all$pip[L6]) %>% rename(L = `"L6"`) 
L6_df = L6_df[,-7]
L_df = rbind(L1_df, L2_df, L3_df, L4_df, L5_df, L6_df) %>% rename(CHR = chr, POS = pos, A1 = alt, A2 = ref) %>% arrange(POS)

GMMAT_ssns = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/GMMAT/ssns_ctr/full_data/glmm.score.bed.all_ssns_ctr_PC_20.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE) 
merge_L = merge(GMMAT_ssns, L_df, by.x = "SNP", by.y = "id", all.y = TRUE) %>% select(CHR = CHR.x, SNP, BP = POS.x, A1 = A1.x, A2 = A2.x, P = PVAL, PIP, L)
write.table(merge_L, "/datacommons/ochoalab/tiffany_data/susie/GMMAT_NS/ssns_ancestry/with_covar/merge_L_ssns.txt", sep = "\t", row.names=FALSE, quote=FALSE)

highlight_all = c(max(L1_df$id), max(L2_df$id), max(L3_df$id), max(L4_df$id), max(L5_df$id), max(L6_df$id))
manhattan(merge_L, annotateTop = FALSE, highlight = highlight_all, main = "SSNS vs Control - susie: plotting p-vals from GMMAT", cex = 0.5, xlim=c(2.9e7,3.4e7))
