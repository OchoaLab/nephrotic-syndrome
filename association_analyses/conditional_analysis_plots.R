library(tidyverse)
# european ancestry
ancestry_w <- read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/GMMAT/ssns_ctr/mac_20/glmm.score.bed.ancestry_white_ssns_ctr_20_PC.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
results_w = ancestry_w %>% select(SNP, CHR, BP = POS, P = PVAL, A1, A2) %>% drop_na(P) %>% filter(CHR == 6 & (BP>29000 | BP < 34000))

#european chr6 SNPs
results_w %>% filter(BP == 32663671) # HLA-DQB1
results_w %>% filter(BP == 32633919) # HLA-DQA1

ancestry_w_1 <- read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/GMMAT/ssns_ctr/mac_20/conditional_chr6/glmm.score.bed.ancestry_white_ssns_ctr_20_PC_chr6.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
results_w_1 = ancestry_w_1 %>% select(SNP, CHR, BP = POS, P = PVAL, A1, A2) %>% drop_na(P) %>% filter(CHR == 6 & (BP>29000 | BP < 34000))
ancestry_w_2 <- read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/GMMAT/ssns_ctr/mac_20/conditional_chr6/conditional_chr6_2/glmm.score.bed.ancestry_white_ssns_ctr_20_PC_chr6_2.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
results_w_2 = ancestry_w_2 %>% select(SNP, CHR, BP = POS, P = PVAL, A1, A2) %>% drop_na(P) %>% filter(CHR == 6 & (BP>29000 | BP < 34000))

png("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/GMMAT/ssns_ctr/mac_20/conditional_chr6/conditional_chr6_2/conditional_chr6_highlight_european.png", 
    width = 350, height = 250, units='mm', res = 300)
par( mfrow = c(3,1) )
manhattan(results_w,
          highlight = c("chr6:32663671:A:G","chr6:32633919:A:G"), 
          #%>% mutate(SNP = ifelse(BP == 32663671,"HLA-DQB1", "")), 
          annotateTop = FALSE, 
          #annotatePval = 0.1,  
          main = "European", xlim=c(3.24e7,3.28e7))
manhattan(results_w_1,
          highlight = "chr6:32633919:A:G",
          #highlight = c("chr6:32663671:A:G","chr6:32633919:A:G"),
          #%>% ,utate(SNP = ifelse(BP == 32597578,"HLA-DRB1", "")), annotatePval = 0.1, 
          annotateTop = FALSE, 
          #highlight = "chr6:32697820:A:G", 
          main = "Condition on HLA-DQB1",  xlim=c(3.24e7,3.28e7))
manhattan(results_w_2, #%>% mutate(SNP = ifelse(BP == 30963282,"NAPGP2", ifelse(BP == 32668717, "HLA-DQB1", ifelse(BP == 32440267, "HLA-DRA", "")))), 
          annotateTop = FALSE, 
          #highlight = "chr6:32633919:A:G",
          #annotatePval = 0.1,
          #highlight = highlight_sa, 
          main = "Condition on HLA-DQA1", xlim=c(3.24e7,3.28e7))
dev.off()


# south asian ancestry
ancestry_sa <- read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/GMMAT/ssns_ctr/mac_20/glmm.score.bed.ancestry_sa_ssns_ctr_20_PC.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
results_sa = ancestry_sa %>% select(SNP, CHR, BP = POS, P = PVAL, A1, A2) %>% drop_na(P) %>% filter(CHR == 6 & (BP>29000 | BP < 34000))

#southasian chr6 SNPs
results_sa %>% filter(BP == 32658788) # HLA DQB1
results_sa %>% filter(BP == 32714675) # AL662789.1

ancestry_sa_1 <- read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/GMMAT/ssns_ctr/mac_20/conditional_chr6/glmm.score.bed.ancestry_sa_ssns_ctr_20_PC_chr6.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
results_sa_1 = ancestry_sa_1 %>% select(SNP, CHR, BP = POS, P = PVAL, A1, A2) %>% drop_na(P) %>% filter(CHR == 6 & (BP>29000 | BP < 34000))
ancestry_sa_2 <- read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/GMMAT/ssns_ctr/mac_20/conditional_chr6/conditional_chr6_2/glmm.score.bed.ancestry_sa_ssns_ctr_20_PC_chr6_2.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
results_sa_2 = ancestry_sa_2 %>% select(SNP, CHR, BP = POS, P = PVAL, A1, A2) %>% drop_na(P) %>% filter(CHR == 6 & (BP>29000 | BP < 34000))

png("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/GMMAT/ssns_ctr/mac_20/conditional_chr6/conditional_chr6_2/conditional_chr6_highlight_southasian.png", 
    width = 350, height = 250, units='mm', res = 300)
par( mfrow = c(3,1) )
manhattan(results_sa,
          highlight = c("chr6:32658788:C:T", "chr6:32714675:G:A"), 
          annotateTop = FALSE, 
          main = "South Asian", xlim=c(3.24e7,3.28e7))
manhattan(results_sa_1,
          highlight = "chr6:32714675:G:A",
          annotateTop = FALSE, 
          main = "Condition on HLA-DQB1",  xlim=c(3.24e7,3.28e7))
manhattan(results_sa_2, 
          annotateTop = FALSE, 
          main = "Condition on AL662789.1", xlim=c(3.24e7,3.28e7))
dev.off()

# african ancestry
ancestry_b <- read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/GMMAT/ssns_ctr/mac_20/glmm.score.bed.ancestry_black_ssns_ctr_mac20_PC.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
results_b = ancestry_b %>% select(SNP, CHR, BP = POS, P = PVAL, A1, A2) %>% drop_na(P) %>% filter(CHR == 6 & (BP>29000 | BP < 34000))

#african chr6 SNPs
results_b %>% filter(BP == 32554331) # HLA DQB1
results_b %>% filter(BP == 32589312) # HLA-DRB6

ancestry_b_1 <- read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/GMMAT/ssns_ctr/mac_20/conditional_chr6/glmm.score.bed.ancestry_black_ssns_ctr_mac20_PC_chr6.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
results_b_1 = ancestry_b_1 %>% select(SNP, CHR, BP = POS, P = PVAL, A1, A2) %>% drop_na(P) %>% filter(CHR == 6 & (BP>29000 | BP < 34000))
ancestry_b_2 <- read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/GMMAT/ssns_ctr/mac_20/conditional_chr6/conditional_chr6_2/glmm.score.bed.ancestry_black_ssns_ctr_mac20_PC_chr6_2.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
results_b_2 = ancestry_b_2 %>% select(SNP, CHR, BP = POS, P = PVAL, A1, A2) %>% drop_na(P) %>% filter(CHR == 6 & (BP>29000 | BP < 34000))

png("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/GMMAT/ssns_ctr/mac_20/conditional_chr6/conditional_chr6_2/conditional_chr6_highlight_african.png", 
    width = 350, height = 250, units='mm', res = 300)
par( mfrow = c(3,1) )
manhattan(results_b,
          highlight = c("chr6:32554331:G:A","chr6:32589312:A:G"), 
          annotateTop = FALSE, 
          main = "African", xlim=c(3.24e7,3.28e7))
manhattan(results_b_1,
          highlight = "chr6:32589312:A:G",
          annotateTop = FALSE, 
          main = "Condition on HLA-DQB1",  xlim=c(3.24e7,3.28e7))
manhattan(results_b_2, 
          annotateTop = FALSE, 
          main = "Condition on HLA-DRB6", xlim=c(3.24e7,3.28e7))
dev.off()

# all ancestry
ancestry_all = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/GMMAT/ssns_ctr/full_data/glmm.score.bed.all_ssns_ctr_PC_20.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
results_all = ancestry_all %>% select(SNP, CHR, BP = POS, P = PVAL, A1, A2) %>% drop_na(P) %>% filter(CHR == 6 & (BP>29000 | BP < 34000))

#all ancestry chr6 SNP
results_all %>% filter(BP == 32684091) # HLA-DQB1
results_all %>% filter(BP == 32618190) # HLA-DQA1

ancestry_all_1 = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/GMMAT/ssns_ctr/full_data/conditional_chr6/glmm.score.bed.all_ssns_ctr_PC_20_chr6.txt", 
                            sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)
results_all_1 = ancestry_all_1 %>% select(SNP, CHR, BP = POS, P = PVAL, A1, A2) %>% drop_na(P) %>% filter(CHR == 6 & (BP>29000 | BP < 34000))
ancestry_all_2 = read.table("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/GMMAT/ssns_ctr/full_data/conditional_chr6/conditional_chr6_2/glmm.score.bed.all_ssns_ctr_PC_20_chr6_2.txt",
                            sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE, fill = TRUE)


x <- count.fields("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/GMMAT/ssns_ctr/full_data/conditional_chr6/conditional_chr6_2/glmm.score.bed.all_ssns_ctr_PC_20_chr6_2.txt", sep="\t")
which(x != 11)
unique(x)
results_all_2 = ancestry_all_2 %>% select(SNP, CHR, BP = POS, P = PVAL, A1, A2) %>% drop_na(P) %>% filter(CHR == 6 & (BP>29000 | BP < 34000))

png("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/GMMAT/ssns_ctr/mac_20/conditional_chr6/conditional_chr6_2/conditional_chr6_highlight_all.png", 
    width = 350, height = 250, units='mm', res = 300)
par( mfrow = c(3,1) )
manhattan(results_all,
          highlight = c("chr6:32684091:A:T", "chr6:32618190:A:G"),
          annotateTop = FALSE, 
          main = "African", xlim=c(3.24e7,3.28e7))
manhattan(results_all_1,
          highlight = "chr6:32618190:A:G",
          annotateTop = FALSE, 
          main = "Condition on HLA-DQB1",  xlim=c(3.24e7,3.28e7))
manhattan(results_all_2, 
          annotateTop = FALSE, 
          main = "Condition on HLA-DQA1", xlim=c(3.24e7,3.28e7))
dev.off()

png("/datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/GMMAT/ssns_ctr/mac_20/conditional_chr6_highlight.png", width = 400, height = 225, units='mm', res = 300)
par( mfrow = c(2,2) )
manhattan(results_all,
          highlight = c("chr6:32684091:A:T", "chr6:32618190:A:G"),
          #%>% mutate(SNP = ifelse(BP == 32684091,"HLA-DQB1", ifelse(BP == 30801705, "LINC00243", ifelse(BP == 30492876, "HLA-E", "")))),
          annotateTop = FALSE,
          #annotatePval = 0.1
          main = "SSNS vs Control: joint analysis",xlim=c(3.22e7,3.3e7))
manhattan(results_w,
          highlight = c("chr6:32663671:A:G","chr6:32633919:A:G"), 
          #%>% mutate(SNP = ifelse(BP == 32663671,"HLA-DQB1", "")), 
          annotateTop = FALSE, 
          #annotatePval = 0.1,  
          main = "European", xlim=c(3.22e7,3.3e7))
manhattan(results_b,
          highlight = c("chr6:32554331:G:A","chr6:32589312:A:G"),
          #%>% mutate(SNP = ifelse(BP == 32597578,"HLA-DRB1", "")), annotatePval = 0.1, 
          annotateTop = FALSE, 
          #highlight = "chr6:32697820:A:G", 
          main = "African",  xlim=c(3.22e7,3.3e7))
manhattan(results_sa, #%>% mutate(SNP = ifelse(BP == 30963282,"NAPGP2", ifelse(BP == 32668717, "HLA-DQB1", ifelse(BP == 32440267, "HLA-DRA", "")))), 
          annotateTop = FALSE, 
          highlight = c("chr6:32658788:C:T", "chr6:32714675:G:A"),
          #annotatePval = 0.1,
          #highlight = highlight_sa, 
          main = "South Asian", xlim=c(3.22e7,3.3e7))
dev.off()