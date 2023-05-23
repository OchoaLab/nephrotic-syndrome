library(tidyverse)
library(rio)
library(plyr)
setwd("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/")
# reading data from all sheets
excel_ssns_ctr <- import_list("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/box_files/ssnsVScontrol_suggestivesig_rsid.xlsx") 
excel_ns_ctr <- import_list("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/box_files/nsVScontrol_suggestivesig_rsid.xlsx") 
excel_srns_ctr <- import_list("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/box_files/srns_suggestivesig_rsid.xlsx") 
excel_meta <- import_list("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/box_files/metaanalysis_suggestivesig_rsid.xlsx") 
excel_bristol <- import_list("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/box_files/Bristol_replication_results.xlsx") 

# read all files in directory/newly calculated case-control AF
files <-  list.files(path = "/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/", pattern = ".cc")
for (i in seq_along(files)) {
    assign(paste(files[i]), read.table(files[i], header = TRUE) %>% 
             separate(SNP, c("chr", "POS", "a1", "a2")) %>% 
  select(-chr, -a1, -a2, CHR, POS, A1, A2, AF_case = MAF_A, AF_control = MAF_U, N_case = NCHROBS_A, N_control = NCHROBS_U))
}

# ssns
new_df_ssns_ctr = merge(excel_ssns_ctr$ssns_ctr, ssns_ctr_casecontrol.frq.cc, by = c("CHR", "POS", "A1", "A2")) %>%
  relocate(c("N_case", "AF_case", "N_control", "AF_control"), .after = AF) %>% arrange(PVAL) %>% select(-N, -AF)
new_df_ssns_ctr_cond1 = merge(excel_ssns_ctr$ssns_ctr_cond1, ssns_ctr_cond1_casecontrol.frq.cc, by = c("CHR", "POS", "A1", "A2")) %>%
  relocate(c("N_case", "AF_case", "N_control", "AF_control"), .after = AF) %>% arrange(PVAL) %>% select(-N, -AF)
new_df_ssns_ctr_cond2 = merge(excel_ssns_ctr$ssns_ctr_cond2, ssns_ctr_cond2_casecontrol.frq.cc, by = c("CHR", "POS", "A1", "A2")) %>%
  relocate(c("N_case", "AF_case", "N_control", "AF_control"), .after = AF) %>% arrange(PVAL) %>% select(-N, -AF)
new_df_ssns_ctr_cond3 = merge(excel_ssns_ctr$ssns_ctr_cond3, ssns_ctr_cond3_casecontrol.frq.cc, by = c("CHR", "POS", "A1", "A2")) %>%
  relocate(c("N_case", "AF_case", "N_control", "AF_control"), .after = AF) %>% arrange(PVAL) %>% select(-N, -AF)
new_df_ssns_ctr_cond4 = merge(excel_ssns_ctr$ssns_ctr_cond4, ssns_ctr_cond4_casecontrol.frq.cc, by = c("CHR", "POS", "A1", "A2")) %>%
  relocate(c("N_case", "AF_case", "N_control", "AF_control"), .after = AF) %>% arrange(PVAL) %>% select(-N, -AF)
# ssns - african
new_df_ssns_afr = merge(excel_ssns_ctr$ssns_ctr_afr, ssns_ctr_afr_casecontrol.frq.cc, by = c("CHR", "POS", "A1", "A2")) %>%
  relocate(c("N_case", "AF_case", "N_control", "AF_control"), .after = AF) %>% arrange(PVAL) %>% select(-N, -AF)
new_df_ssns_afr_cond1 = merge(excel_ssns_ctr$ssns_ctr_afr_cond1, ssns_ctr_afr_cond1_casecontrol.frq.cc,
                              by = c("CHR", "POS", "A1", "A2")) %>%
  relocate(c("N_case", "AF_case", "N_control", "AF_control"), .after = AF) %>% arrange(PVAL) %>% select(-N, -AF)
new_df_ssns_afr_cond2 = merge(excel_ssns_ctr$ssns_ctr_afr_cond2, ssns_ctr_afr_cond2_casecontrol.frq.cc,
                              by = c("CHR", "POS", "A1", "A2")) %>%
  relocate(c("N_case", "AF_case", "N_control", "AF_control"), .after = AF) %>% arrange(PVAL) %>% select(-N, -AF)
# ssns - european
new_df_ssns_eur = merge(excel_ssns_ctr$ssns_ctr_eur, ssns_ctr_eur_casecontrol.frq.cc, by = c("CHR", "POS", "A1", "A2")) %>%
  relocate(c("N_case", "AF_case", "N_control", "AF_control"), .after = AF) %>% arrange(PVAL) %>% select(-N, -AF)
new_df_ssns_eur_cond1 = merge(excel_ssns_ctr$ssns_ctr_eur_cond1, ssns_ctr_eur_cond1_casecontrol.frq.cc, by = c("CHR", "POS", "A1", "A2")) %>%
  relocate(c("N_case", "AF_case", "N_control", "AF_control"), .after = AF) %>% arrange(PVAL) %>% select(-N, -AF)
new_df_ssns_eur_cond2 = merge(excel_ssns_ctr$ssns_ctr_eur_cond2, ssns_ctr_eur_cond2_casecontrol.frq.cc, by = c("CHR", "POS", "A1", "A2")) %>%
  relocate(c("N_case", "AF_case", "N_control", "AF_control"), .after = AF) %>% arrange(PVAL) %>% select(-N, -AF)
# ssns - south asian
new_df_ssns_sas = merge(excel_ssns_ctr$ssns_ctr_sas, ssns_ctr_sas_casecontrol.frq.cc, by = c("CHR", "POS", "A1", "A2"))  %>%
  relocate(c("N_case", "AF_case", "N_control", "AF_control"), .after = AF) %>% arrange(PVAL) %>% select(-N, -AF)
new_df_ssns_sas_cond1 = merge(excel_ssns_ctr$ssns_ctr_sas_cond1, ssns_ctr_sas_cond1_casecontrol.frq.cc, by = c("CHR", "POS", "A1", "A2"))  %>%
  relocate(c("N_case", "AF_case", "N_control", "AF_control"), .after = AF) %>% arrange(PVAL) %>% select(-N, -AF)
new_df_ssns_sas_cond2 = merge(excel_ssns_ctr$ssns_ctr_sas_cond2, ssns_ctr_sas_cond2_casecontrol.frq.cc, by = c("CHR", "POS", "A1", "A2"))  %>%
  relocate(c("N_case", "AF_case", "N_control", "AF_control"), .after = AF) %>% arrange(PVAL) %>% select(-N, -AF)

options(java.parameters = "-Xmx2048m")
library(xlsx)
sheet_list = list(new_df_ssns_ctr, new_df_ssns_ctr_cond1, new_df_ssns_ctr_cond2, new_df_ssns_ctr_cond3, new_df_ssns_ctr_cond4,
                new_df_ssns_afr, new_df_ssns_afr_cond1, new_df_ssns_afr_cond2, new_df_ssns_eur, new_df_ssns_eur_cond1, new_df_ssns_eur_cond2,
                new_df_ssns_sas, new_df_ssns_sas_cond1, new_df_ssns_sas_cond2)

write.xlsx(as.data.frame(sheet_list[[1]]),
           file="/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/excel/ssnsVScontrol.xlsx",
           sheetName=paste(sheet_list[[1]]$source %>% head(1)), row.names=FALSE)
for (x in 2:14) {
  print(x)
  write.xlsx(as.data.frame(sheet_list[[x]]),
             file="/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/excel/ssnsVScontrol.xlsx",
             sheetName=paste(sheet_list[[x]]$source %>% head(1)), append=TRUE, row.names=FALSE)
}



# ns
new_df_ns = merge(excel_ns_ctr$ns_ctr, ns_ctr_casecontrol.frq.cc, by = c("CHR", "POS", "A1", "A2")) %>% 
     relocate(c("N_case", "AF_case", "N_control", "AF_control"), .after = AF) %>% arrange(PVAL) %>% select(-N, -AF)
new_df_ns_cond1 = merge(excel_ns_ctr$ns_ctr_cond1, ns_ctr_cond1_casecontrol.frq.cc, by = c("CHR", "POS", "A1", "A2")) %>% 
  relocate(c("N_case", "AF_case", "N_control", "AF_control"), .after = AF) %>% arrange(PVAL) %>% select(-N, -AF)
new_df_ns_cond2 = merge(excel_ns_ctr$ns_ctr_cond2, ns_ctr_cond2_casecontrol.frq.cc, by = c("CHR", "POS", "A1", "A2")) %>% 
  relocate(c("N_case", "AF_case", "N_control", "AF_control"), .after = AF) %>% arrange(PVAL) %>% select(-N, -AF)
new_df_ns_cond3 = merge(excel_ns_ctr$ns_ctr_cond3, ns_ctr_cond3_casecontrol.frq.cc, by = c("CHR", "POS", "A1", "A2")) %>% 
  relocate(c("N_case", "AF_case", "N_control", "AF_control"), .after = AF) %>% arrange(PVAL) %>% select(-N, -AF) %>% 
  dplyr::mutate(source = "ns_ctr_cond3")

# ns - african
new_df_ns_af = merge(excel_ns_ctr$ns_ctr_afr, ns_ctr_afr_casecontrol.frq.cc, by = c("CHR", "POS", "A1", "A2")) %>% 
  relocate(c("N_case", "AF_case", "N_control", "AF_control"), .after = AF) %>% arrange(PVAL) %>% select(-N, -AF)
new_df_ns_af_cond1 = merge(excel_ns_ctr$ns_ctr_afr_cond1, ns_ctr_afr_cond1_casecontrol.frq.cc, by = c("CHR", "POS", "A1", "A2")) %>% 
  relocate(c("N_case", "AF_case", "N_control", "AF_control"), .after = AF) %>% arrange(PVAL) %>% select(-N, -AF)
# ns - european
new_df_ns_eur = merge(excel_ns_ctr$ns_ctr_eur, ns_ctr_euro_casecontrol.frq.cc, by = c("CHR", "POS", "A1", "A2")) %>% 
  relocate(c("N_case", "AF_case", "N_control", "AF_control"), .after = AF) %>% arrange(PVAL) %>% select(-N, -AF)
new_df_ns_eur_cond1 = merge(excel_ns_ctr$ns_ctr_eur_cond1, ns_ctr_euro_cond1_casecontrol.frq.cc, by = c("CHR", "POS", "A1", "A2")) %>% 
  relocate(c("N_case", "AF_case", "N_control", "AF_control"), .after = AF) %>% arrange(PVAL) %>% select(-N, -AF)
# ns - south asian
new_df_ns_sas = merge(excel_ns_ctr$ns_ctr_sas, ns_ctr_sas_casecontrol.frq.cc, by = c("CHR", "POS", "A1", "A2")) %>% 
  relocate(c("N_case", "AF_case", "N_control", "AF_control"), .after = AF) %>% arrange(PVAL) %>% select(-N, -AF)
new_df_ns_sas_cond1 = merge(excel_ns_ctr$ns_ctr_sas_cond1, ns_ctr_sas_cond1_casecontrol.frq.cc, by = c("CHR", "POS", "A1", "A2")) %>% 
  relocate(c("N_case", "AF_case", "N_control", "AF_control"), .after = AF) %>% arrange(PVAL) %>% select(-N, -AF)


options(java.parameters = "-Xmx2048m")  
library(xlsx)
sheet_list = list(new_df_ns, new_df_ns_cond1, new_df_ns_cond2, new_df_ns_cond3,
                  new_df_ns_af, new_df_ns_af_cond1, new_df_ns_eur, new_df_ns_eur_cond1,
                  new_df_ns_sas, new_df_ns_sas_cond1)

write.xlsx(as.data.frame(sheet_list[[1]]),
           file="/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/excel/nsVScontrol.xlsx",
           sheetName=paste(sheet_list[[1]]$source %>% head(1)), row.names=FALSE)
for (x in 2:10) {
  print(x)
  print(as.character(paste(sheet_list[[x]]$source %>% head(1))))
  write.xlsx(as.data.frame(sheet_list[[x]]),
             file="/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/excel/nsVScontrol.xlsx",
             sheetName=as.character(paste(sheet_list[[x]]$source %>% head(1))), append=TRUE, row.names=FALSE)
}


# srns vs control
new_df_srns_ctr = merge(excel_srns_ctr$srns_ctr, srns_ctr_casecontrol.frq.cc, by = c("CHR", "POS", "A1", "A2")) %>% 
  relocate(c("N_case", "AF_case", "N_control", "AF_control"), .after = AF) %>% arrange(PVAL) %>% select(-N, -AF)
# ssns vs control
new_df_srns_ssns =merge(excel_srns_ctr$srns_ssns, srns_ssns_casecontrol.frq.cc, by = c("CHR", "POS", "A1", "A2")) %>% 
  relocate(c("N_case", "AF_case", "N_control", "AF_control"), .after = AF) %>% arrange(PVAL) %>% select(-N, -AF)

options(java.parameters = "-Xmx2048m")
library(xlsx)
sheet_list = list(new_df_srns_ctr, new_df_srns_ssns)

write.xlsx(as.data.frame(sheet_list[[1]]),
           file="/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/excel/srnsVScontrol.xlsx",
           sheetName=paste(sheet_list[[1]]$source %>% head(1)), row.names=FALSE)
write.xlsx(as.data.frame(sheet_list[[2]]),
             file="/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/excel/srnsVScontrol.xlsx",
             sheetName=paste(sheet_list[[2]]$source %>% head(1)), append=TRUE, row.names=FALSE)

# meta
new_df_ssns_meta = merge(excel_meta$meta_ssns, meta_ssns_ctr_casecontrol.frq.cc, by = c("CHR", "POS", "A1", "A2")) %>% 
  relocate(c("N_case", "AF_case", "N_control", "AF_control"), .after = AF) %>% arrange(PVAL) %>% select(-N, -AF)
# meta ns
new_df_ns_meta = merge(excel_meta$meta_ns, meta_ns_ctr_casecontrol.frq.cc, by = c("CHR", "POS", "A1", "A2")) %>% 
  relocate(c("N_case", "AF_case", "N_control", "AF_control"), .after = AF) %>% arrange(PVAL) %>% select(-N, -AF)

options(java.parameters = "-Xmx2048m")
library(xlsx)
sheet_list = list(new_df_ssns_meta, new_df_ns_meta)

write.xlsx(as.data.frame(sheet_list[[1]]),
           file="/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/excel/meta.xlsx",
           sheetName=paste(sheet_list[[1]]$source %>% head(1)), row.names=FALSE)
write.xlsx(as.data.frame(sheet_list[[2]]),
           file="/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/excel/meta.xlsx",
           sheetName=paste(sheet_list[[2]]$source %>% head(1)), append=TRUE, row.names=FALSE)

# bristol

new_df_bristol_ssns_ctr = merge(excel_bristol$`SSNS vs Control`, 
                                new_df_ssns_ctr %>% select(Chromosome = CHR, Position = POS, A1, A2, 
                                                           original_N_case = N_case, original_AF_case = AF_case, 
                                                           original_N_control = N_control, original_AF_control = AF_control), 
                                by = c("Chromosome", "Position", "A1", "A2")) %>% 
  relocate(c("original_N_case", "original_AF_case", "original_N_control", "original_AF_control"), .after = rsid) %>% 
  arrange(Bristol_GMMAT_PVAL) %>% select(-original_AF) %>% mutate(source = "SSNS vs Control")

new_df_bristol_ssns_ctr_euro = merge(excel_bristol$`SSNS vs Control (euro)`, 
                                new_df_ssns_eur %>% select(Chromosome = CHR, Position = POS, A1, A2, 
                                                           original_N_case = N_case, original_AF_case = AF_case, 
                                                           original_N_control = N_control, original_AF_control = AF_control), 
                                by = c("Chromosome", "Position", "A1", "A2")) %>% 
  relocate(c("original_N_case", "original_AF_case", "original_N_control", "original_AF_control"), .after = rsid) %>% 
  arrange(Bristol_GMMAT_PVAL) %>% select(-original_AF) %>% mutate(source = "SSNS vs Control (euro)")

new_df_bristol_ssns_ctr_sas = merge(excel_bristol$`SSNS vs Control (sas)`, 
                                     new_df_ssns_sas %>% select(Chromosome = CHR, Position = POS, A1, A2, 
                                                                original_N_case = N_case, original_AF_case = AF_case, 
                                                                original_N_control = N_control, original_AF_control = AF_control), 
                                     by = c("Chromosome", "Position", "A1", "A2")) %>% 
  relocate(c("original_N_case", "original_AF_case", "original_N_control", "original_AF_control"), .after = rsid) %>% 
  arrange(Bristol_GMMAT_PVAL) %>% select(-original_AF) %>% mutate(source = "SSNS vs Control (sas)")

new_df_bristol_ns_ctr = merge(excel_bristol$`NS vs Control`, 
                                    new_df_ns %>% select(Chromosome = CHR, Position = POS, A1, A2, 
                                                               original_N_case = N_case, original_AF_case = AF_case, 
                                                               original_N_control = N_control, original_AF_control = AF_control), 
                                    by = c("Chromosome", "Position", "A1", "A2")) %>% 
  relocate(c("original_N_case", "original_AF_case", "original_N_control", "original_AF_control"), .after = rsid) %>% 
  arrange(Bristol_GMMAT_PVAL) %>% select(-original_AF) %>% mutate(source = "NS vs Control")

new_df_bristol_ns_ctr_euro = merge(excel_bristol$`NS vs Control (euro)`, 
                              new_df_ns_eur %>% select(Chromosome = CHR, Position = POS, A1, A2, 
                                                         original_N_case = N_case, original_AF_case = AF_case, 
                                                         original_N_control = N_control, original_AF_control = AF_control), 
                              by = c("Chromosome", "Position", "A1", "A2")) %>% 
  relocate(c("original_N_case", "original_AF_case", "original_N_control", "original_AF_control"), .after = rsid) %>% 
  arrange(Bristol_GMMAT_PVAL) %>% select(-original_AF) %>% mutate(source = "NS vs Control (euro)")

new_df_bristol_ns_ctr_sas = merge(excel_bristol$`NS vs Control (sas)`, 
                                   new_df_ns_sas %>% select(Chromosome = CHR, Position = POS, A1, A2, 
                                                            original_N_case = N_case, original_AF_case = AF_case, 
                                                            original_N_control = N_control, original_AF_control = AF_control), 
                                   by = c("Chromosome", "Position", "A1", "A2")) %>% 
  relocate(c("original_N_case", "original_AF_case", "original_N_control", "original_AF_control"), .after = rsid) %>% 
  arrange(Bristol_GMMAT_PVAL) %>% select(-original_AF) %>% mutate(source = "NS vs Control (sas)")

new_df_bristol_srns_ctr = merge(excel_bristol$`SRNS vs Control`, 
                                  new_df_srns_ctr %>% select(Chromosome = CHR, Position = POS, A1, A2, 
                                                           original_N_case = N_case, original_AF_case = AF_case, 
                                                           original_N_control = N_control, original_AF_control = AF_control), 
                                  by = c("Chromosome", "Position", "A1", "A2")) %>% 
  relocate(c("original_N_case", "original_AF_case", "original_N_control", "original_AF_control"), .after = rsid) %>% 
  arrange(Bristol_GMMAT_PVAL) %>% select(-original_AF) %>% mutate(source = "SRNS vs Control")

new_df_bristol_srns_ssns = merge(excel_bristol$`SRNS vs SSNS`, 
                                new_df_srns_ssns %>% select(Chromosome = CHR, Position = POS, A1, A2, 
                                                           original_N_case = N_case, original_AF_case = AF_case, 
                                                           original_N_control = N_control, original_AF_control = AF_control), 
                                by = c("Chromosome", "Position", "A1", "A2")) %>% 
  relocate(c("original_N_case", "original_AF_case", "original_N_control", "original_AF_control"), .after = rsid) %>% 
  arrange(Bristol_GMMAT_PVAL) %>% select(-original_AF) %>% mutate(source = "SRNS vs SSNS")

options(java.parameters = "-Xmx2048m")
library(xlsx)
sheet_list = list(new_df_bristol_ssns_ctr, new_df_bristol_ssns_ctr_euro, new_df_bristol_ssns_ctr_sas,
                  new_df_bristol_ns_ctr, new_df_bristol_ns_ctr_euro, new_df_bristol_ns_ctr_sas,
                  new_df_bristol_srns_ctr, new_df_bristol_srns_ssns)

write.xlsx(as.data.frame(sheet_list[[1]]),
           file="/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/excel/Bristol.xlsx",
           sheetName=paste(sheet_list[[1]]$source %>% head(1)), row.names=FALSE)
for (x in 2:8) {
  print(x)
  print(as.character(paste(sheet_list[[x]]$source %>% head(1))))
  write.xlsx(as.data.frame(sheet_list[[x]]),
             file="/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/AF_casecontrol/excel/Bristol.xlsx",
             sheetName=as.character(paste(sheet_list[[x]]$source %>% head(1))), append=TRUE, row.names=FALSE)
}