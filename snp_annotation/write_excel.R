library(tidyverse)
#### organize results
merge_rsid = read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/annotate/combined_list_pval_rsid.txt", sep = "\t", header = TRUE) %>% distinct(.keep_all = TRUE)
nexus_neargene_10000 <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/annotate/near_gens_10000.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)  %>% 
  mutate(Chromosome = str_replace(Chromosome, "chr", "")) %>% dplyr::rename(CHR = Chromosome) %>% 
  select(-Position) %>% separate(Variation.ID, c('chr', 'POS', "A1", "A2", "strand")) %>% 
  select(CHR, POS, everything(), -chr, -strand) %>% distinct(.keep_all = TRUE) %>% 
  mutate(SNP = paste0("chr", CHR, ":", POS, ":", A2, ":", A1)) %>% select(-A1, -A2)

nexus_neargene_1102 <- read.table("/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/annotate/near_gens_1102.txt", sep = "\t", stringsAsFactors=FALSE, quote = "", header = TRUE)  %>% 
  mutate(Chromosome = str_replace(Chromosome, "chr", "")) %>% dplyr::rename(CHR = Chromosome) %>% 
  select(-Position) %>% separate(Variation.ID, c('chr', 'POS', "A1", "A2", "strand")) %>% 
  select(CHR, POS, everything(), -chr, -strand) %>% distinct(.keep_all = TRUE) %>% 
  mutate(SNP = paste0("chr", CHR, ":", POS, ":", A2, ":", A1)) %>% select(-A1, -A2)

nexus_neargene_combine = rbind(nexus_neargene_10000, nexus_neargene_1102)
merge_results = merge(nexus_neargene_combine %>% select(-SNP), merge_rsid %>% select(-SNP, -cM), by = c("CHR", "POS"))
colnames(merge_results)

#### split df into separate df's by analysis and write into excel sheets
options(java.parameters = "-Xmx2048m")  
library(xlsx)
merge_results_subset = merge_results %>% arrange(PVAL) %>% group_split(source)
#paste(merge_results_subset[[1]]$source %>% head(1))

write.xlsx(as.data.frame(merge_results_subset[[3]]), file="/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/annotate/filename_new__.xlsx", sheetName=paste(merge_results_subset[[3]]$source %>% head(1)), row.names=FALSE)
#length(merge_results_subset)
for (x in 4:25) {
  #print(as.data.frame(merge_results_subset[[x]]))
  print(x)
  write.xlsx(as.data.frame(merge_results_subset[[x]]), file="/datacommons/ochoalab/ssns_gwas/GMMAT_0418/results_sugg_sig/annotate/filename_new__.xlsx", sheetName=paste(merge_results_subset[[x]]$source %>% head(1)), append=TRUE, row.names=FALSE)
  print(x)
}
