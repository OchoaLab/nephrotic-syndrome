#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=cc-af
#SBATCH --output=cc_af.out
#SBATCH --mem=64G
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=tiffany.tu@duke.edu
#SBATCH --mail-type=END,FAIL

# NOTE: 32G sufficed for HO and HGDP, but increased for TGP

module load Plink/1.90  

# # ssns
# time plink --bfile /datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/ssns_ctr_mac20 --pheno /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/ssns_ctr_pheno.txt --extract snp_id/ssns_ctr_snpid.txt --keep-allele-order --freq case-control --nonfounders --allow-no-sex --out ssns_ctr_casecontrol
# time plink --bfile /datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/ssns_ctr_mac20 --pheno /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/ssns_ctr_pheno.txt --extract snp_id/ssns_ctr_cond1_snpid.txt --keep-allele-order --freq case-control --nonfounders --allow-no-sex --out ssns_ctr_cond1_casecontrol
# time plink --bfile /datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/ssns_ctr_mac20 --pheno /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/ssns_ctr_pheno.txt --extract snp_id/ssns_ctr_cond2_snpid.txt --keep-allele-order --freq case-control --nonfounders --allow-no-sex --out ssns_ctr_cond2_casecontrol
# time plink --bfile /datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/ssns_ctr_mac20 --pheno /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/ssns_ctr_pheno.txt --extract snp_id/ssns_ctr_cond3_snpid.txt --keep-allele-order --freq case-control --nonfounders --allow-no-sex --out ssns_ctr_cond3_casecontrol
# time plink --bfile /datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/ssns_ctr_mac20 --pheno /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/ssns_ctr_pheno.txt --extract snp_id/ssns_ctr_cond4_snpid.txt --keep-allele-order --freq case-control --nonfounders --allow-no-sex --out ssns_ctr_cond4_casecontrol
# # ssns african
# time plink --bfile /datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/afr/ssns_ctr_afr_mac20 --pheno /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/split_filter/ssns_ctr_ancestry/ssns_ctr_pheno_b.txt --keep-allele-order --extract snp_id/ssns_ctr_afr_snpid.txt --freq case-control --nonfounders --allow-no-sex --out ssns_ctr_afr_casecontrol
# time plink --bfile /datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/afr/ssns_ctr_afr_mac20 --pheno /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/split_filter/ssns_ctr_ancestry/ssns_ctr_pheno_b.txt --keep-allele-order --extract snp_id/ssns_ctr_afr_cond1_snpid.txt --freq case-control --nonfounders --allow-no-sex --out ssns_ctr_afr_cond1_casecontrol
# time plink --bfile /datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/afr/ssns_ctr_afr_mac20 --pheno /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/split_filter/ssns_ctr_ancestry/ssns_ctr_pheno_b.txt --keep-allele-order --extract snp_id/ssns_ctr_afr_cond2_snpid.txt --freq case-control --nonfounders --allow-no-sex --out ssns_ctr_afr_cond2_casecontrol
# 
# # ssns european
# time plink --bfile /datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/euro/ssns_ctr_euro_mac20 --pheno /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/split_filter/ssns_ctr_ancestry/ssns_ctr_pheno_w.txt --keep-allele-order --extract snp_id/ssns_ctr_eur_snpid.txt --freq case-control --nonfounders --allow-no-sex --out ssns_ctr_eur_casecontrol
# time plink --bfile /datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/euro/ssns_ctr_euro_mac20 --pheno /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/split_filter/ssns_ctr_ancestry/ssns_ctr_pheno_w.txt --keep-allele-order --extract snp_id/ssns_ctr_eur_cond1_snpid.txt --freq case-control --nonfounders --allow-no-sex --out ssns_ctr_eur_cond1_casecontrol
# time plink --bfile /datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/euro/ssns_ctr_euro_mac20 --pheno /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/split_filter/ssns_ctr_ancestry/ssns_ctr_pheno_w.txt --keep-allele-order --extract snp_id/ssns_ctr_eur_cond2_snpid.txt --freq case-control --nonfounders --allow-no-sex --out ssns_ctr_eur_cond2_casecontrol
# 
# # ssns south asian
# time plink --bfile /datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/sas/ssns_ctr_sas_mac20 --pheno /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/split_filter/ssns_ctr_ancestry/ssns_ctr_pheno_sa.txt --keep-allele-order --extract snp_id/ssns_ctr_sas_snpid.txt --freq case-control --nonfounders --allow-no-sex --out ssns_ctr_sas_casecontrol
# time plink --bfile /datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/sas/ssns_ctr_sas_mac20 --pheno /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/split_filter/ssns_ctr_ancestry/ssns_ctr_pheno_sa.txt --keep-allele-order --extract snp_id/ssns_ctr_sas_cond1_snpid.txt --freq case-control --nonfounders --allow-no-sex --out ssns_ctr_sas_cond1_casecontrol
# time plink --bfile /datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/sas/ssns_ctr_sas_mac20 --pheno /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/split_filter/ssns_ctr_ancestry/ssns_ctr_pheno_sa.txt --keep-allele-order --extract snp_id/ssns_ctr_sas_cond2_snpid.txt --freq case-control --nonfounders --allow-no-sex --out ssns_ctr_sas_cond2_casecontrol
# 
# # ns
# time plink --bfile /datacommons/ochoalab/ssns_gwas/GMMAT_0418/ns_ctr/ssns_tgp_mac_20 --pheno /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/gcta/ssns_tgp_pheno_imp.txt --keep-allele-order --extract snp_id/ns_ctr_snpid.txt --freq case-control --nonfounders --allow-no-sex --out ns_ctr_casecontrol
# time plink --bfile /datacommons/ochoalab/ssns_gwas/GMMAT_0418/ns_ctr/ssns_tgp_mac_20 --pheno /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/gcta/ssns_tgp_pheno_imp.txt --keep-allele-order --extract snp_id/ns_ctr_cond1_snpid.txt --freq case-control --nonfounders --allow-no-sex --out ns_ctr_cond1_casecontrol
# time plink --bfile /datacommons/ochoalab/ssns_gwas/GMMAT_0418/ns_ctr/ssns_tgp_mac_20 --pheno /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/gcta/ssns_tgp_pheno_imp.txt --keep-allele-order --extract snp_id/ns_ctr_cond2_snpid.txt --freq case-control --nonfounders --allow-no-sex --out ns_ctr_cond2_casecontrol
# time plink --bfile /datacommons/ochoalab/ssns_gwas/GMMAT_0418/ns_ctr/ssns_tgp_mac_20 --pheno /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/gcta/ssns_tgp_pheno_imp.txt --keep-allele-order --extract snp_id/ns_ctr_cond3_snpid.txt --freq case-control --nonfounders --allow-no-sex --out ns_ctr_cond3_casecontrol
# 
# 
# time plink --bfile /datacommons/ochoalab/ssns_gwas/GMMAT_0418/ns_ctr/afr/ns_ctr_black_mac20 --pheno /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/ancestry/b_pheno.txt --keep-allele-order --extract snp_id/ns_ctr_afr_snpid.txt --freq case-control --nonfounders --allow-no-sex --out ns_ctr_afr_casecontrol
# time plink --bfile /datacommons/ochoalab/ssns_gwas/GMMAT_0418/ns_ctr/afr/ns_ctr_black_mac20 --pheno /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/ancestry/b_pheno.txt --keep-allele-order --extract snp_id/ns_ctr_afr_cond1_snpid.txt --freq case-control --nonfounders --allow-no-sex --out ns_ctr_afr_cond1_casecontrol
# time plink --bfile /datacommons/ochoalab/ssns_gwas/GMMAT_0418/ns_ctr/euro/ns_ctr_white_mac20 --pheno /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/ancestry/w_pheno.txt --keep-allele-order --extract snp_id/ns_ctr_eur_snpid.txt --freq case-control --nonfounders --allow-no-sex --out ns_ctr_euro_casecontrol
# time plink --bfile /datacommons/ochoalab/ssns_gwas/GMMAT_0418/ns_ctr/euro/ns_ctr_white_mac20 --pheno /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/ancestry/w_pheno.txt --keep-allele-order --extract snp_id/ns_ctr_eur_cond1_snpid.txt --freq case-control --nonfounders --allow-no-sex --out ns_ctr_euro_cond1_casecontrol
# time plink --bfile /datacommons/ochoalab/ssns_gwas/GMMAT_0418/ns_ctr/sas/ns_ctr_sa_mac20 --pheno /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/ancestry/sa_pheno.txt --keep-allele-order --extract snp_id/ns_ctr_sas_snpid.txt --freq case-control --nonfounders --allow-no-sex --out ns_ctr_sas_casecontrol
# time plink --bfile /datacommons/ochoalab/ssns_gwas/GMMAT_0418/ns_ctr/sas/ns_ctr_sa_mac20 --pheno /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/ancestry/sa_pheno.txt --keep-allele-order --extract snp_id/ns_ctr_sas_cond1_snpid.txt --freq case-control --nonfounders --allow-no-sex --out ns_ctr_sas_cond1_casecontrol
# 
# # srns
# time plink --bfile /datacommons/ochoalab/ssns_gwas/GMMAT_0418/srns_ctr/srns_ctr_match_mac20 --pheno /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/split_filter/srns_ctr_ancestry/srns_ctr_match_pheno.txt --keep-allele-order --extract snp_id/srns_ctr_snpid.txt --freq case-control --nonfounders --allow-no-sex --out srns_ctr_casecontrol
# time plink --bfile /datacommons/ochoalab/ssns_gwas/GMMAT_0418/srns_ctr/srns_ssns/srns_ssns_mac20 --pheno /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/ssns_srns_pheno.txt --keep-allele-order --extract snp_id/srns_ssns_snpid.txt --freq case-control --nonfounders --allow-no-sex --out srns_ssns_casecontrol

#meta
time plink --bfile /datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/ssns_ctr_mac20 --pheno /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/case_control/ssns_ctr_pheno.txt --extract snp_id/meta_ssns_ctr_snpid.txt --freq case-control --nonfounders --allow-no-sex --out meta_ssns_ctr_casecontrol
time plink --bfile /datacommons/ochoalab/ssns_gwas/GMMAT_0418/ns_ctr/ssns_tgp_mac_20 --pheno /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205/allele_freq/newfiles_alpha/imputation_10_new/results/gcta/ssns_tgp_pheno_imp.txt --extract snp_id/meta_ns_ctr_snpid.txt --freq case-control --nonfounders --allow-no-sex --out meta_ns_ctr_casecontrol

module unload Plink/1.90  



