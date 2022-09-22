May 2022 (QC, filtering)

* remove bristol samples, create new bim/bed/fam files (remove_B.q)
-> plink2 --bfile ../nephrotic.syndrome.gwas.all.b38.n2656/nephrotic.syndrome.gwas.all.b38.n2656.R2 --remove bristol_remove.txt --out ssns_remove_B --make-bed
-> creates ssns_remove_B.bed/bim/fam 

1. create 'duplicate_related.king.cutoff.out.id (step1.q)
  -> plink2 --bfile ../nephrotic.syndrome.gwas.all.b38.n2656/nephrotic.syndrome.gwas.all.b38.n2656.R2 --king-table-filter 0.354 --make-king-table --out duplicated_related
  -> creates duplicated_related.kin0
  
2. calculate missingness (step2.q)
  -> plink2 --bfile ../nephrotic.syndrome.gwas.all.b38.n2656/nephrotic.syndrome.gwas.all.b38.n2656.R2 --missing --out duplicate_sample_all_missingness
  -> creates duplicate_sample_all_missingness.smiss (columns: FID, ID, MISSING_CT, OBS_CT, FMISS) and duplicate_sample_all_missingness.vmiss (columns: CHROM, ID, MISSING_CT, OBS_CT, FMISS)
  
3. add missingness info to the .kin0 file created in step1
  -> run script /datacommons/ochoalab/tiffany_data/ns_qualitycontrol/missingness_kin.Rmd
  -> output: list of id's to remove in duplicate_remove.txt

4. remove individuals from duplicate_remove.txt and run various filtering steps
  -> run 'step3.q'
  -> output: bim/bed/bam files for 'ssns_gwas_maf'

June 2022 (merge SSNS + TGP) 

5. prepare TGP data
  -> tgp_flip.q
  -> tgp_flip_rm.q
  -> tgp_subset.q (output: 'TGP-subset')

6. prepare SSNS data
  -> recode ssns (recode_ssns.q) -> input: 'ssns_gwas_maf'; output: 'ssns_gwas_maf_recode'
  -> depup ssns (remove duplicate ID's from ssns) -> output: 'ssns_gwas_maf_dedup'

7. merge (merge_ssns_maf_tgp.q)
  -> input: 'ssns_gwas_maf_dedup' + 'TGP-subset', output: 'ssns_tgp_merge'

July/August 2022 (PCA/admixture analysis to identify outliers, allele frequency test)

* Identifying outliers: extract south asian individual id's/remove highly admixed individuals for AF test with TGP
  
* AF tests (likelihoodratiotest_function.R and merge_analysis_AF.Rmd)
  -> identify flip/remove SNPs for SSNS. Split with TGP, flip, merge with TGP, remove. (rm_flip_step.q) 

September 2022 (Imputation on TopMed server)
  -> split (remove TGP), create vcf
  -> split by chromosome
