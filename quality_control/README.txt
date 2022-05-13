README file for May 2022 analysis

* remove bristol samples, create new bim/bed/fam files
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

5. imputation step
  -> run 'ste4.q' in imputation folder