# Alex ran this mostly on interactive shell rather than as submitted jobs, which explains lack of .q files
# start interactive job this way:
srun -p ochoalab --account ochoalab --pty bash -i

# most of the work requires plink2, but we'll load/unload plink1 here too as needed
module load Plink/2.00a2LM
# I used this R version too
module load R/4.0.0

# Run from this location
cd /datacommons/ochoalab/ssns_gwas/nephrotic.syndrome.gwas_proprocessing_202205


# reads 2022-02-25_2PhenotypeDatafirstDataSets.xlsx, PHENOTYPEDATASECONDSETOFPLATES2021_2022.xlsx, and original fam table
# creates patient-data.txt.gz (formerly pheno_excelmerge.csv), includes Bristol!
# creates ids-bristol.txt too
phenotype_excel_merge.Rmd


### REMOVE BRISTOL ###

# remove bristol samples, create new bim/bed/fam files
plink2 --bfile ../nephrotic.syndrome.gwas.all.b38.n2656/nephrotic.syndrome.gwas.all.b38.n2656.R2 --remove ids-bristol.txt --make-bed --out ssns_remove_B
# creates ssns_remove_B.{bed,bim,fam}


### REMOVE DUPLICATE INDIVIDUALS, QC ###

# identify individuals that are likely duplicated (extreme kinship like that of twins)
plink2 --bfile ssns_remove_B --king-table-filter 0.354 --make-king-table --out duplicated_related
# creates duplicated_related.kin0
  
# calculate missingness
plink2 --bfile ssns_remove_B --missing sample-only --out duplicate_sample_all_missingness
# creates duplicate_sample_all_missingness.smiss
  
# combine dups with missingness, demographics, to decide who to remove
missingness_kin.Rmd
# output: duplicate_remove.txt

# remove individuals from duplicate_remove.txt and run various filtering steps
time plink2 \
     --bfile ssns_remove_B \
     --remove duplicate_remove.txt \
     --autosome \
     --snps-only just-acgt \
     --max-alleles 2 --min-alleles 2 \
     --mind 0.1 \
     --geno 0.1 \
     --hwe 0.0001 \
     --maf 0.01 \
     --make-bed --out ssns_gwas_maf

# cleanup: separate files/logs about individual dups
mkdir dedup-inds
mv duplicated_related.* duplicate_sample_all_missingness.* duplicate_remove.txt dedup-inds/

# recode array SNP IDs (needed to dedup positions, and for TGP merge)
# original had rs, JHU, and other IDs; new has chr:pos format that matches TGP
plink2 --bfile ssns_gwas_maf --set-all-var-ids '@:#' --make-just-bim --out ssns_gwas_maf_recode
# replace file while preserving original IDs
mv ssns_gwas_maf.bim ssns_gwas_maf_original.bim
mv ssns_gwas_maf_recode.bim ssns_gwas_maf.bim

# deduplicate array SNPs (based on chr:pos, though inspected and they are equivalent alt/ref alleles too, not multiallelic cases)
# NOTE: this requires `--mem 16G` on DCC!
time Rscript dedup-snps.R
# 1m40.272s DCC
# creates array-clean.{bed,bim,fam}
# 11 SNPs are duplicated, we kept one of each for 10 of them, removed both for one, so 12 SNPs removed total

# this is final clean version of array dataset!
wc -l array-clean.{bim,fam}
# 833047 array-clean.bim
#   1981 array-clean.fam

# remove large intermediate files we can easily reproduce from the raw originals
# left smaller logs and such behind for reference
rm {ssns_remove_B,ssns_gwas_maf}.{bed,bim,fam} ssns_gwas_maf_original.bim


### ADMIXTURE ###

# only run on array data, as needed because TGP merging is complex and many SNPs are misaligned or don't pass af-test, but af-test requires unadmixed subsets (created by this admixture analysis)!

# Run admixture!  this is slow so it's best as unsupervised run
# also compresses and moves outputs
sbatch admixture.q
# 223m13.160s/2141m9.388s DCC 10 threads

# create PCs (only used for one plot in admixture analysis, to validate)
plink2 --bfile array-clean --pca 10 --out array-clean

# make structure and PCA plots, identify highly admixed individuals
cd admixture/
admixture_results.Rmd
cd ..


### MERGE ARRAY AND TGP ###

# create TGP intersected with array SNPs by chr:pos (format of IDs in both cases)
time plink2 --bfile /datacommons/ochoalab/tgp/tgp-nygc-autosomes --extract array-clean.bim --make-bed --out TGP-subset
# 1m55.455s DCC

# reanalize array data for compatibility with tgp
Rscript premerge-tgp-arr.R
# creates premerge-array-exclude.txt (70418 SNPs), premerge-array-flip.txt (50404)

# apply exclude/flip changes to array data
# only plink1 has --flip
module purge # unloads plink2
module load Plink/1.90
plink --keep-allele-order --bfile array-clean --exclude premerge-array-exclude.txt --flip premerge-array-flip.txt --make-bed --out array-clean2

# change array ref to match TGP
module purge
module load Plink/2.00a2LM
plink2 --bfile array-clean2 --ref-allele TGP-subset.bim 6 2 --make-bed --out array-clean3

# TGP still has 81 extra SNPs that were removed in array because, though they match in chr:pos, they weren't a match to ref/alt, remove those now
plink2 --bfile TGP-subset --extract array-clean3.bim --make-bed --out TGP-subset2

# now merge!  only plink1 works for this
module purge
module load Plink/1.90
plink --keep-allele-order --bfile array-clean3 --bmerge TGP-subset2 --out ssns_tgp_merge

# cleanup
rm array-clean{2,3}.* TGP-subset{,2}.??? 

wc -l ssns_tgp_merge.{bim,fam} old/ssns_tpg_merge.{bim,fam}
# 762629 ssns_tgp_merge.bim # one more than expected!
#   4485 ssns_tgp_merge.fam
# 833047 old/ssns_tpg_merge.bim # 70418 more
#   4485 old/ssns_tpg_merge.fam


### ALLELE FREQUENCY TEST ###

# AF tests
# identify flip/remove SNPs for SSNS. Split with TGP, flip, merge with TGP, remove.
# use --mem 16G for this job on slurm!
AF_preprocessing.Rmd

### TODO ###

AF_postprocessing.Rmd

### IMPUTATION ###

# Imputation on TopMed server
# - create vcf
# - split by chromosome

### GWAS PREP ###

# reads phenotype_filtered.txt (its creation is not documented)
# makes merge_pheno.txt, merge_pheno_race_sex.txt
create_phenofile.Rmd

# reads ssns_tgp_pheno_imp.txt (its creation is not documented), merge_pheno_race_sex.txt, patient-data.txt.gz
# create covariate file for merged NS + TGP data `covar_ns_tgp.txt`
# run simple trait ~ sex + race model (confirms need for these covariates)
# obsolete: plot gcta results after incorporating covariates
# obsolete: compare different MAC threshold results
covariate_analysis.Rmd

# reads covar_ns_tgp_new.txt, patient-data.txt.gz
# creates binary case/control phenotype files for subanalyses:
# - ssns vs control: ssns_ctr.txt, ssns_ctr_pheno.txt, ssns_ctr_covar.txt
# - srns vs control: srns_ctr.txt, srns_ctr_pheno.txt, srns_ctr_covar.txt
# - ssns vs srns: ssns_srns.txt, ssns_srns_pheno.txt, ssns_srns_covar.txt
write_txt_files.Rmd


# probably obsolete
hwe_analysis.Rmd
split_allelefreq.Rmd


