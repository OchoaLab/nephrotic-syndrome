# Alex ran this mostly on interactive shell rather than as submitted jobs, which explains lack of .q files
# start interactive job this way:
srun -p ochoalab --account ochoalab --pty bash -i

# most of the work requires plink2, but we'll load/unload plink1 here too as needed
module load Plink/2.00a3LM
# I used this R version too
module load R/4.0.0

# Run from this location
cd /datacommons/ochoalab/ssns_gwas/array


# reads 2022-02-25_2PhenotypeDatafirstDataSets.xlsx, PHENOTYPEDATASECONDSETOFPLATES2021_2022.xlsx, and original fam table
# creates patient-data.txt.gz (formerly pheno_excelmerge.csv), includes Bristol!
# creates ids-bristol.txt too
phenotype_excel_merge.Rmd


### REMOVE BRISTOL ###

# size of rawest data, including Bristol
wc -l ../raw/*.{fam,bim}
#    2656 nephrotic.syndrome.gwas.all.b38.n2656.R2.fam
# 1748250 nephrotic.syndrome.gwas.all.b38.n2656.R2.bim

# and list of individuals to remove
wc -l ids-bristol.txt 
# 590 ids-bristol.txt

# remove bristol samples, create new bim/bed/fam files
plink2 --bfile ../raw/nephrotic.syndrome.gwas.all.b38.n2656.R2 --remove ids-bristol.txt --make-bed --out ssns_remove_B
# creates ssns_remove_B.{bed,bim,fam}

# so output should have had the same number of SNPs, and 2066 individuals!


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

wc -l ssns_tgp_merge.{bim,fam}
# 762629 ssns_tgp_merge.bim # one more than expected!
#   4485 ssns_tgp_merge.fam


### ALLELE FREQUENCY TEST ###

# AF tests
# identify flip/remove SNPs for SSNS
# use --mem 16G for this job on slurm!
AF_preprocessing.Rmd
# creates allele_freq/{remove,flip}.txt

# apply changes to data!
module purge
module load Plink/1.90
time plink \
  --keep-allele-order \
  --bfile ssns_tgp_merge \
  --exclude allele_freq/remove.txt \
  --flip allele_freq/flip.txt \
  --flip-subset array-clean.fam \
  --make-bed --out ssns_tgp_merge_clean
# cleanup
rm ssns_tgp_merge_clean.nosex

wc -l ssns_tgp_merge_clean.{bim,fam}
# 761366 ssns_tgp_merge_clean.bim
#   4485 ssns_tgp_merge_clean.fam

# redo allele frequency calculations to confirm alignment succeeded
mkdir allele_freq2
cd allele_freq2
module purge
module load Plink/2.00a2LM
plink2 --bfile ../ssns_tgp_merge_clean --keep ../admixture/ids_controls_afr.txt --freq --out array_controls_afr
plink2 --bfile ../ssns_tgp_merge_clean --keep ../admixture/ids_controls_eur.txt --freq --out array_controls_eur
plink2 --bfile ../ssns_tgp_merge_clean --keep ../admixture/ids_controls_sas.txt --freq --out array_controls_sas
plink2 --bfile ../ssns_tgp_merge_clean --keep ../allele_freq/tgp_controls_afr.fam --freq --out tgp_controls_afr
plink2 --bfile ../ssns_tgp_merge_clean --keep ../allele_freq/tgp_controls_eur.fam --freq --out tgp_controls_eur
plink2 --bfile ../ssns_tgp_merge_clean --keep ../allele_freq/tgp_controls_sas.fam --freq --out tgp_controls_sas

# make plot that confirms data is now aligned (i.e. SNPs were correctly removed or flipped)
# makes allele_freq2/af-test.pdf
AF_postprocessing.Rmd


### CHECK REF ###

# verify that ref alleles are aligned with reference sequence!
# download hg38 reference sequence suggested in plink2 page
cd /datacommons/ochoalab/
wget https://www.dropbox.com/s/xyggouv3tnamh0j/GRCh38_full_analysis_set_plus_decoy_hla.fa.zst?dl=1
# fix stupid extension issue
mv GRCh38_full_analysis_set_plus_decoy_hla.fa.zst?dl=1 GRCh38_full_analysis_set_plus_decoy_hla.fa.zst
# go back to where the data is
cd /datacommons/ochoalab/ssns_gwas/array/
# perform test!
# I don't really want to make anything, but plink won't let me, specified --make-just-bim to keep it minimal
time plink2 \
     --bfile ssns_tgp_merge_clean \
     --fa /datacommons/ochoalab/GRCh38_full_analysis_set_plus_decoy_hla.fa.zst \
     --ref-from-fa \
     --make-just-bim \
     --out test
# --ref-from-fa: 0 variants changed, 761366 validated.
# 0m8.30s7 DCC
# further confirmation, these are identical!
diff -q test.bim ssns_tgp_merge_clean.bim
# cleanup
rm test.{bim,log}


### IMPUTATION ###

# put imputation server inputs in this directory
mkdir imputation-input/
cd imputation-input/

# create vcfs, split by chromosome
module load bcftools/1.4
for i in {1..22}; do
    plink2 --bfile ../ssns_tgp_merge_clean --chr $i --output-chr chrM --export vcf bgz id-paste=iid --out chr_$i
    bcftools index chr_$i.vcf.gz
done
cd ..

# Upload those to TopMed server!
# - Rsq filter = 0.3

# place raw outputs here: /datacommons/ochoalab/ssns_gwas/imputed/raw/
mkdir ../imputed/
cd ../imputed/
mkdir raw/

# will merge imputation data following steps I previously used for TGP per-chr VCFs
# https://github.com/OchoaLab/data/blob/main/tgp-nygc.bash

# pwd: /datacommons/ochoalab/ssns_gwas/imputed/

# convert each VCF to pgen, while removing SNPs that don't "PASS"
sbatch -a 1-22 vcf-chr-to-pgen.q
# max 18m3.880s DCC

# check logs, which confirm that nothing was actually removed (all "PASS"ed)
grep -P 'variants (loaded|remaining)' chr*.log

# use those new files to merge into combined pgen file
# run interactively (default --mem 1G suffices)
# creates list of files to merge
for chr in {1..22}; do
    echo chr$chr >> files-merge-list.txt
done
# run merge command
# NOTE: for this case only I had to upload a newer plink2 (v2.00a5LM), none of the DCC versions support `--pmerge-list`
time ./plink2 --pmerge-list files-merge-list.txt pfile-vzs --pmerge-output-vzs --out all
# 5m47.951s DCC

# cleanup, leaves raw data but removes new intermediates
for chr in {1..22}; do
    rm chr$chr.{log,pgen,psam,pvar.zst}
done
rm files-merge-list.txt

# filter more and convert to BED (works with `--mem 16G`, requires more than 1G)
time plink2 --pfile all vzs --max-alleles 2 --mac 1 --make-bed --out all
# 3,873,296 variants removed due to allele frequency threshold(s) (--mac 1)
# 3m29.408s

# NOTE: instead of "all", old name was "ssns_tgp_impute"

# data dimensions
zstdcat all.pvar.zst |wc -l
# 83,378,786 # includes header lines
# 83,378,750 # SNPs according to plink's log
wc -l all.{bim,fam}
# 79,505,454 all.bim
#      4,485 all.fam

# cleanup, all pgen intermediates again are unwanted! (gmmat doesn't support them)
rm all.{log,pgen,psam,pvar.zst}
# remove imputation input too, no longer needed (and easy to re-create if needed from rawer data)
rm -r ../array/imputation-input/


### GWAS PREP ###

# merges patient data (subset to remaining genotyped data) with TGP, including full race and sex covariates, binarized traits
# also creates ancestry column, which is like race except it distinguishes SAS from EAS and for array samples uses admixture results to exclude most admixed individuals from 3 largest race groups (and all individuals of minor race groups)
# creates imputed/patient-data.txt.gz, to be used with GMMAT and other analyses
Rscript merge-patient-data-tgp.R

# run simple trait ~ sex + race model (confirms need for these covariates)
# no files are created
Rscript covariate_analysis.R

# make temporary filter files for subanalyses (diagnosis and ancestry)
# creates {ssns_ctrl,srns_ctrl,ssns_srns,afr,eur,sas}/ids.txt
Rscript filter-subanalyses.R

# based on original scripts/data from here
#cd /datacommons/ochoalab/ssns_gwas/GMMAT_0418/ns_ctr

# filter by minor allele count (MAC)
# No other filters needed (there is no missingness, and HWE is not appropriate post-imputation)
# NOTE: requires `--mem 16G` on DCC
time plink2 --bfile all --mac 20 --make-bed --out mac20
# 5m37.128s

# filter subanalyses, always reapply `--mac 20` filter (removes more SNPs as sample sizes decrease)
time plink2 --bfile mac20 --keep ssns_ctrl/ids.txt --mac 20 --make-bed --out ssns_ctrl/mac20
# 2m6.163s DCC
time plink2 --bfile mac20 --keep srns_ctrl/ids.txt --mac 20 --make-bed --out srns_ctrl/mac20
# 1m52.397s DCC
time plink2 --bfile mac20 --keep ssns_srns/ids.txt --mac 20 --make-bed --out ssns_srns/mac20
# 1m20.720s DCC
time plink2 --bfile mac20 --keep afr/ids.txt --mac 20 --make-bed --out afr/mac20
# 1m54.713s DCC
time plink2 --bfile mac20 --keep eur/ids.txt --mac 20 --make-bed --out eur/mac20
# 1m26.657s DCC
time plink2 --bfile mac20 --keep sas/ids.txt --mac 20 --make-bed --out sas/mac20
# 1m27.884s DCC

# also create combined ancestry-subtype subsets!
# ns_ctrl is done above (it is bare {afr,eur,sas}/)
# only other case is ssns_ctrl (srns_ctrl and ssns_srns do not have ancestry subanalyses due to small sample sizes)
# create subdirectories
mkdir ssns_ctrl/afr
mkdir ssns_ctrl/eur
mkdir ssns_ctrl/sas
# create subsets with double filter
time plink2 --bfile ssns_ctrl/mac20 --keep afr/ids.txt --mac 20 --make-bed --out ssns_ctrl/afr/mac20
# 1m46.996s DCC
time plink2 --bfile ssns_ctrl/mac20 --keep eur/ids.txt --mac 20 --make-bed --out ssns_ctrl/eur/mac20
# 1m16.406s DCC
time plink2 --bfile ssns_ctrl/mac20 --keep sas/ids.txt --mac 20 --make-bed --out ssns_ctrl/sas/mac20
# 1m21.433s DCC

wc -l {.,ssns_ctrl,srns_ctrl,ssns_srns,afr,eur,sas,ssns_ctrl/afr,ssns_ctrl/eur,ssns_ctrl/sas}/mac20.{bim,fam}
# 21171018 ./mac20.bim
#     4485 ./mac20.fam
# 20838869 ssns_ctrl/mac20.bim
#     4278 ssns_ctrl/mac20.fam
# 20109279 srns_ctrl/mac20.bim
#     3746 srns_ctrl/mac20.fam
# 12274557 ssns_srns/mac20.bim
#      918 ssns_srns/mac20.fam
# 16098233 afr/mac20.bim
#     1156 afr/mac20.fam
#  8883726 eur/mac20.bim
#      989 eur/mac20.fam
#  9025296 sas/mac20.bim
#     1080 sas/mac20.fam
# 15922513 ssns_ctrl/afr/mac20.bim
#     1110 ssns_ctrl/afr/mac20.fam
#  8746652 ssns_ctrl/eur/mac20.bim
#      913 ssns_ctrl/eur/mac20.fam
#  9002149 ssns_ctrl/sas/mac20.bim
#     1066 ssns_ctrl/sas/mac20.fam

# use gcta (was already present in my path) to calculate GRM and PCs
# used `--mem 16G` on DCC
sbatch grm.q # edit to run each subtype
# time gcta64 --bfile mac20 --make-grm --out mac20
# # 113m14.189s ns_ctrl DCC
# # 105m51.466s ssns_ctrl DCC
# # 79m21.078s srns_ctrl DCC
# # 3m8.423s ssns_srns DCC
# # (other runs time not saved)
# # PCA part runs with default `--mem 1G`
# time gcta64 --grm mac20 --pca 10 --out mac20
# # 0m27.364s ns_ctrl DCC
# # 0m25.056s ssns_ctrl DCC
# # 0m16.534s srns_ctrl DCC
# # 0m0.328s ssns_srns DCC
# remove .out files when done

# this confirms all grm and eigenvec files were created
ls {.,ssns_ctrl,srns_ctrl,ssns_srns,afr,eur,sas,ssns_ctrl/afr,ssns_ctrl/eur,ssns_ctrl/sas}/mac20.{grm.bin,eigenvec}

# TODO: use multiple threads properly!  Manual says
# - only applies to glmm.score and SMMAT, argument ncores (DONE)
# - requires GDS format :( (next thing to do!)
sbatch gmmat.q
# runs: time Rscript gmmat.R 

# (to quickly get back where we left off)
#cd /datacommons/ochoalab/ssns_gwas/imputed
