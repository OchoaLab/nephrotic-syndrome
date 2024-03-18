# raw data was downloaded manually from collaborator to DCC using Globus
# placed here: /datacommons/ochoalab/curegn/raw/

# start interactive shell
# some of these steps need more than the default 1G memory!
srun -p ochoalab --account ochoalab --mem 16G --pty bash -i

# load plink on DCC
#module load Plink/2.00a3LM
# to use --pmerge-list, I needed a newer plink, copied that to my ~/bin/


################
### MANIFEST ###
################

# got CureGN_RNASeq_WGS_ID_Mapping_20230505.xlsx in an email, placed in same location as this script
# outputs indiv-keep.txt
time Rscript 00-remove-withdrawn.R
wc -l indiv-keep.txt # 2012
# copy that to DCC, same location as raw data
scp indiv-keep.txt $dcc:/datacommons/ochoalab/curegn/raw/


#################
### MAKE PGEN ###
#################

# these files are monstrous, so shrink ASAP
# - keep only desired individuals from get go (raw data includes withdrawn individuals)
# - remove non-PASS loci
# - remove fixed loci (since we removed many individuals, this is bound to result in fixed loci which had rare variants present in those individuals)
# - reencode to compact plink2 format (loses likelihoods and other VCF-specific info, keeps hard calls only)
# - ignores X and Y chrs, though we have them

for chr in {1..22}; do
    time plink2 --vcf curegn.wgs.freeze1.chr$chr.vcf.gz --keep indiv-keep.txt --var-filter --mac 1 --make-pgen vzs --out chr$chr --threads 1 --memory 16000
done
# all DCC runs
# files had 1869 individuals, 1855 are kept (14 removed) # OLD, didn't rerun here, but should have removed two more
# about half of SNPs were removed with the above filters!

# chr  real_time    input  output
#   1 20m58.477s 11297391 6472590
#   2 23m33.288s 12382583 7065342
#   3 19m19.350s 10255030 5848506
#   4 19m 3.720s 10020151 5770905
#   5 17m48.995s  9299147 5339134
#   6 16m36.284s  8803598 5095754
#   7 15m47.968s  8326001 4798652
#   8 15m 8.400s  8043369 4602016
#   9 13m 4.067s  6186553 3559603
#  10 13m18.539s  6953171 4026489
#  11 13m40.269s  7067057 4040803
#  12 13m49.597s  6842162 3941262
#  13  9m35.117s  5011805 2886393
#  14  9m32.745s  4624430 2669881
#  15 12m15.745s  4212025 2428264
#  16  9m10.485s  4746449 2716753
#  17  8m 1.261s  4169196 2403738
#  18  7m34.482s  3951228 2276777
#  19  6m38.625s  3386297 1977223
#  20  6m 7.143s  3209704 1844720
#  21  3m49.079s  1905232 1116162
#  22  3m39.276s  1951156 1134295

#################
### RM INDELS ###
#################

# would have done with earlier step, but don't want to repeat slow steps
# indels are a big problem, and we're really only looking at SNPs for discovery data, so let's filter the same here
for chr in {1..22}; do
    time plink2 --pfile chr$chr vzs --snps-only just-acgt --make-pgen vzs --out snps-chr$chr --threads 1 --memory 16000
    # 0m20.588s chr1 DCC
done

# done in part to solve this error in the --pmerge-list step:
# Error: Conflicting REF alleles for variant '.' at 1:30867.
# zstdgrep -P "^1\t30867"  chr1.pvar.zst 
# # 1	30867	.	C	T	616.03	PASS	AC=2;AF=0.0004045;AN=2106;DP=183284;ExcessHet=0;FS=1.612;InbreedingCoeff=0.1983;MQ=9.2;MQRankSum=-0.348;NEGATIVE_TRAIN_SITE;QD=0.29;ReadPosRankSum=0.69;SOR=0.544;VQSLOD=-1.594;culprit=MQRankSum
# # 1	30867	.	CCTCT	C	616.03	PASS	AC=2;AF=0.0002825;AN=2012;DP=183284;ExcessHet=0;FS=1.612;InbreedingCoeff=0.1983;MQ=9.2;MQRankSum=-0.348;NEGATIVE_TRAIN_SITE;QD=0.29;ReadPosRankSum=0.69;SOR=0.544;VQSLOD=-1.594;culprit=MQRankSum

##############
### PMERGE ###
##############

name=curegn-autosomes-snps

# creates list of files to merge
for chr in {1..22}; do
    echo snps-chr$chr >> files-merge-list.txt
done

# NOTES:
# - for `--pmerge-list` I needed a newer plink than the newest on DCC, copied it manually here
# - `--merge-max-allele-ct 2` is lazy hack to solve this error: 
# Error: The biallelic variants with ID '.' at position 1:28588 in chr1.pvar.zst
# appear to be the components of a 'split' multiallelic variant; if so, it must
# be 'joined' (with e.g. "bcftools norm -m") before a correct merge can occur. If
# you are SURE that your data does not contain any same-position same-ID variant
# groups that should be joined, you can suppress this error with
# --multiallelics-already-joined.

# run merge command
time plink2 --pmerge-list files-merge-list.txt pfile-vzs --merge-max-allele-ct 2 --pmerge-output-vzs --out $name
# 3m29.673s DCC
# dims according to report: 1855 x 70,491,341

du -hs $name.*
# 3.5G	curegn-autosomes-snps.pgen
# 96K	curegn-autosomes-snps.psam
# 3.1G	curegn-autosomes-snps.pvar.zst

# cleanup, can toss per-chr copies and other stuff
for chr in {1..22}; do
    rm {,snps-}chr$chr.{log,pgen,psam,pvar.zst}
done
rm files-merge-list.txt
rm $name.log

#########################################
### ASSIGN UNIQUE IDS TO LOCI W/O IDS ###
#########################################

# this data has all IDs missing; regardless plink2 operates poorly without unique IDs

# set missing IDs to unique values to avoid these being detected as repeated IDs
# `--memory 12000` leftover from TGP code it was based from
time plink2 --pfile $name vzs --set-missing-var-ids '@:#' --make-just-pvar zs --out $name-uniq --memory 12000
# 2m23.756s DCC
# replace data after inspection
mv $name-uniq.pvar.zst $name.pvar.zst
# trash
rm $name-uniq.log

################
### MAKE BED ###
################

# minimal filtering to convert to BED
time plink2 --pfile $name vzs --max-alleles 2 --make-bed --out $name
# 1m54.380s DCC

# data dimensions
zstdcat $name.pvar.zst |grep -c -v '^\#'
# 70,491,341 # so multiallelics were handled successfully in earlier step!  Nothing removed here
wc -l $name.{bim,fam}
# 70,491,341 curegn-autosomes-snps.bim
#      1,855 curegn-autosomes-snps.fam

# cleanup
rm $name.log

##########
### MV ###
##########

# separate working data from messier raw data
mv $name.* ..
# rest of the commands happen in that space
cd ..

##############
### COVARS ###
##############

# NOTE: entire original fam file is trivial except for id column

# extracts and cleans up standard covariates
# creates patient-data.txt.gz and indiv-rm-round2.txt
time Rscript 01-covars.R

#####################
### FINAL CLEANUP ###
#####################

# filter data again because 3 individuals don't have phenotypes or covariates
time plink2 --bfile $name --remove indiv-rm-round2.txt --mac 1 --make-bed --out ${name}2
# 2m35.296s DCC

# cleanup
# remove original pvar files, which are no longer filtered correctly (we didn't have immediate plans for them)
rm $name.p???*
# boring log file
rm ${name}2.log
# replace old with new files
rename snps2 snps ${name}2.*
# preserve second removal list under raw
mv indiv-rm-round2.txt  raw/

# add missingness filters.  No samples were removed, but 2.5M variants are removed!
time plink2 --bfile $name --geno 0.1 --mind 0.1 --make-bed --out ${name}2
# 4m37.029s DCC

# cleanup
# remove unfiltered copy
rm $name.{bed,bim,fam}
# boring log file
rm ${name}2.log
# replace old with new files
rename snps2 snps ${name}2.*

# remove two additional individuals, redo SNP filters: 47371 variats were removed due to --mac only
time plink2 --bfile $name --keep raw/indiv-keep.txt --mac 1 --geno 0.1 --make-bed --out ${name}2
# 2m48.763s DCC
# when done and checked, repeat cleanup just above!

# create mac20 filtered version for GWAS
time plink2 --bfile $name --mac 20 --make-bed --out $name-mac20
# 1m45.938s DCC
# cleanup
rm $name-mac20.log


wc -l $name.{bim,fam}
# 67,841,058 curegn-autosomes-snps.bim
#      1,850 curegn-autosomes-snps.fam
wc -l $name-mac20.{bim,fam}
# 12,978,623 curegn-autosomes-snps-mac20.bim
#      1,850 curegn-autosomes-snps-mac20.fam
