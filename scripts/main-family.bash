# steps for processing family dataset, for quick ligera analysis only
# (not what Amika has used in the past)

cd ~/dbs/ssns_wgs_family

# minimal filtering to get plink1 data
time plink2 --vcf genotypes-recal.vcf.gz --var-filter --autosome --snps-only just-acgt --max-alleles 2 --mac 1 --make-bed --out data
# 4m53.799s

wc -l data.{bim,fam}
# 12757747 data.bim
#      102 data.fam

# cleanup
rm data.log

# need to map traits now, which is quite a mess

# start by mapping IDs using a header line from the VCF file, and the actual IDs (data.fam)
zcat genotypes-recal.vcf.gz|head -n 100|grep CombineGVCFs > tmp.txt
#bcftools query -l genotypes-recal.vcf.gz > ids.txt # completely redundant with data.fam
# clean up with perl hacks
# first remove everything up to the *first* listed "--variant" (not greedy)
perl -p -i -e 's/^.*?--variant/--variant/' tmp.txt
# now cut everything after "--reference" (greedy good here)
perl -p -i -e 's/ --reference.*$//' tmp.txt
# now remove all these flags
perl -p -i -e 's/--variant //g' tmp.txt
# remaining spaces are now better as newlines
perl -p -i -e 's/ /\n/g' tmp.txt
# take out shared paths
perl -p -i -e 's/\.\.\/processedData\/d06_gatk\/gvcf_files\///' tmp.txt
# and extensions/suffixes
perl -p -i -e 's/\.g\.vcf\.gz$//' tmp.txt
# there is one more variable prefix, fortunately regular, and not clearly useful for the rest of this work so I'll just toss
perl -p -i -e 's/^R\d+LD0\d-//' tmp.txt
# done!

# create nicer fam table (data2.fam), other demographics (info.txt)
cd ~/docs/ochoalab/nephrotic-syndrome/scripts/
Rscript fam00-demog-pheno.R

# go back to data and clean up
cd ~/dbs/ssns_wgs_family
mv data2.fam data.fam
rm tmp.txt 

