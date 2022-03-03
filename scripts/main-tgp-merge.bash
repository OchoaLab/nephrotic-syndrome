#################
### TGP MERGE ###
#################

# NOTE: both data are assumed to be same build (they are both hg38) and alleles should be aligned to same reference sequence (though I didn't check, Debo said they were and the numbers appear to agree enough)

# run locally (not DCC) from location of SSNS BIM data, place `bim_intersect.pl` there
time perl -w bim_intersect.pl ssns_gwas_hwe.bim ~/dbs/tgp-nygc/tgp-nygc-autosomes.bim
# 2m45.208s viiiaR5

# compare inputs and outputs
wc -l ssns_gwas_hwe.bim ~/dbs/tgp-nygc/tgp-nygc-autosomes.bim
#  1461991 ssns_gwas_hwe.bim
# 91784660 /home/viiia/dbs/tgp-nygc/tgp-nygc-autosomes.bim
wc -l ssns_gwas_hwe_out.bim ~/dbs/tgp-nygc/tgp-nygc-autosomes_out.bim
#  1034655 ssns_gwas_hwe_out.bim    # 70.8% of original data!
#  1034655 /home/viiia/dbs/tgp-nygc/tgp-nygc-autosomes_out.bim

# create TGP filtered to SNPs in common with array data (not strictly needed, but nice to have a small copy and validate/filter more carefully).
cd ~/dbs/tgp-nygc/
time plink2 --bfile tgp-nygc-autosomes --extract tgp-nygc-autosomes_out.bim --make-bed --out tgp-nygc-autosomes_ssns-gwas-intersect
# 5m49.228s viiiaR5

# check sizes, yes as desired!
wc -l tgp-nygc-autosomes_ssns-gwas-intersect.{bim,fam}
# 1034655 tgp-nygc-autosomes_ssns-gwas-intersect.bim
#    2504 tgp-nygc-autosomes_ssns-gwas-intersect.fam

# verified these two are identical!
diff -q tgp-nygc-autosomes_out.bim tgp-nygc-autosomes_ssns-gwas-intersect.bim
# cleanup
rm tgp-nygc-autosomes_ssns-gwas-intersect.log tgp-nygc-autosomes_out.bim

# transfered subset data (tgp-nygc-autosomes_ssns-gwas-intersect.{bed,bim,fam}) and ssns_gwas_hwe_out.bim to DCC

# back in DCC...
cd /datacommons/ochoalab/ssns_gwas

# remove loci also not in the intersection
module load Plink/2.00a3LM
time plink2 --bfile ssns_gwas_hwe --extract ssns_gwas_hwe_out.bim --make-bed --out ssns_gwas_hwe_tgp-intersect
# 0m7.852s DCC

# check that sizes are are expected
wc -l ssns_gwas_hwe{,_out,_tgp-intersect}.bim
# 1461991 ssns_gwas_hwe.bim
# 1034655 ssns_gwas_hwe_out.bim
# 1034656 ssns_gwas_hwe_tgp-intersect.bim # this is too large by one!

# this is the sole difference
diff ssns_gwas_hwe_{out,tgp-intersect}.bim
# 286086a286087
# > 4	401070	103.7638	99412689	C	G
# it's weird because "401070" as an ID doesn't repeat
# nor is there a repeat in that location, these are the lines around it:
# 4       kgp10752409     103.7637        99412455        C       T
# 4       401070  103.7638        99412689        C       G
# 4       kgp3087994      103.7638        99412820        A       T
# text search only finds the number again as a position in chr6:
# 6       rs78073965      0       401070  T       C

# so I have no explanation for this false positive, might be a weird fluke, but let's just remove that one separately
echo 401070 > rm.txt
time plink2 --bfile ssns_gwas_hwe --extract ssns_gwas_hwe_out.bim --exclude rm.txt --make-bed --out ssns_gwas_hwe_tgp-intersect2
# 0m1.965s
wc -l ssns_gwas_hwe_tgp-intersect2.{bim,fam}
 1034655 ssns_gwas_hwe_tgp-intersect2.bim
    1378 ssns_gwas_hwe_tgp-intersect2.fam
# these are indeed identical, as expected
diff -q ssns_gwas_hwe_{out,tgp-intersect2}.bim
# cleanup
rm rm.txt ssns_gwas_hwe_tgp-intersect.??? ssns_gwas_hwe_tgp-intersect2.log ssns_gwas_hwe_out.bim
rename 2 '' ssns_gwas_hwe_tgp-intersect2.???

# and these also agree
wc -l ssns_gwas_hwe_tgp-intersect.bim ../tgp/tgp-nygc-autosomes_ssns-gwas-intersect.bim
# 1034655 ssns_gwas_hwe_tgp-intersect.bim
# 1034655 ../tgp/tgp-nygc-autosomes_ssns-gwas-intersect.bim

# the only differences are the IDs (it's messy but since all else agrees, we'll just have to assume that it's all in agreement now)
# # for easy merge, replace one BIM so the IDs also match!
# mv ssns_gwas_hwe_tgp-intersect.bim{,~} # keep backup just in case
# cp ../tgp/tgp-nygc-autosomes_ssns-gwas-intersect.bim ssns_gwas_hwe_tgp-intersect.bim

# # merge!
# time ../../bin/plink1 --keep-allele-order --bfile ssns_gwas_hwe_tgp-intersect --bmerge ../tgp/tgp-nygc-autosomes_ssns-gwas-intersect --out data-tgp

# # check dimensions
# wc -l data-tgp.{bim,fam}
