# install validator
sudo pip3 install gwas-sumstats-tools

# validate files that Tiffany generated
for file in *.tsv.gz; do gwas-ssf validate $file; done

# checksums required with upload
md5sum *.tsv.gz
# 05de3959a8b49e60be62e341295214a5  ns_ctrl.tsv.gz
# ed17848c0f22fe1c5912d12c954da917  ns_ctrl_afr.tsv.gz
# aa7e740b7ff9c712d99a61ef364b9f8b  ns_ctrl_eur.tsv.gz
# 1b84c282df3d4e21485a0d701a00048d  ns_ctrl_sas.tsv.gz

# 987dae300cc82dc3ea79d229b1759c62  ssns_ctrl.tsv.gz
# 32a5e8ad5431c4abc6ad805e29b54247  ssns_ctrl_afr.tsv.gz
# b648b627356cfb5bbb2da72dc106970d  ssns_ctrl_eur.tsv.gz
# e654593ccdca8e3870010eb5725915f0  ssns_ctrl_sas.tsv.gz

# 823cbcc07aeaf953738cd338db1a1106  srns_ctrl.tsv.gz
# 69513523cec1927d109a6bdcb9915d25  srns_ctrl_afr.tsv.gz
# 38aef6a623412a92d33d4340d6ce4942  srns_ctrl_sas.tsv.gz

# e8b9f89ca8aa0bac5ee361bca0a30e95  ssns_srns.tsv.gz

# varint counts
gzl *.tsv.gz
# 21171019	ns_ctrl.tsv.gz
# 16098234	ns_ctrl_afr.tsv.gz
# 8883727	ns_ctrl_eur.tsv.gz
# 9025297	ns_ctrl_sas.tsv.gz

# 20838870	ssns_ctrl.tsv.gz
# 15922514	ssns_ctrl_afr.tsv.gz
# 8746653	ssns_ctrl_eur.tsv.gz
# 9002150	ssns_ctrl_sas.tsv.gz

# 20109280	srns_ctrl.tsv.gz
# 15504264	srns_ctrl_afr.tsv.gz
# 8404919	srns_ctrl_sas.tsv.gz

# 12274558	ssns_srns.tsv.gz

# calculate MAF thresholds for each dataset (all had fixed MAC>=20, but because sample size varies so does MAF)
time Rscript maf-calc.R
tab maf-cuts.txt
# file                 maf_empir            maf_cut               n    n_case n_ctrl
# ns_ctrl.tsv.gz       0.00222965           0.002229654403567447  4485 932    3553
# ns_ctrl_afr.tsv.gz   0.00865052           0.00865051903114187   1156 194    962
# ns_ctrl_eur.tsv.gz   0.010110999999999981 0.010111223458038422  989  212    777
# ns_ctrl_sas.tsv.gz   0.009259000000000017 0.009259259259259259  1080 343    737

# ssns_ctrl.tsv.gz     0.00233754           0.0023375409069658717 4278 725    3553
# ssns_ctrl_afr.tsv.gz 0.009009000000000045 0.009009009009009009  1110 148    962
# ssns_ctrl_eur.tsv.gz 0.0109529            0.01095290251916758   913  136    777
# ssns_ctrl_sas.tsv.gz 0.00938086           0.009380863039399626  1066 329    737

# srns_ctrl.tsv.gz     0.00266951           0.0026695141484249867 3746 193    3553
# srns_ctrl_afr.tsv.gz 0.00997009           0.009970089730807577  1003 41     962
# srns_ctrl_sas.tsv.gz 0.0133156            0.013315579227696404  751  14     737

# ssns_srns.tsv.gz     0.010893000000000042 0.010893246187363835  918  193    725

# get the precise ancestry breakdown demanded by the gwas catalog
# run on DCC
time Rscript ancestry.R
