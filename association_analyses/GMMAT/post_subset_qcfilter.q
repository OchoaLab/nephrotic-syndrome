#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
#SBATCH --job-name=mac
#SBATCH --output=step3_mac.out
#SBATCH --mem=164G
#SBATCH --ntasks-per-node=1
#SBATCH --mail-user=tiffany.tu@duke.edu
#SBATCH --mail-type=END,FAIL

# NOTE: 32G sufficed for HO and HGDP, but increased for TGP

module load Plink/2.00a2LM

#time plink2 --bfile ssns_ctr --mind 0.1 --make-bed --out ssns_tgp_imp_mind
#time plink2 --bfile ssns_tgp_imp_mind --geno 0.1 --make-bed --out ssns_tgp_imp_geno
#time plink2 --bfile ssns_ctr --hwe 0.0001 --make-bed --out ssns_ctr_hwe
time plink2 --bfile ssns_ctr --mac 20 --make-bed --out /datacommons/ochoalab/ssns_gwas/GMMAT_0418/ssns_ctr/ssns_ctr_mac20

#time plink2 --bfile srns_ctr --hwe 0.0001 --make-bed --out srns_ctr_hwe
time plink2 --bfile srns_ctr --mac 20 --make-bed --out /datacommons/ochoalab/ssns_gwas/GMMAT_0418/srns_ctr/srns_ctr_mac20

#time plink2 --bfile ssns_srns --hwe 0.0001 --make-bed --out ssns_srns_hwe
time plink2 --bfile ssns_srns --mac 20 --make-bed --out /datacommons/ochoalab/ssns_gwas/GMMAT_0418/srns_ctr/srns_ssns_mac20

module unload Plink/2.00a2LM



