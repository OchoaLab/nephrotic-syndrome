#!/bin/bash
#SBATCH -p ochoalab --account=ochoalab
##SBATCH --job-name=grm-ns-ctrl
##SBATCH --output=grm-ns-ctrl.out
##SBATCH --job-name=grm-ssns-ctrl
##SBATCH --output=grm-ssns-ctrl.out
##SBATCH --job-name=grm-srns-ctrl
##SBATCH --output=grm-srns-ctrl.out
##SBATCH --job-name=grm-ssns-srns
##SBATCH --output=grm-ssns-srns.out
##SBATCH --job-name=grm-ns-ctrl-afr
##SBATCH --output=grm-ns-ctrl-afr.out
##SBATCH --job-name=grm-ns-ctrl-eur
##SBATCH --output=grm-ns-ctrl-eur.out
##SBATCH --job-name=grm-ns-ctrl-sas
##SBATCH --output=grm-ns-ctrl-sas.out
#SBATCH --job-name=grm-ssns-ctrl-afr
#SBATCH --output=grm-ssns-ctrl-afr.out
##SBATCH --job-name=grm-ssns-ctrl-eur
##SBATCH --output=grm-ssns-ctrl-eur.out
##SBATCH --job-name=grm-ssns-ctrl-sas
##SBATCH --output=grm-ssns-ctrl-sas.out
#SBATCH --mem=16G
#SBATCH --ntasks-per-node=6
##SBATCH --mail-user=tiffany.tu@duke.edu
##SBATCH --mail-type=END,FAIL

# select subtype
# dir=. # == ns_ctrl
# dir=ssns_ctrl
# dir=srns_ctrl
# dir=ssns_srns
# dir=afr
# dir=eur
# dir=sas
dir=ssns_ctrl/afr
# dir=ssns_ctrl/eur
# dir=ssns_ctrl/sas

# run!
time gcta64 --bfile $dir/mac20 --make-grm --out $dir/mac20
time gcta64 --grm $dir/mac20 --pca 10 --out $dir/mac20
