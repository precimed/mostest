#!/bin/bash

#SBATCH --account=p33_norment
#SBATCH --time=10:00
#SBATCH --mem-per-cpu=4GB
#SBATCH --cpus-per-task=1

source /cluster/bin/jobsetup
module purge

source ${HOME}/.bashrc
source ${HOME}/cluster/bin/activate

module list
which plink

# $1 should be /path/to/file/"chr<N>_*"
SNP_LIST_PATH=$1
SNP_LIST_NAME=$(basename $1 | cut -f1 -d'.')
CHR=$(cut -f1 -d'_' <<< ${SNP_LIST_NAME})
echo "Processing snp list $1 for chromosome ${CHR}"

IN_BFILE=/cluster/projects/p33/projects/mental/geno/generic_qc/ukb_imp_${CHR}_v3_qc
OUT_BFILE=/cluster/projects/p33/projects/mental/geno/generic_qc/snp_chunks/orig.${SNP_LIST_NAME}

plink --bfile ${IN_BFILE} --extract ${SNP_LIST_PATH} --out ${OUT_BFILE} --make-bed --threads 1 --memory 4048

