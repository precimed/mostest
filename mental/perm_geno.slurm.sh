#!/bin/bash

#SBATCH --account=p33_norment
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --cpus-per-task=1

source /cluster/bin/jobsetup
module purge

source ${HOME}/.bashrc
source ${HOME}/cluster/bin/activate

module list
which python

# $1 should be path to snp list /path/to/file/"chr<N>_*.txt"
SNP_LIST_PATH=$1
SNP_LIST_NAME=$(basename $1 | cut -f1 -d'.')
echo "Processing snp list $1"

IN_BFILE=/cluster/projects/p33/projects/mental/geno/generic_qc/snp_chunks/orig.${SNP_LIST_NAME}
OUT_BFILE=/cluster/projects/p33/projects/mental/geno/generic_qc/snp_chunks/perm.${SNP_LIST_NAME}
PY=/cluster/projects/p33/users/alexeas/most_mental/src/permute_bed.py

python ${PY} --bfile ${IN_BFILE} --out ${OUT_BFILE}

