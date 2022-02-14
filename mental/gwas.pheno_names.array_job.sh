#!/bin/bash

#SBATCH --account=p33
#SBATCH --time=86:00:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --cpus-per-task=1

source /cluster/bin/jobsetup
module purge

source ${HOME}/.bashrc

module list
which plink2

BFILE_ARR=( "$@" )
BFILE_PHENO=${BFILE_ARR[${SLURM_ARRAY_TASK_ID}]}
BFILE=${BFILE_PHENO%:*}
PHENO_NAMES=${BFILE_PHENO#*:}

echo "${#BFILE_ARR[@]} bfiles in BFILE_ARR"
echo "GWAS of ${BFILE}"

MENTAL_DIR=/cluster/p/p33/cluster/users/alexeas/most_mental
PHENO=${MENTAL_DIR}/pheno/touchscreen_20_binary.translated.csv
COVAR=${MENTAL_DIR}/covar/mostest_mental_covar_200807.4plink
COVAR_TO_USE="21022_0_0 22001_0_0 22009_0_1 22009_0_2 22009_0_3 22009_0_4 22009_0_5 22009_0_6 22009_0_7 22009_0_8 22009_0_9 22009_0_10"
N_THREADS=1
MEMORY_MB=8192

OUT_NAME=${BFILE##*/}
OUT=${MENTAL_DIR}/new_start/gwas/${OUT_NAME}

echo Input bfile ${IN}
echo Constructing output ${OUT}

plink2 --bfile ${BFILE} --glm hide-covar firth-fallback --covar ${COVAR} --covar-name ${COVAR_TO_USE} --pheno ${PHENO} --threads ${N_THREADS} --memory ${MEMORY_MB} --1 --pheno-name ${PHENO_NAMES} --out ${OUT}

