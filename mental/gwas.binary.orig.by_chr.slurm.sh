#!/bin/bash

#SBATCH --account=p33_norment
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=8GB
#SBATCH --cpus-per-task=1
#SBATCH --array=1-22

source /cluster/bin/jobsetup
module purge

source ${HOME}/.bashrc

module list

PHENO_NAME=$1
ORIG_OR_PERM=$2

echo "GWAS of ${PHENO_NAME} ${ORIG_OR_PERM}"

MENTAL_DIR=/cluster/p/p33/cluster/users/alexeas/most_mental
BFILE_DIR=/cluster/projects/p33/projects/mental/geno_hapmap/by_chr
PHENO=${MENTAL_DIR}/pheno/touchscreen_20_binary.translated.csv
COVAR=${MENTAL_DIR}/covar/mostest_mental_covar_200807.4plink
COVAR_TO_USE="21022_0_0 22001_0_0 22009_0_1 22009_0_2 22009_0_3 22009_0_4 22009_0_5 22009_0_6 22009_0_7 22009_0_8 22009_0_9 22009_0_10"
N_THREADS=1
MEMORY_MB=4000

if [ ${ORIG_OR_PERM} == "perm" ]; then
    IN="${BFILE_DIR}/ukb_imp_chr${SLURM_ARRAY_TASK_ID}_v3_qc.perm"
    OUT="${MENTAL_DIR}/gwas_results_by_chr/binary.${PHENO_NAME}.chr${SLURM_ARRAY_TASK_ID}.perm"
elif [ ${ORIG_OR_PERM} == "orig" ]; then
    IN="${BFILE_DIR}/ukb_imp_chr${SLURM_ARRAY_TASK_ID}_v3_qc"
    OUT="${MENTAL_DIR}/gwas_results_by_chr/binary.${PHENO_NAME}.chr${SLURM_ARRAY_TASK_ID}.orig"
else
    echo "The second argument must be either 'orig' or 'perm'."
    exit 1
fi

echo Input bfile ${IN}
echo Constructing output ${OUT}

plink --bfile ${IN} --logistic 'beta' 'hide-covar' --ci 0.95 --pheno ${PHENO} --pheno-name ${PHENO_NAME} --1 --covar ${COVAR} --covar-name ${COVAR_TO_USE} --out ${OUT} --threads ${N_THREADS} --memory ${MEMORY_MB}

