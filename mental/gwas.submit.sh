#!/bin/bash

BFILE_DIR=/cluster/projects/p33/projects/mental/geno/generic_qc/snp_chunks
N_THROTTLE=600

BFILE_ARR=(${BFILE_DIR}/*.bed)
BFILE_ARR=(${BFILE_ARR[@]%.bed})
N_BFILES=${#BFILE_ARR[@]}
N_ARRAY_JOB_MAX=$(( ${N_BFILES} - 1 ))
echo "${N_BFILES} bfiles"
sbatch --array 0-${N_ARRAY_JOB_MAX}%${N_THROTTLE} gwas.array_job.sh ${BFILE_ARR[@]}

