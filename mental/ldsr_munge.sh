#!/bin/bash

source ${HOME}/ldsr/bin/activate

which python

GWAS_DIR=/cluster/projects/p33/users/alexeas/most_mental/new_start/gwas_merged/touchscreen
LDSR_SRC=${HOME}/src/ldsc
LDSR_OUT=/cluster/projects/p33/users/alexeas/most_mental/new_start/ldsr/munge/touchscreen

for f in $( ls -1 ${GWAS_DIR}/orig.20127*.csv | grep quant ); do
    PHENO=$(cut -f2 -d'.' <<< $f)
    STAT="T_STAT,0"
    #STAT="Z_STAT,0"
    python ${LDSR_SRC}/munge_sumstats.py --sumstats $f --snp ID --signed-sumstats ${STAT} --N-col OBS_CT --a1 A1 --a2 REF --ignore ALT,OR,BETA,ERRCODE --merge-alleles ${LDSR_SRC}/w_hm3.snplist  --out ${LDSR_OUT}/${PHENO}
done
