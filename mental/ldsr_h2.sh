#!/bin/bash

source ${HOME}/ldsr/bin/activate

which python

LDSR_SRC=${HOME}/src/ldsc
LDSR_IN=/cluster/projects/p33/users/alexeas/most_mental/new_start/ldsr/munge/touchscreen
LDSR_OUT=/cluster/projects/p33/users/alexeas/most_mental/new_start/ldsr/h2/touchscreen

for f in ${LDSR_IN}/20127*.sumstats.gz; do
    PHENO=$(cut -f1 -d'.' <<< $(basename $f))
    python ${LDSR_SRC}/ldsc.py --h2 $f --ref-ld-chr ${LDSR_SRC}/eur_w_ld_chr/ --w-ld-chr ${LDSR_SRC}/eur_w_ld_chr/ --out ${LDSR_OUT}/${PHENO}.h2
done
