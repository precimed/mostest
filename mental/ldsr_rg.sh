#!/bin/bash

source ${HOME}/ldsr/bin/activate

which python

LDSR_SRC=${HOME}/src/ldsc
LDSR_IN=/cluster/projects/p33/users/alexeas/most_mental/new_start/ldsr/munge/cognition
LDSR_OUT=/cluster/projects/p33/users/alexeas/most_mental/new_start/ldsr/rg

SS_ARR=(/cluster/projects/p33/users/alexeas/most_mental/new_start/ldsr/munge/cognition/*gz /cluster/projects/p33/users/alexeas/most_mental/new_start/ldsr/munge/touchscreen/*gz)

N=${#SS_ARR[@]}
echo "$N sumstats files in total"
N_MINUS_1=$(( $N - 1 ))
N_MINUS_2=$(( $N - 2 ))

for i in $( seq 0 ${N_MINUS_2} ); do
    ii=$(( $i + 1 ))
    for j in $( seq ${ii} ${N_MINUS_1} ); do
        SSi=${SS_ARR[$i]}
        SSj=${SS_ARR[$j]}
        PHENOi=$(cut -f1 -d'.' <<< $(basename ${SSi}))
        PHENOj=$(cut -f1 -d'.' <<< $(basename ${SSj}))
        if [[ ${PHENOi} == "20127_0_0_quant" || ${PHENOj} == "20127_0_0_quant" ]]; then
            python ${LDSR_SRC}/ldsc.py --rg ${SSi},${SSj} --ref-ld-chr ${LDSR_SRC}/eur_w_ld_chr/ --w-ld-chr ${LDSR_SRC}/eur_w_ld_chr/ --out ${LDSR_OUT}/${PHENOi}_vs_${PHENOj}.rg
        fi
    done
done 

