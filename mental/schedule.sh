#!/bin/bash

genopool=/cluster/projects/p33/projects/mental/geno/generic_qc
chunksize=20000
export phenodata=/cluster/projects/p33/projects/mental/pheno/mostest_mental_pheno_200807.csv
export phenotype=/cluster/projects/p33/projects/mental/pheno/pheno_type_subset.txt
export covardata=/cluster/projects/p33/projects/mental/covar/mostest_mental_covar_200807.csv
export outfolder=/cluster/projects/p33/projects/mental/results/sumstat

mkdir -p $outfolder

n_pheno=`head -n1 $phenodata | awk -F ',' '{print NF-1}'`
for ((j=1; j<=n_pheno; j++)); do
    export phenoname=`head -n1 $phenodata | awk -v n=$((j+1)) -F ',' '{print $n}'`
    if [ `awk -v fieldid=${phenoname%%-*} '$1==fieldid' $phenotype | wc -l` -eq 0 ]; then
        continue
    fi
    for fn in $genopool/*.fam; do
        export genodata=${fn%%.*}
        bn=$(basename $genodata)
        bn=${bn##*chr}
        bn=${bn%%_*}
        n_snps=`cat $genodata.bim | wc -l`
        for ((i=1; i<=n_snps; i=i+chunksize)); do
            from=$i
            to=$((i+chunksize-1))
            if [ $to -gt $n_snps ]; then
                to=$n_snps
            fi
            export from to
            export snplist=$from:$to
            echo $bn $snplist $phenoname
            #sh run_assoc.job
            sbatch --job-name ${phenoname}_${bn}_${from} run_assoc.job
            break
         done
        break
    done
    break
done
