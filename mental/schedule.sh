#!/bin/bash

genopool=/cluster/projects/p33/projects/mental/geno_hapmap/by_chunk_1K
np=1
chr_from=1
chr_to=22
chunksize=1000
export phenodata=/cluster/projects/p33/projects/mental/pheno/mostest_mental_pheno_ordinal_65_200907.negative2na.csv
export phenotype=/cluster/projects/p33/projects/mental/pheno/pheno_type.txt
export covardata=/cluster/projects/p33/projects/mental/covar/mostest_mental_covar_200807.csv
export outfolder=/cluster/projects/p33/projects/mental/results/sumstat3

mkdir -p $outfolder

n_pheno=`head -n1 $phenodata | awk -F '\t' '{print NF-1}'`
for ((j=$np; j<=n_pheno; j++)); do
    export phenoname=`head -n1 $phenodata | awk -v n=$((j+1)) -F '\t' '{print $n}'`
    echo $phenoname
    if [ `awk -v fieldid=${phenoname%%-*} '$1==fieldid' $phenotype | wc -l` -eq 0 ]; then
        continue
    fi
    for chk in $genopool/ukb_imp_v3_qc.chunk_*.bim; do
    #for chk in $genopool/ukb_imp_v3_qc.perm.chunk_*.bim; do
        export genodata=${chk%.*}
        n_snps=`cat $genodata.bim | wc -l`
        for ((i=1; i<=n_snps; i=i+chunksize)); do
            from=$i
            to=$((i+chunksize-1))
            if [ $to -gt $n_snps ]; then
                to=$n_snps
            fi
            export from to
            export snplist=$from:$to
            bn=$(basename $chk)
            bn=${bn%.*}
            bn=${bn##*chunk_}
            echo $chr $bn $snplist $phenoname
            sbatch --job-name ${phenoname}_$bn run_assoc.job
            #sbatch --nice=1000 --job-name ${phenoname}_$bn run_assoc.job
            #break
        done
        #break
    done
    break
done
