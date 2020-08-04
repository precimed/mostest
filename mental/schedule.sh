#!/bin/bash

#genopool=/tsd/p33/data/durable/s3-api/ukbdata/projects/plsa_mixer/ukb_genetics_qc/ukb_bed
genopool=testdata
chunksize=5
export phenodata=testdata/test_25000_pheno.txt
export phenotype=testdata/test_25000_pheno_types.txt
export covardata=testdata/test_25000_covar.txt
export outfolder=results

mkdir -p $outfolder

for fn in $genopool/*.fam; do
    export genodata=${fn%%.*}
    n_snps=`cat $genodata.bim | wc -l`
    n_pheno=`head -n1 $phenodata | awk '{print NF-2}'`
    for ((i=1; i<=n_snps; i=i+chunksize)); do
        export from=$i
        export to=$((i+chunksize-1))
        if [ $to -gt $n_snps ]; then
            export to=$n_snps
        fi
        export snplist=$from:$to
        for ((j=1; j<=n_pheno; j++)); do
            export phenoname=`head -n1 $phenodata | awk -v n=$((j+2)) '{print $n}'`
            echo $(basename $genodata) $snplist $phenoname
            #sh run_assoc.job
            sbatch --job-name ${from}_${to}_$phenoname run_assoc.job
            break
        done
        break
    done
    break
done
