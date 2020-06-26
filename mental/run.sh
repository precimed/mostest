#!/bin/bash

genodata=testdata/test
chunksize=5
phenodata=testdata/test_25000_pheno.txt
phenotype=testdata/test_25000_pheno_types.txt
covardata=testdata/test_25000_covar.txt
outfolder=testdata

n_snps=`cat $genodata.bim | wc -l`
n_pheno=`head -n1 $phenodata | awk '{print NF-2}'`
for ((i=1; i<=n_snps; i=i+chunksize)); do
    from=$i
    to=$((i+chunksize-1))
    if [ $to -gt $n_snps ]; then
        to=$n_snps
    fi
    snplist=$from:$to
    echo $snplist
    for ((j=1; j<=n_pheno; j++)); do
        phenoname=`head -n1 $phenodata | awk -v n=$((j+2)) '{print $n}'`
        echo $phenoname
        Rscript model.R $genodata $snplist $phenodata $phenoname $phenotype $covardata $outfolder/$(basename $genodata)_${from}-${to}_${phenoname}_mostest.txt
    done
    break
done
