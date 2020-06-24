for ((i=1;i<=10;i++)); do
    Rscript model.R testdata/test 1:5 testdata/test_25000_pheno.txt pheno_$i testdata/test_25000_pheno_types.txt testdata/test_25000_covar.txt testdata/test_pheno_${i}_mostest.txt
done
