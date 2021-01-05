folder=$mt3/hapmap-ordinal-64-pheno/raw
outfile=$yc/mental/hapmap-ordinal-64-pheno-raw.csv
n=1
for i in $folder/*.csv; do
    echo $i
    if [ $n -eq 1 ]; then
        cut -f1-3,8 $i | sed 's/Z/Z_1/' > $outfile
    else
        cut -f8 $i | sed "s/Z/Z_$n/" > $outfile.tmp
        paste $outfile $outfile.tmp > $outfile.tmp2
        mv $outfile.tmp2 $outfile
        rm $outfile.tmp
    fi
    n=$((n+1))  
done
