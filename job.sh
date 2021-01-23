#!/bin/bash



simu_dir=/home/aihual/Simu_1
plink_dir=/home/aihual/Plink/plink
most_dir=/home/aihual/mostest_meta


#################### Step 1
## This simulation only used chr21
chr=chr21

## generate effect variance for each trait
no_traits=10
hsq=0.004 ##for nc=100,  0.04 for nc =1000
cd ${simu_dir}
matlab -nodesktop -nodisplay -r "run snp_make.m;exit"



## calculate pheno data
##for res in 1 2 3 4 5 6 7 8 9 10
for res in $(seq 1 $no_traits)
do
	snp_list=${simu_dir}/snp_${res}.txt
	
	./simu_linux --bfile ${simu_dir}/${chr} --qt --causal-variants $snp_list --hsq ${hsq}

	mv simu.pheno pheno_${res}
done

## combine pheno data to a single file, pheno.txt. Alternatively run pheno_combine.m
for res in $(seq 1 $no_traits)
do	
	awk '{print $3}' pheno_$res > p_$res
	sed -i "1s/.*/trait${res}/" p_$res 
done

mkdir phenos
cp p_* $simu_dir/phenos
cd $simu_dir/phenos
paste * | column -s $'\t' -t > pheno.txt

cd ..
################## Step 2: generate n_cohort sub geno data sets randomly, and split pheno.txt and subjects (fam file) into n_cohort cohorts
################## subj_X.txt, cohort_X.txt, snps_X.txt
matlab -nodesktop -nodisplay -r "run splitdata.m; exit"
				
################## Step 3: generate bim, bed, fam files for each cohort
cd ${plink_dir}
n_cohort=3
for res in $(seq 1 $n_cohort)
do
	./plink --bfile ${simu_dir}/${chr} --keep ${simu_dir}/subj_${res}.txt  --make-bed --extract ${simu_dir}/snps_${res}.txt --out ${simu_dir}/cohort${res}
done

################## Step 4: align snps of all cohorts into a comman snps list
for res in $(seq 2 ${n_cohort})
do	
	echo "${simu_dir}/cohort${res}.bed ${simu_dir}/cohort${res}.bim ${simu_dir}/cohort${res}.fam" >> ${simu_dir}/list.txt
done

./plink --bfile ${simu_dir}/cohort1 --merge-list ${simu_dir}/list.txt --make-bed --allow-no-sex --out ${simu_dir}/merged

################## Step 5: extract sub bim,bed,fam files of cohorts without extracting the snps
for res in $(seq 1 $n_cohort)
do
	./plink --bfile ${simu_dir}/merged --keep ${simu_dir}/subj_${res}.txt  --make-bed --out ${simu_dir}/cohort_$res
done


################## Step 6: create new folder input_bfile
cd ${most_dir}
mkdir input_bfile input_pheno
cp ${simu_dir}/cohort_*.bim ${simu_dir}/cohort_*.fam ${simu_dir}/cohort_*.bed ${most_dir}/input_bfile
cp ${simu_dir}/cohort_*.txt ${most_dir}/input_pheno

################## Step 7: run mostest_meta_std.
matlab -nodesktop -nodisplay -r "run mostest_meta_std.m; exit"

