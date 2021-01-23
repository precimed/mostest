Mostest-Meta simulations

mostest_meta_std is a tool to run Mostest on meta analysis based on inverse variance method. 
The following steps is a simulation on how it works. 

Input: pheno data and geno data sets
Output: sumstats from minp and Mostest

Steps
1. Generate pheno.txt (pheno_maga) from geno data (geno_maga) using SIMU
2. Run splitdata.m to generate n_cohort sub geno data sets randomly, and split pheno.txt and subjects (fam file) into n_cohort cohorts
( subj_X.txt,	snps_X.txt and pheno data cohort_X.txt)
note: X=1:n_cohort

3. Using plink to generate bim, bed, fam files for each cohort
./plink --bfile /path/to/chr21 --keep /path/to/subj_X.txt  --make-bed --extract /path/to/snps_X.txt --out /path/to/cohortX


4. Using merge_list to align snps of all cohorts into a common snps list
./plink --bfile /path/to/cohort1 --merge-list /path/to/list.txt --make-bed --allow-no-sex --out /path/to/merged

note: list.txt is a list of .bed .bim .fam files of all cohorts from 2 - no.cohort

5. Extract sub bim,bed,fam files of cohorts without extracting the snps
./plink --bfile /path/to/merged --keep /path/to/subj_X.txt  --make-bed --out /path/to/cohort_X

6. Create new folder input_bfile containing geno data of all cohorts (cohort_X.bim, cohort_X.fam and cohort_X.bed )and folder input_pheno containing pheno data of all cohorts (cohort_X.txt)

It’s important that the corresponding pheno data and geno data have the same prefix in names.
So that the main code mostest_meta_std reads those files in the right order. 

7. Run mostest_meta_std.


Above steps can be run by job.sh




