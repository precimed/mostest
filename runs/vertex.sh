#!/bin/bash

# Name your job. Unless you use the -o and -e options, output will
# go to a unique file name.ojob_id for each job.
#$ -N mostvert

# Execute job in the queue "std.q" unless you have special requirements.
##$ -q std.q
##$ -q all_24.q

# The SGE batch system uses the current directory as working directory.
# Both files (output.dat and error.dat) will be placed in the current
# directory. The batch system assumes to find the executable in this directory.
#$ -cwd

# Redirect output stream to this file.
##$ -o output.dat

# Redirect error stream to this file.
##$ -e error.dat

# Send status information to this email address.
##$ -M oleksandr.frei@gmail.com

# Send an e-mail when the job is done.
##$ -m e

#$ -l h_vmem=120G
#$ -l h_rt=36:00:00
##$ -pe dmp4 16

#$ -l h=!mmil-compute-5-6.local
#$ -l h=!mmil-compute-5-7.local
#$ -l h=!mmil-compute-5-14.local


#tot -keep=      n

#2384-5          2379
#2384-10         2374
#2384-20         2364
#2384-30         2354
#2384-40         2344


#2384-50         2334
#2384-100        2284
#2384-200        2184
#2384-500        1884
#2384-800        1584
#2384-1400       984

bash
cd /home/oleksandr/precimed/mostest
\matlab -nodisplay -nosplash -r "prefix = '/space/gwas-syn1/1/data/GWAS/UKBioBank/projects/mostest_vertex/';pheno = fullfile(prefix, 'correct_resid_format/lh.pial.15.fsaverage4.QCed.resid.csv');bfile = fullfile(prefix, 'UKB39k_QCed_051119_mat0p05');out = fullfile(prefix, 'results3/light_lh.pial.mat0p05_meta50');meta_simu_num_cohors=50;num_eigval_to_regularize=0;mostest_light;exit;"

