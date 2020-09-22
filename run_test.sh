#!/bin/bash

/mnt/seagate10/matlab/bin/matlab -nodisplay -nosplash -r "bfile='chr21'; pheno='pheno.txt'; num_eigenval_to_keep=3; out='chr21.test.stats'; mostest_light_stats; stats_file='chr21.test.stats.mat'; out='chr21.test.pvals'; daletails_quantile=0.99; mostest_light_pvals; exit"

