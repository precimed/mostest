## Contents

* [Prerequisites](#prerequisites)
* [Install MOSTest](#install-mostest)
* [Data preparation](#data-preparation)
* [Run MOSTest](#run-mostest)
* [Other considerations](#other-considerations)

## Introduction

MOSTest is a tool for join genetical analysis of multiple traits, using multivariate analysis to boost the power of discovering associated loci. If you use MOSTest software for your research publication, please cite the following paper(s):

* Dennis van der Meer, Oleksandr Frei, et al. Understanding the genetic determinants of the brain with MOSTest. 
  Nat Commun 11, 3512 (2020). https://doi.org/10.1038/s41467-020-17368-1

The MOSTest software may not be used for commercial purpose or in medical applications.
We encourage all users to familiarize themselves with US patent https://www.google.no/patents/US20150356243 "Systems and methods for identifying polymorphisms".

## Prerequisites

MOSTest is implemented as a script in MATLAB and python. We tested MOSTest with the following software configuration:
* Ubuntu 18.04.
* MATLAB (tested with 2017a)
* Python (tested with 3.7); required libraries: ``pandas, numpy, h5py``

Other versions are likely to work well too. 
Since MATLAB and Python are cross-platform we expect MOSTest to run well on Windows and MacOS too.

## Install MOSTest

To install MOSTest you simply download or git fetch this repository (https://github.com/precimed/mostest) into an empty folder (which we refer to as `<MOST_ROOT>`),

## Data preparation

MOSTest require the following input files

1. plink bfile (.bim, .bed and .fam files) with genotypes to test for association
2. Phenotype file, stored as a tab-separated table. First line must be a header, and subsequent rows corresponding to individuals.
   All columns must indicate phenotype measures to be jointly analyzed by MOSTest.

NB (!) Rows in the phenotype file must correspond to the same set of individuals,
in exactly the same order, as the .fam file of you bfile argument.
NB2 Missing values should be coded as 'NaN' in order to be handled correctly

## Run MOSTest (demo)
Download ``mostest_demo.tar.gz`` file from [here](https://1drv.ms/u/s!Ai1YZmdFa9ati40Inztrv_4erqcdWw?e=ixWDUe)
and extract it into ``<MOST_ROOT>>``, i.e. to the folder that contains ``mostest.m``.
This package contains a small synthetic dataset for chromosome 21, N=10000 subjects, M=149454 SNPs, K=10 phenotypes.

To run MOSTest, set your current folder to ``<MOST_ROOT>``.
Then start matlab, define ``pheno``, ``out``, ``bfile``, ``snps`` and ``nsubj`` variables as shown below,
and execute ``mostest.m`` script:
```
pheno = 'pheno.txt';            % full or relative path to the phenotype file
bfile = 'chr21';                % full or relative path to plink bfile prefix
out = 'results';                % prefix for the output files
mostest                         % starts the analysis
```

Alternatively one may use ``zmat_name`` argument to re-use the original and permuted z-scores from previous MOSTest run:
```
zmat_name='results_zmat.mat'; out = 'results'; mostest
```

To convert MOSTest results from ``.mat`` format to text files, use the following script in python:
```
python process_results.py chr21.bim results
python process_results_ext.py chr21.bim results
```

[process_results.py](process_results.py) and [process_results_ext.py](process_results_ext.py) scripts accept the following arguments:
```
Usage: process_results.py <bim> <fname> [<out>], where
 bim   - path to bim file (reference set of SNPs
 fname - prefix of .mat files output by mostest.m, 
         i.e. fname should be the same as "out" argument of the mostest.m
 out   - optional suffix for output files, by defautl fname will be used
```
The ``process_results.py`` produces text files with MOSTest and MinP p-values, while
``process_results_ext.py`` produce all univariate GWAS results across all input phenotypes.

The demo takes few minutes to complete. Output of the demo is as follows, highlighting
14436 genome-wide significant SNPs discovered by minP, and 20986 by MOSTest.
```
Loading phenotype matrix from /home/oleksanf/vmshare/data/mostest_demo/pheno.txt... Done, 10 phenotypes found
Perform GWAS on /home/oleksanf/vmshare/data/mostest_demo/chr21 (149454 SNPs are expected)...
gwas: loading snps 1 to 10000... processing... done in 6.6 sec, 6.7 % completed
gwas: loading snps 10001 to 20000... processing... done in 4.6 sec, 13.4 % completed
gwas: loading snps 20001 to 30000... processing... done in 4.6 sec, 20.1 % completed
gwas: loading snps 30001 to 40000... processing... done in 4.6 sec, 26.8 % completed
gwas: loading snps 40001 to 50000... processing... done in 4.6 sec, 33.5 % completed
gwas: loading snps 50001 to 60000... processing... done in 4.6 sec, 40.1 % completed
gwas: loading snps 60001 to 70000... processing... done in 4.6 sec, 46.8 % completed
gwas: loading snps 70001 to 80000... processing... done in 4.5 sec, 53.5 % completed
gwas: loading snps 80001 to 90000... processing... done in 4.6 sec, 60.2 % completed
gwas: loading snps 90001 to 100000... processing... done in 4.6 sec, 66.9 % completed
gwas: loading snps 100001 to 110000... processing... done in 4.6 sec, 73.6 % completed
gwas: loading snps 110001 to 120000... processing... done in 4.7 sec, 80.3 % completed
gwas: loading snps 120001 to 130000... processing... done in 4.7 sec, 87.0 % completed
gwas: loading snps 130001 to 140000... processing... done in 4.8 sec, 93.7 % completed
gwas: loading snps 140001 to 149454... processing... done in 4.4 sec, 100.0 % completed
saving /home/oleksanf/vmshare/data/mostest_demo/results_zmat.mat as -v7.3... OK.
running MOSTest analysis... Done.
GWAS yield minP: 14436; MOST: 20986
10	2.77	0.991	9.659	4.993	2.003	
saving /home/oleksanf/vmshare/data/mostest_demo/results.mat... Done.
MOSTest analysis is completed.

(base) oleksanf@mach:~/vmshare/data/mostest_demo$ python ~/github/mostest/process_results.py chr21.bim results 
Load chr21.bim...
Generate results.plot.png...
Generate results.***.sumstats files...
Generate results_***.zmat.csv files...
Done.
```

Note that it is also possible to use ``mostest_light.m`` script, instead of ``mostest.m``.
``mostest_light.m`` uses much less memory by skip saving univariate results.
However it is not possible to run ``process_results_ext.py`` on the results from ``mostest_light.m``.

## MOSTest results

MOSTest results are saved in .mat files (MATLAB format), in files called
```
<out>.mat      - main results, i.e. -log10(pval)
<out>_zmat.mat - gwas Z scores (original and permuted)
```

Output of ``process_results.py``:
* ``<out>.most.sumstats`` - summary statistics file for MOSTest p-values (note that Z values has a fake direction)
* ``<out>.minp.sumstats`` - as above, but based on minP
* ``<out>_most.zmat.csv`` - GWAS z-scores across all measures for genome-wide significant SNPs according to MOSTest
* ``<out>_minp.zmat.csv`` - as above, but for SNPs passing 5e-8 according to minP
* ``<out>.plot.png``      - QQ-like plot highlighting distribution of minP and MOSTest p-values 
* ``<out>.plot.csv``      - data underlying the above plot

## Other considerations

* all phenotypes should be pre-residualized for covariates, for example in R as follows:
  ```
  TempModel <- summary(lm(get(f) ~ Age + Sex + Scanner + Euler + genBatch + C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C10, data=DM)) 
  DM[,which(colnames(DM)==f)] <- TempModel$coefficients[1,1] + TempModel$residuals
  ```

* by default all phenotypes are automatically inverse-rank transformed into a normal distribution, 
  in a way equivalent to the qnorm(ecdf(...)) function in R:
  ```
  DM[,f] <- qnorm(ecdf(DM[,f])(DM[,f]) - 0.5/dim(DM)[1])
  ```
  If you want to disable this step, specify ``apply_int=0; ``  in your matlab script before running ``mostest``.

* we recommend that bfile only include SNPs with MAF above 0.5%

* MOSTest spends most time in matrix operations, which by default in MATLAB utilize all available computation threads
  (unless MATLAB is started with ``-singleCompThread`` flag). Therefore, for large-scale analysis we recommend to run
  MOSTest on a large server with 16 or more CPUs available.
