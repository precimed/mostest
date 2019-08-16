## Contents

* [Install MOSTest](#install-mostest)
* [Data preparation](#data-preparation)
* [Run MOSTest](#run-mostest)

## Install MOSTest

To install MOSTest you simply download this repository, or fetch it via git.
MOSTest is implemented as a script in MATLAB.
It's tested with MATLAB 2017a, but it should work well with other versions as well.

## Data preparation

MOSTest require the following input files

1. Plink bfile (.bim, .bed and .fam files) with genotypes to test for association
2. Phenotype file (tab separated table, first line must be a header, with rows corresponding to individuals, and columns indicating phenotype measures for those individuals)

Rows in the phenotype file must correspond to the same set of individuals, in the same order, as the .fam file.
All phenotypes will be analyzed together.

## Run MOSTest
Set your current folder ``<MOST_ROOT>``, i.e. the folder containing ``mostest.m``.
Open matlab, define ``pheno``, ``out``, ``bfile``, ``snps`` and ``nsubj`` variables, and execute ``mostest.m`` script:
```
pheno = 'SubcorticalVolume.csv';
out = 'SubcorticalVolume';
bfile = 'UKB26502_QCed_230519_maf0p005';
snps = 7428630; nsubj = 26502;
mostest
```
