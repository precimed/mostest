## Contents

* [Install MOSTest](#install-mostest)
* [Data preparation](#data-preparation)
* [Run MOSTest](#run-mostest)
* [Other considerations](#other-considerations)

## Install MOSTest

To install MOSTest you simply download or git fetch this repository into an empty folder (which we refer to as `<MOST_ROOT>`),
MOSTest is implemented as a script in MATLAB.
Tested with MATLAB 2017a, but it should work just as well with other MATLAB versions.

## Data preparation

MOSTest require the following input files

1. plink bfile (.bim, .bed and .fam files) with genotypes to test for association
2. Phenotype file, stored as a tab-separated table. First line must be a header, and subsequent rows corresponding to individuals. All columns must indicate phenotype measures to be jointly analyzed by MOSTest.

NB (!) Rows in the phenotype file must correspond to the same set of individuals, in exactly the same order, as the .fam file of you bfile argument.

## Run MOSTest
Set your current folder to ``<MOST_ROOT>``, i.e. to the folder that contains ``mostest.m``.
Then start matlab, define ``pheno``, ``out``, ``bfile``, ``snps`` and ``nsubj`` variables as shown below,
and execute ``mostest.m`` script:
```
pheno = 'SubcorticalVolume.csv';         % full or relative path to the phenotype file
bfile = 'UKB26502_QCed_230519_maf0p005'; % full or relative path to plink bfile prefix
out = 'SubcorticalVolume';               % prefix for the output files
snps = 7428630; nsubj = 26502;           % number of snps and subjects in bfile
mostest                                  % starts the analysis
```

## MOSTest results

The output is stored in files
```
SubcorticalVolume.mat      - main results, i.e. -log10(pval)
SubcorticalVolume_zmat.mat - gwas Z scores (original and permuted)
```

See [process_results.py](process_results.py) for an example of a post-processing script that coverts results from .mat to .csv files.

## Other considerations

* all phenotypes should be pre-residualized for global measures
* all phenotypes should be inverse-rank transformed into a normal distribution, for example via qnorm(ecdf(...)) in R
* we recommend that bfile only include SNPs with MAF above 0.5%
