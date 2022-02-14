#!/bin/bash

CHROM=$(seq 1 22)
SNP_IN_CHUNK=50000
IN_DIR=/cluster/projects/p33/projects/mental/geno/generic_qc
OUT_DIR=/cluster/projects/p33/projects/mental/geno/generic_qc/snp_chunks

for chr in ${CHROM}; do
    echo ${chr}
    split -l ${SNP_IN_CHUNK} --numeric-suffixes --suffix-length 2 --additional-suffix ".txt" <(cut -f2 ${IN_DIR}/ukb_imp_chr${chr}_v3_qc.bim) ${OUT_DIR}/chr${chr}_snps_
done

