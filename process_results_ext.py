import pandas as pd
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import scipy.stats as stats
import h5py
import sys
import os

if __name__ == '__main__':
    if len(sys.argv) <= 2:
        print('Usage: process_results_ext.py <bim> <fname> [<out>], where')
        print(' bim   - path to bim file (reference set of SNPs')
        print(' fname - prefix of .mat files output by mostest.m, ie. fname should be the same as "out" argument of the mostest.m')
        print(' out   - optional suffix for output files, by defautl fname will be used')
        print(' mode  - "orig" to export non-permuted summary stats; "perm" for permuted')
        sys.exit()

    bim_file = sys.argv[1] #'UKB26502_QCed_230519_maf0p005_chr21.bim'
    fname = sys.argv[2]    # 'all_chr21'
    out = sys.argv[3] if (len(sys.argv) > 3) else sys.argv[2]
    mode = sys.argv[4] if (len(sys.argv) > 4) else 'orig'

    if mode not in ['orig', 'perm']: print('unknown mode, check 4th argument'); sys.exit()

    # read .bim file (reference set of SNPs)
    print('Load {}...'.format(bim_file))
    bim = pd.read_csv(bim_file, sep='\t', header=None, names='CHR SNP GP BP A1 A2'.split())
    del bim['GP']

    mat = sio.loadmat(fname + '.mat')
    bim['N'] = mat['nvec']

    # save matrix of z scores for SNPs passing 5e-08 threshold
    print('Generate {}_***.zmat.csv files...'.format(out))
    with h5py.File(fname + '_zmat.mat', 'r') as h5file:
        # measures = pd.read_csv(pheno, sep='\t').columns
        measures = [u''.join(chr(c) for c in h5file[h5file['measures'][i, 0]]) for i in range(0, h5file['measures'].shape[0])]
        zmat_orig=np.array(h5file['zmat_orig'])
        for test in ['minp_log10pval_orig', 'most_log10pval_orig']:
            pval = np.power(10, -mat[test].flatten())
            df_zmat = pd.DataFrame(np.transpose(zmat_orig[:, pval<5e-08]), columns=measures)
            df_zmat.insert(0, 'SNP', bim.SNP.values[pval<5e-08])
            df_zmat.to_csv(out + "_" + test.replace('_log10pval', '') + '.zmat.csv', index=False, sep='\t')

        # save individual GWAS results ('freqvec' is an indivator that we've saved individual GWAS beta's)
        if 'freqvec' in h5file:
            bim['FRQ'] = np.transpose(np.array(h5file['freqvec']))

            beta_orig = np.array(h5file['beta_orig'])
            se_orig = np.divide(beta_orig, zmat_orig)
            #pval_orig = stats.norm.sf(np.abs(zmat_orig))*2.0

            for measure_index, measure in enumerate(measures):
                if mode != 'orig': break
                fname = '{}.{}.orig.sumstats.gz'.format(out, measure)
                print('Generate {}...'.format(fname))
                #bim['PVAL'] = np.transpose(pval_orig[measure_index, :])
                #bim['Z'] = np.transpose(zmat_orig[measure_index, :])
                bim['BETA'] = np.transpose(beta_orig[measure_index, :])
                bim['SE'] = np.transpose(se_orig[measure_index, :])
                bim["SNP A1 A2 N FRQ BETA SE".split()].to_csv(fname,  sep='\t', index=False)
            del zmat_orig; del beta_orig; del se_orig; #del pval_orig

            beta_perm = np.array(h5file['beta_perm'])
            zmat_perm = np.array(h5file['zmat_perm'])
            se_perm = np.divide(beta_perm, zmat_perm)
            pval_perm = stats.norm.sf(np.abs(zmat_perm))*2.0
            for measure_index, measure in enumerate(measures):
                if mode=='orig': break
                fname = '{}.{}.perm.sumstats.gz'.format(out, measure)
                print('Generate {}...'.format(fname))
                #bim['PVAL'] = np.transpose(pval_perm[measure_index, :])
                #bim['Z'] = np.transpose(zmat_perm[measure_index, :])
                bim['BETA'] = np.transpose(beta_perm[measure_index, :])
                bim['SE'] = np.transpose(se_perm[measure_index, :])
                bim["SNP A1 A2 N FRQ BETA SE".split()].to_csv(fname,  sep='\t', index=False)
    print('Done.')

