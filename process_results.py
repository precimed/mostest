import pandas as pd
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import scipy.stats as stats
import h5py
import sys

if __name__ == '__main__':
    if len(sys.argv) <= 2:
        print('Usage: process_results.py <bim> <fname> [<out>], where')
        print(' bim   - path to bim file (reference set of SNPs')
        print(' fname - prefix of .mat files output by mostest.m, ie. fname should be the same as "out" argument of the mostest.m')
        print(' out   - optional suffix for output files, by defautl fname will be used')
        sys.exit()

    bim_file = sys.argv[1] #'UKB26502_QCed_230519_maf0p005_chr21.bim'
    fname = sys.argv[2]    # 'all_chr21'
    out = sys.argv[3] if (len(sys.argv) > 3) else sys.argv[2]

    # read .bim file (reference set of SNPs)
    print('Load {}...'.format(bim_file))
    bim = pd.read_csv(bim_file, sep='\t', header=None, names='CHR SNP GP BP A1 A2'.split())
    del bim['GP']

    # save figures with  QQ plots MOST and minP
    print('Generate {}.plot.png...'.format(out))
    plt.figure(figsize=(10, 5), dpi=100)
    mat = sio.loadmat(fname + '.mat')
    chc_maxlogpvecs = -np.log10(1-np.cumsum(np.transpose(mat['hc_maxlogpvecs']))/np.sum(mat['hc_maxlogpvecs']))
    df = pd.DataFrame({
    'minp_x':np.transpose(mat['hv_maxlogpvecs'].flatten()),
    'minp_y1':np.transpose(chc_maxlogpvecs).flatten(),
    'minp_y2':np.transpose(-np.log10(1-mat['cdf_minpvecs'])).flatten(),
    'most_x':mat['hv_logpdfvecs'].flatten(),
    'most_y1':np.transpose(-np.log10(1-mat['chc_logpdfvecs'])).flatten(),
    'most_y2':np.transpose(-np.log10(1-mat['cdf_logpdfvecs'])).flatten() })
    df.to_csv(out + '.plot.csv',index=False, sep='\t')
    plt.subplot(2,4,1)
    plt.plot(df['minp_x'], df['minp_y1'])
    plt.plot(df['minp_x'], df['minp_y2'])
    plt.legend(['data (null)', 'minP (model)'])
    plt.title(fname + ' (minP)')
    plt.subplot(2,4,2)
    plt.plot(df['most_x'], df['most_y1'])
    plt.plot(df['most_x'], df['most_y2'])
    plt.title(fname + ' (MOSTest)')
    plt.legend(['data (null)', 'MOSTest (model)'])
    plt.savefig(out + '.plot.png', bbox_inches='tight')

    # generate .sumstats files, compatible with FUMA
    print('Generate {}.***.sumstats files...'.format(out))
    for test in ['minPval', 'mostPval']:   # in fact mat[test] stores -log10(pval)
        bim['PVAL'] = np.power(10, -mat[test].flatten())
        bim['Z'] = -stats.norm.ppf(bim['PVAL'].values*0.5) #*effect_sign.astype(np.float64) - effect size not available from MOSTest and minP
        bim['N'] = mat['nvec']
        bim.to_csv('{}.{}.sumstats'.format(out, test.replace('minPval', 'minp').replace('mostPval', 'most')), sep='\t', index=False)

    # save maxtrix of z scores for SNPs passing 5e-08 threshold
    print('Generate {}_***.zmat.csv files...'.format(out))
    with h5py.File(fname + '_zmat.mat', 'r') as h5file:
        # measures = pd.read_csv(pheno, sep='\t').columns
        measures = [u''.join(chr(c) for c in h5file[h5file['measures'][i, 0]]) for i in range(0, h5file['measures'].shape[0])]
        zmat_orig=np.array(h5file['zmat_orig'])
        for test in ['minPval', 'mostPval']:   # in fact mat[test] stores -log10(pval)
            pval = np.power(10, -mat[test].flatten())
            df_zmat = pd.DataFrame(np.transpose(zmat_orig[:, pval<5e-08]), columns=measures)
            df_zmat.insert(0, 'SNP', bim.SNP.values[pval<5e-08])
            df_zmat.to_csv(out + test.replace('minPval', '_minp').replace('mostPval', '_most') + '.zmat.csv', index=False, sep='\t')

    print('Done.')

