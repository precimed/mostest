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
    mat = sio.loadmat(fname + '.mat')
    with np.errstate(divide='ignore'):
        df = pd.DataFrame({
        'minp_x':np.transpose(mat['hv_maxlogpvecs'].flatten()),
        'minp_y1':np.transpose(-np.log10(1-mat['chc_maxlogpvecs'])).flatten(),
        'minp_y2':np.transpose(-np.log10(1-mat['cdf_minpvecs'])).flatten(),
        'most_x':mat['hv_mostvecs'].flatten(),
        'most_y1':np.transpose(-np.log10(1-mat['chc_mostvecs'])).flatten(),
        'most_y2':np.transpose(-np.log10(1-mat['cdf_mostvecs'])).flatten() })
    df.to_csv(out + '.plot.csv',index=False, sep='\t')
    plt.figure(figsize=(20, 10), dpi=100)
    plt.subplot(2,4,1)
    plt.plot(df['minp_x'], df['minp_y1'])
    plt.plot(df['minp_x'], df['minp_y2'])
    plt.legend(['data (null)', 'minP (model)'])
    plt.title('minP')
    plt.subplot(2,4,2)
    plt.plot(df['most_x'], df['most_y1'])
    plt.plot(df['most_x'], df['most_y2'])
    plt.title('MOSTest')
    plt.legend(['data (null)', 'MOSTest (model)'])
    plt.savefig(out + '.plot.png', bbox_inches='tight')

    # generate .sumstats files, compatible with FUMA
    print('Generate {}.***.sumstats files...'.format(out))
    for test in ['minp_log10pval_orig', 'most_log10pval_orig', 'minp_log10pval_perm', 'most_log10pval_perm']:
        bim['PVAL'] = np.power(10, -mat[test].flatten())
        bim['Z'] = -stats.norm.ppf(bim['PVAL'].values*0.5) #*effect_sign.astype(np.float64) - effect size not available from MOSTest and minP
        bim['N'] = mat['nvec']
        bim.to_csv('{}.{}.sumstats'.format(out, test.replace('_log10pval', '')), sep='\t', na_rep="NA", index=False)

    print('Done.')

