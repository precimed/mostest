import pandas as pd
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import scipy.stats as stats
import h5py

types = ['CorticalArea', 'CorticalThickness', 'SubcorticalVolume']

# read .bim file (reference set of SNPs)
bim = pd.read_csv('UKB26502_QCed_230519_maf0p005.bim', sep='\t', header=None, names='CHR SNP GP BP A1 A2'.split())
del bim['GP']

# save figures with  QQ plots MOST and minP
plt.figure(figsize=(20, 10), dpi=100)
for index, fname in enumerate(types + ['all']):
    mat = sio.loadmat(fname + '.mat')
    chc_maxlogpvecs = -np.log10(1-np.cumsum(np.transpose(mat['hc_maxlogpvecs']))/np.sum(mat['hc_maxlogpvecs']))
    df = pd.DataFrame({
        'minp_x':np.transpose(mat['hv_maxlogpvecs'].flatten()),
        'minp_y1':np.transpose(chc_maxlogpvecs).flatten(),
        'minp_y2':np.transpose(-np.log10(1-mat['cdf_minpvecs'])).flatten(),
        'most_x':mat['hv_logpdfvecs'].flatten(),
        'most_y1':np.transpose(-np.log10(1-mat['chc_logpdfvecs'])).flatten(),
        'most_y2':np.transpose(-np.log10(1-mat['cdf_logpdfvecs'])).flatten() })
    df.to_csv(fname + '.plot.csv',index=False, sep='\t')
    plt.subplot(2,4,1+index)
    plt.plot(df['minp_x'], df['minp_y1'])
    plt.plot(df['minp_x'], df['minp_y2'])
    plt.legend(['data (null)', 'minP (model)'])
    plt.title(fname + ' (minP)')
    plt.subplot(2,4,5+index)
    plt.plot(df['most_x'], df['most_y1'])
    plt.plot(df['most_x'], df['most_y2'])
    plt.title(fname + ' (MOSTest)')
    plt.legend(['data (null)', 'MOSTest (model)'])
plt.savefig('plot.png', bbox_inches='tight')
    
# generate .sumstats files, compatible with FUMA
for index, fname in enumerate(types + ['all']):
    mat = sio.loadmat(fname + '.mat')
    for test in ['minPval', 'mostPval']:   # in fact mat[test] stores -log10(pval)
        bim['PVAL'] = np.power(10, -mat[test].flatten())
        bim['Z'] = -stats.norm.ppf(bim['PVAL'].values*0.5) #*effect_sign.astype(np.float64) - effect size not available from MOSTest and minP
        bim['N'] = mat['nvec']
        bim.to_csv('{}.{}.sumstats'.format(fname, test.replace('minPval', 'minp').replace('mostPval', 'most')), sep='\t', index=False)

# save maxtrix of z scores for SNPs passing 5e-08 threshold
for fname in types + ['all']:
    with h5py.File(fname + '_zmat.mat', 'r') as h5file:
        measures = pd.read_csv(fname+'.csv', sep='\t').columns
        # old way (cell array of strings within .mat file)
        # measures = [u''.join(chr(c) for c in h5file[h5file['measures'][0, i]]) for i in range(0, 250)]
        mat = sio.loadmat(fname + '.mat')
        zmat_orig=np.array(h5file['zmat_orig'])

        for test in ['minPval', 'mostPval']:   # in fact mat[test] stores -log10(pval)
            pval = np.power(10, -mat[test].flatten())
            df_zmat = pd.DataFrame(np.transpose(zmat_orig[:, pval<5e-08]), columns=measures)
            df_zmat.insert(0, 'SNP', bim.SNP.values[pval<5e-08])
            df_zmat.to_csv(fname + test.replace('minPval', '_minp').replace('mostPval', '_most') + '.zmat.csv', index=False, sep='\t')

