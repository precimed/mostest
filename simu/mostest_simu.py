# Good pipeline:
# + is aware of local storage area (/scratch or TMP)
# + can be customezed to use permanent storage  for all temporary files - this helps debugging
# + saves the results to the permanent storage as soon as they became available
# - allow to specify timeout for potentially long system calls
# + produce a detailed log file, containing hostname, and key environmental variables
# + can enable or disable submodules
# + is flexible enough so that one doesn't need to have two slightly updated copies of the same script
# + can continue where it left last time
# + in a loop across files, it should rely on it's knowledge what files to expect - not on "ls" and looking what files are there
# + it should produce files with deterministic names so that they can be treated as a token of success. This help caller to decide whether to re-run the pipeline or not.
# - consider copying input data to scratch area

# Maybe split into more independent units?
# Step 1. Generate phenotypes.

# generate phenotypes - separate
# mphen - separate
# mqfam - separate
# most, minp, abel - separate
# abel is tricky due to sumstats convertion - c
# qq plots - should be working now, corrrect?

# .sumstats.gz files never leave the temporary partition

import numpy as np
import pandas as pd
import os
import random
import scipy
import scipy.linalg
import argparse
import sys, traceback
from datetime import datetime
import getpass
import socket
import time
import six
import glob
import scipy.io as sio
import shutil
import subprocess

#mysystem = 'COMET'  # 'MMIL'
mysystem = 'MMIL'

# globals (not configurable)
if mysystem=='COMET':
    mostest_dir = '/home/oleksanf/github/mostest'
    indep_snps = '/oasis/projects/nsf/csd604/oleksanf/UKBDATA/projects/mostest_ukb/MultiABEL/indep.snps.RData'
    qq_py = '/home/oleksanf/github/python_convert/qq.py'
    bfile21 = '/oasis/projects/nsf/csd604/oleksanf/UKBDATA/projects/mostest_ukb/UKB26502_QCed_230519_maf0p005_chr21'
    bfile2122 = '/oasis/projects/nsf/csd604/oleksanf/UKBDATA/projects/mostest_ukb/UKB26502_QCed_230519_maf0p005_chr21_chr22'
    snps2122 = 203229; nsubj=26502; mostchunk=1000
    mvplink_exec = '/home/oleksanf/github/ofrei_workflows/mostest_simu/plink.multivariate'
    rlibpaths = '.libPaths("/home/oleksanf/R/site-library/");'
    simu_linux = 'simu_linux'
    load_matlab = 'module load matlab'
    load_R = 'module load R'
    load_plink = 'module load plink'
    real_pheno_default = '/oasis/projects/nsf/csd604/oleksanf/UKBDATA/projects/mostest_ukb/SubcorticalVolume_aseg35.csv'
    tmp_prefix_default = '/oasis/projects/nsf/csd604/oleksanf/UKBDATA/projects/mostest_ukb/simu/run17/simu'
    out_prefix_default = '/oasis/projects/nsf/csd604/oleksanf/UKBDATA/projects/mostest_ukb/simu/run17/simu'
    python_cmd = 'python'
    R_cmd = 'R'
    largest_nc = 1000   # largest number of causal variants (used as a bach size in MultiPhen)

if mysystem == 'MMIL':
    mostest_dir = '/home/oleksandr/precimed/mostest'
    indep_snps = '/space/gwas-syn1/1/data/GWAS/UKBioBank/projects/mostest_ukb/MultiABEL/indep.snps.RData'
    qq_py = '/home/oleksandr/precimed/python_convert/qq.py'
    #bfile21 = '/space/gwas-syn1/1/data/GWAS/UKBioBank/projects/mostest_ukb/UKB26502_QCed_230519_maf0p005_chr21'
    #bfile2122 = '/space/gwas-syn1/1/data/GWAS/UKBioBank/projects/mostest_ukb/UKB26502_QCed_230519_maf0p005'; nsubj=26502; mostchunk=2000;   snps2122 = 7428630;   # 36 times slower 
    #bfile2122 = '/space/gwas-syn1/1/data/GWAS/UKBioBank/projects/mostest_ukb/UKB26502_QCed_230519_maf0p005_chr21_chr22';nsubj=26502; mostchunk=2000; snps2122 = 203229;
    #bfile21= '/space/gwas-syn1/1/data/GWAS/UKBioBank/projects/mostest_ukb/UKB26502_QCed_230519_maf0p005_chr21_debug'; 
    #bfile2122 = '/space/gwas-syn1/1/data/GWAS/UKBioBank/projects/mostest_ukb/UKB26502_QCed_230519_maf0p005_chr21_chr22_debug'; snps2122 =1000;nsubj=100; mostchunk=2000;
    mostchunk = 2000
    mvplink_exec = '/home/oleksandr/bin/plink.multivariate'
    rlibpaths = ''  # that's fine to keep empty- all packages installed via anaconda
    simu_linux = '/home/oleksandr/bin/simu_linux'
    load_matlab = 'true'
    load_R = 'true'
    load_plink = 'true'
    real_pheno_default = '/space/gwas-syn1/1/data/GWAS/UKBioBank/projects/mostest_ukb/SubcorticalVolume_aseg35.csv'
    out_prefix_default = '/space/gwas-syn1/1/data/GWAS/UKBioBank/projects/mostest_ukb/simu/run23/simu'
    #out_prefix_default = os.path.join(os.environ['TMP'], 'simu')
    tmp_prefix_default = os.path.join('/scratch/', 'simu')
    python_cmd = '~/miniconda3/bin/python3'
    R_cmd = '~/miniconda3/bin/R'
    largest_nc = 1000   # largest number of causal variants (used as a bach size in MultiPhen)

# Software requirements for this simulation pipeline:
# - parallel (gnu)
# - R with MultiPhen and MultiAbel
# - MATLAB with MOSTest
# - qq.py from python_convert
# Data inputs:
# - bfile

def parse_args(args):
    parser = argparse.ArgumentParser(description="A helper tool to run simulations with MOSTest.")
    parser.add_argument("--real-pheno", type=str, default=real_pheno_default)
    parser.add_argument("--out-prefix", type=str, default=out_prefix_default)
    parser.add_argument("--tmp-prefix", type=str, default=tmp_prefix_default)
    parser.add_argument("--seed", type=int, default=None, help="Random seed")
    parser.add_argument("--T", type=int, default=None, help="number of traits (total)")
    parser.add_argument("--t", type=int, default=None, help="number of traits with an effect")
    parser.add_argument("--nc", type=int, default=None, help="number of causal SNPs per trait")
    parser.add_argument("--dist", type=str, default='norm', choices=['norm', 'cauchy', 'sparse'], help="distribution of genetic effects")
    parser.add_argument("--rg", type=str, default='eye', choices=['eye', 'real'], help="correlation of genetic effects")
    parser.add_argument("--re", type=str, default='eye', choices=['eye', 'real'], help="correlation of environmental effects")
    parser.add_argument("--h2", type=str, default='4e-3', help="heritability")
    parser.add_argument("--comb", type=str, default='none', choices=['none', 'sum', 'prod'], help="combine features (none = keep the original; sum = take all pairwise sums, prod = take all pairwise products)")
    parser.add_argument("--link", type=str, default='id', choices=['id', 'exp'], help="link function (identity or exponent)")
    parser.add_argument("--rep", type=str, default=['0'], nargs='+', help="simulation repeat")
    parser.add_argument("--eig", type=str, default=[], nargs='*', help="a list of num_eigval_to_regularize parameters")
    parser.add_argument("--int", default=False, action="store_true", help='apply rank-based inverse normal transform (as implemented in MOSTest)')
    parser.add_argument("--analysis", type=str, default=['pheno', 'most', 'abel', 'qq', 'mqfam', 'mphen'], nargs='*', choices=['pheno', 'most', 'abel', 'mqfam', 'mphen', 'qq'])
    parser.add_argument("--bfile21", type=str, default=None, help="plink bfile for all chromosomes chr21 (where we draw causal variants)")
    parser.add_argument("--bfile2122", type=str, default=None, help="plink bfile for all chromosomes or chr21+chr22")
    parser.add_argument("--snps2122", type=str, default=None, help="number of SNPs in all chromosomes or chr21+chr22")
    parser.add_argument("--nsubj", type=int, default=None, help="number of subjects in bim file (must be the same in bfile21 and bfile2122")
    
    return parser.parse_args(args)

def process_args(args):
    pass

'''
T - 25, 171 - num traits (total),
t - 1..T - num traits with an effect
nc - 10 | 100 | 1000 - num causals
dist - norm | cauchy - distribution of genetic effects
rg - eye | real - genetic correlation (eye / real)
re - eye | real - environmental noise correlation (eye / real)
h2 - 4e-4 | 4e-3 | 4e-2 - heritability 
comb - none | sum | prod - combine features (none = keep the original; sum = take all pairwise sums, prod = take all pairwise products)
link - id | exp - link function (identity or exponent)
rep - 0..9 - repeat


simu_T=10_t=5_nc=100_dist=norm_rg=eye_re=eye_h2=4e-3_comb=none_link=id_rep=9
'''

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def sec_to_str(t):
    '''Convert seconds to days:hours:minutes:seconds'''
    [d, h, m, s, n] = six.moves.reduce(lambda ll, b : divmod(ll[0], b) + ll[1:], [(t, 1), 60, 60, 24])
    f = ''
    if d > 0:
        f += '{D}d:'.format(D=d)
    if h > 0:
        f += '{H}h:'.format(H=h)
    if m > 0:
        f += '{M}m:'.format(M=m)

    f += '{S}s'.format(S=s)
    return f

class Logger(object):
    '''
    Lightweight logging.
    '''
    def __init__(self, fh, mode):
        self.fh = fh
        self.log_fh = open(fh, mode) if (fh is not None) else None

        # remove error file from previous run if it exists
        try:
            os.remove(fh + '.error')
        except OSError:
            pass

    def system(self, command, timeout=None):  # timeout in seconds, None means no timeout.
        start_time = time.time()
        log.log('>command start time: {T}'.format(T=time.ctime()) )
        self.log(command)        

        os.system(command)
 
        time_elapsed = round(time.time()-start_time,2)
        log.log('=command elapsed time: {T}'.format(T=sec_to_str(time_elapsed)))
        log.log('<command end time: {T}'.format(T=time.ctime()) )        

    def log(self, msg):
        '''
        Print to log file and stdout with a single command.
        '''
        eprint(msg)
        if self.log_fh:
            self.log_fh.write(str(msg).rstrip() + '\n')
            self.log_fh.flush()

    def error(self, msg):
        '''
        Print to log file, error file and stdout with a single command.
        '''
        eprint(msg)
        if self.log_fh:
            self.log_fh.write(str(msg).rstrip() + '\n')
            with open(self.fh + '.error', 'w') as error_fh:
                error_fh.write(str(msg).rstrip() + '\n')

def apply_heritability(pheno, noise, hsq):
    pheno_var = np.sum(np.square(pheno)) / len(pheno)
    noise_var = np.sum(np.square(noise)) / len(noise)
    gen_coef = 0 if (pheno_var < 1e-12) else np.sqrt(hsq / pheno_var)
    env_coef = np.sqrt(1.0 / noise_var) if (pheno_var < 1e-12) else np.sqrt((1.0 - hsq) / noise_var)
    return gen_coef * pheno + env_coef * noise

def generate_pheno(log, T, t, nc, dist, rg, re, h2, comb, link, rep, args, basename_tmp, basename_out):
    if (t != T) and (rg=='real'): raise(ValueError('rg==real require t==T'))
    if (t<=0) or (T<=0) or (T<t): raise(ValueError('invalid option: T={}, t={}'.format(T, t)))
    if dist not in ['norm', 'cauchy', 'sparse']: raise(ValueError('invalid option: dist=={}'.format(dist)))
    if rg not in ['eye', 'real']: raise(ValueError('invalid option: rg=={}'.format(rg)))
    if re not in ['eye', 'real']: raise(ValueError('invalid option: re=={}'.format(re)))
    if comb not in ['none', 'sum', 'prod']: raise(ValueError('invalid option: comb=={}'.format(comb)))
    if link not in ['id', 'exp']: raise(ValueError('invalid option: link=={}'.format(link)))
    if (comb != 'none') and (T > 40): raise(ValueError('invalid option: comb=={}, T=={} (too high number of output phenotypes)'.format(comb, T)))

    pheno_out = basename_out + '_pheno.csv'

    bim = pd.read_csv(args.bfile21+'.bim', sep='\t', header=None, names='CHR SNP GP BP A1 A2'.split())
    fam = pd.read_csv(args.bfile21+'.fam', delim_whitespace=True, header=None, names='FID IID C1 C2 C3 C4'.split())  #1000341 1000341 0 0 1 -9
    log.log('read bim file {}'.format(args.bfile21+'.bim'))
    log.log('read fam file {}'.format(args.bfile21+'.fam'))

    log.log(fam.head())

    if (rg=='real') or (re=='real'):
        real_corr = np.corrcoef(np.transpose(pd.read_csv(args.real_pheno, sep='\t').values))
        real_corr = real_corr[:T, :T]
        real_corr_sqrt = scipy.linalg.sqrtm(real_corr)
        real_dim = real_corr.shape[0]
        log.log('read real correlation matrix from {}, size {}'.format(args.real_pheno, real_dim))
        if (real_dim != T): raise(ValueError('the number of real phenotypes does not match --T'))

    envr_corr = real_corr if (re=='real') else np.identity(T)

    Teff = T if (dist != 'sparse') else t

    NC = int(np.ceil(nc*Teff/t))   
    causal_snps = list(bim.SNP.sample(NC).values)
    if dist in ['norm', 'sparse']: beta = np.random.normal(size=(nc, Teff))
    elif dist=='cauchy': beta = np.random.standard_cauchy(size=(nc, Teff))
    
    df=pd.DataFrame(index=causal_snps)
    for trait_index in range(0, T):
        causal_snps_trait = random.sample(causal_snps, nc)
        trait_name = 'trait{}'.format(trait_index+1)
        df[trait_name] = 0
        if trait_index >= Teff: continue
        df.loc[causal_snps_trait, trait_name] = beta[:, trait_index]

    if rg=='real': df.values[:] = np.transpose(np.dot(real_corr_sqrt, np.transpose(df.values)))

    df['SNP']=df.index

    df[['SNP'] + ['trait{}'.format(trait_index+1) for trait_index in range(0, T)]].to_csv(basename_out + '_beta.csv', header=None, sep='\t',index=False)
    for trait_index in range(0, T):
        df[['SNP', 'trait{}'.format(trait_index+1)]].to_csv(basename_tmp + '_tmp_trait{}_beta.csv'.format(trait_index+1), header=None, sep='\t',index=False)

    log.system('for PHENO in {1..'+str(T)+'}; do echo "' + simu_linux + ' --bfile '+args.bfile21+' --causal-variants ' + basename_tmp + '_tmp_trait${PHENO}_beta.csv --out ' + basename_tmp + '_tmp_trait${PHENO} --qt --hsq 1"; done | /home/oleksandr/bin/parallel -k --lb -j8')

    df=pd.read_csv(basename_tmp + '_tmp_trait{}.pheno'.format('1'), sep='\t'); del df['FID']; del df['IID']
    for i in range(2, T+1):
        df2=pd.read_csv(basename_tmp + '_tmp_trait{}.pheno'.format(i),sep='\t')
        df['trait{}'.format(i)] = df2['trait1']

    noise = np.random.multivariate_normal(np.zeros((T, 1)).flatten(), envr_corr, len(df))

    for trait_index in range(T):
        df.values[:, trait_index] = apply_heritability(df.values[:, trait_index], noise[:, trait_index], float(h2))

    if (comb=='sum') or (comb=='prod'):
        func = np.add if comb=='sum' else np.multiply
        for i in range(T):
            for j in range(i+1, T):
                df['trait_{}x{}'.format(i+1, j+1)] = func(df['trait{}'.format(i+1)].values, df['trait{}'.format(j+1)].values)
        for i in range(T):
            del df['trait{}'.format(i+1)]

    if link=='exp': df.values[:] = np.exp(df.values)

    df.columns = ['Var{}'.format(i+1) for i,c in enumerate(df.columns)]
    df.to_csv(pheno_out, sep='\t', index=False)
    df.insert(0, "IID", fam['IID'].values)
    df.to_csv(pheno_out + '.iid', sep='\t', index=False)
    df.insert(1, "FID", fam['FID'].values)
    df.to_csv(pheno_out + '.iid.fid', sep='\t', index=False)
    log.system('rm {}_tmp_*'.format(basename_tmp))

def aggregate_results(basename_tmp, basename_out, sumstats_file, aparam, bparam, C0_cond, test, mode, test2):
    ts=[0.05, 0.01, 1e-4, 1e-6, 5e-7, 5e-8]

    beta = pd.read_csv(basename_out + '_beta.csv',sep='\t',header=None,usecols=[0]).rename(columns={0:'SNP'})
    if test in ['most', 'minp', 'abel']: sumstats = pd.read_csv(basename_tmp + sumstats_file, delim_whitespace=True)
    if test in ['mqfam']: sumstats = pd.read_csv(basename_tmp + sumstats_file, delim_whitespace=True).rename(columns={'P':'PVAL'}) 
    if test in ['mphen']: sumstats = pd.read_csv(basename_tmp + sumstats_file, delim_whitespace=True, header=None,usecols=[1]).rename(columns={1:'PVAL'}) 
    if test=='mphen': sumstats['SNP'] = beta
    
    df_data={}
    for chri in ['chr21', 'chr22']:
        if (test in ['mqfam', 'mphen']) and ((chri!='chr21') or (mode != 'orig')): continue

        if (test in ['mqfam', 'mphen']): pvalues = sumstats['PVAL']
        elif (mode == 'orig') and (chri == 'chr21'): pvalues = pd.merge(sumstats, beta)['PVAL']
        elif (mode == 'orig') and (chri != 'chr21'): pvalues = sumstats[sumstats.CHR!=21]['PVAL'] 
        elif (mode != 'orig') and (chri == 'chr21'): pvalues = sumstats[sumstats.CHR==21]['PVAL']
        elif (mode != 'orig') and (chri != 'chr21'): pvalues = sumstats[sumstats.CHR!=21]['PVAL'] 

        for thresh in ts:
            for key,value in [x.split('=') for x in basename_out.split('/')[-1].split('_')[1:]]:
                insert_key_to_dictionary_as_list(key, value, df_data)

            insert_key_to_dictionary_as_list('test', test, df_data)
            insert_key_to_dictionary_as_list('test2', test2, df_data)
            insert_key_to_dictionary_as_list('matfile', "", df_data)
            insert_key_to_dictionary_as_list('basename', basename_out, df_data)
            insert_key_to_dictionary_as_list('valA', aparam, df_data)
            insert_key_to_dictionary_as_list('valB', bparam, df_data)
            insert_key_to_dictionary_as_list('valRCOND', C0_cond, df_data)
            insert_key_to_dictionary_as_list('mode', mode, df_data)
            insert_key_to_dictionary_as_list('chr', chri, df_data)
            insert_key_to_dictionary_as_list('thresh', thresh, df_data)
            insert_key_to_dictionary_as_list('num_sig', np.sum(pvalues<=thresh), df_data)
            insert_key_to_dictionary_as_list('num_val', len(pvalues), df_data)
            insert_key_to_dictionary_as_list('frac_sig', np.sum(pvalues<=thresh)/len(pvalues), df_data)

    pd.DataFrame(df_data).to_csv(basename_out + sumstats_file + '.table.csv', sep='\t', index=False)

def mostest(args, log, basename_tmp, basename_out):

    pheno_tmp = '{}_pheno.csv'.format(basename_tmp)
    pheno_out = '{}_pheno.csv'.format(basename_out)
    zmat_out = '{}_zmat.mat'.format(pheno_out)
    zmat_tmp = '{}_zmat.mat'.format(pheno_tmp)

    if os.path.exists(zmat_out):
        log.log("{} exists, skip MOSTest analysis".format(zmat_out))
        return
    if os.path.exists(zmat_tmp):
        log.log("{} exists, skip MOSTest analysis".format(zmat_tmp))
        return

    matlab_command = "pheno='{pheno_out}'; out = '{pheno_tmp}'; bfile='{bfile}';chunk={chunk};snps={snps2122};nsubj={nsubj}; apply_int={applyint}; mostest; exit;".format(pheno_out=pheno_out, pheno_tmp=pheno_tmp, snps2122=args.snps2122, bfile=args.bfile2122, nsubj=args.nsubj, chunk=mostchunk, applyint=('true' if args.int else 'false'))
    command = '{load_matlab} && cd {mostest} && matlab -nodisplay -nosplash -r "'.format(mostest=mostest_dir, load_matlab=load_matlab) + matlab_command + '"'
    log.system(command)
    log.log('Done mostest')

def mostest_rerun(args, log, num_eigval_to_regularize, use_pheno_corr, basename_tmp, basename_out):
    label='and'.join(['PCOR' if use_pheno_corr else 'GCOR', 'EIG{}'.format(num_eigval_to_regularize)])

    pheno_out = '{}_pheno.csv'.format(basename_out)
    pheno_tmp = '{}_pheno.csv'.format(basename_tmp)
    zmat_tmp = '{}_zmat.mat'.format(pheno_tmp)
    zmat_out = '{}_zmat.mat'.format(pheno_out)

    if not os.path.exists(zmat_tmp): log.system("cp {} {}".format(zmat_out, zmat_tmp))

    matlab_command = "zmat_name='{pheno_tmp}_zmat.mat'; out = '{pheno_tmp}_{label}'; num_eigval_to_regularize={eig}; use_pheno_corr={pcr}; mostest; exit;".format(pheno_tmp=pheno_tmp, label=label, eig=num_eigval_to_regularize, pcr=('true' if use_pheno_corr else 'false'))
    command = '{load_matlab} && cd {mostest} && matlab -nodisplay -nosplash -r "'.format(mostest=mostest_dir, load_matlab=load_matlab) + matlab_command + '"'
    log.system(command)
    log.system('cd {mostest} && {python_cmd} process_results.py {bfile}.bim {pheno_tmp}_{label}'.format(mostest=mostest_dir, bfile=args.bfile2122, pheno_tmp=pheno_tmp, label=label, python_cmd=python_cmd))
    log.system("gzip  -f {pheno_tmp}_{label}*sumstats".format(pheno_tmp=pheno_tmp, label=label))
    log.log('Done mostest re_run')

    mat=sio.loadmat('{pheno_tmp}_{label}.mat'.format(pheno_tmp=pheno_tmp, label=label))
    gamma_a,gamma_b=list(mat['gamma_params'][0])
    beta_a,beta_b=list(mat['beta_params'][0])
    C0_cond = np.linalg.cond(mat['C0'])

    for test in 'minp most'.split():
        for mode in 'orig perm'.split():
            sumstats_file = '_pheno.csv_{}.{}_{}.sumstats.gz'.format(label, test, mode)
            aparam, bparam = (gamma_a, gamma_b) if (test=='most') else (beta_a, beta_b)
            aggregate_results(basename_tmp, basename_out, sumstats_file, aparam, bparam, C0_cond, test, mode, label)
      

def multiabel(args, log, perm_or_orig, basename_tmp, basename_out):
    pheno_tmp = '{}_pheno.csv'.format(basename_tmp)
    pheno_out = '{}_pheno.csv'.format(basename_out)
    zmat_tmp = '{}_zmat.mat'.format(pheno_tmp)
    zmat_out = '{}_zmat.mat'.format(pheno_out)
    if not os.path.exists(zmat_tmp): log.system("cp {} {}".format(zmat_out, zmat_tmp))

    # ensure we have .mat file that go together with zmat_file
    if not os.path.exists('{}.mat'.format(pheno_tmp)):
        log.log("re-generate {}.mat".format(pheno_tmp))
        matlab_command = "zmat_name='{pheno_tmp}_zmat.mat'; out = '{pheno_tmp}'; mostest; exit;".format(pheno_tmp=pheno_tmp)
        log.system('{load_matlab} && cd {mostest} && matlab -nodisplay -nosplash -r "'.format(mostest=mostest_dir, load_matlab=load_matlab) + matlab_command + '"')

    varpheno=pheno_tmp  # ok as long as we use /scratch
    log.system('cd {mostest} && {python_cmd} process_results_ext.py {bfile}.bim {pheno_tmp} {varpheno} {what}'.format(mostest=mostest_dir, 
               bfile=args.bfile2122, varpheno=varpheno, pheno_tmp=pheno_tmp, python_cmd=python_cmd, what=perm_or_orig))

    # 1894/64969 =  2.9% , 203229/7428630 = 2.7%
    # indep.snps.RData contain 64969, out of which 1894 snps match our 21_22 reference . This is 2.9% of their SNPs
    # our chr21_chr22 reference is 2.7% of the total set of autosomes - so this is a good match, indicating that MultiABEL SNPs match very well.
    r_command = '''{rlibpaths}
require(MultiABEL);
load("{indep}");
pattern="{varpheno}.Var*.{what}.sumstats.gz";
data <- load.summary(Sys.glob(pattern), columnNames = c("SNP", "A1", "A2", "FRQ", "BETA", "SE", "N"), indep.snps = indep.snps);
result <- MultiSummary(data);
write.table(result$scan, "{pheno_tmp}.MultiABEL.{what}.sumstats", sep="\\t", row.names=FALSE, quote=FALSE);
'''.format(rlibpaths=rlibpaths, indep=indep_snps, pheno_tmp=pheno_tmp, what=perm_or_orig, varpheno=varpheno)
    command = "{load_R} && {R_cmd} -e '{}'".format(r_command, load_R=load_R, R_cmd=R_cmd).replace('\n', ' ')
    log.system(command)
    try:
        df=pd.read_csv("{pheno_tmp}.MultiABEL.{what}.sumstats".format(pheno_tmp=pheno_tmp,what=perm_or_orig), delim_whitespace=True)
    except:
        log.log('multiabel failed to produce .sumstats file')
        log.system("rm {varpheno}*Var*".format(varpheno=varpheno))
        return
    df.rename(columns={'marker':'SNP',  'freq':'FRQ',     'n':'N',       'p':'PVAL'}, inplace=True)
    bim2122 = pd.read_csv(args.bfile2122+'.bim', sep='\t', header=None, names='CHR SNP GP BP A1 A2'.split())
    df=pd.merge(bim2122[['SNP', 'CHR', 'BP', 'A1', 'A2']], df, on='SNP', how='left')
    df.to_csv("{pheno_tmp}.MultiABEL.{what}.sumstats".format(pheno_tmp=pheno_tmp, what=perm_or_orig), sep='\t', index=False)
    log.system("gzip  -f {pheno_tmp}.MultiABEL.{what}.sumstats".format(pheno_tmp=pheno_tmp, what=perm_or_orig))

    sumstats_file = '_pheno.csv.MultiABEL.{what}.sumstats.gz'.format(what=perm_or_orig)
    aparam, bparam, C0_cond, test = 0, 0, 0, 'abel',
    aggregate_results(basename_tmp, basename_out, sumstats_file, aparam, bparam, C0_cond, test, perm_or_orig, test)

    log.system("rm {varpheno}*Var*".format(varpheno=varpheno))
    log.log('Done multiabel')
    
def mvplink(args, log, basename_tmp, basename_out):
    command = """
cut {basename}_beta.csv -f1 > {basename}_causals.csv
{load_plink} && plink --bfile {bfile} --extract {basename}_causals.csv --out {basename}_causals_bfile --make-bed
{mvplink} --bfile {basename}_causals_bfile --mult-pheno {basename}_pheno.csv.iid.fid  --noweb --mqfam --out {basename}_pheno.csv_mvplink_orig
""".format(basename=basename_out, bfile=args.bfile21, mvplink=mvplink_exec, load_plink=load_plink)
    log.system(command)

    aparam, bparam, C0_cond, test = 0, 0, 0, 'mqfam'
    sumstats_file = '_pheno.csv_mvplink_orig.mqfam.total'
    aggregate_results(basename_out, basename_out, sumstats_file, aparam, bparam, C0_cond, test, 'orig', test)

    log.log('Done mvplink')

def multiphen(args, log, basename_tmp, basename_out):
    r_command =  """{rlibpaths}
library(MultiPhen);
pheno.opts = mPhen.options("pheno.input");
pheno = mPhen.readPhenoFiles("{basename}_pheno.csv.iid", opts = pheno.opts);
opts = mPhen.options("regression");
phenoObject = mPhen.preparePheno(pheno,opts = opts);
numPhenos = length(phenoObject$phenN);
geno.opts = mPhen.options("geno.input");
geno.opts$mPhen.batch = {largest_nc};
geno.opts$mPhen.format = "GT";
genoConnection <-mPhen.readGenotypes("{basename}_causals_bfile", opts = geno.opts, indiv = rownames(pheno$pheno));
geno <-genoConnection$genoData;
resultsJoint = mPhen.assoc(geno, phenoObject, opts = opts);
write.table(resultsJoint$Res[,,numPhenos+1,2], "{basename}_pheno.csv_mPhen_orig", sep="\\t", row.names=TRUE, quote=FALSE, col.names=FALSE);
""".format(rlibpaths=rlibpaths, basename=basename_out, largest_nc=largest_nc)
    command = "{load_R} && {R_cmd} -e '{rcmd}'".format(load_R=load_R, R_cmd=R_cmd, basename=basename_out, bfile=args.bfile21, rcmd=r_command).replace('\n', '')
    log.system(command)
    sumstats_file = '_pheno.csv_mPhen_orig'
    aparam, bparam, C0_cond, test = 0, 0, 0, 'mphen'
    aggregate_results(basename_out, basename_out, sumstats_file, aparam, bparam, C0_cond, test, 'orig', test)
    log.log('Done multiphen')

def qqplots(log, basename_tmp, basename_out):
    pheno_tmp = '{}_pheno.csv'.format(basename_tmp)
    log.system('ls {pheno_tmp}*sumstats.gz | /home/oleksandr/bin/parallel -j1 echo "{python_cmd} {qq} {{}} --top-as-dot 100 --strata {{}} --strata-cat CHR --strata-cat-ids chr21=21,chrNOT21=1:2:3:4:5:6:7:8:9:10:11:12:13:14:15:16:17:18:19:20:22 --out {{}}.png" | bash'.format(qq=qq_py, pheno_tmp=pheno_tmp, python_cmd=python_cmd))
    log.system('cp {}*json {}'.format(pheno_tmp, os.path.dirname(basename_out)))

def run_analysis(args, log, basename_tmp, basename_out):
    pheno_out = '{}_pheno.csv'.format(basename_out)
    pheno_tmp = '{}_pheno.csv'.format(basename_tmp)
    zmat_out = '{}_zmat.mat'.format(pheno_out)

    if os.path.exists(pheno_out):
        log.log('skip making {}, file exists'.format(pheno_out))
        args.analysis = [x for x in args.analysis if x != 'pheno']

    #if os.path.exists(zmat_out):
    #    log.log('skip making {}, file exists'.format(zmat_out))
    #    args.analysis = [x for x in args.analysis if x != 'most']

    if 'pheno' in args.analysis:
        # saves output to out folder
        generate_pheno(log=log, T=args.T, t=args.t, nc=args.nc, dist=args.dist, rg=args.rg, re=args.re, h2=args.h2, comb=args.comb, link=args.link, rep=repeat, args=args, basename_tmp=basename_tmp, basename_out=basename_out)

    if 'most' in args.analysis:
        mostest(args, log, basename_tmp, basename_out)
        log.system("cp {}_zmat.mat {}_zmat.mat".format(pheno_tmp, pheno_out))

        for eig in args.eig:
            for use_pheno_corr in [False, True]:
                # generates .table.csv and saves results directly to basename_out
                # keeps .mat and .sumstats.gz files in /tmp so we could also make QQ plots
                mostest_rerun(args, log, num_eigval_to_regularize=eig, use_pheno_corr=use_pheno_corr, basename_tmp=basename_tmp, basename_out=basename_out)

    if 'mqfam' in args.analysis:
        mvplink(args, log, basename_tmp, basename_out)

    if 'mphen' in args.analysis:
        multiphen(args, log, basename_tmp, basename_out)

    if ('abel' in args.analysis) and (args.T >= 2):
        multiabel(args, log, 'perm', basename_tmp, basename_out )
        multiabel(args, log, 'orig', basename_tmp, basename_out)

    if 'qq' in args.analysis:
        qqplots(log, basename_tmp, basename_out)

def insert_key_to_dictionary_as_list(key, value, df_data):
    if key not in df_data:
        df_data[key] = []
    df_data[key].append(value)

if __name__ == "__main__":
    args = parse_args(sys.argv[1:])
    process_args(args)
    if args.seed != None: np.random.seed(args.seed)
   
    for repeat in args.rep:
        try:
            basename_tmp = '{tmp_prefix}_T={T}_t={t}_nc={nc}_dist={dist}_rg={rg}_re={re}_h2={h2}_comb={comb}_link={link}_int={int}_rep={repeat}'.format(repeat=repeat, **vars(args))
            basename_out = '{out_prefix}_T={T}_t={t}_nc={nc}_dist={dist}_rg={rg}_re={re}_h2={h2}_comb={comb}_link={link}_int={int}_rep={repeat}'.format(repeat=repeat, **vars(args))

            log = Logger(basename_out + '.log', 'a')
            start_time = time.time()

            defaults = vars(parse_args([]))
            opts = vars(args)
            non_defaults = [x for x in opts.keys() if opts[x] != defaults[x]]
            header = "Call: \n"
            header += './mostest_simu.py \\\n'
            options = ['\t--'+x.replace('_','-')+' '+str(opts[x]).replace('\t', '\\t')+' \\' for x in non_defaults]
            header += '\n'.join(options).replace('True','').replace('False','')
            header = header[0:-1]+'\n'
            log.log(header)
            log.log('Beginning analysis at {T} by {U}, host {H}'.format(T=time.ctime(), U=getpass.getuser(), H=socket.gethostname()))

            for envvar in 'JOB_ID HOSTNAME QUEUE JOB_NAME JOB_SCRIPT TMP TMPDIR'.split():
                if envvar not in os.environ: continue
                log.log('os.environ["{}"] = {}'.format(envvar, os.environ[envvar]))

            run_analysis(args, log, basename_tmp, basename_out)
            log.system('rm {}*'.format(basename_tmp))

        except Exception:
            ex_type, ex, tb = sys.exc_info()
            log.error( traceback.format_exc(ex) )
            raise

        finally:
            log.log('Analysis finished at {T}'.format(T=time.ctime()) )
            time_elapsed = round(time.time()-start_time,2)
            log.log('Total time elapsed: {T}'.format(T=sec_to_str(time_elapsed)))
