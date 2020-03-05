plots
- nc x h2            (obviously, the effect size drives power; also - PCOR MOSTest behaves as GCOR )
- dist x t           (most is always better than minp; minp power is growing when t increases; cauchy distribution => minp has thesame power as mostest)
- dist x T           (more traits = better power; and too many traits = problem for MultiABEL)
- link_int x rg_re   (less obviously, the effect size is best in the presence of phenotypic correlation AND in the absence of genetic correlation)
                     (also, it is not obvious - but known from before - that INT reduces power)
- comb x link_int    (type-I error analysis)
                     (INT is absolutely essential in the presence of non-linear effects)
                     (regularization is important - sometimes it improves power; in terms of improving correctness of the test, this needs better visualization - currently unclear)

ideas
- show QQ plots on the same FacetGrids
- show how a&b&rcond depends on regularization? (can put labels on the QQ plots)
- draw causals on chr22, evaluate on all other chromosomes to have lower error bars. chr21 is enough to evaluate power, but for the null we may need a lot of observation.
- really nice to include MultiPheno
+ show that PCOR is the same as GCOR

findings (also described in the "plots" above)
- enabling INT may reduce the power (e.g comb=none & link=id) but it may improve it also (comb=none & link=exp)
- doing residualization does improve the power (e.g. comb=sum)





import pandas as pd
import glob
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

files = glob.glob('/home/oleksanf/vmshare/data/mostest_simu/run21/*table*')
print(len(files))
#files = ['/home/oleksanf/vmshare/data/mostest_simu/run21/simu_T=25_t=5_nc=100_dist=norm_rg=eye_re=eye_h2=0.004_comb=none_link=id_int=True_rep=4.table.csv']
df=pd.concat([pd.read_csv(file, sep='\t') for file in files])
df.loc[df['test']=='abel', 'test2']='abel'

lens= [len(pd.read_csv(file, sep='\t')) for file in files]
print(df.columns)

df_true = df.copy()
df_true['test'] = 'true'
df_true['test2'] = 'true'
df_true['valA'] = 0
df_true['valB'] = 0
df_true['valRCOND'] = 0
df_true['basename'] = ''
df_true['matfile'] = ''
df_true['num_sig'] = 1
df_true['num_val'] = np.divide(1, df_true['thresh'].values)
df_true['frac_sig'].values[:] = df_true['thresh'].values
print(len(df_true))
df_true.drop_duplicates(inplace=True)
print(len(df_true))

df = pd.concat([df, df_true])

df['rg_re'] = ['{}_{}'.format(x, y) for x, y in zip(df['rg'].values, df['re'].values)]
df['link_int'] = ['{}_{}'.format(x, y) for x, y in zip(df['link'].values, df['int'].values)]
df['test_test2'] = ['{}_{}'.format(x, y) for x, y in zip(df['test'].values, df['test2'].values)]
df['comb_test_test2'] = ['{}_{}_{}'.format(x, y, z) for x, y, z in zip(df['comb'], df['test'].values, df['test2'].values)]

df=df[df.thresh != 5e-7].copy()

#df.to_excel('/home/oleksanf/vmshare/data/mostest_simu/run21.xlsx', index=False)

if 0:
    df_plot = df.reset_index().copy()

    groups=['T', 't', 'nc', 'dist', 'rg', 're', 'h2', 'comb', 'link', 'int', 'test', 'test2', 'mode', 'thresh', 'chr']
    params=['num_sig', 'num_val', 'valA', 'valB', 'valRCOND']
    params_agg={'num_sig':['sum', 'count'], 'num_val':['sum'], 'valA':['mean', 'std'], 'valB':['mean', 'std'], 'valRCOND':['mean', 'std']}

    df_agg = df_plot[groups+params].groupby(groups).agg(params_agg)
    df_agg.columns = ['_'.join(col).strip() for col in df_agg.columns.values]
    df_agg = df_agg.reset_index()
    df_agg['frac_sig'] = np.divide(df_agg['num_sig_sum'].values, df_agg['num_val_sum'].values)
    del df_agg['num_sig_sum'] 
    del df_agg['num_val_sum'] 
    df_agg=df_agg.rename(columns={'num_sig_count':'num_rep'}).copy()  #.to_excel('/home/oleksanf/vmshare/data/mostest_simu/run21.xlsx',index=False)




def filter(df, keep=[], excl=[], scope=[], verbose=False):
    mask_keep = df.index.notnull()
    for key, val in keep:
        if not isinstance(val, list): val=[val]
        mask_keep = mask_keep & (df[key].isin(val))
        if verbose: print('keep', key, val, sum(mask_keep))

    mask_excl = ~df.index.notnull()
    for key, val in excl:
        if not isinstance(val, list): val=[val]
        mask_excl = mask_excl | (df[key].isin(val))
        if verbose: print('excl', key, val, sum(mask_excl))

    mask_scope = df.index.notnull()
    for key, val in scope:
        if not isinstance(val, list): val=[val]
        mask_scope = mask_scope & (df[key].isin(val))
        if verbose: print('scope', key, val, sum(mask_scope))

    return df[~mask_scope | (mask_keep & ~mask_excl)].copy()







for mode in ['orig', 'perm', 'null'] :
    mode_filt = [('chr', 'chr21' if (mode=='orig') else 'chr22') , ('mode', 'orig')] if (mode != 'perm') else [('mode', 'perm')]
    data = filter(df, keep=[('T', 25), ('t', 25), ('dist', 'norm'), ('rg', 'eye'), ('re', 'eye'), ('comb', 'none'), ('link', 'id'), ('int', True)] + mode_filt)
    data = filter(data, keep=[('test2', ['GCORandEIG0', 'PCORandEIG0'])], scope=[('test', ['most', 'minp'])])
    data = filter(data, excl=[('test2', ['PCORandEIG0'])], scope=[('test', ['minp'])])
    data['test'] = [('_'.join([a, b]) if (a != b) else a) for a, b in zip(data['test'], data['test2'])]
    #data = filter(data, keep=[('thresh', '0.01'), ('h2', '0.004'), ('nc', 100)])   # for debugging, to keep just one plot

    g = sns.catplot(x="thresh", y="frac_sig", hue="test", row="nc", col="h2", data=data, height=5, aspect=1.5, kind='bar', ci='sd',
                    row_order=[1000, 100, 10], col_order = [0.0004, 0.004, 0.04, 0.4],
                    hue_order=[ 'true', 'most_GCORandEIG0', 'minp_GCORandEIG0', 'abel', 'most_PCORandEIG0'])  # kind - “point”, “bar”, “strip”, “swarm”, “box”, “violin”, or “boxen”
    if mode!='orig': g.set(yscale="log")
    g.savefig("/home/oleksanf/vmshare/data/mostest_simu/nc_h2_{}.png".format(mode))

    print('max per group', data.groupby(['thresh', 'test', 'nc', 'h2']).agg({'frac_sig':'count'}).reset_index()['frac_sig'].max())







for mode in ['orig', 'perm', 'null'] :
    mode_filt = [('chr', 'chr21' if (mode=='orig') else 'chr22') , ('mode', 'orig')] if (mode != 'perm') else [('mode', 'perm')]
    data = filter(df, keep=[('T', 25), ('t', 25), ('nc', 100),  ('h2', 0.004), ('dist', 'norm'),  ('rg', 'eye'), ('re', 'eye')] + mode_filt)
    data = filter(data, excl=[('test2', ['_pheno.csv', 'PCORandEIG275', 'PCORandEIG100', 'PCORandEIG20', 'PCORandEIG10', 'PCORandEIG0'])])
    data = filter(data, excl=[('test2', ['GCORandEIG100', 'GCORandEIG10', 'GCORandEIG20'])], scope=[('test', ['most', 'minp'])])
    data['test'] = [('_'.join([a, b]) if (a != b) else a) for a, b in zip(data['test'], data['test2'])]
    #data = filter(data, keep=[('thresh', '0.01'), ('h2', '0.004'), ('nc', 100)])   # for debugging, to keep just one plot

    g = sns.catplot(x="thresh", y="frac_sig", hue="test", row="comb", col="link_int", data=data, height=5, aspect=1.5, kind='bar', ci='sd',
                    row_order=['none', 'sum', 'prod'], col_order = ['id_False', 'exp_False', 'id_True', 'exp_True'],
                    hue_order=[ 'true', 'most_GCORandEIG0', 'minp_GCORandEIG0', 'abel', 'most_GCORandEIG275', 'minp_GCORandEIG275'])
    if mode!='orig': g.set(yscale="log")
    g.savefig("/home/oleksanf/vmshare/data/mostest_simu/comb_link_int_{}.png".format(mode))

    print(mode, 'max per group', data.groupby(['thresh', 'test', 'comb', 'link_int']).agg({'frac_sig':'count'}).reset_index()['frac_sig'].max())




for mode in ['orig', 'perm', 'null'] :
    mode_filt = [('chr', 'chr21' if (mode=='orig') else 'chr22') , ('mode', 'orig')] if (mode != 'perm') else [('mode', 'perm')]
    data = filter(df, keep=[('T', 25), ('t', 25), ('nc', 100),  ('h2', 0.004), ('dist', 'norm'), ('comb', 'none')] + mode_filt, verbose=False)
    data = filter(data, keep=[('test2', 'GCORandEIG0')], scope=[('test', ['most', 'minp'])])
    data = filter(data, excl=[('link', 'exp')], scope=[('int', False)])
    
    data['test'] = [('_'.join([a, b]) if (a != b) else a) for a, b in zip(data['test'], data['test2'])]
    #data = filter(data, keep=[('thresh', '0.01')])   # for debugging, to keep just one plot

    g = sns.catplot(x="thresh", y="frac_sig", hue="test", row="link_int", col="rg_re", data=data, height=5, aspect=1.5, kind='bar', ci='sd',
                    row_order=['id_False', 'id_True', 'exp_True'], col_order = ['eye_eye', 'eye_real', 'real_eye', 'real_real'],
                    hue_order=[ 'true', 'most_GCORandEIG0', 'minp_GCORandEIG0', 'abel'])
    if mode!='orig': g.set(yscale="log")
    g.savefig("/home/oleksanf/vmshare/data/mostest_simu/link_int_rg_re_{}.png".format(mode))

    print(mode, 'max per group', data.groupby(['thresh', 'test', 'link_int', 'rg_re']).agg({'frac_sig':'count'}).reset_index()['frac_sig'].max())




for mode in ['orig', 'perm', 'null'] :
    mode_filt = [('chr', 'chr21' if (mode=='orig') else 'chr22') , ('mode', 'orig')] if (mode != 'perm') else [('mode', 'perm')]
    data = filter(df, keep=[('T', 25), ('nc', 100),  ('h2', 0.004),  ('rg', 'eye'), ('re', 'eye'), ('link', 'id'),  ('int', True), ('comb', 'none')] + mode_filt, verbose=False)
    data = filter(data, keep=[('test2', 'GCORandEIG0')], scope=[('test', ['most', 'minp'])])
        
    data['test'] = [('_'.join([a, b]) if (a != b) else a) for a, b in zip(data['test'], data['test2'])]
    #data = filter(data, keep=[('thresh', '0.01')])   # for debugging, to keep just one plot

    g = sns.catplot(x="thresh", y="frac_sig", hue="test", row="dist", col="t", data=data, height=4, aspect=1.3, kind='bar', ci='sd',
                    row_order=['norm', 'cauchy'], hue_order=[ 'true', 'most_GCORandEIG0', 'minp_GCORandEIG0', 'abel'])
    if mode!='orig': g.set(yscale="log")
    g.savefig("/home/oleksanf/vmshare/data/mostest_simu/dist_t_{}.png".format(mode))

    print(mode, 'max per group', data.groupby(['thresh', 'test', 'dist', 't']).agg({'frac_sig':'count'}).reset_index()['frac_sig'].max())







for mode in ['orig', 'perm', 'null'] :
    mode_filt = [('chr', 'chr21' if (mode=='orig') else 'chr22') , ('mode', 'orig')] if (mode != 'perm') else [('mode', 'perm')]
    data = filter(df, keep=[('nc', 100), ('h2', [0.004, 0.0004]), ('dist', 'norm') ,('rg', 'eye'), ('re', 'eye'), ('link', 'id'),  ('int', True), ('comb', 'none')] + mode_filt, verbose=False)
    data = filter(data, keep=[('test2', 'GCORandEIG0')], scope=[('test', ['most', 'minp'])])
        
    data['test'] = [('_'.join([a, b]) if (a != b) else a) for a, b in zip(data['test'], data['test2'])]
    #data = filter(data, keep=[('thresh', '0.01')])   # for debugging, to keep just one plot

    g = sns.catplot(x="thresh", y="frac_sig", hue="test", row="h2", col="T", data=data, height=4, aspect=1.3, kind='bar', ci='sd',
                    hue_order=[ 'true', 'most_GCORandEIG0', 'minp_GCORandEIG0', 'abel'])
    if mode!='orig': g.set(yscale="log")
    g.savefig("/home/oleksanf/vmshare/data/mostest_simu/dist_T_{}.png".format(mode))

    print(mode, 'max per group', data.groupby(['thresh', 'test', 'dist', 't']).agg({'frac_sig':'count'}).reset_index()['frac_sig'].max())

