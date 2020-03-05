import os.path
import sys
import glob

#design = 'debug'
#design = 'chr2122'
#design = 'chrALL'
design = 'real'

if design == 'debug':
    bfile21= '/space/gwas-syn1/1/data/GWAS/UKBioBank/projects/mostest_ukb/UKB26502_QCed_230519_maf0p005_chr21_debug'; 
    bfile2122 = '/space/gwas-syn1/1/data/GWAS/UKBioBank/projects/mostest_ukb/UKB26502_QCed_230519_maf0p005_chr21_chr22_debug'; snps2122 =1000;nsubj=100; mostchunk=2000;
    out_prefix = '/space/gwas-syn1/1/data/GWAS/UKBioBank/projects/mostest_ukb/simu/run31/simuV2'
    tmp_prefix = '/space/gwas-syn1/1/data/GWAS/UKBioBank/projects/mostest_ukb/simu/run31/simuV2_tmp'

if design == 'chr2122':
    bfile21 = '/space/gwas-syn1/1/data/GWAS/UKBioBank/projects/mostest_ukb/UKB26502_QCed_230519_maf0p005_chr21'
    bfile2122 = '/space/gwas-syn1/1/data/GWAS/UKBioBank/projects/mostest_ukb/UKB26502_QCed_230519_maf0p005_chr21_chr22';nsubj=26502; mostchunk=2000; snps2122 = 203229;
    out_prefix = '/space/gwas-syn1/1/data/GWAS/UKBioBank/projects/mostest_ukb/simu/run30/simuV2'
    tmp_prefix = '/scratch/simuV2'

if design == 'chrALL':
    bfile21 = '/space/gwas-syn1/1/data/GWAS/UKBioBank/projects/mostest_ukb/UKB26502_QCed_230519_maf0p005_chr21'
    bfile2122 = '/space/gwas-syn1/1/data/GWAS/UKBioBank/projects/mostest_ukb/UKB26502_QCed_230519_maf0p005'; nsubj=26502; mostchunk=2000;   snps2122 = 7428630;   # 36 times slower 
    out_prefix = '/space/gwas-syn1/1/data/GWAS/UKBioBank/projects/mostest_ukb/simu/run32/simuV3'
    tmp_prefix = '/scratch/simuV3'

if design == 'real':
    bfile21 = '/space/gwas-syn1/1/data/GWAS/UKBioBank/projects/mostest_ukb/UKB26502_QCed_230519_maf0p005_chr21'
    bfile2122 = '/space/gwas-syn1/1/data/GWAS/UKBioBank/projects/mostest_ukb/UKB26502_QCed_230519_maf0p005'; nsubj=26502; mostchunk=2000;   snps2122 = 7428630;   # 36 times slower 
    out_prefix = '/space/gwas-syn1/1/data/GWAS/UKBioBank/projects/mostest_ukb/simu/run35/simuV5'
    tmp_prefix = '/scratch/simuV5'


nrep=10
dry_run=False
group_size = 1



num_submit = int(sys.argv[1])

mostest_simu = '/home/oleksandr/precimed/ofrei_workflows/mostest_simu/mostest_simu.py'
#real_pheno = '/space/gwas-syn1/1/data/GWAS/UKBioBank/projects/mostest_ukb/SubcorticalVolume_aseg35.csv'

basename_pattern = '{out_prefix}_T={T}_t={t}_nc={nc}_dist={dist}_rg={rg}_re={re}_h2={h2}_comb={comb}_link={link}_int={int}_rep={rep}'
out_pattern = '{out_prefix}_T={T}_t={t}_nc={nc}_dist={dist}_rg={rg}_re={re}_h2={h2}_comb={comb}_link={link}_int={int}_rep={rep}_pheno.csv_PCORandEIG0.most_orig.sumstats.gz.table.csv'
cmd_pattern = 'bash && ~/miniconda3/bin/python3 {mostest_simu} --T {T} --t {t} --nc {nc} --dist {dist} --rg {rg} --re {re} --h2 {h2} --comb {comb} --link {link} --rep {rep} --int {int} --eig {eig} --tmp-prefix {tmp_prefix} --out-prefix {out_prefix} --bfile21 {bfile21} --bfile2122 {bfile2122} --snps2122 {snps2122} --nsubj {nsubj} --real-pheno /space/gwas-syn1/1/data/GWAS/UKBioBank/projects/mostest_ukb/{rph}.csv'

# vals = [ [('T', 25), ('t', '25'), ('nc', 100), ('dist', 'norm'), ('rg', 'eye'), ('re', 'eye'), ('h2', '0.004'), ('comb', 'none'), ('link', 'id'), ('eig', '0'), ('int', 'True') ] ]

all_lists = []
real_lists = []

vals = [ [('nc', 100), ('T', 171), ('t', 17), ('h2', '0.004'), ('dist', 'norm'), ('rg', 'eye'), ('re', 'real'), ('comb', 'none'), ('link', 'id'), ('eig', '0'), ('int', 'True'), ('rph', 'all_aseg35') ] ]
vals = [ (v + [('rep', str(x))]) for v in vals for x in range(10) ]
real_lists.extend(vals)
vals = [ [('nc', 100), ('T', 68), ('t', 7), ('h2', '0.004'), ('dist', 'norm'), ('rg', 'eye'), ('re', 'real'), ('comb', 'none'), ('link', 'id'), ('eig', '0'), ('int', 'True'), ('rph', 'CorticalArea') ] ]
vals = [ (v + [('rep', str(x))]) for v in vals for x in range(10) ]
#real_lists.extend(vals)
vals = [ [('nc', 100), ('T', 68), ('t', 6), ('h2', '0.004'), ('dist', 'norm'), ('rg', 'eye'), ('re', 'real'), ('comb', 'none'), ('link', 'id'), ('eig', '0'), ('int', 'True'), ('rph', 'CorticalThickness') ] ]
vals = [ (v + [('rep', str(x))]) for v in vals for x in range(10) ]
#real_lists.extend(vals)
vals = [ [('nc', 100), ('T', 35), ('t', 3), ('h2', '0.004'), ('dist', 'norm'), ('rg', 'eye'), ('re', 'real'), ('comb', 'none'), ('link', 'id'), ('eig', '0'), ('int', 'True'), ('rph', 'SubcorticalVolume_aseg35') ] ]
vals = [ (v + [('rep', str(x))]) for v in vals for x in range(10) ]
#real_lists.extend(vals)

# hard-code nrep=10 for the example run
vals = [ [('nc', 100), ('T', 25), ('t', 25), ('h2', '0.004'), ('dist', 'norm'), ('rg', 'eye'), ('re', 'eye'), ('comb', 'none'), ('link', 'id'), ('eig', '0'), ('int', 'True'), ('rph', 'all_aseg35') ] ]
vals = [ (v + [('rep', str(x))]) for v in vals for x in range(10) ]
all_lists.extend(vals)

vals = [ [('nc', 100), ('dist', 'norm'), ('rg', 'eye'), ('re', 'eye'), ('comb', 'none'), ('link', 'id'), ('eig', '0'), ('int', 'True'), ('rph', 'all_aseg35') ] ]
vals = [ (v + [('t', str(x)), ('T', str(x))]) for v in vals for x in [ 2, 10, 25, 50, 100, 171 ] ]
vals = [ (v + [('h2', str(x))]) for v in vals for x in ['0.04', '0.004', '0.0004' ] ]
#vals = [ (v + [('rg', str(x))]) for v in vals for x in 'real eye'.split() ]
#vals = [ (v + [('re', str(x))]) for v in vals for x in 'real eye'.split() ]
vals = [ (v + [('rep', str(x))]) for v in vals for x in range(nrep) ]
all_lists.extend(vals)
#print('after adding T-t-h2-rg-re: {} '.format(len(all_lists)))

# play with t, dist
vals = [ [('T', 25), ('nc', 100), ('rg', 'eye'), ('re', 'eye'), ('h2', '0.04'), ('comb', 'none'), ('link', 'id'), ('eig', '0'), ('int', 'True'), ('rph', 'all_aseg35') ] ]
vals = [ (v + [('t', str(x))]) for v in vals for x in [ 1, 2, 3, 4, 5 ] ]
vals = [ (v + [('dist', str(x))]) for v in vals for x in ['norm', 'cauchy', 'sparse'] ]
vals = [ (v + [('rep', str(x))]) for v in vals for x in range(nrep) ]
all_lists.extend(vals)
#print('after adding t-dist: {} '.format(len(all_lists)))

vals = [ [('dist', 'sparse'), ('nc', 100), ('rg', 'eye'), ('re', 'eye'), ('h2', '0.04'), ('comb', 'none'), ('link', 'id'), ('eig', '0'), ('int', 'True'), ('rph', 'all_aseg35') ] ]
vals = [ (v + [('t', str(x))]) for v in vals for x in [ 1, 2, 4, 10 ] ]
vals = [ (v + [('T', str(x))]) for v in vals for x in [ 1, 2, 4, 10, 20, 50, 100 ] ]
vals = [ (v + [('rep', str(x))]) for v in vals for x in range(nrep) ]
for val in vals:
    if int(dict(val)['t']) <= int(dict(val)['T']):
        all_lists.append(val)   # add only if t<=T

# play with rg, re, comb, link
vals = [ [('T', 25), ('t', '25'), ('nc', 100), ('dist', 'norm'), ('h2', '0.004'), ('eig', '0'), ('rph', 'all_aseg35') ]]
vals = [ (v + [('rg', str(x))]) for v in vals for x in 'real eye'.split() ]
vals = [ (v + [('re', str(x))]) for v in vals for x in 'real eye'.split() ]
vals = [ (v + [('comb', str(x))]) for v in vals for x in 'none sum'.split() ]   # prod
vals = [ (v + [('link', str(x))]) for v in vals for x in 'exp id'.split() ]
vals = [ (v + [('int', str(x))]) for v in vals for x in 'True False'.split() ]
vals = [ (v + [('rep', str(x))]) for v in vals for x in range(nrep) ]
for val in vals:
    dval = dict(val)
    if ('real' in [dval['re'], dval['rg']]) and (dval['comb']!='none'): continue
    if dval['comb']=='none': all_lists.append(val+[('eig', '0')])
    else: all_lists.append(val+[('eig', '0 275')])
#print('after adding rg-re-comblink: {} '.format(len(all_lists)))


# play with nc and h2
vals = [ [('T', 25), ('t', '25'), ('dist', 'norm'), ('rg', 'eye'), ('re', 'eye'), ('comb', 'none'), ('link', 'id'), ('eig', '0'), ('int', 'True'), ('rph', 'all_aseg35') ]]
vals = [ (v + [('nc', str(x))]) for v in vals for x in [ 10, 100, 1000 ] ]
vals = [ (v + [('h2', str(x))]) for v in vals for x in [ '0.04', '0.004', '0.0004' ] ]
vals = [ (v + [('rep', str(x))]) for v in vals for x in range(nrep) ]
all_lists.extend(vals)
#print('after adding nc-h2: {} '.format(len(all_lists)))

eig_lists = []
vals = [ [('T', 25), ('eig', '0 1 2 3 4 5 10 15 20 21 22 23 24'), ('t', '25'), ('link', 'exp'), ('int', 'True'),  ('comb', 'none'), ('nc', 100), ('dist', 'norm'), ('h2', '0.004'), ('rph', 'all_aseg35') ]]
vals = [ (v + [('rg', str(x))]) for v in vals for x in 'real eye'.split() ]
vals = [ (v + [('re', str(x))]) for v in vals for x in 'real eye'.split() ]
vals = [ (v + [('rep', str(x))]) for v in vals for x in range(nrep) ]
eig_lists.extend(vals)

#all_lists = eig_lists
all_lists = real_lists

skipped = 0; duplicated = 0;
processed = set()  # exclude duplicates

final_list = []
table_list = []
for val in all_lists:
    #table_file = out_pattern.format(**locals(), **dict(val))
    basename = basename_pattern.format(**locals(), **dict(val))
    if basename in processed:
        duplicated+=1
        continue

    processed.add(basename)
    cmd = cmd_pattern.format(**locals(), **dict(val)).replace('--int True', '--int').replace('--int False', '')

    analyses = ['abel', 'mphen', 'mqfam']
    '''
    has_pheno = os.path.exists(basename + "_pheno.csv")
    has_abel = (len(glob.glob(basename + "_pheno.csv.MultiABEL.*.sumstats.gz.table.csv")) == 2)
    
    #has_most = os.path.exists(basename + "_pheno.csv_zmat.mat")
    has_most = False

    has_mphen = (dict(val)['comb'] != 'none') or os.path.exists(basename + "_pheno.csv_mPhen_orig.table.csv")
    has_mqfam = (dict(val)['comb'] != 'none') or os.path.exists(basename + "_pheno.csv_mvplink_orig.mqfam.total.table.csv")

    #eig_corrected = []
    #for eig_val in dict(val)['eig'].split():
    #   if len(glob.glob("_pheno.csv_*andEIG{}.*_*.sumstats.gz.table.csv".format(eig_val))) == 8: continue
    #   eig_corrected.append(eig_val)
    #cmd = cmd.replace('--eig {}'.format(dict(val)['eig']), '--eig {}'.format(' '.join(eig_corrected)))

    #if not has_pheno: analyses.append('pheno')   # safery precaution  never re-generate phenotpe (we've got them by now)
    if not has_most: analyses.append('most')
    if not has_abel: analyses.append('abel')
    if not has_mqfam: analyses.append('mqfam')
    if not has_mphen: analyses.append('mphen')
    
    if len(analyses) == 0:  continue
    '''

    cmd += ' --analysis qq ' + ' '.join(analyses)
    final_list.append(cmd)
    table_list.append(basename+'.touch')

template='''
# Execute job in the queue "std.q" unless you have special requirements.
##$ -q std.q
##$ -q all_24.q

# The SGE batch system uses the current directory as working directory.
# Both files (output.dat and error.dat) will be placed in the current
# directory. The batch system assumes to find the executable in this directory.
#$ -cwd

# Redirect output stream to this file.
##$ -o output.dat

# Redirect error stream to this file.
##$ -e error.dat

# Send status information to this email address.
##$ -M oleksandr.frei@gmail.com

# Send an e-mail when the job is done.
##$ -m e

#$ -l h_vmem=120G
#$ -l h_rt=36:00:00
##$ -pe dmp4 16

##$ -l h=!(mmil-compute-5-7.local|mmil-compute-5-14.local|mmil-compute-5-5.local|mmil-compute-5-2.local|mmil-compute-5-10.local|mmil-compute-8-1|mmil-compute-8-2|mmil-compute-8-3|mmil-compute-8-4|mmil-compute-7-0|mmil-compute-7-1)

## https://unix.stackexchange.com/questions/277981/gnu-parallel-immediately-display-job-stderr-stdout-one-at-a-time-by-jobs-order
bash
printf "{}" | /home/oleksandr/bin/parallel -j''' +str(group_size)+ ''' -k --lb  {{}}
'''


#print(len(final_list))
#print(final_list[0], table_list[0])
groups = []
for rep in range(10): #[0,2,3,4]: #range(5):
    for x, y in zip(final_list, table_list):
        if '--rep {}'.format(rep) in x:
            if len(groups) >= ((num_submit+1) * group_size): break
            if os.path.exists(y): print('skip {}'.format(y)); continue
            groups.append((x, y))

groupby = lambda l, n: [tuple(l[i:i+n]) for i in range(0, len(l), n)] 

#print(len(groups))

for values in groupby(groups, group_size):
    values = list(values)
    #print(values)
    with open('run_script.sh', 'w') as f: f.write(template.format('\\n'.join([x for x, y in values])))
    if num_submit > 0:
        if not dry_run: os.system('qsub run_script.sh')
        for x, y in values:
            if not dry_run: os.system('touch {}'.format(y))
            print('qsub {}'.format(y))
            print('qsub {}'.format(x))
            print('')
    num_submit -= 1
    if num_submit <= 0: exit()


