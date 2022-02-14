import sys
import pandas as pd
import numpy as np
import scipy.io
import scipy.stats
import argparse
import os


def parse_args(args):
    parser = argparse.ArgumentParser(
                description="Produce csv files for FUMA from mostest_light_stats results.")
    parser.add_argument("--mat", help="Path to mostest_light_pvals result mat file.")
    parser.add_argument("--bim", help="Path to corresponding bim file.")
    parser.add_argument("--n", default=None, help="Path to file with sample size.")
    parser.add_argument("--out", default=None, help="Output prefix.")
    return parser.parse_args(args)


args = parse_args(sys.argv[1:])
if args.out is None:
    args.out = os.path.splitext(args.mat)[0]

bim = pd.read_csv(args.bim, sep='\t', usecols=(1,0,3,4,5), header=None,
                  names='CHR SNP BP A1 A2'.split(), na_filter=False)

mat = scipy.io.loadmat(args.mat)
bim = bim[mat['ivec_snp_good'].flatten().astype('bool')]
if args.n:
    print(f"Reading N from {args.n}")
    n_df = pd.read_csv(args.n, sep='\t')
    bim = bim.merge(n_df, on="SNP", how='left')
else:
    bim['N'] = mat['nvec'].flatten().astype('int')

for m in ('most','minp'):
    for s in ('orig','perm'):
        print(f'processing {m} {s} ... ', end='', flush=True)
        pval = 10**(-mat[f'{m}_log10pval_{s}'].flatten())
        bim['PVAL'] = pval
        bim['Z'] = -scipy.stats.norm.ppf(0.5*pval)
        out_path = f'{args.out}.{m}.{s}.csv.gz'
        bim.to_csv(out_path, sep='\t', float_format='%g', index=False)
        print(out_path)

