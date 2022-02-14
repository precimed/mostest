import pandas as pd
import scipy.stats as ss
import numpy as np
import os

mental_dir = "/mnt/seagate10/projects/mostest_mental"
zmat_file = os.path.join(mental_dir,"mental_export/zmat/ordilal_binary_continuous_multresp_zmat_201218.orig.csv")
bim_file = os.path.join(mental_dir,"ukb_imp_all_chr_v3_qc.bim")
out_dir = os.path.join(mental_dir,"mental_export/zmat/sumstats4ldsc")

zmat_df = pd.read_csv(zmat_file,sep='\t')
bim_df = pd.read_csv(bim_file,sep='\t',header=None,names='CHR SNP GP BP A1 A2'.split(),
        usecols='CHR SNP BP A1 A2'.split())

print(f'{zmat_df.shape[1]} columns in zmat file')
for i,c in enumerate(zmat_df.columns):
    print(f"processing columns {c} [{i+1}]")
    z = zmat_df[c].values
    p = 2*ss.norm.sf(np.abs(z))
    bim_df['Z'] = z
    bim_df['PVAL'] = p
    out_f = os.path.join(out_dir,f"{c}.sumstats.csv.gz")
    bim_df.to_csv(out_f,sep='\t',index=False)

