import pandas as pd
import numpy as np

pheno_type = "multresp"
zmat_file_orig = f"zmat/{pheno_type}_orig_zmat_201217.csv.gz"
zmat_file_perm = f"zmat/{pheno_type}_perm_zmat_201217.csv.gz"
out_file_orig = f"zmat/{pheno_type}_orig_zmat_201217.filtered.csv.gz"
out_file_perm = f"zmat/{pheno_type}_perm_zmat_201217.filtered.csv.gz"

z_threshold = 30
pheno_max_na_content = 0.05

df_orig = pd.read_csv(zmat_file_orig, sep='\t')
df_perm = pd.read_csv(zmat_file_perm, sep='\t')
assert (df_orig.columns == df_perm.columns).all() and (df_orig.index == df_perm.index).all()
df_orig.values[np.abs(df_orig.values)>z_threshold] = np.nan
df_perm.values[np.abs(df_perm.values)>z_threshold] = np.nan

cols2use = []
N = df_orig.shape[0]
for col in df_orig.columns:
    frac_good_orig = np.isfinite(df_orig[col]).sum()/N
    frac_good_perm = np.isfinite(df_perm[col]).sum()/N
    if (frac_good_orig > (1 - pheno_max_na_content)) and (frac_good_perm > (1 - pheno_max_na_content)):
        cols2use.append(col)

cols2drop = set(df_orig.columns) - set(cols2use)
print(f"{len(cols2drop)} columns with NA content > {pheno_max_na_content}:\n\t{', '.join(cols2drop)}")

if len(cols2drop) > 0:
    for df, outf in zip((df_orig, df_perm), (out_file_orig,out_file_perm)):
        df = df[cols2use]
        df.to_csv(outf,sep='\t',index=False,float_format='%.5f')
        print(f"saved to {outf}")

print("Done")
