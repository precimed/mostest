import os
import pandas as pd

def get_z_df(gwas_dir, bim_ref_df, traits, orig_perm):
    # orig_perm = 'orig' or 'perm'
    z_series_list = []
    for t in traits:
        print(f"Processing {t} {orig_perm}")
        stat_col = "T_STAT" if "quant" in t else "Z_STAT"
        gwas_file = os.path.join(gwas_dir, f"{orig_perm}.{t}.csv")
        # all SNP ids in these files are unique
        gwas_df = pd.read_csv(gwas_file, sep='\t', usecols=["ID", "A1", stat_col], index_col="ID")
        gwas_df = gwas_df.loc[bim_ref_df.SNP,:] # reorder based on the SNP order in the reference bim
        i_swap = (gwas_df.A1.values != bim_ref_df.A1.values)
        print(f"    {i_swap.sum()} signes swapped")
        gwas_df.loc[i_swap,stat_col] *= -1. # swap the sign of test statistics where allele order is different
        z_series = gwas_df[stat_col].reset_index(drop=True)
        i_na = z_series.isna().sum()
        print(f"    {i_na} NA values")
        z_series = z_series.fillna(0.).astype(float) # replace NA with 0 and convert to float
        z_series.name = f"p{t}" # change the name to pheno name, this will be column name in concatenated data frame
        z_series_list.append(z_series)
    z_df = pd.concat(z_series_list, axis=1)
    return z_df

# parameters ------------------
gwas_dir = "/cluster/projects/p33/users/alexeas/most_mental/new_start/gwas_merged/touchscreen"
out_dir = "/cluster/projects/p33/users/alexeas/most_mental/new_start/zmat"
out_suffix = "touchscreen"
bim_dir = "/cluster/projects/p33/projects/mental/geno/generic_qc"
# -----------------------------

orig_z_out_fname = os.path.join(out_dir, f"zmat.orig.{out_suffix}.csv.gz")
perm_z_out_fname = os.path.join(out_dir, f"zmat.perm.{out_suffix}.csv.gz")

orig_traits = [de.name.split('.')[1] for de in os.scandir(gwas_dir) if de.name.startswith('orig.')]
perm_traits = [de.name.split('.')[1] for de in os.scandir(gwas_dir) if de.name.startswith('perm.')]
orig_traits.sort()
perm_traits.sort()

assert orig_traits == perm_traits

print(f"{len(orig_traits)} orig and perm traits found")

# Read, concat and sort  reference bim files
print("Reading bim ref files")
bim_dfs = [pd.read_csv(de.path, na_filter=False, sep='\t', header=None, usecols=[0,1,3,4,5], names=["CHR","SNP","BP","A1", "A2"])
        for de in os.scandir(bim_dir) if de.name.endswith('.bim')]
print(f"{len(bim_dfs)} bim files found")
bim_ref_df = pd.concat(bim_dfs, axis=0, ignore_index=True)
bim_ref_df.sort_values(by=["CHR","BP","SNP"], inplace=True)

# Make merged bim file with ordering of SNPs same as in zmat
#bim_ref_df.to_csv("/cluster/projects/p33/users/alexeas/most_mental/new_start/mostest/ref.bim", sep='\t', index=False)

if False:
    # combine and save z data
    orig_z_df = get_z_df(gwas_dir, bim_ref_df, orig_traits, 'orig')
    orig_z_df.to_csv(orig_z_out_fname, sep='\t', index=False)
    print(f"Orig z mat saved to: {orig_z_out_fname}")

if True:
    perm_z_df = get_z_df(gwas_dir, bim_ref_df, orig_traits, 'perm') # orig_traits = perm_traits
    perm_z_df.to_csv(perm_z_out_fname, sep='\t', index=False)
    print(f"Perm z mat saved to: {perm_z_out_fname}")

