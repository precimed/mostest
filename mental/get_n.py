import os
import pandas as pd

gwas_dir = "/cluster/projects/p33/users/alexeas/most_mental/new_start/gwas_merged/touchscreen"
outf = "gwas_n.touchscreen.csv"
outf_max = outf.replace("gwas_n", "n_max")
series_list = []
for de in os.scandir(gwas_dir):
    if de.name.startswith("orig"):
        trait = de.name.split('.')[1]
        print(f"Processing {de.path}")
        s = pd.read_csv(de.path, sep='\t', usecols=["ID", "OBS_CT"], index_col="ID", squeeze=True)
        s.name = trait
        print(f"    {s.isna().sum()} NA values")
        series_list.append(s)
df = pd.concat(series_list, axis=1)
df.fillna(0, inplace=True)
df.to_csv(outf, sep='\t')
print(f"Saved to: {outf}")

n_max = df.max(axis=1)
n_max.name = "N"
n_max.to_csv(outf_max, sep='\t')

