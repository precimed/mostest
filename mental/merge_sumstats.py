import re
import pandas as pd
import os
from collections import defaultdict

gwas_dir = "/cluster/projects/p33/users/alexeas/most_mental/new_start/gwas/touchscreen_20127"
out_dir = "/cluster/projects/p33/users/alexeas/most_mental/new_start/gwas_merged/touchscreen"

# {pheno_fid: {orig_perm: [file_names] }, ... }
pheno_file_dict = defaultdict(lambda: defaultdict(list))
for de in os.scandir(gwas_dir):
    if ".glm." in de.name: # sumstats file
        splt_name = de.name.split('.') # perm.chr19_snps_00.1940_0_0.glm.logistic.hybrid | orig.chr4_snps_07.4609_0_0_quant.glm.linear
        pheno_file_dict[splt_name[2]][splt_name[0]].append(de.name)

print(f"{len(pheno_file_dict)} phenotypes detected")

for pheno_fid in pheno_file_dict:
    for orig_perm, file_names in pheno_file_dict[pheno_fid].items():
        file_names.sort(key=lambda fname: fname.split('.')[1].split('_snps_'))
        print(f"Processing {pheno_fid} {orig_perm} with {len(file_names)} files")
        dfs = [pd.read_table(os.path.join(gwas_dir,fname)) for fname in file_names]
        df = pd.concat(dfs, axis=0, sort=False)
        outf = os.path.join(out_dir, f"{orig_perm}.{pheno_fid}.csv")
        df.to_csv(outf, sep='\t', index=False)
        print(f"    {len(df)} SNPs saved to {outf}")

print(f"{len(pheno_file_dict)} different phenotypes")

