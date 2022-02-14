import os
import re

def get_bfile(slurm_lines):
    bfile = None
    for l in slurm_lines:
        l = l.strip()
        if l.startswith("--bfile"):
            bfile = l.replace('--bfile','').strip()
            break
    return bfile

def get_complete_pheno(slurm_lines):
    complete_pheno = []
    for l in slurm_lines:
        if "CANCELLED AT" in l:
            break
        else:
            l = l.strip()
            if l.startswith('--glm'):
                p = re.search("'([0-9_]+)'",l)
                if p:
                    complete_pheno.append(p.group(1))
    return complete_pheno


def get_pheno(slurm_lines):
    pheno_file = None
    for l in slurm_lines:
        l = l.strip()
        if l.startswith("--pheno"):
            pheno_file = l.replace('--pheno','').strip()
            break
    if pheno_file is None:
        raise Exception("No --pheno file is found in slurm file.")
    else:
        with open(pheno_file) as f:
            header = f.readline().strip().split("\t")
            pheno = [n for n in header if n!="IID" and n!="FID"]
    return pheno

def parse_slurm(slurm_file):
    # Return: string = "BFILE_NAME:pheno1,pheno2", only phenotypes which are not completed are listed after BFILE_NAME
    with open(slurm_file) as f:
        slurm_lines = f.readlines()
    is_cancelled = any("CANCELLED AT" in l for l in slurm_lines)
    if is_cancelled:
        pheno = get_pheno(slurm_lines)
        bfile = get_bfile(slurm_lines)
        complete_pheno = get_complete_pheno(slurm_lines)
        cancelled_pheno = [p for p in pheno if not p in complete_pheno]
        if cancelled_pheno:
            bfile_pheno = f"{bfile}:{','.join(cancelled_pheno)}"
            return bfile_pheno
        else:
            return None
    else:
        return None


def main():
    slurm_dir = "/cluster/projects/p33/users/alexeas/most_mental/new_start/src" # folder with slurm*.out files

    bfile_pheno_arr = []
    for de in os.scandir(slurm_dir):
        if de.is_file() and de.name.startswith("slurm-") and de.name.endswith(".out"):
            bfile_pheno = parse_slurm(de.path)
            if bfile_pheno:
                bfile_pheno_arr.append(bfile_pheno)
    print(f"{len(bfile_pheno_arr)} cancelled jobs will be restarted.")
    bfile_pheno_str = ' '.join(bfile_pheno_arr)
    command = f"sbatch --array 0-{len(bfile_pheno_arr)-1}%600 gwas.pheno_names.array_job.sh {bfile_pheno_str}"
    #os.system(command)

main()

