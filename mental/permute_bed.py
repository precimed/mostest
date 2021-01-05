import numpy as np
import pandas as pd
import sys
import argparse

def parse_args(args):
    parser = argparse.ArgumentParser(description="Permute genotypes for each SNP in bed file.")
    parser.add_argument("--bfile", help="Plink bfile prefix to permute.")
    parser.add_argument("--out", help="Permuted output bfile prefix.")
    return parser.parse_args(args)


def get_byte_geno_map():
    """
    Construct mapping between bytes 0..255 and 4-element arrays of a1 genotypes
    from plink bed file.
    Return 256 x 4 array A, where A[i] = [a1, a2, a3, a4], each ai from {2, -1, 1, 0}.
    """
    genotype_codes = np.array([2, -1, 1, 0],dtype=np.int8)
    byte_map = np.empty((256,4), dtype=np.int8)
    for b in range(256):
        for a in range(4):
            byte_map[b,a] = genotype_codes[(b >> 2*a) & 3]
    return byte_map


def get_geno_byte_map():
    genotype_codes = np.array([2, -1, 1, 0],dtype=np.int8)
    geno_map = {}
    for b in range(256):
        four_geno = tuple(genotype_codes[(b >> 2*a) & 3] for a in range(4))
        geno_map[four_geno] = b
    return geno_map


def get_snp_geno(bed, i_snp, n_samples, byte_map=None):
    """
    Get {2, -1, 1, 0} genotype for i-th SNP from bed file.
    """
    if byte_map is None:
        byte_map = get_byte_geno_map()
    i_snp_geno = np.empty(bed.shape[1]*4, dtype=np.int8)
    j = 0
    for b in bed[i_snp]:
        i_snp_geno[j:j+4] = byte_map[b]
        j += 4
    return i_snp_geno[:n_samples]


def get_snp_bytes(snp_geno, geno_map=None):
    if geno_map is None:
        geno_map = get_geno_byte_map()
    n_samples = len(snp_geno)
    n_samples_mod4 = n_samples%4
    i_end = n_samples - n_samples_mod4
    bb = bytearray(geno_map[tuple(snp_geno[i:i+4])] for i in range(0,i_end,4))
    if n_samples_mod4 != 0:
        four_geno = tuple(list(snp_geno[-n_samples_mod4:]) + [0]*(4 - n_samples_mod4))
        bb.append(geno_map[four_geno])
    return bytes(bb)
           
 
def permute_bed(bfile_orig, bfile_perm):
    print(f'Permuting {bfile_orig}')
    # Read original bfile
    bim_path = f'{bfile_orig}.bim'
    bim = pd.read_csv(bim_path, sep='\t', header=None)
    fam_path = f'{bfile_orig}.fam'
    fam = pd.read_csv(fam_path, delim_whitespace=True, header=None)
    n_snps = bim.shape[0]
    n_samples = fam.shape[0]
    print(f'    {n_snps} SNPs')
    print(f'    {n_samples} samples')
    
    # Wright bim and fam
    bim_perm_path = f'{bfile_perm}.bim'
    bim.to_csv(bim_perm_path, sep='\t', index=False, header=False)
    fam_perm_path = f'{bfile_perm}.fam'
    fam.to_csv(fam_perm_path, sep='\t', index=False, header=False)
    
    # Read bed one SNP at a time, permute, and wright to permuted bed
    n_cols = n_samples//4
    if 4*n_cols != n_samples:
        n_cols += 1
    bed_path = f'{bfile_orig}.bed'
    bed = np.memmap(bed_path, dtype=np.uint8, offset=3, mode='r', shape=(n_snps,n_cols))
    
    byte_map = get_byte_geno_map()
    geno_map = get_geno_byte_map()
    magic = bytes([0x6c, 0x1b, 0x01])
    bed_perm_path = f'{bfile_perm}.bed'
    with open(bed_perm_path, 'wb') as bed_perm:
        bed_perm.write(magic)
        for i_snp in range(n_snps):
            snp_geno = get_snp_geno(bed, i_snp, n_samples, byte_map)
            np.random.shuffle(snp_geno)
            snp_bytes = get_snp_bytes(snp_geno, geno_map)
            bed_perm.write(snp_bytes)
            if i_snp%1000 == 0:
                print(f'    {i_snp+1} SNPs processed')
    

if __name__ == '__main__':
    args = parse_args(sys.argv[1:])
    permute_bed(args.bfile, args.out)
    print('Done')

