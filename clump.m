function survive = clump(bfile, pval, pval_thresh, r2_thresh)

% Args:
%   bfile:       plink bfile prefix (N SNPs and M samples).
%   pval:        array of p-values, size(pval) = N, the order of SNPs must correspond to the order in bfile.
%   pval_thresh: p-value threshold for clumping.
%   r2_thresh:   r2 threshold for clumping
% Return:
%   survive:     logical array, size(survive) = N, survive(i) == 1 if i-th SNP
%                survives clumping, else survive(i) == 0.

fam_file_id = fopen(sprintf('%s.fam', bfile));
fam_ids = textscan(fam_file_id,'%s %*[^\n]');
fclose(fam_file_id);
n_samples=length(fam_ids{1});

survive = pval < pval_thresh; % == 0 for NaN values
idx_to_clump = find(survive);

geno_int8 = PlinkRead_binary2(n_samples, idx_to_clump, bfile);
geno = single(geno_int8);
geno(geno < 0) = nan;
geno_r2_mat = nancorr(geno, geno).^2;

[sorted_pval, sorted_i] = sort(pval(idx_to_clump));
clump_order = idx_to_clump(sorted_i);

for i = 1:numel(clump_order)
    i_snp = clump_order(i);
    j_snp = sorted_i(i);
    if survive(i_snp)
        i_in_ld = idx_to_clump(geno_r2_mat(:,j_snp) > r2_thresh);
        survive(i_in_ld) = 0;
        survive(i_snp) = 1;
    end
end