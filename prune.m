function survive = prune(bfile, pval_vec, pval_thresh, r2_thresh)

% Args:
%   bfile:       plink bfile prefix (N SNPs and M samples). All chromosomes in corresponding bim file must be int.
%   pval_vec:    array of p-values, size(pval_vec) = N, the order of SNPs must correspond to the order in bfile.
%   pval_thresh: p-value threshold for pruning.
%   r2_thresh:   r2 threshold for pruning.
% Return:
%   survive:     logical array, size(survive) = N, survive(i) == 1 if i-th SNP
%                survives pruning, else survive(i) == 0.

verbose = false;

fam_file_id = fopen(sprintf('%s.fam', bfile));
fam_ids = textscan(fam_file_id,'%s %*[^\n]');
fclose(fam_file_id);
n_samples=length(fam_ids{1});

bim_file_id = fopen(sprintf('%s.bim', bfile));
chrs = textscan(bim_file_id,'%d %*[^\n]');
chrs = chrs{1}; % int32 vector
chrs = reshape(chrs, size(pval_vec));
fclose(bim_file_id);

survive = pval_vec < pval_thresh; % == 0 for NaN values

unique_chrs = unique(chrs);
for chr_i = 1:numel(unique_chrs)
    chr = unique_chrs(chr_i);
    if verbose fprintf('Processing chr %d\n', chr); end
    on_chr = chrs == chr; % on_chr(i) == 1 if i-th SNP is on chromosome chr, else on_chr(i) == 0
    if verbose fprintf('%d SNPs on chr\n', sum(on_chr)); end
    survive_on_chr = survive & on_chr;
    if verbose fprintf('%d SNPs have p < %f\n', sum(survive_on_chr), pval_thresh); end
    idx_to_prune = find(survive_on_chr);

    geno_int8 = PlinkRead_binary2(n_samples, idx_to_prune, bfile);
    geno = single(geno_int8);
    geno(geno < 0) = nan;
    geno_r2_mat = nancorr(geno, geno).^2;

    [sorted_pval, sorted_i] = sort(pval_vec(idx_to_prune));
    prune_order = idx_to_prune(sorted_i);

    if verbose fprintf('Start pruning ... '); end;
    for i = 1:numel(prune_order)
        i_snp = prune_order(i);
        j_snp = sorted_i(i);
        if survive_on_chr(i_snp)
            i_in_ld = idx_to_prune(geno_r2_mat(:,j_snp) > r2_thresh);
            survive_on_chr(i_in_ld) = 0;
            survive_on_chr(i_snp) = 1;
        end
    end
    if verbose fprintf('Done\n'); end
    if verbose fprintf('%d SNPs survive pruning\n', sum(survive_on_chr)); end
    survive(on_chr) = survive_on_chr(on_chr);
end
