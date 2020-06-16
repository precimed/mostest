% =============== parameters section =============== 

% Required args:
%   zmat_names_array: array of file names with univariate GWAS results from the MOSTest analysis
%   out:              output file prefix
if ~exist('zmat_names_array', 'var'), error('zmat files array is required'); end              
if ~exist('out', 'var'), error('out file prefix is required'); end

% Optional args:
% num_eigval_to_keep:   how many largest eigenvalues of C0 matrix (z score correlation)
%                       to keep, the remaining will be assigned to the num_eigval_to_keep-th
%                       eigenvalue, num_eigval_to_keep = 0 - keep all
% use_paretotails:      fit tail of the mostest and minp statistics with pareto
% paretotails_quantile: a number close to 1.0, used as a second argument in MATLAB's paretotails
if ~exist('num_eigval_to_keep', 'var'), num_eigval_to_keep = 0; end
if ~exist('use_paretotails', 'var'), use_paretotails = false; end
if ~exist('paretotails_quantile', 'var'), paretotails_quantile = 0.99; end
% =============== end of parameters section ===============

tic

n_zmat = numel(zmat_names_array);
fprintf('Estimating MOSTest and MinP statistics\n')
fprintf('%d zmat files will be processed\n', n_zmat);

% combine z-score mtrixes
mostvecs = NaN(2,0); minpvecs = NaN(2,0); maxlogpvecs = NaN(2,0);
combined_nvec = zeros(0); combined_freqvec = zeros(0);
ivec_snp_good = false(0);
for i_zmat = 1:n_zmat
    zmat_name = zmat_names_array(i_zmat);
    fprintf('loading %s... ', zmat_name);
    load(zmat_name);
    fprintf('OK.\n')
    nsnps = size(zmat_orig, 1);
    npheno = size(zmat_orig, 2);

    combined_nvec = [combined_nvec; nvec];
    combined_freqvec = [combined_freqvec; freqvec];

    ivec_snp_good_i = all(isfinite(zmat_orig) & isfinite(zmat_perm), 2);
    % estimate correlation matrix and regularize only for the first iteration
    % This is not 100% accurate but should be accurate enough.
    if i_zmat == 1
        
        % use correlation structure of the z scores, calculated under permutation
        % we don't weight SNPs by LD because the permutation scheme breaks the LD structure
        % correlation structure of the null z scores
        C0 = corr(zmat_perm(ivec_snp_good_i, :));
        C1 = corr(zmat_orig(ivec_snp_good_i, :));
        [U S]  = svd(C0); s = diag(S);

        if (num_eigval_to_keep == 0) max_lambda = 0; else max_lambda = s(num_eigval_to_keep); end
        C0_reg = U*diag(max(max_lambda,s))*U'; % Good gamma fit
    end

    mostvecs_i = NaN(2,nsnps); minpvecs_i = NaN(2,nsnps); maxlogpvecs_i = NaN(2,nsnps);
    for i  = 1:2
      if i==1, zmat=zmat_orig; else zmat=zmat_perm; end;
      mostvecs_i(i,:) = dot(inv(C0_reg)*zmat', zmat');
      minpvecs_i(i,:) = 2*normcdf(-max(abs(zmat), [], 2));
      maxlogpvecs_i(i, :) = -log10(minpvecs_i(i, :));
    end
    mostvecs = [mostvecs, mostvecs_i];
    minpvecs = [minpvecs, minpvecs_i];
    maxlogpvecs = [maxlogpvecs, maxlogpvecs_i];
    ivec_snp_good = [ivec_snp_good; ivec_snp_good_i];
end

nvec = combined_nvec;
freqvec = combined_freqvec;

[hc_maxlogpvecs hv_maxlogpvecs] = hist(maxlogpvecs(2,ivec_snp_good),1000); chc_maxlogpvecs = cumsum(hc_maxlogpvecs)/sum(hc_maxlogpvecs);
[hc_mostvecs hv_mostvecs] = hist(mostvecs(2,ivec_snp_good),1000); chc_mostvecs = cumsum(hc_mostvecs)/sum(hc_mostvecs);

if use_paretotails
  pd_maxlogpvecs = paretotails(maxlogpvecs(2,ivec_snp_good), 0.0, paretotails_quantile);
  pd_minpvecs_params = upperparams(pd_maxlogpvecs);
  cdf_minpvecs = 1.0 - fixed_paretotails_cdf(pd_maxlogpvecs,hv_maxlogpvecs);
  maxlogpvecs_corr = -log10(fixed_paretotails_cdf(pd_maxlogpvecs, maxlogpvecs));

  pd_mostvecs = paretotails(mostvecs(2,ivec_snp_good),  0.0, paretotails_quantile);
  pd_mostvecs_params = upperparams(pd_mostvecs);
else
  pd_minpvecs = fitdist(colvec(minpvecs(2,ivec_snp_good)),'beta'); % Not a great fit
  pd_minpvecs_params = [pd_minpvecs.a, pd_minpvecs.b];
  cdf_minpvecs=cdf(pd_minpvecs,10.^-hv_maxlogpvecs,'upper');
  maxlogpvecs_corr = -log10(cdf(pd_minpvecs,minpvecs));

  pd_mostvecs = fitdist(colvec(mostvecs(2,ivec_snp_good)),'gamma'); % Seems to work -- beta and wbl  do not
  pd_mostvecs_params = [pd_mostvecs.a, pd_mostvecs.b];
end

if use_paretotails
    cdf_mostvecs = 1.0 - fixed_paretotails_cdf(pd_mostvecs,hv_mostvecs);
    mostvecs_corr = -log10(fixed_paretotails_cdf(pd_mostvecs,mostvecs));
else
    cdf_mostvecs = pd_mostvecs.cdf(hv_mostvecs);
    mostvecs_corr = -log10(cdf(pd_mostvecs,mostvecs,'upper'));
end

fprintf('Done.\n')

fprintf('GWAS yield minP: %d; MOST: %d\n',sum(maxlogpvecs_corr(1,ivec_snp_good)>-log10(5e-8)),sum(mostvecs_corr(1,ivec_snp_good)>-log10(5e-8)));
fprintf('%i\t%.2f\t%.3f\t%.3f\t%.3f\t%.3f\t\n', npheno, cond(C0), pd_minpvecs_params(1), pd_minpvecs_params(2), pd_mostvecs_params(1), pd_mostvecs_params(2)) 

most_time_sec = toc;

minp_log10pval_orig = maxlogpvecs_corr(1, :);
most_log10pval_orig = mostvecs_corr(1, :);
minp_log10pval_perm = maxlogpvecs_corr(2, :);
most_log10pval_perm = mostvecs_corr(2, :);
fname=sprintf('%s.mat', out);
fprintf('saving %s... ', fname);
save(fname, '-v7', ...
 'most_log10pval_orig', 'minp_log10pval_orig', ...
 'most_log10pval_perm', 'minp_log10pval_perm', ...
 'nvec', 'freqvec', 'ivec_snp_good', ...
 'measures', 'ymat_corr', 'C0', 'C1', ...
 'minpvecs', 'mostvecs', ...
 'hv_maxlogpvecs', 'hc_maxlogpvecs', 'chc_maxlogpvecs', 'cdf_minpvecs', ...
 'hv_mostvecs', 'hc_mostvecs', 'chc_mostvecs', 'cdf_mostvecs', ...
 'pd_minpvecs_params', 'pd_mostvecs_params', 'most_time_sec');
fprintf('Done.\n')
fprintf('MOSTest analysis is completed in %.2f sec.\n', most_time_sec)
