% parameters

if ~exist('zmat_orig_file','var'), zmat_orig_file = '/cluster/projects/p33/users/alexeas/most_mental/new_start/zmat/zmat.orig.personality12_cognition21.no_sumscores.csv'; end
if ~exist('zmat_perm_file','var'), zmat_perm_file = '/cluster/projects/p33/users/alexeas/most_mental/new_start/zmat/zmat.perm.personality12_cognition21.no_sumscores.csv'; end
if ~exist('num_eigval_to_keep', 'var'), num_eigval_to_keep = 10; end
if ~exist('out', 'var'), out=sprintf('results/personality12_cognition21_no_sumscores/personality12_cognition21_no_sumscores_eig%d',num_eigval_to_keep); end

if ~exist('daletails_quantile', 'var') daletails_quantile = 0.9999; end
if ~exist('daletails_quantile_up', 'var') daletails_quantile_up = nan; end
if ~exist('daletails_distrib_most', 'var') daletails_distrib_most = 'gamma'; end
if ~exist('daletails_distrib_minp', 'var') daletails_distrib_minp = 'beta'; end
if ~exist('maf_threshold', 'var'), maf_threshold = 0.005; end;

tic

if ~exist('zmat_orig_tab', 'var')
    fprintf('reading %s ... ', zmat_orig_file);
    zmat_orig_tab = readtable(zmat_orig_file, 'Delimiter', 'tab');
    fprintf('Done\n');
    fprintf('reading %s ... ', zmat_perm_file);
    zmat_perm_tab = readtable(zmat_perm_file, 'Delimiter', 'tab');
    fprintf('Done\n');
end

zmat_orig = table2array(zmat_orig_tab);
zmat_perm = table2array(zmat_perm_tab);

[n_snps, n_pheno] = size(zmat_orig);
[n_snps_perm, n_pheno_perm] = size(zmat_perm);
fprintf('%d SNPs and %d pheno in orig file\n', n_snps, n_pheno);
fprintf('%d SNPs and %d pheno in perm file\n', n_snps_perm, n_pheno_perm);

assert((n_snps == n_snps_perm) && (n_pheno == n_pheno_perm));

ivec_snp_good = all(isfinite(zmat_orig) & isfinite(zmat_perm), 2);
n_good = sum(ivec_snp_good);
fprintf('%d SNPs have valid zscores in all phenotypes\n', n_good);

C0 = corr(zmat_perm(ivec_snp_good, :));
[U S]  = svd(C0); s = diag(S);
if (num_eigval_to_keep == 0) max_lambda = 0; else max_lambda = s(num_eigval_to_keep); end
C0_reg = U*diag(max(max_lambda,s))*U';

mostvecs_orig = dot(inv(C0_reg)*zmat_orig', zmat_orig')';
minpvecs_orig = 2*normcdf(-max(abs(zmat_orig), [], 2));
mostvecs_perm = dot(inv(C0_reg)*zmat_perm', zmat_perm')';
minpvecs_perm = 2*normcdf(-max(abs(zmat_perm), [], 2));


mostvecs_orig = mostvecs_orig(ivec_snp_good);
minpvecs_orig = minpvecs_orig(ivec_snp_good);
mostvecs_perm = mostvecs_perm(ivec_snp_good);
minpvecs_perm = minpvecs_perm(ivec_snp_good);

log_minpvecs_orig = -log10(minpvecs_orig);
log_minpvecs_perm = -log10(minpvecs_perm);


if isnan(daletails_quantile_up) daletails_quantile_up = 1 - 10/n_good; end

fprintf('Estimating probability density and p-values for MinP ... ');
pd_log_minpvecs_perm = daletails(log_minpvecs_perm, daletails_quantile, daletails_quantile_up, daletails_distrib_minp);
% for permuted statistics take p-values only for the first permutation 
minp_log10pval_perm = -log10(pd_log_minpvecs_perm.cdf(log_minpvecs_perm(1:n_good),'upper'));
minp_log10pval_orig = -log10(pd_log_minpvecs_perm.cdf(log_minpvecs_orig,'upper'));
fprintf('Done.\n');


fprintf('Estimating probability density and p-values for MOSTest ... ');
pd_mostvecs_perm = daletails(mostvecs_perm, daletails_quantile, daletails_quantile_up, daletails_distrib_most);
most_log10pval_perm = -log10(pd_mostvecs_perm.cdf(mostvecs_perm(1:n_good),'upper'));
most_log10pval_orig = -log10(pd_mostvecs_perm.cdf(mostvecs_orig,'upper'));
fprintf('Done.\n')

%fprintf('Estimating probability density and p-values for MinP ... ');
%pd_log_minpvecs_perm = daletails(log_minpvecs_perm, 0.999, 1, 'gamma');
%minp_log10pval_perm = -log10(pd_log_minpvecs_perm.cdf(log_minpvecs_perm,'upper')); % -log10(fixed_paretotails_cdf(pd_log_minpvecs_perm, log_minpvecs_perm(ind_in_original_template)));
%minp_log10pval_orig = -log10(pd_log_minpvecs_perm.cdf(log_minpvecs_orig,'upper')); % -log10(fixed_paretotails_cdf(pd_log_minpvecs_perm, log_minpvecs_orig));
%fprintf('Done.\n');

%fprintf('Estimating probability density and p-values for MOSTest ... ');
%pd_mostvecs_perm = daletails(mostvecs_perm, 0.999, 1, 'gamma'); % paretotails(mostvecs_perm,  0.0, paretotails_quantile);
%most_log10pval_perm = -log10(pd_mostvecs_perm.cdf(mostvecs_perm,'upper')); % -log10(fixed_paretotails_cdf(pd_mostvecs_perm, mostvecs_perm(ind_in_original_template)));
%most_log10pval_orig = -log10(pd_mostvecs_perm.cdf(mostvecs_orig,'upper')); % -log10(fixed_paretotails_cdf(pd_mostvecs_perm, mostvecs_orig));
%fprintf('Done.\n')

% fix potential log10(0) = Inf and log(NaN) = NaN issues for valid SNPs
minp_log10pval_orig(~isfinite(minp_log10pval_orig)) = -log10(eps(0));
minp_log10pval_perm(~isfinite(minp_log10pval_perm)) = -log10(eps(0));
most_log10pval_orig(~isfinite(most_log10pval_orig)) = -log10(eps(0));
most_log10pval_perm(~isfinite(most_log10pval_perm)) = -log10(eps(0));

fprintf('GWAS yield minP: %d; MOST: %d\n',sum(minp_log10pval_orig > -log10(5e-8)), sum(most_log10pval_orig > -log10(5e-8)));

% nvec and freqvec arrays should be replaced with the actual values
nvec = 300000*ones(n_good,1);
freqvec = 0.5*ones(n_good,1);

most_time_sec = toc;

fname=sprintf('%s.mat', out);
fprintf('saving %s... ', fname);

save(fname, '-v7', ...
    'most_log10pval_orig', 'minp_log10pval_orig', ...
    'most_log10pval_perm', 'minp_log10pval_perm', ...
    'nvec', 'freqvec', 'ivec_snp_good', 'most_time_sec');
fprintf('Done.\n')
fprintf('MOSTest analysis is completed in %.2f sec.\n', most_time_sec)

