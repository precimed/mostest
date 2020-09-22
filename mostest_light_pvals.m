% =============== parameters section =============== 

% Required args:
%   stats_file: mat file with mostest_light_stats.m results. Must contain variables:
%               nvec, freqvec, mostvecs_orig, minpvecs_orig, mostvecs_perm, minpvecs_perm.
%   out:        output file prefix
if ~exist('stats_file', 'var'), error('stats_file file is required'); end              
if ~exist('out', 'var'), error('out file prefix'); end

% Optional args:
% daletails_quantile: a number slightly below 1.0 (argument of daletails), ecdf is used to approximate distribution of most/minp statistics below this quantile.
% daletails_quantile_up: parameters of distribution are fitted based on most/minp statistics values between daletails_quantile and daletails_quantile_up. Nan results in automatic selection.
% daletails_distrib_most/daletails_distrib_minp: distribution type  to uses for extrapolation of the upper tail (above daletails_quantile) of most/minp statistics values (argument of daletails).
% maf_threshold: ignore all variants with maf < maf_threshold in MOSTest analysis
if ~exist('daletails_quantile', 'var') daletails_quantile = 0.9999; end
if ~exist('daletails_quantile_up', 'var') daletails_quantile_up = nan; end
if ~exist('daletails_distrib_most', 'var') daletails_distrib_most = 'gamma'; end
if ~exist('daletails_distrib_minp', 'var') daletails_distrib_minp = 'beta'; end
if ~exist('maf_threshold', 'var'), maf_threshold = 0.005; end;
% =============== end of parameters section ===============

tic

fprintf('loading %s... ', stats_file);
load(stats_file);
fprintf('Done.\n');
ivec_snp_good = isfinite(mostvecs_orig) & all(isfinite(mostvecs_perm),2) & isfinite(minpvecs_orig) & all(isfinite(minpvecs_perm),2);
ivec_snp_good = ivec_snp_good & (freqvec > maf_threshold); % ignore all SNPs with maf < maf_threshold
n_good = sum(ivec_snp_good);

% take further only valid SNPs
nvec = nvec(ivec_snp_good);
freqvec = freqvec(ivec_snp_good);
mostvecs_orig = mostvecs_orig(ivec_snp_good);
mostvecs_perm = mostvecs_perm(ivec_snp_good,:);
minpvecs_orig = minpvecs_orig(ivec_snp_good);
minpvecs_perm = minpvecs_perm(ivec_snp_good,:);

fprintf('%d/%d valid elements in original/permuted statistics\n', numel(mostvecs_orig), numel(mostvecs_perm));

% flatten permuted statistics
mostvecs_perm = reshape(mostvecs_perm,numel(mostvecs_perm),1);
minpvecs_perm = reshape(minpvecs_perm,numel(minpvecs_perm),1);
log_minpvecs_orig = -log10(minpvecs_orig);
log_minpvecs_perm = -log10(minpvecs_perm);

if isnan(daletails_quantile_up) daletails_quantile_up = 10/n_good; end

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

% fix potential log10(0) = Inf and log(NaN) = NaN issues for valid SNPs
minp_log10pval_orig(~isfinite(minp_log10pval_orig)) = -log10(eps(0));
minp_log10pval_perm(~isfinite(minp_log10pval_perm)) = -log10(eps(0));
most_log10pval_orig(~isfinite(most_log10pval_orig)) = -log10(eps(0));
most_log10pval_perm(~isfinite(most_log10pval_perm)) = -log10(eps(0));

fprintf('GWAS yield minP: %d; MOST: %d\n',sum(minp_log10pval_orig > -log10(5e-8)), sum(most_log10pval_orig > -log10(5e-8)));

most_time_sec = toc;

fname=sprintf('%s.mat', out);
fprintf('saving %s... ', fname);

save(fname, '-v7', ...
    'most_log10pval_orig', 'minp_log10pval_orig', ...
    'most_log10pval_perm', 'minp_log10pval_perm', ...
    'nvec', 'freqvec', 'ivec_snp_good', 'most_time_sec');
fprintf('Done.\n')
fprintf('MOSTest analysis is completed in %.2f sec.\n', most_time_sec)
