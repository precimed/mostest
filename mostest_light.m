% optional arguments
if ~exist('chunk', 'var'), chunk = 10000; end;
if ~exist('lam_reg', 'var'), lam_reg = 1.0; end;  %  default is to disable pre-whitening filter
if ~exist('num_eigval_to_regularize', 'var'), num_eigval_to_regularize = 0; end;  %  how many smallest eigenvalues of C0 matrix (z score correlation) to regularize
if ~exist('apply_int', 'var'), apply_int = true; end;       % apply rank-based inverse normal transform

% required input
if ~exist('out', 'var'),   error('out file prefix is required'); end
if ~exist('pheno', 'var'), error('pheno file is required'); end
if ~exist('bfile', 'var'), error('bfile is required'); end
if ~exist('snps', 'var'), snps=nan; end; 
if ~exist('nsubj', 'var'), nsubj=nan; end;

if exist('Shuffle') ~= 3, mex 'Shuffle.c'; end;   % ensure Shuffle is compiled

tic

fileID = fopen(sprintf('%s.bim', bfile));
bim_file = textscan(fileID,'%s %s %s %s %s %s');
fclose(fileID);
if isfinite(snps) && (snps ~= length(bim_file{1})), error('snps=%i is incompatible with .bim file; please check your snps parameter (or remove it to auto-detect #snps)', snps);end
snps=length(bim_file{1});

fileID = fopen(sprintf('%s.fam', bfile));
fam_file = textscan(fileID,'%s %s %s %s %s %s');
fclose(fileID);
if isfinite(nsubj) && (nsubj ~= length(fam_file{1})), error('nsubj=%i is incompatible with .fam file; please check your snps parameter (or remove it to auto-detect nsubj)', nsubj);end
nsubj=length(fam_file{1});

fprintf('%i snps and %i subjects detected in bfile\n', snps, nsubj);

fprintf('Loading phenotype matrix from %s... ', pheno);
if 1 
    ymat_df = readtable(pheno, 'Delimiter', 'tab','HeaderLines',0);
    measures = ymat_df.Properties.VariableNames;
    ymat_orig = table2array(ymat_df);
else
    ymat_orig=dlmread(pheno); ymat_orig=ymat_orig(:, 2:end);
    measures = cell(size(ymat_orig, 2), 1);
    for i=1:length(measures), measures{i} = sprintf('V%i', i); end;
end

npheno=size(ymat_orig, 2);
fprintf('Done, %i phenotypes found\n', npheno);
if size(ymat_orig, 1) ~= nsubj, error('roi matrix has info for %i subjects, while nsubj argument is specified as %i. These must be consistent.', size(ymat_orig, 1), nsubj); end;

keep = (min(ymat_orig)~=max(ymat_orig));
fprintf('Remove %i phenotypes (no variation)\n', length(keep) - sum(keep));
ymat_orig = ymat_orig(:, keep);
measures = measures(keep);
npheno=size(ymat_orig, 2);

% perform rank-based inverse normal transform, equivalently to the following R code:
% DM[,f] <- qnorm(ecdf(DM[,f])(DM[,f]) - 0.5/dim(DM)[1])
kurt = nan(npheno, 2);
for pheno_index=1:npheno
  vals = ymat_orig(:, pheno_index); kurt(pheno_index, 1) = kurtosis(vals);
  if apply_int
    [F, X] = ecdf(vals); F=transpose(F(2:end)); X=transpose(X(2:end));
    vals = norminv(interp1(X,F,vals,'nearest') - 0.5 / length(vals));
  end
  ymat_orig(:, pheno_index) = vals; kurt(pheno_index, 2) = kurtosis(vals);
end
fprintf('kurtosis before INT - %.2f %.2f (mean, max)\n', mean(kurt(:, 1)), max(kurt(:, 1)))
if apply_int, fprintf('kurtosis after  INT - %.2f %.2f (mean, max)\n', mean(kurt(:, 2)), max(kurt(:, 2))); end;

C = corr(ymat_orig);
C_reg = (1-lam_reg)*C + lam_reg*diag(max(0.01,diag(C))); % Ridge regularized covariance matrix
C_inv = inv(C_reg);
W_wht = chol(C_inv); % Whitening filter
ymat = ymat_orig*W_wht'; % Whitened residualized data

ymat_corr = corr(ymat);

% use correlation structure of the phenotypes, instead of genotype permutation scheme
C0 = ymat_corr;
C1 = ymat_corr;

[U S]  = svd(C0); s = diag(S);

% C0_reg = diag(diag(C0)); % Complete regularization -- results in imperfect gamma fit
% C0_reg = eye(size(C0)); % Complete regularization -- results in imperfect gamma fit
% max_lambda = s(min(10, length(s)));
% max_lambda = min(0.1, s(min(10, length(s)))); % oleksanf: don't regularize unless it's too bad

if (num_eigval_to_regularize > 0), max_lambda=s(end-num_eigval_to_regularize); else max_lambda = 0; end;
C0_reg = U*diag(max(max_lambda,s))*U'; % Good gamma fit

% C0_reg = U*diag(max(s(40),s))*U';
% C0_reg = C0;  % no regularization

fprintf('Perform GWAS on %s (%i SNPs are expected)...\n', bfile, snps)
nvec=zeros(snps, 1, 'single');
freqvec=zeros(snps, 1, 'single');

logpdfvecs = NaN(2,snps);  % logpdfvec
minpvecs = NaN(2,snps);    % minpvecs

for i=1:chunk:snps
  j=min(i+chunk-1, snps);
  fprintf('gwas: loading snps %i to %i... ', i, j);    tic;
  geno_int8 = PlinkRead_binary2(nsubj, i:j, bfile);
  fprintf('processing... ', i, j);   
  geno = nan(size(geno_int8), 'single'); for code = int8([0,1,2]), geno(geno_int8==code) = single(code); end;

  shuffle_geno = Shuffle(geno);
  [~, zmat_orig_chunk] = nancorr(ymat, geno);
  [~, zmat_perm_chunk] = nancorr(ymat, shuffle_geno);

  for orig_or_perm  = 1:2
    if orig_or_perm==1, zmat=zmat_orig_chunk'; else zmat=zmat_perm_chunk'; end;
    logpdfvecs(orig_or_perm,i:j) = dot(inv(C0_reg)*zmat', zmat');
    minpvecs(orig_or_perm,i:j) = 2*normcdf(-max(abs(zmat), [], 2));
  end
  
  nvec(i:j) = sum(isfinite(geno))';
  freqvec(i:j) = (1*sum(geno==1) + 2*sum(geno==2))' ./ (2*nvec(i:j));
  fprintf('done in %.1f sec, %.1f %% completed\n', toc, 100*(j+1)/snps);
end

gwas_time_sec = toc; tic

fprintf('running MOSTest analysis...')

maxlogpvecs = -log10(minpvecs);
ivec_snp_good = all(isfinite(logpdfvecs+minpvecs+maxlogpvecs));

[hc_maxlogpvecs hv_maxlogpvecs] = hist(maxlogpvecs(2,ivec_snp_good),1000); chc_maxlogpvecs = cumsum(hc_maxlogpvecs)/sum(hc_maxlogpvecs);
pd_minpvecs = fitdist(colvec(minpvecs(2,ivec_snp_good)),'beta'); % Not a great fit
%  pd_minpvecs.a = 1; % Hard-code known parameter (see http://www.di.fc.ul.pt/~jpn/r/prob/range.html)
[hc_logpdfvecs hv_logpdfvecs] = hist(logpdfvecs(2,ivec_snp_good),1000); chc_logpdfvecs = cumsum(hc_logpdfvecs)/sum(hc_logpdfvecs);
pd_logpdfvecs = fitdist(colvec(logpdfvecs(2,ivec_snp_good)),'gamma'); % Seems to work -- beta and wbl  do not

cdf_minpvecs=cdf(pd_minpvecs,10.^-hv_maxlogpvecs,'upper');
cdf_logpdfvecs = pd_logpdfvecs.cdf(hv_logpdfvecs);

maxlogpvecs_corr = -log10(cdf(pd_minpvecs,minpvecs));
logpdfvecs_corr = -log10(cdf(pd_logpdfvecs,logpdfvecs,'upper'));
fprintf('Done.\n')

fprintf('GWAS yield minP: %d; MOST: %d\n',sum(maxlogpvecs_corr(1,ivec_snp_good)>-log10(5e-8)),sum(logpdfvecs_corr(1,ivec_snp_good)>-log10(5e-8)));
fprintf('%i\t%.2f\t%.3f\t%.3f\t%.3f\t%.3f\t\n', npheno, cond(C0), pd_minpvecs.a, pd_minpvecs.b, pd_logpdfvecs.a, pd_logpdfvecs.b) 

most_time_sec = toc;

beta_params = [pd_minpvecs.a, pd_minpvecs.b]
gamma_params = [pd_logpdfvecs.a, pd_logpdfvecs.b]

minp_log10pval_orig = maxlogpvecs_corr(1, :);
most_log10pval_orig = logpdfvecs_corr(1, :);
minp_log10pval_perm = maxlogpvecs_corr(2, :);
most_log10pval_perm = logpdfvecs_corr(2, :);
fname=sprintf('%s.mat', out);
fprintf('saving %s... ', fname);
save(fname, '-v7', ...
 'most_log10pval_orig', 'minp_log10pval_orig', ...
 'most_log10pval_perm', 'minp_log10pval_perm', ...
 'nvec', 'freqvec', 'ivec_snp_good', ...
 'measures', 'ymat_corr', 'C0', 'C1', ...
 'hv_maxlogpvecs', 'hv_logpdfvecs', 'hc_maxlogpvecs', ...
 'chc_logpdfvecs', 'cdf_minpvecs', 'cdf_logpdfvecs', ...
 'beta_params', 'gamma_params', 'gwas_time_sec', 'most_time_sec');
fprintf('Done.\n')

fprintf('MOSTest analysis is completed.\n')
