% example commands:
if 0
clear; pheno = 'SubcorticalVolume.csv'; out = 'SubcorticalVolume'; bfile = 'UKB26502_QCed_230519_maf0p005'; snps = 7428630; nsubj = 26502; mostest;
clear; pheno = 'CorticalArea.csv'; out = 'CorticalArea';           bfile = 'UKB26502_QCed_230519_maf0p005'; snps = 7428630; nsubj = 26502; mostest;
clear; pheno = 'CorticalThickness.csv'; out = 'CorticalThickness'; bfile = 'UKB26502_QCed_230519_maf0p005'; snps = 7428630; nsubj = 26502; mostest;
clear; pheno = 'all.csv'; out = 'all';                             bfile = 'UKB26502_QCed_230519_maf0p005'; snps = 7428630; nsubj = 26502; mostest;
clear; pheno = 'all.csv'; out = 'all_chr21'; chunk=200;            bfile = 'UKB26502_QCed_230519_maf0p005_chr21'; snps = 102079; nsubj = 26502; mostest;
clear; pheno = 'simu/rep1.csv'; out = 'simu/rep1';                 bfile = 'UKB26502_QCed_230519_maf0p005'; snps = 7428630; nsubj = 26502; mostest;
clear; pheno = 'simu/rep2.csv'; out = 'simu/rep2';                 bfile = 'UKB26502_QCed_230519_maf0p005'; snps = 7428630; nsubj = 26502; mostest;
clear; pheno = 'simu/rep3.csv'; out = 'simu/rep3';                 bfile = 'UKB26502_QCed_230519_maf0p005'; snps = 7428630; nsubj = 26502; mostest;
clear; zmat_name='all.nonorm_zmat.mat'; out = 'all.nonorm';  mostest;
python process_results.py UKB26502_QCed_230519_maf0p005.bim SubcorticalVolume.nonorm


clear; pheno = '/oasis/projects/nsf/csd604/oleksanf/MOSTest/all.csv'; out = '/oasis/projects/nsf/csd604/oleksanf/MOSTest/all';   chunk=200;   perform_cca=1;                       bfile = '/oasis/projects/nsf/csd604/oleksanf/MOSTest/UKB26502_QCed_230519_maf0p005_chr21'; snps = 102079; nsubj = 26502; mostest;

end

% optional arguments
if ~exist('chunk', 'var'), chunk = 10000; end;
if ~exist('lam_reg', 'var'), lam_reg = 1.0; end;  %  default is to disable pre-whitening filter
if ~exist('zmat_name', 'var'), zmat_name = ''; end;
if ~exist('perform_cca', 'var'), perform_cca = false; end;  % perform canonical correlation analysis
% required input
if ~exist('out', 'var'),   error('out file prefix is required'); end
if isempty(zmat_name)
  if ~exist('pheno', 'var'), error('pheno file is required'); end
  if ~exist('bfile', 'var'), error('bfile is required'); end
  if ~exist('snps', 'var'),  error('number of SNPs in --bfile must be specified as "snps" parameter'); end
  if ~exist('nsubj', 'var'), error('number of subjects in --bfile and pheno files must be specified as "nsubj" parameter'); end
end

if exist('Shuffle') ~= 3, mex 'Shuffle.c'; end;   % ensure Shuffle is compiled

if isempty(zmat_name)
  fprintf('Loading phenotype matrix from %s... ', pheno);
  ymat_df = readtable(pheno, 'Delimiter', 'tab');
  measures = ymat_df.Properties.VariableNames;
  ymat_orig = table2array(ymat_df);
  npheno=size(ymat_orig, 2);
  fprintf('Done, %i phenotypes found\n', npheno);
  if size(ymat_orig, 1) ~= nsubj, error('roi matrix has info for %i subjects, while nsubj argument is specified as %i. These must be consistent.', size(ymat_orig, 1), nsubj); end;

  C = corr(ymat_orig);
  C_reg = (1-lam_reg)*C + lam_reg*diag(max(0.01,diag(C))); % Ridge regularized covariance matrix
  C_inv = inv(C_reg);
  W_wht = chol(C_inv); % Whitening filter
  ymat = ymat_orig*W_wht'; % Whitened residualized data

  fprintf('Perform GWAS on %s (%i SNPs are expected)...\n', bfile, snps)
  zmat=zeros(snps, npheno, 2, 'single'); 
  nvec=zeros(snps, 1, 'single');
  zvec_cca=nan(snps, 2);

  for i=1:chunk:snps
    j=min(i+chunk-1, snps);
    fprintf('gwas: loading snps %i to %i... ', i, j);    tic;
    geno_int8 = PlinkRead_binary2(nsubj, i:j, bfile);
    fprintf('processing... ', i, j);   
    geno = nan(size(geno_int8), 'single'); for code = int8([0,1,2]), geno(geno_int8==code) = single(code); end;

    shuffle_geno = Shuffle(geno);
    [~, zmat_orig] = nancorr(ymat, geno);
    [~, zmat_perm] = nancorr(ymat, shuffle_geno);

    if perform_cca
      fprintf('cca... ');   
      ymat1 = [ymat, ones(size(ymat, 1), 1)];
      for k=i:j
        % These two are equivalent:
        % [b, bint, r, rint, stats] = regress(y, [X ones(n, 1)]);  stats(3)      % based on F-test
        % [A, B, r, U, V, statsCCA] = canoncorr(X, y);             statsCCA.p  
        [b, bint, r, rint, stats] = regress(geno(:,         k-i+1), ymat1); zvec_cca(k, 1) = stats(3);
        [b, bint, r, rint, stats] = regress(shuffle_geno(:, k-i+1), ymat1); zvec_cca(k, 2) = stats(3);
      end
    end

    zmat(i:j, :, 1) = zmat_orig';
    zmat(i:j, :, 2) = zmat_perm';
    nvec(i:j) = sum(isfinite(geno))';
    fprintf('done in %.1f sec, %.1f %% completed\n', toc, 100*(j+1)/snps);
  end

  zmat_orig = zmat(:, :, 1); zmat_perm = zmat(:, :, 2);
  fname = sprintf('%s_zmat.mat', out);
  fprintf('saving %s as -v7.3... ', fname);
  save(fname, '-v7.3', 'zmat_orig', 'zmat_perm', 'measures', 'nvec', 'zvec_cca');
  fprintf('OK.\n')
else
  fprintf('loading %s... ', zmat_name);
  load(zmat_name);
  fprintf('OK.\n')
  snps=size(zmat_orig, 1);
  npheno=size(zmat_orig, 2);
  zmat = zeros(snps, npheno, 2, 'single'); 
  zmat(:, :, 1) = zmat_orig;
  zmat(:, :, 2) = zmat_perm;
end

fprintf('running MOSTest analysis...')
ivec_snp_good = all(all(isfinite(zmat), 2), 3);

% correlation structure of the null z scores
C0 = corr(zmat(ivec_snp_good, :, 2));

% correlation structure of the real z scores
C1 = corr(zmat(ivec_snp_good, :, 1)); % & Hvec>0.1 & CRvec>0.95 & max(abs(zmat(:,:,1)),[],1)>abs(norminv(1e-5))),1)');

[U S]  = svd(C0); s = diag(S);

%  C0_reg = diag(diag(C0)); % Complete regularization -- results in imperfect gamma fit
%  C0_reg = eye(size(C0)); % Complete regularization -- results in imperfect gamma fit
%  max_lambda = s(min(10, length(s)));
%  max_lambda = min(0.1, s(min(10, length(s)))); % oleksanf: don't regularize unless it's too bad

% C0_reg = U*diag(max(max_lambda,s))*U'; % Good gamma fit
%  C0_reg = U*diag(max(s(40),s))*U';
C0_reg = C0;  % no regularization

logpdfvecs = NaN(size(zmat,3),size(zmat,1)); minpvecs = NaN(size(zmat,3),size(zmat,1)); maxlogpvecs = NaN(size(zmat,3),size(zmat,1));
for i  = 1:size(zmat,3)
  logpdfvecs(i,:) = -(mvnpdfln(zmat(:,:,i),0,C0_reg)-mvnpdfln(zeros(size(C0,1),1),0,C0_reg))/log(10); % Should be using mvncdf instead?
  minpvecs(i,:) = 2*normcdf(-max(abs(zmat(:, :, i)), [], 2));
  maxlogpvecs(i, :) = -log10(minpvecs(i, :));
end

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

minPval = maxlogpvecs_corr(1, :);
mostPval = logpdfvecs_corr(1, :);
fname=sprintf('%s.mat', out);
fprintf('saving %s... ', fname);
save(fname, '-v7', ...
 'mostPval', 'minPval', 'nvec', 'C0', 'C1', 'ivec_snp_good', ...
 'hv_maxlogpvecs', 'hv_logpdfvecs', 'hc_maxlogpvecs', ...
 'chc_logpdfvecs', 'cdf_minpvecs', 'cdf_logpdfvecs');
fprintf('Done.\n')

fprintf('MOSTest analysis is completed.\n')
