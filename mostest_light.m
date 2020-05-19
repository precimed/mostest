% =============== parameters section =============== 

% required input
if ~exist('out', 'var'),   error('out file prefix is required'); end
if ~exist('mat_name', 'var'), mat_name = ''; end;                                 % re-use minpvecs and mostvecs (i.e. minP and MOSTest test statistics) from previously saved <out>.mat file
if isempty(mat_name)  
  if ~exist('pheno', 'var'), error('pheno file is required'); end
  if ~exist('bfile', 'var'), error('bfile is required'); end
end

% optional arguments
if ~exist('chunk', 'var'), chunk = 10000; end;                                    % chunk size (how many SNPs to read at a time)
if ~exist('num_eigval_to_keep', 'var'), num_eigval_to_keep = 0; end;              % how many largest eigenvalues of C0 matrix (z score correlation) to keep, the remaining will be assigned to the num_eigval_to_keep-th eigenvalue, num_eigval_to_keep = 0 - keep all
if ~exist('apply_int', 'var'), apply_int = true; end;                             % apply rank-based inverse normal transform
if ~exist('auto_compile_shuffle', 'var'), auto_compile_shuffle = 1; end;          % automatically compile shuffle.mex
if ~exist('use_paretotails', 'var'), use_paretotails = false; end;                % use paretotails instead of the gamma and beta functions to fit the distribution of the MOSTest & minP test statistic under null
  
% debug features - internal use only
if ~exist('lam_reg', 'var'), lam_reg = 1.0; end;  %  default is to disable pre-whitening filter
if ~exist('meta_simu_num_cohors', 'var'), meta_simu_num_cohors = 1; end;  % number of cohorts to simulate meta-analysis
if ~exist('snps', 'var'), snps=nan; end;                                          % number of SNPs in the analysis
if ~exist('nsubj', 'var'), nsubj=nan; end;                                        % number of subjects in the analysis
if ~exist('paretotails_quantile', 'var'), paretotails_quantile = 0.99; end;       % a number close to 1.0, used as a second argument in MATLAB's paretotails

% =============== end of parameters section =============== 

if auto_compile_shuffle && (exist('Shuffle') ~= 3), mex 'Shuffle.c'; end;   % ensure Shuffle is compiled

tic

if isempty(mat_name)
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
      ymat_df = readtable(pheno, 'Delimiter', 'tab');
      measures = ymat_df.Properties.VariableNames;
      ymat_orig = table2array(ymat_df);
  else
      % an alternative helper code that reads phenotype matrix without a header
      ymat_orig=dlmread(pheno); ymat_orig=ymat_orig(:, 2:end);
      measures = cell(size(ymat_orig, 2), 1);
      for i=1:length(measures), measures{i} = sprintf('V%i', i); end;
  end

  npheno=size(ymat_orig, 2);
  fprintf('Done, %i phenotypes found\n', npheno);
  if size(ymat_orig, 1) ~= nsubj, error('roi matrix has info for %i subjects, while nsubj argument is specified as %i. These must be consistent.', size(ymat_orig, 1), nsubj); end;

  if (meta_simu_num_cohors > 1)
    p=rand(meta_simu_num_cohors, 1); p=p/sum(p); 
    nsubj_per_cohort=mnrnd(nsubj, p, 1);
    nsubj_per_cohort = nsubj_per_cohort(nsubj_per_cohort>=2);
    meta_simu_num_cohors = length(nsubj_per_cohort);
    fprintf('simulated distribution of subjects across %i meta-cohorts:\n', meta_simu_num_cohors);
    disp(nsubj_per_cohort);
  end

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

  if (num_eigval_to_keep > 0), max_lambda=s(num_eigval_to_keep); else max_lambda = 0; end;
  C0_reg = U*diag(max(max_lambda,s))*U'; % Good gamma fit

  % C0_reg = U*diag(max(s(40),s))*U';
  % C0_reg = C0;  % no regularization

  fprintf('Perform GWAS on %s (%i SNPs are expected)...\n', bfile, snps)
  nvec=zeros(snps, 1, 'single');
  freqvec=zeros(snps, 1, 'single');

  mostvecs = NaN(2,snps);
  minpvecs = NaN(2,snps);

  for i=1:chunk:snps
    j=min(i+chunk-1, snps);
    fprintf('gwas: loading snps %i to %i... ', i, j);    tic;
    geno_int8 = PlinkRead_binary2(nsubj, i:j, bfile);
    fprintf('processing... ', i, j);   
    geno = nan(size(geno_int8), 'single'); for code = int8([0,1,2]), geno(geno_int8==code) = single(code); end;

    if meta_simu_num_cohors==1
      shuffle_geno = Shuffle(geno);
      [~, zmat_orig_chunk] = nancorr(ymat, geno);
      [~, zmat_perm_chunk] = nancorr(ymat, shuffle_geno);
    else
      % experimental "what-if" scenario to validate MOSTest in a meta-analytical setting
      zmat_orig_chunk = zeros(npheno, j-i+1);
      zmat_perm_chunk = zeros(npheno, j-i+1);

      subj_per_cohort_start = [1, 1+cumsum(nsubj_per_cohort(1:end-1))];
      subj_per_cohort_end = cumsum(nsubj_per_cohort);
      weight_per_cohort=sqrt(nsubj_per_cohort)/sqrt(sum(nsubj_per_cohort));


      for meta_cohort = 1:meta_simu_num_cohors
        subj_per_cohort_idx = subj_per_cohort_start(meta_cohort):subj_per_cohort_end(meta_cohort);
        ymat_cohort = ymat(subj_per_cohort_idx, :);
        geno_cohort = geno(subj_per_cohort_idx, :);
        shuffle_geno_cohort = Shuffle(geno_cohort);


        [~, zmat_orig_chunk_cohort] = nancorr(ymat_cohort, geno_cohort);
        [~, zmat_perm_chunk_cohort] = nancorr(ymat_cohort, shuffle_geno_cohort);
        
        zmat_orig_chunk = zmat_orig_chunk + weight_per_cohort(meta_cohort) * zmat_orig_chunk_cohort;
        zmat_perm_chunk = zmat_perm_chunk + weight_per_cohort(meta_cohort) * zmat_perm_chunk_cohort;
      end
    end
    
    for orig_or_perm  = 1:2
      if orig_or_perm==1, zmat=zmat_orig_chunk'; else zmat=zmat_perm_chunk'; end;
      mostvecs(orig_or_perm,i:j) = dot(inv(C0_reg)*zmat', zmat');
      minpvecs(orig_or_perm,i:j) = 2*normcdf(-max(abs(zmat), [], 2));
    end
    
    nvec(i:j) = sum(isfinite(geno))';
    freqvec(i:j) = (1*sum(geno==1) + 2*sum(geno==2))' ./ (2*nvec(i:j));
    fprintf('done in %.1f sec, %.1f %% completed\n', toc, 100*(j+1)/snps);
  end
else
  fprintf('loading %s... ', mat_name);
  load(mat_name);
  fprintf('OK.\n')
  npheno=size(C0, 1);
end

gwas_time_sec = toc; tic

fprintf('running MOSTest analysis...')

maxlogpvecs = -log10(minpvecs);
ivec_snp_good = all(isfinite(mostvecs+minpvecs+maxlogpvecs));

[hc_maxlogpvecs hv_maxlogpvecs] = hist(maxlogpvecs(2,ivec_snp_good),1000); chc_maxlogpvecs = cumsum(hc_maxlogpvecs)/sum(hc_maxlogpvecs);
[hc_mostvecs hv_mostvecs] = hist(mostvecs(2,ivec_snp_good),1000); chc_mostvecs = cumsum(hc_mostvecs)/sum(hc_mostvecs);

if use_paretotails
  pd_maxlogpvecs = paretotails(maxlogpvecs(2,ivec_snp_good), 0.0, paretotails_quantile);
  pd_minpvecs_params = upperparams(pd_maxlogpvecs);
  cdf_minpvecs = pd_maxlogpvecs.cdf(hv_maxlogpvecs);
  maxlogpvecs_corr = -log10(pd_maxlogpvecs.cdf(maxlogpvecs,'upper'));

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

cdf_mostvecs = pd_mostvecs.cdf(hv_mostvecs);
mostvecs_corr = -log10(cdf(pd_mostvecs,mostvecs,'upper'));
fprintf('Done.\n')

fprintf('GWAS yield minP: %d; MOST: %d\n',sum(maxlogpvecs_corr(1,ivec_snp_good)>-log10(5e-8)),sum(mostvecs_corr(1,ivec_snp_good)>-log10(5e-8)));
fprintf('%i\t%.2f\t%.3f\t%.3f\t%.3f\t%.3f\t\n', npheno, cond(C0), pd_minpvecs_params(1), pd_minpvecs_params(2), pd_mostvecs_params(1), pd_mostvecs_params(2)) 

most_time_sec = toc;

minp_log10pval_orig = maxlogpvecs_corr(1, :);
most_log10pval_orig = mostvecs_corr(1, :);
minp_log10pval_perm = maxlogpvecs_corr(2, :);
most_log10pval_perm = mostvecs_corr(2, :);
fname=sprintf('%s.mat', out);
hv_logpdfvecs=hv_mostvecs; cdf_logpdfvecs=cdf_mostvecs; chc_logpdfvecs=chc_mostvecs;
fprintf('saving %s... ', fname);
save(fname, '-v7', ...
 'most_log10pval_orig', 'minp_log10pval_orig', ...
 'most_log10pval_perm', 'minp_log10pval_perm', ...
 'nvec', 'freqvec', 'ivec_snp_good', ...
 'measures', 'ymat_corr', 'C0', 'C1', ...
 'minpvecs', 'mostvecs', ...
 'hv_maxlogpvecs', 'hc_maxlogpvecs', 'chc_maxlogpvecs', 'cdf_minpvecs', ...
 'hv_mostvecs', 'hc_mostvecs', 'chc_mostvecs', 'cdf_mostvecs', ...
 'pd_minpvecs_params', 'pd_mostvecs_params', 'gwas_time_sec', 'most_time_sec');
fprintf('Done.\n')

fprintf('MOSTest analysis is completed.\n')
