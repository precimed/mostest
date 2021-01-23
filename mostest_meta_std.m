%% =============== parameters section =============== 
% Inverse variance based meta
clear all;
close all ;
clc;

% required input
%if ~exist('out', 'var'),   error('out file prefix is required'); end
%if ~exist('pheno', 'var'), error('pheno file is required'); end
%if ~exist('bfile', 'var'), error('bfile is required'); end

% optional arguments
if ~exist('chunk', 'var'), chunk = 10000; end;                                    % chunk size (how many SNPs to read at a time)
if ~exist('num_eigval_to_regularize', 'var'), num_eigval_to_regularize = 0; end;  % how many smallest eigenvalues of C0 matrix (z score correlation) to regularize
if ~exist('apply_int', 'var'), apply_int = true; end;                             % apply rank-based inverse normal transform
if ~exist('auto_compile_shuffle', 'var'), auto_compile_shuffle = 1; end;          % automatically compile shuffle.mex
  
% debug features - internal use only
if ~exist('lam_reg', 'var'), lam_reg = 1.0; end;  %  default is to disable pre-whitening filter
%if ~exist('meta_simu_num_cohors', 'var'), meta_simu_num_cohors = 1; end;  % number of cohorts to simulate meta-analysis
if ~exist('snps', 'var'), snps=nan; end;                                          % number of SNPs in the analysis
if ~exist('nsubj', 'var'), nsubj=nan; end;                                        % number of subjects in the analysis
  

if auto_compile_shuffle && (exist('Shuffle') ~= 3), mex 'Shuffle.c'; end;   % ensure Shuffle is compiled
  
tic
%% Input files

input_bfile_folder = '/home/aihual/mostest_meta/input_bfile';%uigetdir('Select the bfile folder')
output_folder = '/home/aihual/mostest_meta/output';%uigetdir(input_folder,'Select the folder to save output');
input_pheno_folder = '/home/aihual/mostest_meta/input_pheno';%uigetdir('Select the pheno folder')

bim_files = dir(fullfile(input_bfile_folder,'*.bim'));
no_bfiles=length(bim_files);

%% Calculate covariance matrix C0 from pheno data out of chunk loop

nsubj=[];
for i =1:no_bfiles   
    [~,file_name] = fileparts(bim_files(i).name);
    fileID = fopen(strcat(input_bfile_folder,'/',file_name,'.bim'));
    geno_cohort = textscan(fileID,'%s %s %s %s %s %s');
    snps=length(geno_cohort{1});
    pheno_cohort=strcat(input_pheno_folder,'/',file_name,'.txt');
    ymat_cohort=dlmread(pheno_cohort);
    nsubj_cohort=size(ymat_cohort,1);
    nsubj=[nsubj nsubj_cohort];   
    % we suppose there is always variablity on pheno data, otherwise need to implement across all cohorts
    %measures=pheno.names, no pheno names in simulated pheno data yet, can
    %be added if needed
    % calculate covariance matrix of each cohort C0_cohort from pheno data
    
    % perform rank-based inverse normal transform, equivalently to the following R code:
    % DM[,f] <- qnorm(ecdf(DM[,f])(DM[,f]) - 0.5/dim(DM)[1])
    % corr(ymat)=cov(ymat) after performing rank-based inverse normal
    % transform?
    npheno=size(ymat_cohort,2);
    kurt = nan(npheno, 2);
    for pheno_index=1:npheno
     vals = ymat_cohort(:, pheno_index); kurt(pheno_index, 1) = kurtosis(vals);
    if apply_int
        [F, X] = ecdf(vals); F=transpose(F(2:end)); X=transpose(X(2:end));
        vals = norminv(interp1(X,F,vals,'nearest') - 0.5 / length(vals));
    end
    ymat_cohort(:, pheno_index) = vals; kurt(pheno_index, 2) = kurtosis(vals);
    end
    fprintf('kurtosis before INT - %.2f %.2f (mean, max)\n', mean(kurt(:, 1)), max(kurt(:, 1)))
    if apply_int, fprintf('kurtosis after  INT - %.2f %.2f (mean, max)\n', mean(kurt(:, 2)), max(kurt(:, 2))); end;
   
    ymat = ymat_cohort;
    ymat_corr = corr(ymat);
    % weight C0_cohort by sample size
    if i==1
        C0 = ymat_corr*nsubj_cohort;
        C1 = ymat_corr*nsubj_cohort;
    else
        C0 = C0+ymat_corr*nsubj_cohort;
        C1 = C1+ymat_corr*nsubj_cohort;
    end
end
C0=C0/sum(nsubj);
C1=C1/sum(nsubj);

[U S]  = svd(C0); s = diag(S);
if (num_eigval_to_regularize > 0), max_lambda=s(end-num_eigval_to_regularize); else max_lambda = 0; end;
C0_reg = U*diag(max(max_lambda,s))*U'; % Good gamma fit

%nvec=zeros(snps, 1, 'single');
%freqvec=zeros(snps, 1, 'single');

mostvecs = NaN(2,snps);  % logpdfvec
minpvecs = NaN(2,snps);    % minpvecs
chunk=10000;

for i=1:chunk:snps %all cohorts have been aligned to the same snps list
  j=min(i+chunk-1, snps);
  fprintf('gwas: loading snps %i to %i... \n', i, j);    tic;
  
  beta_orig_chunk = zeros(npheno, j-i+1);
  beta_perm_chunk = zeros(npheno, j-i+1);
  
  w_orig_chunk = zeros(npheno, j-i+1);
  w_perm_chunk = zeros(npheno, j-i+1);
  
  
  for no_cohort = 1:no_bfiles
      [~,file_name] = fileparts(bim_files(no_cohort).name);
      fileID = fopen(strcat(input_bfile_folder,'/',file_name,'.bim'));
      geno_cohort = textscan(fileID,'%s %s %s %s %s %s');
      bfile_cohort=strcat(input_bfile_folder,'/',file_name);
      geno_int8 = PlinkRead_binary2(nsubj(no_cohort), i:j, bfile_cohort);
      geno_cohort = nan(size(geno_int8), 'single');   
      for code = int8([0,1,2]), geno_cohort(geno_int8==code) = single(code); end;  
      pheno_cohort=strcat(input_pheno_folder,'/',file_name,'.txt');
      ymat_cohort=dlmread(pheno_cohort);
      
      % calculate beta and w 
      shuffle_geno_cohort = Shuffle(geno_cohort);
      
      [beta_chunk_cohort, w_chunk_cohort] = nancorr_std(ymat_cohort, geno_cohort);
      [beta_perm_chunk_cohort, w_perm_chunk_cohort] = nancorr_std(ymat_cohort, shuffle_geno_cohort);
  

      beta_chunk_cohort(~isfinite(w_chunk_cohort)) = 0;
      beta_perm_chunk_cohort(~isfinite(w_perm_chunk_cohort)) = 0;
      
      w_chunk_cohort(~isfinite(w_chunk_cohort)) = 0;
      w_perm_chunk_cohort(~isfinite(w_perm_chunk_cohort)) = 0;
      
      % sum(beta* w)
      beta_orig_chunk = beta_orig_chunk + beta_chunk_cohort .* w_chunk_cohort;
      beta_perm_chunk = beta_perm_chunk + beta_perm_chunk_cohort .* w_perm_chunk_cohort;
      
      % sum(w)
      w_orig_chunk = w_orig_chunk + w_chunk_cohort;
      w_perm_chunk = w_perm_chunk + w_perm_chunk_cohort;
  end
 
  % z=beta/SE
  SE_orig_chunk=sqrt(1.0./w_orig_chunk);
  beta_orig_chunk=beta_orig_chunk./w_orig_chunk;
  zmat_orig_chunk=beta_orig_chunk./SE_orig_chunk;
  
  SE_perm_chunk=sqrt(1.0./w_perm_chunk);
  beta_perm_chunk=beta_perm_chunk./w_perm_chunk;
  zmat_perm_chunk=beta_perm_chunk./SE_perm_chunk;
  
  for orig_or_perm  = 1:2
    if orig_or_perm==1, zmat=zmat_orig_chunk'; else zmat=zmat_perm_chunk'; end;
    mostvecs(orig_or_perm,i:j) = dot(inv(C0_reg)*zmat', zmat');    
    minpvecs(orig_or_perm,i:j) = 2*normcdf(-max(abs(zmat), [], 2));
  end
  
  %% no geno for final meta?
  %nvec(i:j) = sum(isfinite(geno))';
  %freqvec(i:j) = (1*sum(geno==1) + 2*sum(geno==2))' ./ (2*nvec(i:j));
  fprintf('done in %.1f sec, %.1f %% completed\n', toc, 100*(j+1)/snps);
end

gwas_time_sec = toc; tic

%% Mostest analysis

fprintf('running MOSTest analysis...')

maxlogpvecs = -log10(minpvecs);
ivec_snp_good = all(isfinite(mostvecs+minpvecs+maxlogpvecs));

[hc_maxlogpvecs hv_maxlogpvecs] = hist(maxlogpvecs(2,ivec_snp_good),1000); chc_maxlogpvecs = cumsum(hc_maxlogpvecs)/sum(hc_maxlogpvecs);
pd_minpvecs = fitdist(colvec(minpvecs(2,ivec_snp_good)),'beta'); % Not a great fit
%  pd_minpvecs.a = 1; % Hard-code known parameter (see http://www.di.fc.ul.pt/~jpn/r/prob/range.html)
[hc_logpdfvecs hv_logpdfvecs] = hist(mostvecs(2,ivec_snp_good),1000); chc_logpdfvecs = cumsum(hc_logpdfvecs)/sum(hc_logpdfvecs);
pd_logpdfvecs = fitdist(colvec(mostvecs(2,ivec_snp_good)),'gamma'); % Seems to work -- beta and wbl  do not

cdf_minpvecs=cdf(pd_minpvecs,10.^-hv_maxlogpvecs,'upper');
cdf_logpdfvecs = pd_logpdfvecs.cdf(hv_logpdfvecs);

maxlogpvecs_corr = -log10(cdf(pd_minpvecs,minpvecs));


logpdfvecs_corr = -log10(cdf(pd_logpdfvecs,mostvecs,'upper'));
fprintf('Done.\n')

%% Saving results

fprintf('GWAS yield minP: %d; MOST: %d\n',sum(maxlogpvecs_corr(1,ivec_snp_good)>-log10(5e-8)),sum(logpdfvecs_corr(1,ivec_snp_good)>-log10(5e-8)));
fprintf('%i\t%.2f\t%.3f\t%.3f\t%.3f\t%.3f\t\n', npheno, cond(C0), pd_minpvecs.a, pd_minpvecs.b, pd_logpdfvecs.a, pd_logpdfvecs.b) 

most_time_sec = toc;

beta_params = [pd_minpvecs.a, pd_minpvecs.b];
gamma_params = [pd_logpdfvecs.a, pd_logpdfvecs.b];

minp_log10pval_orig = maxlogpvecs_corr(1, :);

most_log10pval_orig = logpdfvecs_corr(1, :);
minp_log10pval_perm = maxlogpvecs_corr(2, :);
most_log10pval_perm = logpdfvecs_corr(2, :);

fname=strcat(output_folder,'/','result.mat');

fprintf('saving %s... ', fname);
save(fname, '-v7', ...
 'most_log10pval_orig', 'minp_log10pval_orig', ...
 'most_log10pval_perm', 'minp_log10pval_perm', ...
 'ivec_snp_good', ...
 'ymat_corr', 'C0', 'C1', ...
 'hv_maxlogpvecs', 'hv_logpdfvecs', 'hc_maxlogpvecs', ...
 'chc_logpdfvecs', 'cdf_minpvecs', 'cdf_logpdfvecs', ...
 'beta_params', 'gamma_params', 'gwas_time_sec', 'most_time_sec');
fprintf('Done.\n')

fprintf('MOSTest analysis is completed.\n')
