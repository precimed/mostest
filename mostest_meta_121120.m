%% =============== parameters section =============== 

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
  
%% Load geno data
if auto_compile_shuffle && (exist('Shuffle') ~= 3), mex 'Shuffle.c'; end;   % ensure Shuffle is compiled
  
tic
 

%% Input files

input_bfile_folder = '/home/aihual/meta_29112020/input_bfile';%uigetdir('Select the bfile folder')
output_folder = '/home/aihual/meta_29112020/output';%uigetdir(input_folder,'Select the folder to save output');
input_pheno_folder = '/home/aihual/meta_29112020/input_pheno';%uigetdir('Select the pheno folder')

bim_files = dir(fullfile(input_bfile_folder,'*.bim'));
no_bfiles=length(bim_files);

%% Calculate weight_cohort, nsnps_cohort, npheno, C0_reg
nsubj_cohort=[];
%nsnps_cohort=[];

for i =1:no_bfiles   
    [~,file_name] = fileparts(bim_files(i).name);
    fileID = fopen(strcat(input_bfile_folder,'/',file_name,'.bim'));
    geno_cohort = textscan(fileID,'%s %s %s %s %s %s');
    snps=length(geno_cohort{1});
    %nsnps_cohort=[nsnps_cohort,snps];% nsnps_cohort was created for cases when cohorts have different snps lists. Now the plan is to align snps, so no use anymore. 
    
    pheno_cohort=strcat(input_pheno_folder,'/',file_name,'.txt');
    ymat_cohort=dlmread(pheno_cohort);
    nsubj=size(ymat_cohort,1);
    nsubj_cohort=[nsubj_cohort,nsubj]; 
    
    % check variability of pheno data, 
    %cohort_keep = (min(ymat_cohort)~=max(ymat_cohort));
    %fprintf('Remove %i phenotypes (no variation) from %s\n', length(cohort_keep) - sum(cohort_keep),file_name);
    %ymat_cohort = ymat_cohort(:, cohort_keep);
    
    % take the overlap of ymat_cohort outside of loop?
    
    % What's the strategy for the cases that there is no variability for one/multiple phenos in one/multiple cohorts?
    % Should we only keep the overlapping phenos for all cohorts?
    % Or should we filter the results after calculations as we did at mostest_light? By setting the zscore=0 for the cohorts for phenos that have no variability. 
    
    
    if i==1
        npheno=size(ymat_cohort,2);
    end
    fprintf('%i snps and %i subjects detected in bfile %s\n', snps, nsubj,file_name);
    fprintf('Final %i phenotypes found in pheno file %s\n', npheno,file_name);
    
    % perform rank-based inverse normal transform, equivalently to the following R code:
    % DM[,f] <- qnorm(ecdf(DM[,f])(DM[,f]) - 0.5/dim(DM)[1])
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

   % use correlation structure of the phenotypes, instead of genotype permutation scheme
   % here use C0_reg from cohort 1, will update to use the biggest cohort
   if i==1
     C0 = ymat_corr;
     C1 = ymat_corr;
     [U S]  = svd(C0); s = diag(S);

     if (num_eigval_to_regularize > 0), max_lambda=s(end-num_eigval_to_regularize); else max_lambda = 0; end;
     C0_reg = U*diag(max(max_lambda,s))*U'; % Good gamma fit
   end 
   
   % test if C0_reg is from the a different cohort
   %if i==2
   %  C0 = ymat_corr;
   %  C1 = ymat_corr;

   %  [U S]  = svd(C0); s = diag(S);
   %  if (num_eigval_to_regularize > 0), max_lambda=s(end-num_eigval_to_regularize); else max_lambda = 0; end;
   %  C0_reg = U*diag(max(max_lambda,s))*U'; % Good gamma fit
   %end 

end

weight_cohort=sqrt(nsubj_cohort)/sqrt(sum(nsubj_cohort))


%% Perform GWAS by chunk

fprintf('Perform GWAS on %i cohorts\n', no_bfiles)
nvec=zeros(snps, 1, 'single');
freqvec=zeros(snps, 1, 'single');

logpdfvecs = NaN(2,snps);  % logpdfvec
minpvecs = NaN(2,snps);    % minpvecs

chunk=10000;

for i=1:chunk:snps %need to change for different cohort with different nsnps
  j=min(i+chunk-1, snps);
  fprintf('gwas: loading snps %i to %i... \n', i, j);    tic;
  
  zmat_orig_chunk = zeros(npheno, j-i+1);
  zmat_perm_chunk = zeros(npheno, j-i+1);

  for no_cohort = 1:no_bfiles
      [~,file_name] = fileparts(bim_files(no_cohort).name);
      fileID = fopen(strcat(input_bfile_folder,'/',file_name,'.bim'));
      geno_cohort = textscan(fileID,'%s %s %s %s %s %s');
      
      bfile_cohort=strcat(input_bfile_folder,'/',file_name);
      
      geno_int8 = PlinkRead_binary2(nsubj_cohort(no_cohort), i:j, bfile_cohort);
   
      geno_cohort = nan(size(geno_int8), 'single'); 
      for code = int8([0,1,2]), geno_cohort(geno_int8==code) = single(code); end;
      
      pheno_cohort=strcat(input_pheno_folder,'/',file_name,'.txt');
      ymat_cohort=dlmread(pheno_cohort);
      
      shuffle_geno_cohort = Shuffle(geno_cohort);
      [~, zmat_orig_chunk_cohort] = nancorr(ymat_cohort, geno_cohort);
      [~, zmat_perm_chunk_cohort] = nancorr(ymat_cohort, shuffle_geno_cohort);
  


      zmat_orig_chunk_cohort(~isfinite(zmat_orig_chunk_cohort)) = 0;
      zmat_perm_chunk_cohort(~isfinite(zmat_perm_chunk_cohort)) = 0;
      
      zmat_orig_chunk = zmat_orig_chunk + weight_cohort(no_cohort) * zmat_orig_chunk_cohort;
      zmat_perm_chunk = zmat_perm_chunk + weight_cohort(no_cohort) * zmat_perm_chunk_cohort;

  end
  
  for orig_or_perm  = 1:2
    if orig_or_perm==1, zmat=zmat_orig_chunk'; else zmat=zmat_perm_chunk'; end;
    logpdfvecs(orig_or_perm,i:j) = dot(inv(C0_reg)*zmat', zmat');    % calculate MOSTest test statistic (ToDo: rename logpdfvecs -> mostestvec)
    minpvecs(orig_or_perm,i:j) = 2*normcdf(-max(abs(zmat), [], 2));
  end
  
  %% nvec and freqvec are not calculated yet, need to know what they are
  %nvec(i:j) = sum(isfinite(geno))';
  %freqvec(i:j) = (1*sum(geno==1) + 2*sum(geno==2))' ./ (2*nvec(i:j));
  fprintf('done in %.1f sec, %.1f %% completed\n', toc, 100*(j+1)/snps);
end

return
gwas_time_sec = toc; tic

%% Mostest analysis

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

%need to update after figure out what 'measure' is
measures=1;

%% Saving results

fprintf('GWAS yield minP: %d; MOST: %d\n',sum(maxlogpvecs_corr(1,ivec_snp_good)>-log10(5e-8)),sum(logpdfvecs_corr(1,ivec_snp_good)>-log10(5e-8)));
fprintf('%i\t%.2f\t%.3f\t%.3f\t%.3f\t%.3f\t\n', npheno, cond(C0), pd_minpvecs.a, pd_minpvecs.b, pd_logpdfvecs.a, pd_logpdfvecs.b) 

most_time_sec = toc;

beta_params = [pd_minpvecs.a, pd_minpvecs.b]
gamma_params = [pd_logpdfvecs.a, pd_logpdfvecs.b]

minp_log10pval_orig = maxlogpvecs_corr(1, :);

most_log10pval_orig = logpdfvecs_corr(1, :);
minp_log10pval_perm = maxlogpvecs_corr(2, :);
most_log10pval_perm = logpdfvecs_corr(2, :);

fname=strcat(output_folder,'/','result.mat')

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
