% =============== parameters section =============== 

% Required args:
%   out:   output file prefix
%   pheno: file with phenotyes
%   bfile: plink bfile prefix 
if ~exist('out', 'var'),   error('out file prefix is required'); end
if ~exist('pheno', 'var'), error('pheno file is required'); end
if ~exist('bfile', 'var'), error('bfile is required'); end

% Optional args:
%   chunk:                chunk size (how many SNPs to read at a time)
%   apply_int:            apply rank-based inverse normal transform
%   auto_compile_shuffle: automatically compile shuffle.mex
%   perform_cca:          perform canonical correlation analysis
%   lam_reg:              strength of pre-whitening filter, default is to disable
if ~exist('chunk', 'var'), chunk = 10000; end
if ~exist('apply_int', 'var'), apply_int = true; end
if ~exist('auto_compile_shuffle', 'var'), auto_compile_shuffle = 1; end
if ~exist('perform_cca', 'var'), perform_cca = false; end
if ~exist('lam_reg', 'var'), lam_reg = 1.0; end
      
% =============== end of parameters section =============== 

if auto_compile_shuffle && (exist('Shuffle') ~= 3), mex 'Shuffle.c'; end;   % ensure Shuffle is compiled

tic

fileID = fopen(sprintf('%s.bim', bfile));
bim_file = textscan(fileID,'%s %s %s %s %s %s');
fclose(fileID);
snps=length(bim_file{1});

fileID = fopen(sprintf('%s.fam', bfile));
fam_file = textscan(fileID,'%s %s %s %s %s %s');
fclose(fileID);
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

fprintf('Perform GWAS on %s (%i SNPs are expected)...\n', bfile, snps)
zmat_orig=zeros(snps, npheno, 'single');
zmat_perm=zeros(snps, npheno, 'single');
beta_orig=zeros(snps, npheno, 'single');  % skip saving p-values and standard errors (SE)
beta_perm=zeros(snps, npheno, 'single');  % (can be derived from Z and BETA)
nvec=zeros(snps, 1, 'single');
freqvec=zeros(snps, 1, 'single');
zvec_cca=nan(snps, 2);

for i=1:chunk:snps
  j=min(i+chunk-1, snps);
  fprintf('gwas: loading snps %i to %i... ', i, j);    tic;
  geno_int8 = PlinkRead_binary2(nsubj, i:j, bfile);
  fprintf('processing... ', i, j);   
  geno = nan(size(geno_int8), 'single'); for code = int8([0,1,2]), geno(geno_int8==code) = single(code); end;

  shuffle_geno = Shuffle(geno);
  [rmat_orig_chunk, zmat_orig_chunk] = nancorr(ymat, geno);
  [rmat_perm_chunk, zmat_perm_chunk] = nancorr(ymat, shuffle_geno);

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

  zmat_orig(i:j, :) = zmat_orig_chunk';
  zmat_perm(i:j, :) = zmat_perm_chunk';
  beta_factor = std(ymat)' * (1./std(geno, 'omitnan'));
  beta_orig(i:j, :) = transpose(rmat_orig_chunk .* beta_factor);
  beta_perm(i:j, :) = transpose(rmat_perm_chunk .* beta_factor);
  nvec(i:j) = sum(isfinite(geno))';
  freqvec(i:j) = (1*sum(geno==1) + 2*sum(geno==2))' ./ (2*nvec(i:j));
  fprintf('done in %.1f sec, %.1f %% completed\n', toc, 100*(j+1)/snps);
end

gwas_time_sec = toc;

fname = sprintf('%s_zmat.mat', out);
fprintf('saving %s as -v7.3... ', fname);
save(fname, '-v7.3', 'zmat_orig', 'zmat_perm', 'beta_orig', 'beta_perm', 'measures', 'nvec', 'zvec_cca', 'freqvec', 'ymat_corr', 'gwas_time_sec');
fprintf('OK.\n')
