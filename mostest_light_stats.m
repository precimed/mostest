% =============== parameters section =============== 

% Required args:
%   out:   output file prefix
%   pheno: file with phenotyes
%   bfile: plink bfile prefix 
if ~exist('out', 'var'),   error('out file prefix is required'); end
if ~exist('pheno', 'var'), error('pheno file is required'); end
if ~exist('bfile', 'var'), error('bfile is required'); end

% Optional args:
% num_eigval_to_keep:     how many largest eigenvalues of C0 matrix (z score correlation)
%                         to keep, the remaining will be assigned to the num_eigval_to_keep-th
%                         eigenvalue, num_eigval_to_keep = 0 - keep all
%   num_perm:             number of genotype permutations
%   pheno_file_sep:       separator in the pheno file  
%   chunk:                chunk size (how many SNPs to read at a time)
%   apply_int:            apply rank-based inverse normal transform
%   auto_compile_shuffle: automatically compile shuffle.mex
if ~exist('num_eigval_to_keep', 'var'), num_eigval_to_keep = 0; end;
if ~exist('num_perm', 'var'), num_perm = 1; end;
if ~exist('pheno_file_sep', 'var'), pheno_file_sep = 'tab'; end;
if ~exist('chunk', 'var'), chunk = 10000; end
if ~exist('apply_int', 'var'), apply_int = true; end
if ~exist('auto_compile_shuffle', 'var'), auto_compile_shuffle = 1; end
      
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
[filepath,name,ext] = fileparts(pheno);
if strcmp(ext,'.mat')
    load(pheno);
    ymat_orig = pheno_mat;
else
    % load from text file
    pheno_tab = readtable(fname_in, 'Delimiter', pheno_file_sep);
    ymat_orig = table2array(pheno_tab);
end

npheno=size(ymat_orig, 2);
fprintf('Done, %i phenotypes found\n', npheno);
if size(ymat_orig, 1) ~= nsubj, error('roi matrix has info for %i subjects, while nsubj argument is specified as %i. These must be consistent.', size(ymat_orig, 1), nsubj); end;

keep = (min(ymat_orig)~=max(ymat_orig));
fprintf('Remove %i phenotypes (no variation)\n', length(keep) - sum(keep));
ymat_orig = ymat_orig(:, keep);
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

ymat = ymat_orig;
ymat_corr = corr(ymat);

% use correlation structure of the phenotypes, instead of genotype permutation scheme
C0 = ymat_corr;
C1 = ymat_corr;

[U S]  = svd(C0); s = diag(S);

if (num_eigval_to_keep > 0), max_lambda=s(num_eigval_to_keep); else max_lambda = 0; end;
C0_reg = U*diag(max(max_lambda,s))*U'; % Good gamma fit


fprintf('Perform GWAS on %s (%i SNPs are expected)...\n', bfile, snps)
nvec=zeros(snps, 1, 'single');
freqvec=zeros(snps, 1, 'single');
mostvecs_orig = NaN(snps, 1);
minpvecs_orig = NaN(snps, 1);
mostvecs_perm = NaN(snps, num_perm);
minpvecs_perm = NaN(snps, num_perm);

Shuffle(rand(4,1),'seed');
inv_C0_reg = inv(C0_reg);

for i=1:chunk:snps
  j=min(i+chunk-1, snps);
  fprintf('gwas: loading snps %i to %i... ', i, j);    tic;
  geno_int8 = PlinkRead_binary2(nsubj, i:j, bfile);
  fprintf('processing... ', i, j);   
  geno = nan(size(geno_int8), 'single'); for code = int8([0,1,2]), geno(geno_int8==code) = single(code); end;

  [rmat_chunk, zmat_chunk] = nancorr(ymat, geno);
  mostvecs_orig(i:j) = dot(inv_C0_reg*zmat_chunk, zmat_chunk);
  minpvecs_orig(i:j) = 2*normcdf(-max(abs(zmat_chunk), [], 1));

  % process permuted genotypes
  for i_perm = 1:num_perm
      shuffle_geno = Shuffle(geno);
      [rmat_chunk, zmat_chunk] = nancorr(ymat, shuffle_geno);
      mostvecs_perm(i:j, i_perm) = dot(inv_C0_reg*zmat_chunk, zmat_chunk); 
      minpvecs_perm(i:j, i_perm) = 2*normcdf(-max(abs(zmat_chunk), [], 1));
  end

  nvec(i:j) = sum(isfinite(geno))';
  freqvec(i:j) = (1*sum(geno==1) + 2*sum(geno==2))' ./ (2*nvec(i:j));
  fprintf('done in %.1f sec, %.1f %% completed\n', toc, 100*(j+1)/snps);
end

% ensure that freqvec contains frequency of minor allele
i_major = freqvec > 0.5;
freqvec(i_major) = 1.0 - freqvec(i_major);

gwas_time_sec = toc;

fname = sprintf('%s.mat', out);
fprintf('saving %s as -v7.3... ', fname);
save(fname, '-v7.3', 'nvec', 'freqvec', 'gwas_time_sec', 'mostvecs_orig', 'minpvecs_orig', 'mostvecs_perm', 'minpvecs_perm');
fprintf('OK.\n')
