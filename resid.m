function pheno_resid = resid(pheno_file, ids_file, covar_file)

% Args:
%   pheno_file: a file with phenotypes to residualize (lh.pial.0.fsaverage3.csv)
%                  - no header
%                  - columns contain phenotypes (no other columns)
%   ids_file:   a file with imaging ids (vertex_subjectlist)
%                  - no header
%                  - first column must contain imaging ids
%                  - the number of rows should be the same as in phenotype file
%   covar_file: a file with covariates to use for pre-residualization (demographics.mostest.csv)
%                  - has header
%                  - must contain column "MRID" with imaging ids
%                  - order of ids in "FID" and "IID" columns must correspond to the order in plink fam
%
% Return:
%   pheno_resid: a matrix of residualized phenotypes [n_sampels x n_phenotypes]


% read phenotypes and remove constant
pheno_tab = readtable(pheno_file, 'Delimiter', ';');
pheno = table2array(pheno_tab);
[n_pheno_samples, n_pheno] = size(pheno);
fprintf('%d sampels, %d phenotypes in %s\n', n_pheno_samples, n_pheno, pheno_file);
i_const_pheno = all(pheno(1,:) == pheno, 1);
fprintf('%d constant phenotypes will be dropped\n', nnz(i_const_pheno));
pheno_tab = pheno_tab(:,~i_const_pheno);
n_nonconst_pheno = numel(pheno_tab.Properties.VariableNames);

% read imaging ids
fid = fopen(ids_file);
im_ids = textscan(fid,'%s');
fclose(fid);
pheno_tab.MRID = im_ids{1}; % add ids as 'MRID' column

% read covariates
opts = detectImportOptions(covar_file);
opts.VariableTypes{1} = 'char'; % read FID column as char
opts.VariableTypes{2} = 'char'; % read IID column as char
opts.VariableTypes{3} = 'char'; % read MRID column as char
if contains(pheno_file,'thickness','IgnoreCase', true)
    meanCovar2use = 'MeanThickness';
else
    meanCovar2use = 'MeanArea';
end
% select columns to read from the covar_file
id_cols = {'FID','IID','MRID'};
covar_cols = {'Age','Sex','Scanner','C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14','C15','C16','C17','C18','C19','C20','Euler',meanCovar2use};
opts.SelectedVariableNames = [id_cols covar_cols];
covar_tab = readtable(covar_file, opts);
% set categorical covars
covar_tab.Sex = categorical(covar_tab.Sex);
covar_tab.Scanner = categorical(covar_tab.Scanner);

% join covar and pheno on MRID
merged_tab = join(covar_tab, pheno_tab, 'Key', 'MRID');
[n_pheno_samples, n_pheno] = size(merged_tab);
fprintf('%d sampels (raws), %d columns in covar merged with pheno table\n', n_pheno_samples, n_pheno);

modelspec_covars = strjoin(covar_cols,' + ');
fprintf('Covariates model:\n%s\n', modelspec_covars);
pheno_resid = zeros(n_pheno_samples, n_nonconst_pheno);
fprintf('Residualizing %d phenotypes (columns) for %d samples (raws)\n', n_nonconst_pheno, n_pheno_samples);
for i=1:n_nonconst_pheno
    pheno_id = pheno_tab.Properties.VariableNames{i};
    modelspec = sprintf('%s ~ %s', pheno_id, modelspec_covars);
    mdl = fitlm(merged_tab, modelspec);
    pheno_resid(:,i) = mdl.Residuals.Raw;
    if mod(i,50) == 0
        fprintf('%d out of %d phenotypes processed\n', i, n_nonconst_pheno);
    end
end

