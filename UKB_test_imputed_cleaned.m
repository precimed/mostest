% download and compile https://se.mathworks.com/matlabcentral/fileexchange/27076-shuffle  (just once)
mex Shuffle.c 

snps = 10408776; nsubj = 26502;
df = readtable('/home/oleksanf/vmshare/data/UKBDATA/phenotypes/FS/UKB26502_GWAS_qnorm_ecdf_ordered_230519.txt');
fname = '/home/oleksanf/vmshare/data/UKBDATA/genotypes/imputed/UKB26502_QCed_230519';
ymat = table2array(df);  %d.Properties.VariableNames - column names
ymat_orig = ymat;
rois = df.Properties.VariableNames;

% in the data dictionary, preppend 'x' to 3rd_Ventricle, 4th_Ventricle, 5th_Ventricle, because that's what MATLAB does to column names in df table
df_dict = readtable('/home/oleksanf/vmshare/data/UKBDATA/phenotypes/FS/FSdict_ver53_230519.txt', 'Delimiter', 'tab');
for i=1:length(df_dict.ROI_name)
  if ismember(df_dict.ROI_name{i}(1), '0123456789'), df_dict.ROI_name{i} = sprintf('x%s', df_dict.ROI_name{i}); end
end

% Sparse binary matrix for random pruning, LD r2 of .1
load('/home/oleksanf/vmshare/data/UKBDATA/genotypes/imputed/UKB26502_QCed_230519_r2p1.mat'); % LDmat, chrnumvec, posvec, mafvec
defvec = randn(size(LDmat, 1), 1); defvec(mafvec < 0.05) = nan; defvec = FastPrune(defvec, LDmat);
snp_idx = find(isfinite(defvec));

mask = false(length(rois), 1);
table = struct('is_DK_AREA', mask, 'is_DK_VOLUME', mask, 'is_DK_THICKNESS', mask, 'is_ASEG', mask);
table.measure = {};
for i=1:length(rois), 
  subset = df_dict(strcmp(df_dict.ROI_name, rois{i}), :);
  if isempty(subset), error('%s\n', rois{i}); end;
  table.is_DK_AREA(i, 1)      = strcmp('CorticalArea', subset.ROI_type{1}); 
  table.is_DK_VOLUME(i, 1)    = strcmp('CorticalVolume', subset.ROI_type{1}); 
  table.is_DK_THICKNESS(i, 1) = strcmp('CorticalThickness', subset.ROI_type{1}); 
  table.is_ASEG(i, 1)         = strcmp('SubcorticalVolume', subset.ROI_type{1}); 
  table.measure{end+1, 1} = rois{i};
end

data = {};

if 0
  npheno=size(ymat, 2);
  zmat=zeros(length(snp_idx), npheno, 2, 'single'); 

  lamvec = 0:0.2:1; 
  for lami=1:length(lamvec)

  lam_reg = lamvec(lami);
  data{lami}.lam_reg = lam_reg;

  chunk = 300000;
  for i=1:chunk:length(snp_idx)
    j=min(i+chunk-1, length(snp_idx));
    fprintf('loading snps %i to %i... ', i, j);    tic;
    %geno_int8 = PlinkRead_binary2(nsubj, snp_idx(i:j), fname);
    %fprintf('processing... ', i, j);   
    %geno = nan(size(geno_int8), 'single'); for code = int8([0,1,2]), geno(geno_int8==code) = single(code); end;

    if 0
      C = corr(ymat_orig);
    
      C_reg = (1-lam_reg)*C + lam_reg*diag(max(0.01,diag(C))); % Ridge regularized covariance matrix
      C_inv = inv(C_reg);
      W_wht = chol(C_inv); % Whitening filter
      ymat = ymat_orig*W_wht'; % Whitened residualized data

    end

    tic
    [~, zmat_orig] = nancorr(ymat, geno);
    [~, zmat_perm] = nancorr(ymat, Shuffle(geno));
    zmat(i:j, :, 1) = zmat_orig';
    zmat(i:j, :, 2) = zmat_perm';
    fprintf('done in %.1f sec, %.1f %% completed\n', toc, 100*(j+1)/length(snp_idx));
  end
  data{lami}.zmat = zmat;
  data{lami}.ymat = ymat;

  %zmat_bak = zmat;
  end

  save('gwas2.UKB26502_GWAS_qnorm_ecdf_ordered_230519.mat', '-v7.3', 'zmat');
else
  load('gwas.UKB26502_GWAS_qnorm_ecdf_ordered_230519.mat');
end


ivec_list = {{'dk_area',       find(table.is_DK_AREA)}, ...
             {'dk_thick',      find(table.is_DK_THICKNESS)}, ...
             {'dk_vol',        find(table.is_DK_VOLUME)}, ...
             {'dk_area_thick', find(table.is_DK_AREA | table.is_DK_THICKNESS)}, ...
             {'dk_all',        find(table.is_DK_AREA | table.is_DK_THICKNESS | table.is_DK_VOLUME)}, ...
             {'aseg_vol',      find(table.is_ASEG)}, ...
             {'all',           find(table.is_DK_AREA | table.is_DK_THICKNESS | table.is_ASEG)}};
%ivec_list = {{'aseg_vol',      find(table.is_ASEG)}};
for iveci = 1:length(ivec_list)
  ivec = ivec_list{iveci}{2}; ivec_name = ivec_list{iveci}{1};
  fprintf('%s: ', ivec_name);
  for j=1:min(length(ivec), 4)
      fprintf('%s, ', table.measure{ivec(j)});
  end    
  fprintf('... (%i in total)\n', length(ivec));
end

fh1 = sfigure(1001); clf; 
fh2 = sfigure(1002); clf; 
fh3 = sfigure(1003); clf; 
for iveci = 1:length(ivec_list)

  ivec = ivec_list{iveci}{2}; ivec_name = ivec_list{iveci}{1};
  zmat = zmat_bak(:, ivec, :);
  % logpmat = logpmat_bak(:, ivec, :);

  ivec_snp_good = all(all(isfinite(zmat), 2), 3);
  
  C0 = corr(zmat(ivec_snp_good, :, 2));
  C1 = corr(zmat(ivec_snp_good, :, 1)); % & Hvec>0.1 & CRvec>0.95 & max(abs(zmat(:,:,1)),[],1)>abs(norminv(1e-5))),1)');
  sfigure(fh1); 
  subplot(2,length(ivec_list),0*length(ivec_list)+iveci); imagesc(C0,1*[-1 1]); colormap(blueblackred); axis tight equal; title(ivec_name); % Check if this is identical to cov(ymat_resid)
  if iveci == 1, ylabel('perm'); end;
  subplot(2,length(ivec_list),1*length(ivec_list)+iveci); imagesc(C1,1*[-1 1]); colormap(blueblackred); axis tight equal; % Check if this is identical to cov(ymat_resid)
  if iveci == 1, ylabel('true'); end;
  
  [U S]  = svd(C0); s = diag(S);
%  C0_reg = diag(diag(C0)); % Complete regularization -- results in imperfect gamma fit
%  C0_reg = eye(size(C0)); % Complete regularization -- results in imperfect gamma fit
%  max_lambda = s(min(10, length(s)));
  max_lambda = min(0.1, s(min(10, length(s)))); % oleksanf: don't regularize unless it's too bad

  C0_reg = U*diag(max(max_lambda,s))*U'; % Good gamma fit
%  C0_reg = U*diag(max(s(40),s))*U';

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

  sfigure(fh2); 
  subplot(2,length(ivec_list),0*length(ivec_list)+iveci); plot(hv_maxlogpvecs,-log10(1-chc_maxlogpvecs),hv_maxlogpvecs,-log10(1-cdf(pd_minpvecs,10.^-hv_maxlogpvecs,'upper')), '.', 'LineWidth',2);title(ivec_name); 
  if iveci == 1, ylabel('minP'); end;
  subplot(2,length(ivec_list),1*length(ivec_list)+iveci); plot(hv_logpdfvecs,-log10(1-chc_logpdfvecs),hv_logpdfvecs,-log10(1-pd_logpdfvecs.cdf(hv_logpdfvecs)),'.', 'LineWidth',2)
  if iveci == 1, ylabel('MOST'); end;

% Plot FWER-corrected GWAS stats

  maxlogpvecs_corr = -log10(cdf(pd_minpvecs,minpvecs));
  logpdfvecs_corr = -log10(cdf(pd_logpdfvecs,logpdfvecs,'upper'));

  sfigure(fh3); 
  subplot(4,length(ivec_list),0*length(ivec_list)+iveci); plot(maxlogpvecs_corr(1,ivec_snp_good),'*'); axis tight;title(ivec_name); 
  if iveci == 1, ylabel('minP, true'); end;
  subplot(4,length(ivec_list),1*length(ivec_list)+iveci); plot(maxlogpvecs_corr(2,ivec_snp_good),'*'); axis tight;
  if iveci == 1, ylabel('minP, perm'); end;
  subplot(4,length(ivec_list),2*length(ivec_list)+iveci); plot(logpdfvecs_corr(1,ivec_snp_good),'*'); axis tight;
  if iveci == 1, ylabel('MOST, true'); end;
  subplot(4,length(ivec_list),3*length(ivec_list)+iveci); plot(logpdfvecs_corr(2,ivec_snp_good),'*'); axis tight;
  if iveci == 1, ylabel('MOST, perm'); end;

%  Compare power of logpdfvecs vs max(logpmat) -- much greater power!!
  is_low_maf = (mafvec < 0.05) | (mafvec > 0.95);
  minp_prunevec = maxlogpvecs_corr(1,ivec_snp_good); minp_prunevec(is_low_maf) = nan; minp_prunevec(minp_prunevec < -log10(5e-8)) = nan; 
  most_prunevec = logpdfvecs_corr(1,ivec_snp_good); most_prunevec(is_low_maf) = nan; most_prunevec(most_prunevec < -log10(5e-8)) = nan; 
  
  fprintf('\n');
  fprintf(1,'%s: GWAS yield                       minP: %d; MOST: %d\n',ivec_name,sum(maxlogpvecs_corr(1,ivec_snp_good)>-log10(5e-8)),sum(logpdfvecs_corr(1,ivec_snp_good)>-log10(5e-8)));
  fprintf(1,'%s: GWAS yield, maf > 0.05           minP: %d; MOST: %d\n',ivec_name,sum(isfinite(minp_prunevec)),sum(isfinite(most_prunevec)));

  minp_prunevec = FastPrune(minp_prunevec, LDmat); most_prunevec = FastPrune(most_prunevec, LDmat);
  fprintf(1,'%s: GWAS yield, maf > 0.05, prunning minP: %d; MOST: %d\n',ivec_name,sum(isfinite(minp_prunevec)),sum(isfinite(most_prunevec)));
  fprintf(1,'%s: #measures: %i; #eff(minP): %.1f; #eff(MOST)=%.1f;\n',ivec_name, length(ivec), pd_minpvecs.b, pd_logpdfvecs.a);

  minPval = maxlogpvecs_corr(1, :);
  mostPval = logpdfvecs_corr(1, :);
  save(sprintf('MOST_minP_%s.mat', ivec_name), '-v7', 'mostPval', 'minPval');
  drawnow;
end

saveas(fh1,'fh1.pdf')
saveas(fh2,'fh2.pdf')
saveas(fh3,'fh3.pdf')
