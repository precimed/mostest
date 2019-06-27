% addpath(genpath('/home/cfan/codes/matlab')); % get Chun's /home/cfan/codes/matlab/PlinkRead_binary.m

% Freesurfer ROI spreadsheet:i
% /space/gwas-syn1/1/data/GWAS/UKBioBank/Phenotypes/NORMENT_T1_FS.ROI.RData

% Raw genotype data:
% /space/gwas-ds2/1/data/GWAS/UKBio/genotype_calls/ukb_cal_chr*_v2.bed
% /space/gwas-ds2/1/data/GWAS/UKBio/genotype_calls/ukb2741_cal_chr*_v2_s488366.fam
% /space/gwas-ds2/1/data/GWAS/UKBio/genotype_calls/ukb_snp_chr*_v2.bim

% Imputed genotype data:
% /space/gwas-syn1/1/data/GWAS/UKBioBank/genotypes/imputed_plink_app27412

% Look in ~/Dropbox/data/UKB for
%   KeepCauc.txt : IDs of Caukacians
%   UK500k_BasicCovs.txt : standard covariates

data_folder = 'Z:\vmshare\data\UKB';
%data_folder = '/space/syn03/1/data/oleksandr/UKB'p

data_covars = readtext_amd(fullfile(data_folder, 'UK500k_BasicCovs.txt'),char(9));
IDvec_covars = data_covars(2:end,1);
data_covars(strcmp('NA',data_covars(:))) = {'NaN'};
agevec_tmp = cell2mat_amd(data_covars(2:end,2));
sexvec_tmp = strcmp('Male',data_covars(2:end,3));
scannervec_tmp = data_covars(2:end,4); 
scannerlist = setdiff(unique(scannervec_tmp),'NaN');
scannermat_tmp = zeros(length(agevec_tmp),length(scannerlist)-1);
for sci = 2:length(scannerlist)
  scannermat_tmp(strcmp(scannerlist{sci},scannervec_tmp)) = 1;
end
scannermat_tmp(~ismember(scannervec_tmp,scannerlist),:) = NaN;
datamat_covars = cat(2,agevec_tmp,sexvec_tmp,scannermat_tmp,cell2mat_amd(data_covars(2:end,5:end)));

data_keep = readtext_amd(fullfile(data_folder, 'KeepCauc.txt'),char(9));
IDvec_keep = data_keep(:,1);

% Limit to Caucasian subjects -- should perhap;s avoid re-using variables?
ivec_tmp = ismember(IDvec_covars,IDvec_keep);
IDvec_covars = IDvec_covars(ivec_tmp);
datamat_covars = datamat_covars(ivec_tmp,:);

data = readtext_amd(fullfile(data_folder, 'NORMENT_T1_FS.ROI.csv'));
colnames = StripQuotes_amd(data(1,2:end));
IDvec_FS = strrep(data(2:end,1),'"','');
data = StripQuotes_amd(data(2:end,2:end));

chrlist = cat(2,cellfun(@(x)num2str(x),num2cell([1:22]),'UniformOutput',false),'X','Y','XY','MT');

% Should load in datamat_covars, etc. from file

genomat_cat = []; chrvec_cat = []; bpvec_cat = []; snplist_cat = {};
for chri = 1:length(chrlist)
  chr = chrlist{chri};
  fprintf(1,'%d: chr%s (now=%s)\n',chri,chr,datestr(now));
  load(fullfile(data_folder, sprintf('UKB_snap%02d.mat',chri)));
  [dummy IA IB] = intersect(IID,IDvec_covars,'stable');
  genomat_cat = cat(2,genomat_cat,genomat(IA,:));
  chrvec_cat = cat(1,chrvec_cat,chrvec0);
  bpvec_cat = cat(1,bpvec_cat,bpvec0);
  snplist_cat = cat(1,snplist_cat,snplist);
end
datamat_covars_merged = datamat_covars(IB,:); IDvec = IDvec_covars(IB);

nsnp = size(genomat_cat,2);
Hvec = NaN(1,nsnp); CRvec = NaN(1,nsnp);
for snpi = 1:nsnp
  if mod(snpi,1000)==0
    fprintf(1,'snpi=%d of %d (now=%s)\n',snpi,nsnp,datestr(now));
  end
  gvec = double(genomat_cat(:,snpi));
  gvec(gvec<0) = NaN;
  Hvec(snpi) = var(gvec(isfinite(gvec))); % Note that genotype variance (H) can be as high as 1! -- should compute 2*maf*(1-maf) instead?
  CRvec(snpi) = mean(isfinite(gvec));
end

[dummy IA2 IB2] = intersect(IDvec,IDvec_FS,'stable');
datamat = cell2mat_amd(data(IB2,:));

ivec_pheno = var(datamat,[],1)>0; % All measures with non-zero variance
%ivec_pheno = find_cell(cellfun(@(x)regexp(x,'.vol$|^ICV$','match'),colnames,'UniformOutput',false))&var(datamat,[],1)>0; % Subcortical volumes
%ivec_pheno = find_cell(cellfun(@(x)regexp(x,'.area.|^area','match'),colnames,'UniformOutput',false))&var(datamat,[],1)>0; % Cortical area
%ivec_pheno = find_cell(cellfun(@(x)regexp(x,'.thick.|^thick','match'),colnames,'UniformOutput',false))&var(datamat,[],1)>0; % Cortical thickness
%ivec_pheno = find_cell(cellfun(@(x)regexp(x,'.vol.','match'),colnames,'UniformOutput',false)); % Cortical volumes

% residualize phenotypes on covariates
ymat = datamat(:,ivec_pheno); phenonames = colnames(ivec_pheno);
X_tmp = cat(2,datamat_covars_merged,ones(size(ymat,1),1));
defvec =  isfinite(sum(ymat,2)+sum(X_tmp,2));
ymat_resid = ymat - X_tmp*(X_tmp(defvec,:)\ymat(defvec,:));
ymat_resid = ymat_resid ./ repmat(max(0.01,std(ymat_resid,[],1)),[size(ymat_resid,1) 1]);

ivec_area = find(find_cell(regexp(phenonames,'.area.','match'))); 
ivec_thick = find(find_cell(regexp(phenonames,'.thick.','match'))); 
ivec_cortvol = find(find_cell(regexp(phenonames,'.vol.','match'))); 
ivec_vol = find(find_cell(regexp(phenonames,'.vol$','match'))); 
ivec_all = [ivec_area,ivec_thick,ivec_cortvol,ivec_vol]; 

% residualize phenotypes on corresponding global measures
valvec_area_norm = (ymat_resid(:,strcmp('area.left.total',phenonames))+ymat_resid(:,strcmp('area.right.total',phenonames)));
valvec_thick_norm = (ymat_resid(:,strcmp('thick.left.mean',phenonames))+ymat_resid(:,strcmp('thick.right.mean',phenonames)))/2;
valvec_cortvol_norm = sum(ymat_resid(:,ivec_cortvol),2);
valvec_vol_norm = sum(ymat_resid(:,ivec_vol),2);
ymat_resid(:,ivec_area) = ymat_resid(:,ivec_area) - valvec_area_norm*(valvec_area_norm\ymat_resid(:,ivec_area));
ymat_resid(:,ivec_thick) = ymat_resid(:,ivec_thick) - valvec_thick_norm*(valvec_thick_norm\ymat_resid(:,ivec_thick));
ymat_resid(:,ivec_cortvol) = ymat_resid(:,ivec_cortvol) - valvec_cortvol_norm*(valvec_cortvol_norm\ymat_resid(:,ivec_cortvol));
ymat_resid(:,ivec_vol) = ymat_resid(:,ivec_vol) - valvec_vol_norm*(valvec_vol_norm\ymat_resid(:,ivec_vol));

% transform measures into gaussian shape. Filter out measures that has too
% large kurtosis after transformation (indicate large probability mass at certain points)
% => won't work for case/control phenotypes?
% ??? what if there is inflation in summary statistics? Should we
% incorporate intergenic GC ? Or is it taken care by permutation?
ymat_resid_trans = NaN(size(ymat_resid));
nroi = size(ymat,2);
hv = linspace(-10,10,201); 
for roii = 1:nroi
  yvec = ymat_resid(:,roii);  
  hc = hist(yvec,hv); pdfhat = hc/sum(hc)/(hv(2)-hv(1));
  pd = fitdist(yvec,'kernel','Kernel','normal');
  yvec_trans = norminv(pd.cdf(yvec));
  hc_trans = hist(yvec_trans,hv); pdfhat_trans = hc_trans/sum(hc_trans)/(hv(2)-hv(1));
  ymat_resid_trans(:,roii) = yvec_trans;
  fprintf(1,'roii=%d (%s): kurtosis= %f / %f\n',roii,phenonames{roii},kurtosis(yvec),kurtosis(yvec_trans));
  if abs(kurtosis(yvec_trans)-3)>1
    sfigure(1000+roii); plot(hv,pdfhat,hv,pd.pdf(hv),'LineWidth',3);
    sfigure(2000+roii); plot(hv,pdfhat_trans,hv,normpdf(hv),'LineWidth',3);
    drawnow;
  end
end

ymat_resid = ymat_resid_trans;
 
load(fullfile(data_folder , 'UKB_zmat_oleksanf.mat'));
 
if 0
    
% What does Voxelwise_GLM?
X = cat(2,NaN(size(ymat,1),1),ones(size(ymat,1),1)); contrastvec = [1 0];
nsnp = size(genomat_cat,2); npheno = sum(ivec_pheno);
bmat = NaN(npheno,nsnp,2,'single'); zmat = NaN(npheno,nsnp,2,'single'); logpmat = NaN(npheno,nsnp,2,'single');
for snpi = 1:nsnp
  gvec = double(genomat_cat(:,snpi));
  gvec(gvec<0) = NaN;
  X(:,1) = gvec;
  defvec = isfinite(sum(X,2)+sum(ymat_resid,2)); ivec_def = find(defvec); 
  for permflag = 0:1
    if permflag
      ivec_perm = randperm(length(ivec_def));
      % what's the second component in betavols?
      % important: permute all phenotypes in a consistent fasion, to keep
      % correlation structure among phenotypes.
      % here permutation is applied to ymat_resid, not to the genotype matrix; that's because we've already residualized the phenotypes?
      [pvols,betavols,beta_cov,sigvol] = Voxelwise_GLM(X(ivec_def,:),ymat_resid(ivec_perm,:),contrastvec); betavec = betavols(1,:); zvec = betavec./(sqrt(beta_cov(1,1))*sigvol);
    else
      [pvols,betavols,beta_cov,sigvol] = Voxelwise_GLM(X(ivec_def,:),ymat_resid(ivec_def,:),contrastvec); betavec = betavols(1,:); zvec = betavec./(sqrt(beta_cov(1,1))*sigvol);
    end
    bmat(:,snpi,permflag+1) = betavec; zmat(:,snpi,permflag+1) = zvec; logpmat(:,snpi,permflag+1) = -log10(normcdf(-abs(zvec))*2);
  end
  if mod(snpi,1000)==0
    fprintf(1,'snpi=%d of %d (now=%s)\n',snpi,nsnp,datestr(now));
    tmp = zmat(:,1:snpi,:); tmp(~isfinite(tmp)) = 0;
    ivec_snp = find(max(max(abs(tmp),[],3),[],1)>abs(norminv(1e-6)));
    tmp = tmp(:,ivec_snp,:);
    if ~isempty(tmp)
%      [mv mi] = max(abs(tmp),[],1);
%      for j = 1:length(mi), tmp(:,j) = tmp(:,j)*sign(tmp(mi(j),j)); end
      for j = 1:size(tmp,2)
        v = tmp(:,j,1); [mv mi] = max(abs(v)); tmp(:,j,1) = v/sign(v(mi));
        v = tmp(:,j,2); [mv mi] = max(abs(v)); tmp(:,j,2) = v/sign(v(mi));
      end
      sfigure(101); clf; imagesc(tmp(:,:,1),max(abs(colvec(tmp(:,:,1))))*[-1 1]); colormap(blueblackred); colorbar;
      sfigure(102); clf; imagesc(tmp(:,:,2),max(abs(colvec(tmp(:,:,2))))*[-1 1]); colormap(blueblackred); colorbar;
      sfigure(201); clf; plot(max(logpmat(:,1:snpi,1),[],1),'*'); axis tight;
      sfigure(202); clf; plot(max(logpmat(:,1:snpi,2),[],1),'*'); axis tight;
      drawnow;
    end
  end
end

%    fname_snap = '~/UKB/UKB_zmat.mat';
%   save(fullfile(data_folder ,'UKB_zmat_oleksanf.mat'),'zmat','bmat','logpmat'); 
end % if 0

ivec_snp_good = find(isfinite(sum(zmat(:,:,2),1)) & Hvec>0.1 & CRvec>0.95);

zmat_bak = zmat; logpmat_bak = logpmat; ymat_resid_bak = ymat_resid;

if 0  % unused?
% Compute covariance
C = cov(ymat_resid); % Should residualize for all covariates!
sfigure(1); clf; imagesc(C,[-1 1]); colormap(blueblackred); axis equal tight;
lam_reg = 1.0; % No whitening
%lam_reg = 0.1;
%lam_reg = 0.5;
C_reg = (1-lam_reg)*C + lam_reg*diag(max(0.01,diag(C))); % Ridge regularized covariance matrix
C_inv = inv(C_reg);
%W_wht = chol(C_inv);
%ymat_resid_wht = ymat_resid*W_wht'; % Whitened data
%ymat_resid_wht(abs(ymat_resid_wht)>5) = 0; % Censor outliers -- should set to NaN instead?
end

ivec_list = {{'area' ivec_area} {'thick' ivec_thick} {'cortvol' ivec_cortvol} {'vol' ivec_vol} {'all' ivec_all}};

for iveci = 1:length(ivec_list)
  ivec = ivec_list{iveci}{2}; ivec_name = ivec_list{iveci}{1};
  fprintf('%s: ', ivec_name);
  for j=1:min(length(ivec), 3)
      fprintf('%s, ', colnames{ivec(j)});
  end    
  fprintf('and %i others\n', length(ivec)-j);
end
    
fh1 = sfigure(1001); clf; 
fh2 = sfigure(1002); clf; 
fh3 = sfigure(1003); clf; 
for iveci = 1:length(ivec_list)

  ivec = ivec_list{iveci}{2}; ivec_name = ivec_list{iveci}{1};
  zmat = zmat_bak(ivec,:,:); logpmat = logpmat_bak(ivec,:,:); ymat_resid = ymat_resid_bak(:,ivec);

  C0 = corr(zmat(:,ivec_snp_good,2)');
  C1 = corr(zmat(:,find(isfinite(sum(zmat(:,:,2),1)) & Hvec>0.1 & CRvec>0.95 & max(abs(zmat(:,:,1)),[],1)>abs(norminv(1e-5))),1)');
  sfigure(fh1); 
  subplot(2,length(ivec_list),0*length(ivec_list)+iveci); imagesc(C0,1*[-1 1]); colormap(blueblackred); axis tight equal; title(ivec_name); % Check if this is identical to cov(ymat_resid)
  if iveci == 1, ylabel('perm'); end;
  subplot(2,length(ivec_list),1*length(ivec_list)+iveci); imagesc(C1,1*[-1 1]); colormap(blueblackred); axis tight equal; % Check if this is identical to cov(ymat_resid)
  if iveci == 1, ylabel('true'); end;
  
  [U S]  = svd(C0); s = diag(S);
%  C0_reg = diag(diag(C0)); % Complete regularization -- results in imperfect gamma fit
%  C0_reg = eye(size(C0)); % Complete regularization -- results in imperfect gamma fit
  C0_reg = U*diag(max(s(10),s))*U'; % Good gamma fit
%  C0_reg = U*diag(max(s(40),s))*U';

  logpdfvecs = NaN(size(zmat,3),size(zmat,2)); maxlogpvecs = NaN(size(zmat,3),size(zmat,2));
  for i  = 1:size(zmat,3)
    logpdfvecs(i,:) = -(mvnpdfln(zmat(:,:,i)',0,C0_reg)-mvnpdfln(zeros(size(C0,1),1),0,C0_reg))/log(10); % Should be using mvncdf instead?
    maxlogpvecs(i,:) = max(logpmat(:,:,i),[],1);
  end
  minpvecs = 10.^-maxlogpvecs;

  [hc_maxlogpvecs hv_maxlogpvecs] = hist(maxlogpvecs(2,ivec_snp_good),1000); chc_maxlogpvecs = cumsum(hc_maxlogpvecs)/sum(hc_maxlogpvecs);
  pd_minpvecs = fitdist(colvec(minpvecs(2,ivec_snp_good)),'beta'); % Not a great fit
%  pd_minpvecs.a = 1; % Hard-code known parameter (see http://www.di.fc.ul.pt/~jpn/r/prob/range.html)
  [hc_logpdfvecs hv_logpdfvecs] = hist(logpdfvecs(2,ivec_snp_good),1000); chc_logpdfvecs = cumsum(hc_logpdfvecs)/sum(hc_logpdfvecs);
  pd_logpdfvecs = fitdist(colvec(logpdfvecs(2,ivec_snp_good)),'gamma'); % Seems to work -- beta and wbl  do not

  sfigure(fh2); 
  subplot(2,length(ivec_list),0*length(ivec_list)+iveci); plot(hv_maxlogpvecs,-log10(1-chc_maxlogpvecs),hv_maxlogpvecs,-log10(1-cdf(pd_minpvecs,10.^-hv_maxlogpvecs,'upper')),'LineWidth',2);title(ivec_name); 
  if iveci == 1, ylabel('minP'); end;
  subplot(2,length(ivec_list),1*length(ivec_list)+iveci); plot(hv_logpdfvecs,-log10(1-chc_logpdfvecs),hv_logpdfvecs,-log10(1-pd_logpdfvecs.cdf(hv_logpdfvecs)),'LineWidth',2)
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
  fprintf(1,'%s: GWAS yield minP: %d; MOST: %d\n',ivec_name,sum(maxlogpvecs_corr(1,ivec_snp_good)>-log10(5e-8)),sum(logpdfvecs_corr(1,ivec_snp_good)>-log10(5e-8)));
  fprintf(1,'%s: #measures: %i; #eff(minP): %.1f; #eff(MOST)=%.1f;\n',ivec_name, length(ivec), pd_minpvecs.b, pd_logpdfvecs.a);

  drawnow;

end

















%stuff below is not relevant
return

% Estimate components
if 1
  C0 = cov(zmat(:,ivec_snp_good,2)');
  C1 = cov(zmat(:,find(isfinite(sum(zmat(:,:,2),1)) & Hvec>0.1 & CRvec>0.95 & max(abs(zmat(:,:,1)),[],1)>abs(norminv(1e-5))),1)');

  ncomp_svd = 9;
  [U_svd S_svd]  = svd(C1);
  sfigure(1000); imagesc(U_svd(:,1:ncomp_svd),max(abs(U_svd(:)))*0.5*[-1 1]); colormap(blueblackred);

  ncomp_geig = 9;
  [U_geig S_geig]  = eig(C1,C0);
  sfigure(1001); imagesc(U_geig(:,1:ncomp_geig),max(abs(U_geig(:)))*0.5*[-1 1]); colormap(blueblackred);

%  ncomp = ncomp_svd; U = U_svd(:,1:ncomp);
  ncomp = ncomp_geig; U = U_geig(:,1:ncomp);

  bvecs = NaN(ncomp,nsnp,2); zvecs = NaN(ncomp,nsnp,2); logpvecs = NaN(ncomp,nsnp,2);
  fh1=sfigure(11201); clf; fh2=sfigure(11202); clf;
  nr = round(sqrt(ncomp)); nc = ceil(ncomp/nr);
  C_inv_tmp = inv(C_reg(ivec,ivec));
  for compi = 1:ncomp
    fprintf(1,'compi=%d of %d (now=%s)\n',compi,ncomp,datestr(now));
    compvec = U(:,compi);
    yvec = ymat_resid*compvec;
    for snpi = 1:nsnp
      gvec = double(genomat_cat(:,snpi));
      gvec(gvec<0) = NaN;
      X(:,1) = gvec;
      defvec = isfinite(sum(X,2)+sum(yvec,2)); ivec_def = find(defvec);
      for permflag = 0:1
        if permflag
          ivec_perm = randperm(length(ivec_def));
          [pvols,betavols,beta_cov,sigvol] = Voxelwise_GLM(X(ivec_def,:),yvec(ivec_perm,:),contrastvec); b = betavols(1,:); se = sqrt(beta_cov(1,1))*sigvol; z = b/se;
        else
          [pvols,betavols,beta_cov,sigvol] = Voxelwise_GLM(X(ivec_def,:),yvec(ivec_def,:),contrastvec); b = betavols(1,:); se = sqrt(beta_cov(1,1))*sigvol; z = b/se;
        end
        bvecs(compi,snpi,permflag+1) = b; zvecs(compi,snpi,permflag+1) = z; logpvecs(compi,snpi,permflag+1) = -log10(normcdf(-abs(z))*2);
        if mod(snpi,1000)==0
          if ~permflag
            sfigure(fh1); subplot(nc,nr,compi); plot(logpvecs(compi,:,1),'*'); axis tight;
          else
            sfigure(fh2); subplot(nc,nr,compi); plot(logpvecs(compi,:,2),'*'); axis tight;
          end
          drawnow;
        end
      end
    end
  end

  save(sprintf('~/UKB/UKB_zvecs_%s.mat',ivec_name),'ncomp','U','bvecs','zvecs','logpvecs');
end

% ToDo
%  Update with 25k subjects (Chun)
%  Implement multivariate F-test detection
%     Apply to subsets of phenotypes
%  Define optimal projections eig(C1,C0) -- where C1 and C0 are based on covariance


% Extras

sfigure(1); clf;
nr = round(sqrt(nroi)); nc = ceil(nroi/nr);
for roii = 1:nroi
  fprintf(1,'roii=%d of %d (now=%s)\n',roii,nroi,datestr(now));
  [logqmat hv_logp] = GenStats_QQ_plot(squeeze(logpmat(roii,ivec_snp_good,:)));
  sfigure(1); subplot(nr,nc,roii); plot([0 7],[0 7],'k:',logqmat,hv_logp,'LineWidth',2); xlim([0 7]); ylim([0 10]);
%  drawnow
end

custpdf = @(data,Neff) Neff*(1-data).^(Neff-1);
custcdf = @(data,Neff) (1-data).^(Neff);
phat = mle(10.^-maxlogpvecs(2,ivec_snp),'pdf',custpdf,'cdf',custcdf,'start',size(C0,1));

%pd = fitdist(colvec(logpdfvecs(2,ivec_snp)),'kernel','Kernel','normal');

% Compute covariance
C = cov(ymat_resid); % Should residualize for all covariates!
sfigure(1); clf; imagesc(C,[-1 1]); colormap(blueblackred); axis equal tight;
lam_reg = 0.1;
%lam_reg = 0.5;
C_reg = (1-lam_reg)*C + lam_reg*diag(max(0.01,diag(C))); % Ridge regularized covariance matrix
C_inv = inv(C_reg);
W_wht = chol(C_inv);
ymat_resid_wht = ymat_resid*W_wht'; % Whitened data
ymat_resid_wht(abs(ymat_resid_wht)>5) = 0; % Censor outliers -- should set to NaN instead?


X = cat(2,NaN(size(ymat,1),1),ones(size(ymat,1),1)); contrastvec = [1 0];
nsnp = size(genomat_cat,2); npheno = sum(ivec_pheno);
bmat = NaN(npheno,nsnp,2,'single'); zmat = NaN(npheno,nsnp,2,'single'); logpmat = NaN(npheno,nsnp,2,'single');
for snpi = 1:nsnp
  gvec = double(genomat_cat(:,snpi));
  gvec(gvec<0) = NaN;
  X(:,1) = gvec;
  for whitenflag = 0:1
    defvec = isfinite(sum(X,2)+sum(ymat_resid,2));
    if whitenflag
      [pvols,betavols,beta_cov,sigvol] = Voxelwise_GLM(X(defvec,:),ymat_resid_wht(defvec,:),contrastvec); betavec = betavols(1,:); zvec = betavec./(sqrt(beta_cov(1,1))*sigvol);
    else
      [pvols,betavols,beta_cov,sigvol] = Voxelwise_GLM(X(defvec,:),ymat_resid(defvec,:),contrastvec); betavec = betavols(1,:); zvec = betavec./(sqrt(beta_cov(1,1))*sigvol);
    end
    bmat(:,snpi,whitenflag+1) = betavec; zmat(:,snpi,whitenflag+1) = zvec; logpmat(:,snpi,whitenflag+1) = -log10(normcdf(-abs(zvec))*2);
  end
  if mod(snpi,1000)==0
    fprintf(1,'snpi=%d of %d (now=%s)\n',snpi,nsnp,datestr(now));
    tmp = zmat(:,1:snpi,:); tmp(~isfinite(tmp)) = 0; 
    ivec_snp = find(max(max(abs(tmp),[],3),[],1)>abs(norminv(1e-6)));
    tmp = tmp(:,ivec_snp,:);
    if ~isempty(tmp)
%      [mv mi] = max(abs(tmp),[],1); 
%      for j = 1:length(mi), tmp(:,j) = tmp(:,j)*sign(tmp(mi(j),j)); end
      for j = 1:size(tmp,2)
        v = tmp(:,j,1); [mv mi] = max(abs(v)); tmp(:,j,1) = v/sign(v(mi));
        v = tmp(:,j,2); [mv mi] = max(abs(v)); tmp(:,j,2) = v/sign(v(mi));
      end
      sfigure(101); clf; imagesc(tmp(:,:,1),max(abs(colvec(tmp(:,:,1))))*[-1 1]); colormap(blueblackred); colorbar;
      sfigure(102); clf; imagesc(tmp(:,:,2),max(abs(colvec(tmp(:,:,2))))*[-1 1]); colormap(blueblackred); colorbar;
      sfigure(201); clf; plot(max(logpmat(:,1:snpi,1),[],1),'*'); axis tight;
      sfigure(202); clf; plot(max(logpmat(:,1:snpi,2),[],1),'*'); axis tight;
      drawnow;
    end
  end
end

fname_snap = '~/UKB/UKB_zmat.mat';
save(fname_snap,'zmat','bmat','logpmat');

sfigure(1); clf;
nr = round(sqrt(nroi)); nc = ceil(nroi/nr);
for roii = 1:nroi
  fprintf(1,'roii=%d of %d (now=%-35s)\n',roii,nroi,datestr(now));
  [logqmat hv_logp] = GenStats_QQ_plot(squeeze(logpmat(roii,:,:)));
  sfigure(1); subplot(nr,nc,roii); plot([0 7],[0 7],'k:',logqmat,hv_logp,'LineWidth',2); xlim([0 7]); ylim([0 10]);
%  drawnow
end


% ToDo
%   Identify SNPs with the most significant bvec or zvec, before and after whitening (use respective phenotypic covariance matrices, mapped to respecive estimation spaces))

%   Perform residualization for all local / regional measures

% Estimate componnets

ivec = ivec_area;
%ivec = ivec_thick;
%ivec = ivec_cortvol;
%ivec = ivec_vol;

tmp = zmat(ivec,:,:); tmp(~isfinite(tmp)) = 0;
ivec_snp = find(max(max(abs(tmp),[],3),[],1)>abs(norminv(1e-5)));
tmp = tmp(:,ivec_snp,:);
for j = 1:size(tmp,2)
  v = tmp(:,j,1); [mv mi] = max(abs(v)); tmp(:,j,1) = v/v(mi);
  v = tmp(:,j,2); [mv mi] = max(abs(v)); tmp(:,j,2) = v/v(mi);
end

sfigure(1101); clf;
subplot(2,2,1); imagesc(tmp(:,:,1),max(abs(colvec(tmp(:,:,1))))*[-1 1]); 
subplot(2,2,2); imagesc(tmp(:,:,2),max(abs(colvec(tmp(:,:,2))))*[-1 1]);
subplot(2,2,3); imagesc(corr(tmp(:,:,1)'),[-1 1]); axis equal tight;
subplot(2,2,4); imagesc(corr(tmp(:,:,2)'),[-1 1]); axis equal tight;


ncomp_svd = 4;
[U_svd1 S_svd1]  = svd(tmp(:,:,1),'econ');
sfigure(1014); imagesc(U_svd1(:,1:ncomp_svd),max(abs(U_svd1(:)))*0.5*[-1 1]); colormap(blueblackred); 
%sfigure(1015); plot(diag(S_svd1),'*-','lineWidth',2); xlim([1 4*ncomp_svd]);
[U_svd2 S_svd2]  = svd(tmp(:,:,2),'econ');
sfigure(1024); imagesc(U_svd2(:,1:ncomp_svd),max(abs(U_svd2(:)))*0.5*[-1 1]); colormap(blueblackred); 
%sfigure(1025); plot(diag(S_svd2),'*-','lineWidth',2); xlim([1 4*ncomp_svd]);


% k-means clustering
ncomp_kmeans = 9;
[T_kmeans1 U_kmeans1] = kmeans(tmp(:,:,1)',ncomp_kmeans,'Distance','cosine','Replicates',100); U_kmeans1 = U_kmeans1';
sfigure(2014); imagesc(U_kmeans1,max(abs(U_kmeans1(:)))*0.5*[-1 1]); colormap(blueblackred);
[T_kmeans2 U_kmeans2] = kmeans(tmp(:,:,2)',ncomp_kmeans,'Distance','cosine','Replicates',100); U_kmeans2 = U_kmeans2';
sfigure(2024); imagesc(U_kmeans2,max(abs(U_kmeans2(:)))*0.5*[-1 1]); colormap(blueblackred);


ncomp = ncomp_svd; U1 = U_svd1(:,1:ncomp); U2 = U_svd2(:,1:ncomp);
%ncomp = ncomp_kmeans; U1 = U_kmeans1(:,1:ncomp); U2 = U_kmeans2(:,1:ncomp);

bvecs = NaN(ncomp,nsnp,2); zvecs = NaN(ncomp,nsnp,2); logpvecs = NaN(ncomp,nsnp,2);
fh1=sfigure(11201); clf; fh2=sfigure(11202); clf;
nr = round(sqrt(ncomp)); nc = ceil(ncomp/nr);
C_inv_tmp = inv(C_reg(ivec,ivec));
for compi = 1:ncomp
  fprintf(1,'compi=%d of %d (now=%s)\n',compi,ncomp,datestr(now));
%  compvec = U1(:,compi); 
  compvec = U2(:,compi); 
  for whitenflag = [1 0]
    if whitenflag
%      yvec = ymat_resid_wht*U2(:,compi); % Use projections based on whitened stats
%      yvec = ymat_resid_wht*U1(:,compi); % Use projections based on non-whitened stats
%      yvec = ymat_resid*C_inv*U1(:,compi); % Use projections based on non-whitened stats
      yvec = ymat_resid(:,ivec)*C_inv_tmp*compvec;
    else
%      yvec = ymat_resid_wht*U2(:,compi); % Use projections based on whitened stats
%      yvec = ymat_resid*U1(:,compi); % Use projections based on non-whitened stats
      yvec = ymat_resid(:,ivec)*compvec;
    end
    includevec = (abs(yvec-nanmean(yvec))/nanstd(yvec))<4; % Exclude outliers
    for snpi = 1:nsnp
      gvec = double(genomat_cat(:,snpi));
      gvec(gvec<0) = NaN;
      X(:,1) = gvec;
      defvec = isfinite(sum(X,2)+sum(yvec,2))&includevec;
      [pvols,betavols,beta_cov,sigvol] = Voxelwise_GLM(X(defvec,:),yvec(defvec,:),contrastvec); b = betavols(1,:); se = sqrt(beta_cov(1,1))*sigvol; z = b/se;
      bvecs(compi,snpi,whitenflag+1) = b; zvecs(compi,snpi,whitenflag+1) = z; logpvecs(compi,snpi,whitenflag+1) = -log10(normcdf(-abs(z))*2);
      if mod(snpi,1000)==0
        if ~whitenflag
          sfigure(fh1); subplot(nc,nr,compi); plot(logpvecs(compi,:,1),'*'); axis tight;
        else
          sfigure(fh2); subplot(nc,nr,compi); plot(logpvecs(compi,:,2),'*'); axis tight;
        end
        drawnow;
      end
    end
  end
end

fname_snap = '~/UKB/UKB_zvecs.mat';
save(fname_snap,'zvecs','bvecs','logpvecs');


sfigure(666); clf; 


% figure; GenStats_QQ_plot(squeeze(logpvecs(1,:,:))); 


% ToDo

%   Inspect zmat covariance matrix for specific sets of measures

%   Compute omnibut stat based bvec and standard error covariance matrix   

%   Inspect ymat_resid_wht for outliers
%     Look for SNPs with multiple ROI hits for ymat_resid_wht -- should be ones that would benefit most from "projection" (e.g., snpi=614489)

%   Normalize zmat columns as input to svd

%   Re-map measures to fit Gaussian distribution, rather than eliminating measures

%   Split measures into different "types", and residualize regional  measures by total/global

%   Test entire procedure under Null / w. permutation

%   Re-run all with imputed genotypes (Rob?)

%   Generate "de-noised" ymat_resid 
%      Use beta-hat for each snp projected onto signal subspace?
%        Does this improve replication performance across independent samples? -- test as function of sample seize, w. eandomly selected training and validation set, for specific confidently non-null SNPs 

%   Find classification method invariant to sign of projection (e.g., cos^2 rather than cos)
%     Alternatively, eliminate reduntant classes

%  Check on precision of inputs to Voxelwise_GLM -- should be double


% Extras

ncomp_kmeans = 5;
% k-means clustering
[T_kmeans U_kmeans] = kmeans(tmp(:,:,1)',ncomp_kmeans,'Distance','cosine','Replicates',100); U_kmeans = U_kmeans';
sfigure(1001); imagesc(U_kmeans,max(abs(U_kmeans(:)))*0.5*[-1 1]); colormap(blueblackred);
[T_kmeans U_kmeans] = kmeans(tmp(:,:,2)',ncomp_kmeans,'Distance','cosine','Replicates',100); U_kmeans = U_kmeans';
sfigure(1002); imagesc(U_kmeans,max(abs(U_kmeans(:)))*0.5*[-1 1]); colormap(blueblackred);

% Compute histograms
nr = round(sqrt(ncomp)); nc = ceil(ncomp/nr);
hv = linspace(-7,7,101); hcvecs = NaN(ncomp,length(hv),2);
for compi = 1:ncomp
  for whitenflag = 0:1
    compvec = U(:,compi);
    if whitenflag
      yvec = ymat_resid*C_inv*compvec;
    else
      yvec = ymat_resid*compvec;
    end
    defvec = isfinite(sum(yvec,2));
    yvec = yvec./std(yvec(defvec));
    hcvecs(compi,:,whitenflag+1) = hist(yvec(defvec),hv);
  end
end

sfigure(21); clf; sfigure(22); clf;
for compi = 1:ncomp
  sfigure(21); subplot(nc,nr,compi); bar(hv,hcvecs(compi,:,1))
  sfigure(22); subplot(nc,nr,compi); bar(hv,hcvecs(compi,:,2))
  drawnow;
end


bvecs = NaN(ncomp,nsnp,2); zvecs = NaN(ncomp,nsnp,2);
%for compi = 1:ncomp
for compi = 2:ncomp
  fprintf(1,'compi=%d of %d (now=%s)\n',compi,ncomp,datestr(now));
  for whitenflag = [0 1]
    compvec = U(:,compi);
    if whitenflag
      yvec = ymat_resid*C_inv*compvec;
    else
      yvec = ymat_resid*compvec;
    end
    for snpi = 1:nsnp
      gvec = double(genomat_cat(:,snpi));
      gvec(gvec<0) = NaN;
      X(:,1) = gvec;
      defvec = isfinite(sum(X,2)+sum(yvec,2));
      [pvols,betavols,beta_cov,sigvol] = Voxelwise_GLM(X(defvec,:),yvec(defvec,:),contrastvec); b = betavols(1,:); se = sqrt(beta_cov(1,1))*sigvol; z = b/se;
      bvecs(compi,snpi,whitenflag+1) = b; zvecs(compi,snpi,whitenflag+1) = z;
    end
  end
  logpvecs = -log10(normcdf(-abs(zvecs))*2);
  sfigure(11); subplot(nc,nr,compi); plot(logpvecs(compi,:,1),'*'); %ylim([0 8]); xlim(1e5+[-1e3 1e3]);
  sfigure(12); subplot(nc,nr,compi); plot(logpvecs(compi,:,2),'*'); %ylim([0 8]); xlim(1e5+[-1e3 1e3]);
  drawnow;
end

nr = round(sqrt(ncomp)); nc = ceil(ncomp/nr);
sfigure(11); clf; sfigure(12); clf; 
for compi = 1:ncomp
  sfigure(11); subplot(nc,nr,compi); plot(logpvecs(compi,:,1),'*'); %ylim([0 8]); xlim(1e5+[-1e3 1e3]);
  sfigure(12); subplot(nc,nr,compi); plot(logpvecs(compi,:,2),'*'); %ylim([0 8]); xlim(1e5+[-1e3 1e3]);
end

save('~/UKB_zvecs_kmeans.mat','zvecs','bvecs','logpvecs','U_kmeans','ncomp');


%sfigure(666); clf; imagesc(corr(

% ToDo
%   Identify SNPs with combination of associations (betas) with largest Malanobis distance / lowest log-likelihood (e.g., strong asymetry?)
%   Compute GWAS on both unwhitened and whitened ymat_resid
%   Derive optimal projections / filters based on zmat and/or bmat
%     Construct ymat_resid based on estimated signal and noise covariance matrices?

%   Generate separate clusters / components for cortical area & thickness, subcortical volumes, after regressing out respective "global" measures
%      Perform all steps above for each set of measures

%   Get ico surface maps from Chun




% This is not working
ymat_proj_wht = ymat_resid*C_inv*U_kmeans(:,1:ncomp);
bmat2 = NaN(ncomp,nsnp); zmat2 = NaN(ncomp,nsnp);
for snpi = 1:nsnp
  gvec = double(genomat_cat(:,snpi));
  gvec(gvec<0) = NaN;
  X(:,1) = gvec;
  defvec = isfinite(sum(X,2)+sum(ymat_proj_wht,2));
  [pvols,betavols,beta_cov,sigvol] = Voxelwise_GLM(X(defvec,:),ymat_proj_wht(defvec,:),[1 0 0]); b = betavols(1,:); z = betavec./(sqrt(beta_cov(1,1))*sigvol);
  bmat2(:,snpi) = b; zmat2(:,snpi) = z;
  if mod(snpi,1000)==0
    fprintf(1,'snpi=%d of %d (now=%s)\n',snpi,nsnp,datestr(now));
    tmp = zmat2(:,1:snpi); tmp(~isfinite(tmp)) = 0;
    tmp = tmp(:,max(abs(tmp),[],1)>7);
    [mv mi] = max(abs(tmp),[],1);
    for j = 1:length(mi), tmp(:,j) = tmp(:,j)*sign(tmp(mi(j),j)); end
%    tmp = tmp.*repmat(sign(mean(tmp,1)),[size(tmp,1) 1]);
    sfigure(1); imagesc(tmp,max(abs(colvec(tmp)))*[-1 1]); colormap(blueblackred); colorbar; drawnow;
  end
end



% Examine stats for individual measures, with and without whitening -- does not always improve power
bvecs = NaN(2,nsnp); zvecs = NaN(2,nsnp);
for pheni = 1:length(phenonames)
  for whitenflag = [0 1]
    compvec = zeros(size(ymat,2),1); compvec(pheni) = 1;
    if whitenflag
      yvec = ymat_resid*C_inv*compvec;
    else
      yvec = ymat_resid*compvec;
    end
    for snpi = ivec_snp
      gvec = double(genomat_cat(:,snpi));
      gvec(gvec<0) = NaN;
      X(:,1) = gvec;
      defvec = isfinite(sum(X,2)+sum(yvec,2));
      [pvols,betavols,beta_cov,sigvol] = Voxelwise_GLM(X(defvec,:),yvec(defvec,:),[1 0 0]); b = betavols(1,:); se = sqrt(beta_cov(1,1))*sigvol; z = b/se;
      bvecs(whitenflag+1,snpi) = b; zvecs(whitenflag+1,snpi) = z;
    end 
  end
  fprintf(1,'%03d (%-20s): %6.2f %6.2f\n',pheni,phenonames{pheni},max(abs(zvecs),[],2));
end



% Extras

% Hierarchical clustering -- need to check with Chun on this
Y = pdist(tmp','cosine');
Z = linkage(Y,'average');
sfigure(666); clf; dendrogram(Z)
T_pdist = cluster(Z,'maxclust',ncomp);


% Old Code

%fname_snap2 = '/usr/tmp/UKB_snap2.mat';
fname_snap2 = '~/Downloads/UKB_snap2.mat';
load(fname_snap2);

ntrait = size(zmat,2);


M = [ones(length(Hvec),1) colvec(Hvec.*TLDvec)]; W = pinv(M);
C0 = NaN(ntrait); Cg = NaN(ntrait);
for i = 1:ntrait
  fprintf(1,'i = %d of %d (now=%s)\n',i,ntrait,datestr(now,'HH:MM:SS'));
  for j = 1:i
    b = W*(zmat(:,i).*zmat(:,j));
    C0(i,j) = b(1); C0(j,i) = b(1);
    Cg(i,j) = b(2); Cg(j,i) = b(2);
  end
end

figure(1); clf; imagesc(C0,[-1 1]); colormap(blueblackred); axis equal xy tight;
figure(2); clf; imagesc(Cg,max(abs(Cg(:)))*[-1 1]); colormap(blueblackred); axis equal xy tight; % Not sure why Cg doesn't look like a covariance matrix

[mv coli] = max(mean(zmat.^2,1));
figure; plot(Hvec.*TLDvec,zmat(:,coli),'*')

xvec = Hvec.*TLDvec; yvec = zmat(:,coli).^2; xvals = linspace(min(xvec),max(xvec),11);
[yvals_mean yvals_count] = histmean(xvec,yvec,xvals);
figure; plot(xvals,yvals_mean,'*')
figure; plot(xvals,yvals_count,'*')

C0 = corr(zmat((xvec<prctile(xvec,1))&(min(pmat,[],2)>1e-4),:)); % Correlation matrix for depleted, insignificant SNPs
Cg = corr(zmat((xvec>prctile(xvec,50))&(min(pmat,[],2)<1e-4),:)); % Correlation matrix for enriched, significant SNPs

figure(1); clf; imagesc(C0,[-1 1]); colormap(blueblackred); axis equal xy tight;
figure(2); clf; imagesc(Cg,[-1 1]); colormap(blueblackred); axis equal xy tight;

% Should compare C0 to correlation matrix for phenotypes


% Load in summary stats from Chun
UKBstruct = load('data/UKB/UKB10336_cal_v1_orig_cont_summary.mat');

% Load in annotation data, and merge with summary stats
annostruct = load('~/GWAS_Annot/9m/annomat.mat');
infostruct = load('~/GWAS_Annot/9m/infomat.mat');
if 0 % Merge based on SNP_ID
  [dummy IA IB] = intersect(UKBstruct.snplist,infostruct.snpidlist,'stable'); % 580703 of 594503 matches
else % Merge based on chromosomal position
  [dummy IA IB] = intersect(cat(2,UKBstruct.chrvec,UKBstruct.posvec),cat(2,infostruct.chrnumvec,infostruct.posvec),'rows','stable'); % 589963 of 594503 matches
end

logpmat = UKBstruct.logpmat(IA,:);
betamat = UKBstruct.betamat(IA,:);
semat = UKBstruct.betamat(IA,:);
annomat = annostruct.annomat(IB,:);
TLDvec = annomat(:,20);
mafvec = infostruct.mafvec(IB);
Hvec = 2*mafvec.*(1-mafvec);
pnamesvec = UKBstruct.pnamesvec;

ivec0 = sum(annomat(:,1:end-1),2) < 1.0; % Should update this to include TLD and H
ivec1 = sum(annomat(:,[4 8 9 10]),2) > 2.0; % Identify enriched SNPs -- should use optimal weighting

zmat = -norminv(10.^-logpmat/2).*sign(betamat);
zmat(find(~isfinite(zmat))) = NaN;

z2meanvec = colvec(nanmean(zmat.^2,1));
z4meanvec = colvec(nanmean(zmat.^4,1));
kurtvec = z4meanvec./z2meanvec.^2;
%z4meanvec = colvec(mean(zmat.^4,1));

figure(1); clf; plot([z2meanvec-1 kurtvec-3],'*');

%[sv  si] = sort(kurtvec,'descend'); % Sort based on kurtosis
[sv  si] = sort(z2meanvec,'descend'); % Sort based on variance (~heritability)

for sii = 1:length(si)
  fprintf(1,'%3d: exess kurtosis = %4.1f excess varance = %4.1f %3d %s\n',sii,kurtvec(si(sii))-3,z2meanvec(si(sii))-1,si(sii),pnamesvec{si(sii)});
end

%phenogroup = 'brainvols';
phenogroup = 'cognitive';
%phenogroup = 'anthropomorphic';
switch lower(phenogroup)
  case {'brainvols'}
    collist1 = find(find_cell(strfind(pnamesvec,'Volume of ')));
    collist2 = find(find_cell(strfind(pnamesvec,'Volume of grey matter in ')));
%    collist = setdiff(collist1,collist2); % Exclude 'Volume of grey matter in '
    collist = collist1;
  case {'cognitive'}
    collist = [167 407 169 168 409 85 86 410 408 183 411 134 88 140 145 139 84 141 135 143 142];
  case {'anthropomorphic'}
    collist = [26 40 43 42 25 24 132 133 81 41];
end 

lamvec = ones(1,length(collist)); % Vector of Lambdas (sig0^2) per trait

zmat_tmp = zmat(:,collist);
C0 = nancov(zmat_tmp(ivec0,:)); 
M = diag(sqrt(lamvec./diag(C0))); C0 = M*C0*M; 
C = nancov(zmat_tmp(ivec1,:)); % Should look at the most enriched SNPs rather than all

zmat_whitened = zmat_tmp*pinv(chol(C0));
C0_w = nancov(zmat_whitened(ivec0,:));
M = diag(sqrt(lamvec./diag(C0_w))); C0_w = M*C0_w*M; 
C_w = nancov(zmat_whitened(ivec1,:));

figure(2); clf; crange = [-1 1];
subplot(2,2,1); imagesc(C0,crange); colormap(blueblackred); axis xy equal tight;
subplot(2,2,2); imagesc(C,crange); colormap(blueblackred); axis xy equal tight;
subplot(2,2,3); imagesc(C0_w,0.1*crange); colormap(blueblackred); axis xy equal tight;
subplot(2,2,4); imagesc(C_w,0.1*crange); colormap(blueblackred); axis xy equal tight;

sort(diag(C),'descend')'
sort(diag(C_w),'descend')'
svd(C_w)'

% Should look for projections that maximize kurtosis (try kICA.m, fastICA.m -- look at demo_ICA.m)


% ToDo
%   Re-scale diagonals of C0, C0_w based on sig0^2 estimate for each trait (start by assuming unity)
%   Extract an d save genotypes and phenotypes (esp. imaging)
%   Optimize linear projection(s), test on raw data, try to predict based on summary stats
%   Should correct for sample size (Neff) per phenotype

%   Identify sets of variables: 
%     'Volume of grey matter in '
%     'Volume of '
%     ' T2star'
%     ' BOLD'
%     ' fluid intelligence '
%     ' diabetes '
%     ' circumference '
%     ' height '
%     '   intake '
%     '  blood pressure'
%     'Waist circumference '
%     'Body mass index '
%     'Weight '
%     'Sleep '
%     'Pulse rate'


% Go over UKB variables w. Olav -- identify interesting ones

