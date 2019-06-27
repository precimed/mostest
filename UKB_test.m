addpath(genpath('/home/cfan/codes/matlab')); % get Chun's /home/cfan/codes/matlab/PlinkRead_binary.m

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

data_covars = readtext_amd('~dale/Dropbox/data/UKB/UK500k_BasicCovs.txt',char(9));
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

data_keep = readtext_amd('~dale/Dropbox/data/UKB/KeepCauc.txt',char(9));
IDvec_keep = data_keep(:,1);

% Limit to Caucasian subjects -- should perhap;s avoid re-using variables?
ivec_tmp = ismember(IDvec_covars,IDvec_keep);
IDvec_covars = IDvec_covars(ivec_tmp);
datamat_covars = datamat_covars(ivec_tmp,:);

data = readtext_amd('~dale/NORMENT_T1_FS.ROI.csv');
colnames = StripQuotes_amd(data(1,2:end));
IDvec_FS = strrep(data(2:end,1),'"','');
data = StripQuotes_amd(data(2:end,2:end));

chrlist = cat(2,cellfun(@(x)num2str(x),num2cell([1:22]),'UniformOutput',false),'X','Y','XY','MT');

for chri = 1:length(chrlist)
%for chri = 24:length(chrlist)
  chr = chrlist{chri};
  fname_fam = sprintf('/space/gwas-ds2/1/data/GWAS/UKBio/genotype_calls/ukb2741_cal_chr%s_v2_s488366.fam',chr);
  fname_bed = sprintf('/space/gwas-ds2/1/data/GWAS/UKBio/genotype_calls/ukb_cal_chr%s_v2.bed',chr);
  fname_bim = sprintf('/space/gwas-ds2/1/data/GWAS/UKBio/genotype_calls/ukb_snp_chr%s_v2.bim',chr);
  tic; [chrvec0 snplist A1vec A2vec cMvec bpvec0 FID IID sexvec phenovec genomat] = PlinkRead_binary_separate(fname_bed, fname_fam, fname_bim); toc
  IDvec_plink = IID;
  [dummy IA IB] = intersect(IDvec_FS,IDvec_plink,'stable');
  genomat = genomat(IB,:); FID = FID(IB); IID = IID(IB);
  save(sprintf('~/UKB/UKB_snap%02d.mat',chri),'genomat','chrvec0','bpvec0','snplist','FID','IID','sexvec','phenovec');
end

genomat_cat = []; chrvec_cat = []; bpvec_cat = []; snplist_cat = {};
for chri = 1:length(chrlist)
  chr = chrlist{chri};
  fprintf(1,'%d: chr%s (now=%s)\n',chri,chr,datestr(now));
  load(sprintf('~/UKB_snap%02d.mat',chri));
  [dummy IA IB] = intersect(IID,IDvec_covars,'stable');
  genomat_cat = cat(2,genomat_cat,genomat(IA,:));
  chrvec_cat = cat(1,chrvec_cat,chrvec0);
  bpvec_cat = cat(1,bpvec_cat,bpvec0);
  snplist_cat = cat(1,snplist_cat,snplist);
end
datamat_covars_merged = datamat_covars(IB,:); IDvec = IDvec_covars(IB);

[dummy IA2 IB2] = intersect(IDvec,IDvec_FS,'stable');
datamat = cell2mat_amd(data(IB2,:));

ivec_pheno = var(datamat,[],1)>0; % All measures with non-zero variance
%ivec_pheno = find_cell(cellfun(@(x)regexp(x,'.vol$|^ICV$','match'),colnames,'UniformOutput',false))&var(datamat,[],1)>0; % Subcortical volumes
%ivec_pheno = find_cell(cellfun(@(x)regexp(x,'.area.|^area','match'),colnames,'UniformOutput',false))&var(datamat,[],1)>0; % Cortical area
%ivec_pheno = find_cell(cellfun(@(x)regexp(x,'.thick.|^thick','match'),colnames,'UniformOutput',false))&var(datamat,[],1)>0; % Cortical thickness
%ivec_pheno = find_cell(cellfun(@(x)regexp(x,'.vol.','match'),colnames,'UniformOutput',false)); % Cortical volumes

ymat = datamat(:,ivec_pheno); phenonames = colnames(ivec_pheno);

nroi = size(ymat,2);
skewvec = NaN(1,nroi); kurtvec = NaN(1,nroi);
for roii = 1:nroi
  yvec = ymat(:,roii);
  defvec = isfinite(yvec);
  skewvec(roii) = skewness(yvec(defvec)); kurtvec(roii) = kurtosis(yvec(defvec));
  fprintf(1,'%3d %-40s: skewness=%5.2f kurtosis=%5.2f\n',roii,phenonames{roii},skewvec(roii),kurtvec(roii));
end

ymat(:,kurtvec>6) = 0; % Eliminate "bad" measures
X_tmp = cat(2,datamat_covars_merged,ones(size(ymat,1),1));
defvec =  isfinite(sum(ymat,2)+sum(X_tmp,2));
ymat_resid = ymat - X_tmp*(X_tmp(defvec,:)\ymat(defvec,:));
ymat_resid = ymat_resid ./ repmat(max(0.01,std(ymat_resid,[],1)),[size(ymat_resid,1) 1]);
ymat_resid(abs(ymat_resid)>5) = 0; % Censor outliers -- should set to NaN instead?

% Compute covariance
lam_reg = 0.1;
%lam_reg = 0.5;
C = cov(ymat_resid); % Should residualize for all covariates!
C_reg = (1-lam_reg)*C + lam_reg*diag(max(0.01,diag(C))); % Ridge regularized covariance matrix
C_inv = inv(C_reg);
W_wht = chol(C_inv);
ymat_resid_wht = ymat_resid*W_wht'; % Whitened data
ymat_resid_wht(abs(ymat_resid_wht)>5) = 0; % Censor outliers -- should set to NaN instead?
sfigure(1); clf; imagesc(C,[-1 1]); colormap(blueblackred); axis equal tight;

%C_wht = cov(ymat_resid_wht);
%sfigure(2); clf; imagesc(C_wht,[-1 1]); colormap(blueblackred); axis equal tight;


X = cat(2,NaN(size(ymat,1),1),ones(size(ymat,1),1)); contrastvec = [0 1];
nsnp = size(genomat_cat,2); npheno = sum(ivec_pheno);
bmat = NaN(npheno,nsnp,2,'single'); zmat = NaN(npheno,nsnp,2,'single'); logpmat = NaN(npheno,nsnp,2,'single');
for snpi = 1:nsnp
%for snpi = 263000:nsnp
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


% Estimate componnets

tmp = zmat; tmp(~isfinite(tmp)) = 0;
ivec_snp = find(max(max(abs(tmp),[],3),[],1)>abs(norminv(1e-5)));
tmp = tmp(:,ivec_snp,:);
for j = 1:size(tmp,2)
  v = tmp(:,j,1); [mv mi] = max(abs(v)); tmp(:,j,1) = v/v(mi);
  v = tmp(:,j,2); [mv mi] = max(abs(v)); tmp(:,j,2) = v/v(mi);
end
sfigure(1101); imagesc(tmp(:,:,1),max(abs(colvec(tmp(:,:,1))))*[-1 1]); colormap(blueblackred); colorbar; drawnow;
sfigure(1102); imagesc(tmp(:,:,2),max(abs(colvec(tmp(:,:,2))))*[-1 1]); colormap(blueblackred); colorbar; drawnow;


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


%ncomp = ncomp_svd; U1 = U_svd1(:,1:ncomp); U2 = U_svd2(:,1:ncomp);
ncomp = ncomp_kmeans; U1 = U_kmeans1(:,1:ncomp); U2 = U_kmeans2(:,1:ncomp);

bvecs = NaN(ncomp,nsnp,2); zvecs = NaN(ncomp,nsnp,2); logpvecs = NaN(ncomp,nsnp,2);
sfigure(1201); clf; sfigure(1202); clf;
nr = round(sqrt(ncomp)); nc = ceil(ncomp/nr);
for compi = 1:ncomp
  fprintf(1,'compi=%d of %d (now=%s)\n',compi,ncomp,datestr(now));
  compvec = U1(:,compi); 
  for whitenflag = [1 0]
    if whitenflag
%      yvec = ymat_resid_wht*U2(:,compi); % Use projections based on whitened stats
%      yvec = ymat_resid_wht*U1(:,compi); % Use projections based on non-whitened stats
%      yvec = ymat_resid*C_inv*U1(:,compi); % Use projections based on non-whitened stats
      yvec = ymat_resid*C_inv*compvec;
    else
%      yvec = ymat_resid_wht*U2(:,compi); % Use projections based on whitened stats
%      yvec = ymat_resid*U1(:,compi); % Use projections based on non-whitened stats
      yvec = ymat_resid*compvec;
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
          sfigure(1201); subplot(nc,nr,compi); plot(logpvecs(compi,:,1),'*'); axis tight;
        else
          sfigure(1202); subplot(nc,nr,compi); plot(logpvecs(compi,:,2),'*'); axis tight;
        end
        drawnow;
      end
    end
  end
end

fname_snap = '~/UKB/UKB_zvecs.mat';
save(fname_snap,'zvecs','bvecs','logpvecs');

% figure; GenStats_QQ_plot(squeeze(logpvecs(1,:,:))); 


% ToDo

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


% Eliminate ROIs with non-Gaussian distributions

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



% ToDo

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

