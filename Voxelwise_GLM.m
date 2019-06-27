function [pvols,betavols,beta_cov,sig_vol] = Voxelwise_GLM(X,datamat,R,options)

if ~exist('options','var')
  options = struct();
end

if isfield(options,'synthflag'), synthflag = options.synthflag; else synthflag = false; end
if isfield(options,'plotflag'), plotflag = options.plotflag; else plotflag = false; end
if isfield(options,'pvolflag'), pvolflag = options.pvolflag; else pvolflag = true; end

Xi = pinv(X);

XiX = Xi*X; % Crosstalk matrix
xix = diag(XiX);
Xi(find(xix<0.9),:) = 0; % Eliminate non-identifiable params

C = Xi*Xi';
nvols = size(X,1);
nparams = size(X,2);
nvox = size(datamat,1);

if synthflag % Synthesize data?
  datamat = 2.1234*randn(size(datamat));
end

betavols = Xi*datamat;
datamat_hat = X*betavols;
vol_sum2 = sum((datamat_hat-datamat).^2,1);
sig_vol = sqrt(vol_sum2/(nvols-nparams));

beta_cov = C;

if isempty(R)
  pvols = [];
elseif size(R,1) == 1
  contrast_vec = R;
  sigma_contrast = sqrt(contrast_vec*C*contrast_vec');
  contrast_vol = contrast_vec*betavols;
  contrast_std = sigma_contrast*sig_vol;
  if pvolflag
    dof = nvols - nparams;
    contrast_vol(find(abs(contrast_vol)<10*eps)) = 0;
    tstats = contrast_vol./(sigma_contrast*sig_vol+eps);
    pvols = (tstats<0).*log10(2*tcdf(tstats,dof)) + (tstats>0).*(-log10(2*tcdf(-tstats,dof)));
    if plotflag
      [hc,hv] = hist(tstats(:),1000); hc = hc / sum(hc) / (hv(2)-hv(1));
      figure(1); plot(hv,hc,hv,tpdf(hv,dof));
    end
  else
    pvols = struct();
    pvols.contrast_mean = contrast_vol;
    pvols.contrast_std = contrast_std;
  end
else
  betavols2 = R*betavols;
  if pvolflag
    Ci = inv(R*inv(X'*X)*R');
    dof_numer = size(R,1);
    dof_denom = nvols - nparams;
    sumsq_vol = 0;
    for i = 1:size(betavols2,1)
      for j = 1:size(betavols2,1)
        sumsq_vol = sumsq_vol + betavols2(i,:).*betavols2(j,:)./(sig_vol.^2)*Ci(i,j);
      end
    end
    Fstats = sumsq_vol/dof_numer;
    pvols = -log10(max(eps,(1-fcdf(Fstats,dof_numer,dof_denom))));
    if plotflag
      [hc,hv] = hist(Fstats(:),1000); hc = hc / sum(hc) / (hv(2)-hv(1));
      figure(2); plot(hv,hc,hv,fpdf(hv,dof_numer,dof_denom));
    end
  else
    pvols = [];
  end
end
