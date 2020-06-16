function pval_vec = fixed_paretotails_cdf(orig_paretotails, x)

% WARNING: this function is used due to potential bug in MATLAB paretotails.
%          The lowest value produced by paretotails.cdf is ~ 1E-17 (regardless whether lower or
%          upper tail is taken). All lower values are truncated to 0. Directly using
%          GeneralizedPareto with corresponding parameters allows to handle values down to ~ 1E-300.

% Args:
%   orig_paretotails: an object created with MATLAB's paretotails(data,0,pu) function, i.e. fitting only the upper tail.
%   x:                matrix of points to use for cdf evaluation.

% Make GeneralizedPareto distribution using shape, scale and location paraeters from the orig_paretotails.
shape_scale = upperparams(orig_paretotails);
[p,loc] = boundary(orig_paretotails);
gpd = makedist('GeneralizedPareto','k',shape_scale(1),'sigma',shape_scale(2),'theta',loc);

pval_vec = cdf(orig_paretotails, x, 'upper');
% Estimate pvalues for the tail (which has probability p) with GeneralizedPareto.
i_tail = x > loc;
pval_vec(i_tail) = (1.0 - p)*cdf(gpd, x(i_tail), 'upper');

