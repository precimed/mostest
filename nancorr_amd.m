function [coef, t, n] = nancorr_amd(A, B)
  %NANCORR - Pearson correlation coefficient
  %
  %  coef = NANCORR(A, B) is equivalent to 
  %  coef = corr(A, B, 'rows','pairwise'),
  %  but NANCORR works much faster.
  %  
  %  [coef, t, n] = NANCORR(A, B)
  %  INPUT:
  %    A, B - input matrices, single or double, with equal number of rows
  %
  %  OUTPUT:
  %    coef - matrix of Pearson correlation coefficients
  %    t    - matrix of t-statistics
  %    n    - matrix containing the number of defined values
  %
  % NOTES  
  %  pvalue can be calculated as 2*tcdf(-abs(t), n - 2)

  %Am, Ap: m=missing elements, p=present elements

  gpudev = [];
  try
    gpudev = gpuDevice; % Get GPU info
  end
  usegpuflag = ~isempty(gpudev); % Should check type of GPU, datatype(s) and perhaps load vs. CPU?

  Am=~isfinite(A); Bm=~isfinite(B); 

  if strcmp(class(A), 'single')
    Ap=single(~Am); Bp=single(~Bm);    
    if usegpuflag
       Am=gpuArray(Am); Bm=gpuArray(Bm);
       Ap=gpuArray(Ap); Bp=gpuArray(Bp);
       A=gpuArray(A); B=gpuArray(B); % Should perhaps check if sufficient GPU  memory is available?
    end
  else
    Ap=double(~Am); Bp=double(~Bm);
  end

  % zero out nan elements
  A(Am)=0; B(Bm)=0;

  % code one of the formulas from https://en.wikipedia.org/wiki/Pearson_correlation_coefficient
  % this procedure might be numericaly unstable for large values,
  % it might be reasonable to center each column before calling nancorr.

  xy = A' * B;          % sum x_i y_i
  n  = Ap' * Bp;        % number of items defined both in x and y
  
  mx = A' * Bp ./ n;    % mean values in x, calculated across items defined both in x and y
  my = Ap' * B ./ n;    % mean values in y, calculated across items defined both in x and y
  
  x2 = (A.*A)' * Bp;    % sum x^2_i, calculated across items defined both in x and y
  y2 = Ap' * (B.*B);    % sum y^2_i, calculated across items defined both in x and y
  
  sx   = sqrt(x2 - n .* (mx.^2));  % sx, sy - standard deviations 
  sy   = sqrt(y2 - n .* (my.^2));
  
  coef = (xy - n .* mx .* my) ./ (sx .* sy);      % correlation coefficient
  t    = coef .* sqrt((n - 2) ./ (1 - coef.^2));  % t-test statistic

  if strcmp(class(coef),'gpuArray'), coef = single(gather(coef)); end
  if strcmp(class(t),'gpuArray'), t = single(gather(t)); end
  if strcmp(class(n),'gpuArray'), n = single(gather(n)); end

% keyboard

