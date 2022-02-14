function pd = weibulltails(x,pl,pu,distrib)

if ~exist('distrib','var') | isempty(distrib)
  distrib = 'weibull';
end

[yvals xvals ylvals yuvals] = ecdf(x); 
yvals = 1-yvals; ylvals = 1-ylvals; yuvals = 1-yuvals; % Use upper tail

xvals = (xvals(1:end-1)+xvals(2:end))/2; % Piecewise linear interpolation, to enble interp1
yvals = (yvals(1:end-1)+yvals(2:end))/2; 
ylvals = (ylvals(1:end-1)+ylvals(2:end))/2; 
yuvals = (yuvals(1:end-1)+yuvals(2:end))/2;

%figure(666); plot(xvals,-log10(yvals),xvals,-log10(ylvals),xvals,-log10(yuvals),'LineWidth',2);

xvec = colvec(linspace(0,5*max(x),1001));
yvec = interp1(xvals,yvals,xvec,'linear');
ylvec = interp1(xvals,ylvals,xvec,'linear');
yuvec = interp1(xvals,yuvals,xvec,'linear');

%figure(667); plot(xvec,-log10(yvec),xvec,-log10(ylvec),xvec,-log10(yuvec),'LineWidth',2);


%figure(668); plot(log(xvec(ivec_fit)),log(-log(yvec(ivec_fit))),log(xvec(ivec_fit)),log(-log(ylvec(ivec_fit))),log(xvec(ivec_fit)),log(-log(yuvec(ivec_fit))),'LineWidth',2);

xvec_tmp = colvec(log(xvec)); 
switch lower(distrib)
  case {'weibull' 'wbl'}
    wvec = colvec((log(-log(yuvec))-log(-log(ylvec))).^2);
    ivec_fit = find(xvec>=prctile(x,100*pl)&xvec<=prctile(x,100*pu)&isfinite(wvec));
    yvec_tmp = colvec(log(-log(yvec))); 
    M = cat(2,xvec_tmp,ones(size(xvec_tmp))); % Linear -- liberal?
    betahat = lscov(M(ivec_fit,:),yvec_tmp(ivec_fit),wvec(ivec_fit));
    yvec_fit = exp(-exp(M*betahat)); 
  case {'beta' 'logbeta'}
    wvec = colvec((-log(yuvec)-(-log(ylvec))).^2);
    ivec_fit = find(xvec>=prctile(x,100*pl)&xvec<=prctile(x,100*pu)&isfinite(wvec));
    yvec_tmp = colvec(-log(yvec));
    parvec0 = betafit(10.^(-x));
    costfun = @(x)sum(wvec(ivec_fit).*(-log(betacdf(10.^(-xvec(ivec_fit)),x(1),x(2)))-yvec_tmp(ivec_fit)).^2);
    parvec = fminsearch(costfun,parvec0,optimset('MaxFunEvals',1000));
    yvec_fit = betacdf(10.^(-xvec),parvec(1),parvec(2));
  case {'sidak' 'minp'}
    wvec = colvec((-log(yuvec)-(-log(ylvec))).^2);
    ivec_fit = find(xvec>=prctile(x,100*pl)&xvec<=prctile(x,100*pu)&isfinite(wvec));
    yvec_tmp = colvec(-log(yvec));
    parvec0 = betafit(10.^(-x)); parvec0 = parvec0(2); 
    costfun = @(x)sum(wvec(ivec_fit).*(-log(betacdf(10.^(-xvec(ivec_fit)),1,x(1)))-yvec_tmp(ivec_fit)).^2);
    parvec = fminsearch(costfun,parvec0,optimset('MaxFunEvals',1000));
    yvec_fit = betacdf(10.^(-xvec),1,parvec(1));
  case {'gamma'}
    wvec = colvec((-log(yuvec)-(-log(ylvec))).^2);
    ivec_fit = find(xvec>=prctile(x,100*pl)&xvec<=prctile(x,100*pu)&isfinite(wvec));
    yvec_tmp = colvec(-log(yvec)); 
    parvec0 = gamfit(x);
    costfun = @(x)sum(wvec(ivec_fit).*(-log(gamcdf(xvec(ivec_fit),x(1),x(2),'upper'))-yvec_tmp(ivec_fit)).^2);
    parvec = fminsearch(costfun,parvec0,optimset('MaxFunEvals',1000));
    yvec_fit = gamcdf(xvec,parvec(1),parvec(2),'upper');
  case {'chi2'}
    wvec = colvec((-log(yuvec)-(-log(ylvec))).^2);
    ivec_fit = find(xvec>=prctile(x,100*pl)&xvec<=prctile(x,100*pu)&isfinite(wvec));
    yvec_tmp = colvec(-log(yvec));
    parvec0 = gamfit(x); parvec0 = parvec0(1);
    costfun = @(x)sum(wvec(ivec_fit).*(-log(gamcdf(xvec(ivec_fit),x(1),2,'upper'))-yvec_tmp(ivec_fit)).^2);
    parvec = fminsearch(costfun,parvec0,optimset('MaxFunEvals',1000));
    yvec_fit = gamcdf(xvec,parvec(1),2,'upper');
end

%figure(669); plot(log(xvec(ivec_fit)),log(-log(yvec(ivec_fit))),log(xvec(ivec_fit)),log(-log(yvec_fit(ivec_fit))),'LineWidth',2);

yvec_pred = yvec_fit;
blendvec = interp1([xvec(1) xvec(ivec_fit(1)) xvec(ivec_fit(end)) xvec(end)],[0 0 1 1],xvec,'linear','extrap');
ivec_blend = find(blendvec<1);
yvec_pred(ivec_blend) = exp(-exp((1-blendvec(ivec_blend)).*log(-log(yvec(ivec_blend)))+blendvec(ivec_blend).*log(-log(yvec_fit(ivec_blend)))));
yvec_pred(1:(min(find(isfinite(yvec_pred)))-1)) = 1;

%figure(670); plot(xvec,-log10(yvec),xvec,-log10(yvec_pred),'LineWidth',2);

pd = struct;
pd.nlogcdf = @(x)exp(interp1(xvec,log(-log(yvec_pred)),x,'linear'));
pd.cdf = @(x,y)condexp(exist('y','var')&strcmp(y,'upper'),exp(-pd.nlogcdf(x)),1-exp(-pd.nlogcdf(x)));

%figure(671); plot(xvec,-log10(yvec),xvec,-log10(pd.cdf(xvec,'upper')),'LineWidth',2);

% keyboard

return


