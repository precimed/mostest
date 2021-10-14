clear all;
close all;
clc;

n_cohort=11;
n_sample=randi([250 270],1,n_cohort);
n_sample_t=sum(n_sample);
index_cohort=cumsum(n_sample);

mu=0;
sigma=10;
x_t=normrnd(mu,sigma,[n_sample_t,1]);
beta=0.06;
eps=rand(n_sample_t,1);
y_t=x_t*beta+eps;
rat=var(x_t*beta)/var(y_t);

cov=0:500:index_cohort(end);% number of covariates
n_cov=size(cov,2)% loop for different number of covariates

for i=1:n_cohort 
    x_i=x_t(1:index_cohort(i));
    y_i=y_t(1:index_cohort(i));
    y_i=y_i/norm(y_i);% normalize
    for j=1:n_cov
        mc=normrnd(mu,sigma,[index_cohort(i),cov(j)]);
        % generate effect size of covariates
        r=linspace(0.1,0.2,cov(j))';
        beta_es=mc\y_i;%(y-eps);
        y_res=y_i-mc*beta_es; % residulize
        
        [coef, t(i,j)] = nancorr(y_res, x_i);
    end
end
z2=t.*t;
  
figure
for j1=1:n_cov
    plot(index_cohort,z2(:,j1),'LineWidth',2)
    %plot(index_cohort,z2(:,j1),'LineWidth',2)
    legendInfo{j1} = ['No.cov= ' num2str(cov(j1))]; 
    hold on
end
xlabel('Sample size')
ylabel('Z^2')
xlim([-index_cohort(end) index_cohort(end)])
xlim([index_cohort(1) index_cohort(end)])
legend(legendInfo,'Position',[0.7 0.2 0.1 0.2])

%z_eff=z2/coe;
z_eff=index_cohort(end)*z2./z2(:,1)


%% effective sample size
figure
X = cov

Y = floor(z_eff(n_cohort,:))

plot(X,Y,'k-s')
xlim([0 cov(end)*1.1])
ylim([0 Y(1)*1.1])
strValues = strtrim(cellstr(num2str([X(:) Y(:)],'(%d,%d)')));
text(X,Y,strValues,'VerticalAlignment','bottom');
xlabel('Number of covariates')
ylabel('Effective sample size')


