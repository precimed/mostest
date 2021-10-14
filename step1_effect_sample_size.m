clear all;
close all;
clc;

n_cohort=100;
n_sample=randi([200 260],1,n_cohort);
n_sample_t=sum(n_sample);
index_cohort=cumsum(n_sample);

mu=0;
sigma=10;
x_t=normrnd(mu,sigma,[n_sample_t,1]);
beta=0.06;
eps=rand(n_sample_t,1);
y_t=x_t*beta+eps;
rat=var(x_t*beta)/var(y_t);

for i=1:n_cohort
    x_i=x_t(1:index_cohort(i));
    y_i=y_t(1:index_cohort(i)); 
    [coef, t(i)] = nancorr(y_i, x_i);
end

z2=t.*t;

figure
plot(index_cohort,z2,'b','LineWidth',2)
xlabel('Sample size')
ylabel('Z^2')
xlim([index_cohort(1) index_cohort(end)])
return
