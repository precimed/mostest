
mu = 0;
sigma = 1;
beta = 0.4;

n_sample_t = 5002;   % _t means 'total'
n_cov_t = 300;
cov_size=0:100:n_cov_t;              % number of covariates
cohort_size=2:1000:n_sample_t;

x_t = normrnd(mu, sigma, [n_sample_t, 1]);
c_t = normrnd(mu, sigma, [n_sample_t, n_cov_t]);

eps = rand(n_sample_t, 1);
y_t = x_t*beta + eps;
rat = var(x_t*beta)/var(y_t);fprintf('variance: %f\n',rat)
y_t=y_t/norm(y_t);                   % normalize

t_preres = nan(length(cohort_size), length(cov_size));
t_multreg = nan(length(cohort_size), length(cov_size));

effN_preres = nan(length(cohort_size), length(cov_size));
effN_multreg = nan(length(cohort_size), length(cov_size));

%rand('seed', 123)

for i=1:length(cohort_size) 
    fprintf('subjects: %f\n',cohort_size(i))
    x = x_t(1:cohort_size(i));
    y = y_t(1:cohort_size(i));
    for j = 1:length(cov_size)
        fprintf('covariates: %f\n',cov_size(j))
        c = c_t(1:cohort_size(i), 1:cov_size(j));
        y_res   = y - c*(c\y); % pre-residulize
        [coef, t_preres(i,j)] = nancorr(y_res, x);
    
        if (cov_size(j) ~= 0) && (cov_size(j)+2 >= cohort_size(i)), continue; end
       
        % try multiple regression (instead of pre-residualization)
        % https://stats.stackexchange.com/questions/27916/standard-errors-for-multiple-regression-coefficients
        X = [x c];
        hat_beta = X \ y;
        s2=var(y - X * hat_beta);
        %Sigma = s2 * inv(X' * X); se = sqrt(Sigma(1,1)); 
        b1 = zeros(size(X, 2), 1); b1(1)=1; a=(X' * X \ b1); se = sqrt(s2 * a(1)); 
        t_multreg(i, j) = hat_beta(1) / se;
    end
end
z2_preres=t_preres.*t_preres;
z2_multreg=t_multreg.*t_multreg;

for j = 1:length(cov_size)
    effN_multreg(:, j) =  z2_multreg(:, j) ./ z2_multreg(:, 1) .* cohort_size';
    effN_preres(:, j) =  z2_preres(:, j) ./ z2_preres(:, 1) .* cohort_size';
end

figure(1); clf; hold on;legendInfo={};
for j=1:length(cov_size)
    plot(cohort_size, z2_preres(:,j), 'LineWidth', 2)
    legendInfo{j} = ['No.cov= ' num2str(cov_size(j))]; 
end
set(gca,'ColorOrderIndex',1)
for j=1:length(cov_size)
    plot(cohort_size, z2_multreg(:,j), '--', 'LineWidth', 2)
end
xlabel('Sample size')
ylabel('Z^2')
legend(legendInfo,'Position',[0.7 0.2 0.1 0.2])
title({'solid lines: pre-regularization' 'dashed lines: multiple regression'})

figure(2); clf; hold on; legendInfo={};
for i=1:length(cohort_size)
    plot(cov_size, effN_preres(i, 1) - effN_preres(i,:), 'LineWidth', 2)
    legendInfo{i} = ['No.subj= ' num2str(cohort_size(i))]; 
end
set(gca,'ColorOrderIndex',1)
for i=1:length(cohort_size)
    plot(cov_size, effN_multreg(i,1) - effN_multreg(i,:), '--', 'LineWidth', 2)
end
xlabel('No.cov')
ylabel('totalN - effN (delta)')
legend(legendInfo, 'location', 'northwest')


