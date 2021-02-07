
function [index_cohort] = divide(n_subj,n_cohort)

x = randperm(n_subj);

x_s=x(2:end-1);
idx=randperm(numel(x_s),n_cohort-1);

midpoint=x_s(idx);
[~,index] = intersect(x,midpoint);
index=sort(index);

index_cohort{1}=sort(x(1:index(1)));
for i=2:n_cohort-1
    index_cohort{i}=sort(x(index(i-1)+1:index(i)));
end
index_cohort{n_cohort}=sort(x(index(n_cohort-1)+1:end));
clear x;

