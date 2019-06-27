function R = corr_reg(M)

C = cov(M);
d = sqrt(max(eps,diag(C)).^-1);
R = diag(d)*C*diag(d); R = (R+R')/2;
ivec = find(diag(R)==0);
R(ivec,ivec) = eye(length(ivec));

if max(colvec(abs(R-R')))>0
  keyboard
end

