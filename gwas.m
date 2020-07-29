function [z,beta,beta_se] = gwas(P, G)
    % Args:
    %   P: nsampels x npheno phenotype matrix.
    %   G: nsampels x nsnps genotype matrix (e.g. produced by bed2mat).
    %
    % Return:
    %   z:       npheno x nsnps matrix of z-scores = beta/beta_se.
    %   beta:    npheno x nsnps matrix with GWAS beta coefficients.
    %   beta_se: npheno x nsnps matrix with standard errors of beta estimates.

    iP = isfinite(P);
    iG = G ~= -1;
    
    if strcmp(class(P), 'single')
        vP = single(iP);
        vG = single(iG);
        G = single(G);
    else
        vP = double(iP);
        vG = double(iG);
        G = double(G);
    end
    
    P(~iP) = 0;
    G(~iG) = 0;

    yx = P' * G;
    n = vP' * vG;
    y = P' * vG;
    x = vP' * G;
    yy = (P.*P)' * vG;
    xx = vP' * (G.*G);

    beta = (n.*yx - y.*x) ./ (n.*xx - x.*x);
    beta_se = sqrt( ((n.*yy - y.*y) ./ (n.*xx - x.*x) - beta.^2) ./ (n - 2) );
    z = beta ./ beta_se;

    %cov2 = (n.*yx - y.*x).^2;
    %z = sqrt( cov2.*(n-2) ./ ((n.*yy - y.*y).*(n.*xx - x.*x) - cov2) );
end
