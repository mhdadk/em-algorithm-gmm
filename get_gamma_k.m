function gamma_k = get_gamma_k(k,X,comp_PMF,mu,Sigma)

% compute all the numerators for all gamma_jk

gamma_numerators = comp_PMF(k) * mvnpdf(X,mu(k,:),Sigma(:,:,k));

% initialize array to store all gamma denominators

gamma_denominators = zeros(length(X),1);

for i = 1:length(X)
    
    % compute all the denominators for all gamma_jk
    
    gamma_denominators(i) = dot(comp_PMF,...
                                mvnpdf(X(i,:),mu,Sigma));
    
end

% compute gamma_jk for all j values

gamma_k = rdivide(gamma_numerators,gamma_denominators);
    
end