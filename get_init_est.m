function [comp_PMF_init,mu_init,Sigma_init] = get_init_est()

% initial guesses for component probabilities

comp_PMF_init = [0.2,0.5,0.3];

% initial guesses for mean vectors

mu_init(1,:) = [1,3];
mu_init(2,:) = [5,2];
mu_init(3,:) = [4,6];

% initial guesses for covariance matrices

Sigma_init(:,:,1) = 0.1*eye(2);
Sigma_init(:,:,2) = 0.2*eye(2);
Sigma_init(:,:,3) = 0.3*eye(2);

end