function [N,comp_PMF,mu,Sigma] = get_GMM_parameters()

% original component probabilities used to generate the GMM

comp_PMF = [0.3,0.6,0.1];

% number of random vectors to sample from GMM

N = 500;

% original mean vectors used to generate GMM

mu(1,:) = [1,1];
mu(2,:) = [10,1];
mu(3,:) = [1,10];

% original covariance matrices used to generate GMM

Sigma(:,:,1) = [2,-0.5;-0.5,1];
Sigma(:,:,2) = [2,0.8;0.8,4];
Sigma(:,:,3) = [1,0.9;0.9,3];

end