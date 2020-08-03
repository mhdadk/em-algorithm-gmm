function [X,comp_idx] = sample_GMM(N,comp_PMF,mu,Sigma)

% for reproducibility

rng('default')  

% create a gmdistribution object to sample the random vectors from

gm = gmdistribution(mu,Sigma,comp_PMF);

% sample N random vectors from the GMM

[X,comp_idx] = random(gm,N);

end