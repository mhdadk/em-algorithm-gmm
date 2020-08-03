%{

get the original GMM parameters.

Inputs: None

Outputs:

N - The number of vectors to sample from the GMM. This is a scalar.

comp_PMF_true - A 1 x M vector of probabilities representing the
component probabilities of the GMM. M is the number of components in the
GMM

mu_true - A M x d array containing the invidual mean vectors for each
mixture component. d is the dimensionality of each mean vector.
Each row of this array is a mean vector.

Sigma_true - A d x d x M array containing the covariance matrices for each
mixture component. Each page (3rd dimension) of this array contains an
individual covariance matrix.

%}

[N,comp_PMF_true,mu_true,Sigma_true] = get_GMM_parameters();

%{

sample N random vectors from the GMM.

Inputs:

N - The number of vectors to sample from the GMM. This is a scalar.

comp_PMF_true - A 1 x M vector of probabilities representing the
component probabilities of the GMM. M is the number of components in the
GMM

mu_true - A M x d array containing the invidual mean vectors for each
mixture component. Each row of this array is a mean vector.

Sigma_true - A d x d x M array containing the covariance matrices for each
mixture component. Each page (3rd dimension) of this array contains an
individual covariance matrix.

Outputs:

X - A N x d array containing the random vectors that were sampled from the
GMM. N is the number of vectors to be sampled and d is the dimensionality
of each random vector.

comp_idx - A N x 1 array containing the index of the component that the
corresponding random vector was sampled from. This can be used to estimate
the component probabilities so that they are compared to comp_PMF_true.
For example: sum(comp_idx == 1)/N

Note that it is assumed that beyond this point, only the samples in X
are available, and that the true parameters above are not.

%}

[X,comp_idx] = sample_GMM(N,comp_PMF_true,mu_true,Sigma_true);

%{
 
get the initial estimates of the GMM parameters.

Inputs: None

Outputs:

comp_PMF - A 1 x M vector of probabilities representing the initial
guesses for the component probabilities of the GMM. M is the number of
components in the GMM.

mu - A M x d array containing the initial guesses for the invidual mean
vectors for each mixture component. d is the dimensionality of each mean
vector. Each row of this array is a mean vector.

Sigma - A d x d x M array containing the initial guesses for the
covariance matrices for each mixture component. Each page (3rd dimension)
of this array contains an individual covariance matrix.

%}

[comp_PMF,mu,Sigma] = get_init_est();

% store all parameter estimates for future plotting

comp_PMF_save = [];
mu_save = [];
Sigma_save = [];

%{

initialize the value that will be used to detect the convergence of the EM
algorithm. This value was chosen to start the while loop later in the
program. epsilon will then be computed using the old and new parameter
estimates.

Note that it is assumed here that the true parameter values are not
available. If they were available, then epsilon could have been initialized
as follows:

epsilon = sum(abs(comp_PMF_true - comp_PMF),'all') + ...
          sum(abs(mu_true - mu),'all') + ...
          sum(abs(Sigma_true - Sigma),'all');

%}

epsilon = 10;

% record epsilons for future plotting

epsilon_save = [];
    
% count the number of iterations required for convergence    

num_iter = 0;
    
%{ 

while the sum of absolute differences between the old parameters and the
new parameters is greater than a certain threshold
        
%}

while epsilon > 1e-5
    
    num_iter = num_iter + 1;
    
    epsilon_save = [epsilon_save,epsilon];
    
    comp_PMF_save = cat(3,comp_PMF_save,comp_PMF);
    
    mu_save = cat(3,mu_save,mu);
    
    Sigma_save = cat(4,Sigma_save,Sigma);
    
    % store the old parameters for future comparison
    
    comp_PMF_old = comp_PMF;
    
    mu_old = mu;
    
    Sigma_old = Sigma;
    
    %{
    
    initialize the array that will store each of the gamma_k vectors for
    each component. The gamma_k vector contains values for gamma_jk for
    all j
    
    %}
    
    gamma_k = zeros(length(X),length(comp_PMF));
    
    %{
    
    initialize the array that will store each of the N_k values for each
    component
    
    %}
    
    N_k = zeros(1,length(comp_PMF));
    
    %{
    
    compute parameter estimates for each component in the GMM
    
    %}
    
    for k = 1:length(comp_PMF)
        
        %{
        
        get all gamma_jk values for all j and the k^th component.
        This is a N x 1 array, where N is the number of sampled vectors
        from the GMM. Store this column vector in the k^th column of the
        gamma_k array.
        
        %}

        gamma_k(:,k) = get_gamma_k(k,X,comp_PMF,mu,Sigma);
        
        %{
        
        compute N_k for the k^th component by summing up all gamma_jk
        
        %}
              
        N_k(k) = sum(gamma_k(:,k));

        % compute the mean vector estimate for the k^th component

        mu(k,:) = (1/N_k(k))*sum(X.*gamma_k(:,k),1);
        
        %{
        
        compute the covariance matrix estimate  for the k^th component of
        the GMM as an inner product (sum) of outer products, which is
        simply matrix multiplication.
        
        transpose((X - mu(k,:)).*gamma_k(:,k)) is a d x N matrix, where d
        is the dimensionality of the sampled vectors and N is the number
        of sampled vectors, which means that X - mu(k,:) is a N x d
        matrix. This means that the multiplication of these two matrices
        yields a d x d matrix. This is equivalent to the sum of the
        individual outer products resulting from the multiplication of the
        column vectors in the d x N matrix by the row vectors in the N x d
        matrix.
        
        %}
        
        Sigma(:,:,k) = (1/N_k(k)) *...
                       (transpose((X - mu(k,:)).*gamma_k(:,k)) *...
                       (X - mu(k,:)));

        % compute the component probability estimates

        comp_PMF(k) = N_k(k)/N;
        
    end
    
    % check for convergence by comparing the old parameter values to the
    % new parameter values
    
    epsilon = sum(abs(comp_PMF_old - comp_PMF),'all') + ...
              sum(abs(mu_old - mu),'all') + ...
              sum(abs(Sigma_old - Sigma),'all');
    
end

% plot epsilon against iteration number

figure
plot(epsilon_save)
set(gcf,'color','w')
xlim([1,inf])
xlabel('Iteration Number')
ylabel('Sum of absolute differences')

%{

create an animation showing the shape of the PDF for different GMM's for
each iteration of the EM algorithm. This shows how the EM algorithm
converges. Additionally, save this animation as 'GMM.gif' in your working
directory.

%}

frames = save_GIF(X,num_iter,mu_save,Sigma_save,comp_PMF_save);

% play this animation 3 times at 2 frames per second

figure
movie(gcf,frames,3,2)

figure
plot(X(:,1),X(:,2),'+')
set(gcf,'color','w')
lim = axis;
axis([lim(1)-1,lim(2)+1,lim(3)-1,lim(4)+1])
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')

% create a gmdistribution object and get its PDF

gm = gmdistribution(mu,Sigma,comp_PMF);

gm_PDF = @(x,y)reshape(pdf(gm,[x(:) y(:)]),size(x));

% plot the contours for the PDF of the GMM on the same figure

hold on
fc = fcontour(gm_PDF);
hold off
