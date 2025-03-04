function varargout = EM_gamma(D,y,options)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Expectation-Maximization method for training of BMA model if predictive %
% PDF of ensemble members is a gamma distribution with constant or        %
% non-constant (heteroscedastic) group or individual variances            %
% This function is OBSOLETE, try EM_bma or MODELAVG toolbox instead       %
%                                                                         %
% SYNOPSIS: [beta,sigma,loglik,t] = EM_gamma(D,y);                        %
%           [beta,sigma,loglik,t] = EM_gamma(D,Y,options);                %
%  where                                                                  %
%   D         [input] nxK matrix with forecasts of ensemble members       %
%   y         [input] nx1 vector with verifying observations              %
%   options   [input] OPT: structure with BMA algorithmic variables       %
%   beta      [outpt] 1xK vector with maximum likelihood BMA weights      %
%   sigma     [outpt] 1xK vector of BMA standard deviations - or          %
%   c         [outpt] (1 x K)-vector of multipliers c: sigma = c*|D|      %
%   loglik    [outpt] BMA log-likelihood                                  %
%   it        [outpt] number of iterations to reach convergence           %
%                                                                         %
% (c) Written by Jasper A. Vrugt, Feb 2006                                %
% University of California Irvine                                         %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

if nargin < 2  
    error(['EM ERROR: TooFewInputs - Function EM_gamma requires at ' ...
        'least two input arguments.']); 
end
if nargin < 3
    warning(['EM WARNING: Variance option not specified - resort to ' ...
        'options.VAR = ''1''']);
    options.VAR = '1'; 
end

[n,K] = size(D);                              % Matrix of ensemble forecasts
beta = ones(1,K)/K;                           % Initial values weights
mu = abs(D);
z_t = zeros(n,K);                             % Initial values latent variables
loglik_t = -Inf; err = 1; t = 0; max_t = 1e4; % Settings/constraints while loop
% 
switch options.VAR % Initial values standard devs 
    case '1'
        sigma = std(y)*rand;
    case '2'
        sigma = std(y)*rand(1,8);
    case '3'
        c = 0.5*rand;
    case '4'
        c = 0.2*rand(1,K);
end

while ( max(err) > 1e-4 ) && ( t < max_t )       % Until ... do    
    loglik = loglik_t; z = z_t;                           % Copy z and loglik 
    % Compute A and B
    switch options.VAR
        case {'1','2'}
            % do nothing
        case {'3','4'}
            sigma = bsxfun(@times,c,abs(D));
    end
    A = mu.^2./sigma.^2; B = sigma.^2./mu;
    for k = 1:K                                           % EXPECTATION STEP
        z_t(:,k) = beta(k)*gampdf(y,A(1:n,k),B(1:n,k));   % Update latent variables
    end
    loglik_t = sum(log(sum(z_t,2)));                      % Log-likelihood BMA model
    z_t = bsxfun(@rdivide,z_t,sum(z_t,2));                % Normalize latent variables   
    
    beta_t = sum(z_t)/n;                                  % MAXIMIZATION STEP
    % Determine sigma using Nelder-Mead simplex algorithm
    switch options.VAR
        case {'1','2'}
            sigma2_t = fminsearch(@(sigma) ell_gamma(sigma,beta_t,mu,y,D,n,K,options),sigma).^2;
        case {'3','4'}
            c_t = fminsearch(@(c) ell_gamma(c,beta_t,mu,y,D,n,K,options),c);
            sigma2_t = bsxfun(@times,c_t,abs(D)).^2;      % matrix: works if c_t is scalar or 1xK vector
            c = c_t;                                      % make sure we use last value(s) of c_t 
    end
    err(1) = max(abs(beta_t - beta));                     % Convergence: weights
    err(2) = mean(max(abs(log(sigma2_t./sigma.^2))));     % Convergence: variance(s)
    err(3) = max(max(abs(z - z_t)));                      % Convergence: latent variables
    err(4) = max(abs(loglik - loglik_t));                 % Convergence: log-likelihood
    loglik_t
    beta = beta_t; sigma = sqrt(sigma2_t);                % Update BMA weights and variances
    t = t + 1;                                            % Iteration counter
end
% End while loop

% Populate return argument
switch options.VAR
    case {'1','2'}
        varargout = {beta , sigma , loglik , t};
    case {'3','4'}
        varargout = {beta , c , loglik , t};
end        

end

% Secondary function
function ell = ell_gamma(x,beta,mu,Y,D,n,K,options)

if any(x < 0)
    ell = inf; return
end

switch options.VAR
    case '1' % constant: sigma = x (scalar)
        sigma = x*ones(n,K);
    case '2' % constant: sigma = x (vector)
        sigma = repmat(x,n,1);
    case '3' % nonconstant: c = x (scalar) 
        sigma = x*abs(D);
    case '4' % nonconstant: c = x (vector c1...cK)
        sigma = repmat(x,n,1) .* abs(D);
end
A = mu.^2./sigma.^2; B = sigma.^2./mu;
ell = - sum(log(gampdf(repmat(Y,1,K),A,B) * beta'));

end
