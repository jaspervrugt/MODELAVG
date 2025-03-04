function [beta,sigma,loglik,t] = EM_normal(D,y,options)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% Expectation-Maximization method for training of BMA model if predictive %
% PDF of ensemble members is a normal distribution with constant          %
% individual variances                                                    %
% This function is OBSOLETE, try EM_bma or MODELAVG toolbox instead       %
%                                                                         %
% SYNOPSIS: [beta,sigma,loglik,t] = EM_normal(D,y);                       %
%           [beta,sigma,loglik,t] = EM_normal(D,Y,options);               %
%  where                                                                  %
%   D         [input] nxK matrix with forecasts of ensemble members       %
%   y         [input] nx1 vector with verifying observations              %
%   options   [input] OPT: structure with BMA algorithmic variables       %
%   beta      [outpt] 1xK vector with maximum likelihood BMA weights      %
%   sigma     [outpt] 1xK vector of BMA standard deviations               %
%   loglik    [outpt] BMA log-likelihood                                  %
%   it        [outpt] number of iterations to reach convergence           %
%                                                                         %
% (c) Written by Jasper A. Vrugt, Feb 2006                                %
% University of California Irvine                                         %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

if nargin < 2  
    error(['EM ERROR: TooFewInputs - Function EM_normal requires at ' ...
        'least two input arguments.']); 
end
if nargin < 3
    warning(['EM WARNING: Variance option not specified - resort to ' ...
        'options.VAR = ''1''']);
    options.VAR = '1'; 
end

[n,K] = size(D);                              % Matrix of ensemble forecasts
beta = ones(1,K)/K;                           % Initial values weights
sigma = std(y)*rand(1,K);                     % Initial values standard devs 
z_t = zeros(n,K);                             % Initial values latent variables
loglik_t = -Inf; err = 1; t = 0; max_t = 1e4; % Settings/constraints while loop

while ( max(err) > 1e-6 ) && ( t < max_t )       % Until ... do    
    loglik = loglik_t; z = z_t;                           % Copy z and loglik 
    for k = 1:K                                           % EXPECTATION STEP
        z_t(:,k) = beta(k)*normpdf(y,D(:,k),sigma(k));    % Update latent variables
    end
    loglik_t = sum(log(sum(z_t,2)));                      % Log-likelihood BMA model
    z_t = bsxfun(@rdivide,z_t,sum(z_t,2));                % Normalize latent variables   
    
    beta_t = sum(z_t)/n;                                  % MAXIMIZATION STEP
    sigma2_t = sum(z_t.*bsxfun(@minus,D,y).^2)./sum(z_t);   
    if strcmp(options.VAR,'1')                            % If common constant variance
        sigma2_t = mean(sigma2_t)*ones(1,K);              % Use mean value   
    end
    
    err(1) = max(abs(beta_t - beta));                     % Convergence: weights
    err(2) = max(abs(log(sigma2_t./sigma.^2)));           % Convergence: variance(s)
    err(3) = max(max(abs(z - z_t)));                      % Convergence: latent variables
    err(4) = max(abs(loglik - loglik_t));                 % Convergence: log-likelihood
    
    beta = beta_t; sigma = sqrt(sigma2_t);                % Update BMA weights and variances
    t = t + 1;                                            % Iteration counter
end                                             % End while loop

end
