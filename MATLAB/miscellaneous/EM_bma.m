function varargout = EM_bma(D,y,options,errtol,maxit)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Expectation-Maximization method for training of BMA model when the    %%
%% conditional PDF is a normal, lognormal, truncated normal or gamma     %%
%% distribution with constant or nonconstant variance                    %%
%%                                                                       %%
%% SYNOPSIS: [beta,sigma,loglik,it] = EM_bma(D,y,options)                %%
%%           [beta,sigma,loglik,it] = EM_bma(D,y,options,errtol)         %%
%%           [beta,sigma,loglik,it] = EM_bma(D,y,options,errtol,maxit)   %%
%%  where                                                                %%
%%   D         [input]  REQUIRED: n x K matrix of ensemble forecasts     %%
%%   y         [input]  REQUIRED: n x 1 vector of verifying observations %%
%%   options   [input]  REQUIRED: structure of BMA algorithmic variables %%
%%   errtol    [input]  OPTIONAL: error tolerance: default = 1e-4        %%
%%   maxit     [input]  OPTIONAL: maximum # iterations: default = 10000  %%
%%   beta      [output] (1 x K)-vector of BMA weights                    %%
%%   sigma     [output] (1 x K)-vector of BMA standard deviations - or   %%
%%   c         [output] (1 x K)-vector of multipliers c: sigma = c*|D|   %%
%%   loglik    [output] BMA log-likelihood                               %%
%%   it        [output] number of iterations to reach convergence        %%
%%                                                                       %%
%% (c) Written by Jasper A. Vrugt, Feb 2006                              %%
%% Los Alamos National Laboratory                                        %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

if nargin < 5, maxit = 1e4; end
if nargin < 4, errtol = 1e-4; end
if nargin < 3
    error(['EM_bma ERROR: TooFewInputs - Function EM_bma requires at ' ...
        'least three input arguments.']);
end
if ~sum(strcmp(options.PDF,{'normal','lognormal','tnormal','gamma'}))
    error(['EM_bma ERROR: Code only works if conditional PDF is a ' ...
        'normal, lognormal, truncated normal or gamma distribution']);
end
if ~sum(strcmp(options.VAR,{'1','2','3','4'}))
    error(['EM_bma ERROR: Code only works if variance treatment VAR is' ...
        'equal to ''1'', ''2'', ''3'' or ''4'' ']);
end
warning off                             % for fminsearch
[n,K] = size(D);                        % Ensemble forecasts
beta = ones(1,K)/K;                     % Initial weights
Y = repmat(y,1,K);                      % Duplicate y model K times
z_it = zeros(n,K);                      % Initial latent variables
ell_it = -inf; err = 1; it = 0;         % Constraints while loop

switch options.VAR                      % Variance treatment forecast PDF?
    case '1', s = std(y)*rand;
    case '2', s = std(y)*rand(1,K);
    case '3', c = 0.5*rand;
    case '4', c = 0.2*rand(1,K);
end

while (max(err) > errtol) && (it < maxit)   % Until ... do
    ell = ell_it; z = z_it;                     % Copy loglik and z
    switch options.VAR                          % Compute sigma
        case {'3','4'}, s = bsxfun(@times,c,abs(D));
    end
    %% EXPECTATION STEP
    [~,L] = BMA_loglik(s,beta,Y,D,n,K,options,0);   % Compute likelihood
    z_it = bsxfun(@times,beta,L);                   % Latent variable
    ell_it = sum(log(sum(z_it,2)));                 % Log-likelihood BMA model
    z_it = bsxfun(@rdivide,z_it,sum(z_it,2));       % Norm. latent variables
    %% MAXIMIZATION STEP
    beta_it = sum(z_it)/n;                          % New weights
    switch options.VAR                              % Compute new sigma2's
        case {'1','2'}
            switch options.PDF
                case 'normal'
                    s2_it = sum(z_it.*bsxfun(...        % sigma2 estimate
                        @minus,D,Y).^2)./sum(z_it);
                    if strcmp(options.VAR,'1')          % Comm. const. var.
                        s2_it = mean(s2_it)*ones(1,K);  % Copy mean value
                    end
                otherwise % Nelder-Mead simplex
                    s2_it = fminsearch(@(s) ...         % s2_it from minimization
                        BMA_loglik(s,beta_it,Y,...
                        D,n,K,options,1),s).^2;
            end
        case {'3','4'}
            c_it = fminsearch(@(c) BMA_loglik(c,...     % c_it from minimization     
                beta_it,Y,D,n,K,options,1),c);
            s2_it = bsxfun(@times,c_it,abs(D)).^2;      % matrix: c_t scalar or 1xK vector
            c = c_it;                                   % use last value(s) of c_t
    end
    %% CONVERGENCE DIAGNOSTICS
    err(1) = max(abs(beta_it - beta));              % Conv. weights
    err(2) = mean(max(abs(log(s2_it./s.^2))));      % Conv. variance(s)
    err(3) = max(max(abs(z - z_it)));               % Conv. latent variables
    err(4) = max(abs(ell - ell_it));                % Conv. log-likelihood
    beta = beta_it; s = sqrt(s2_it);                % update weights/vars
    it = it + 1;                                    % Iteration counter
end                                             % End while loop

if strcmp(options.VAR,'1')  % Duplicate K times sigma or multiplier c
    s = s*ones(1,K);
elseif strcmp(options.VAR,'3')
    c = c*ones(1,K);
end

switch options.VAR          % Group return argument
    case {'1','2'}, varargout = {beta , s , ell , it};
    case {'3','4'}, varargout = {beta , c , ell , it};
end

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local subroutines
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ell,L] = BMA_loglik(x,beta,Y,D,n,K,options,flag)
% Make sure that weights and variances/multipliers do not go negative
if sum(any(x < 0))
    ell = inf; return
end
% Check implementation
if flag == 1 % fminsearch implementation
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
else % sigma is already known
    switch options.VAR
        case '1' % constant: sigma = x (scalar)
            sigma = x * ones(n,K);
        case '2' % constant: sigma = x (vector)
            sigma = repmat(x,n,1); 
        otherwise
            sigma = x;
    end
end

% Compute likelihood
switch options.PDF
    case 'normal'
        L = exp(-1/2*((Y-D)./sigma).^2)./(sqrt(2*pi).*sigma);
    case 'lognormal'
        A = sigma.^2./(D.^2);
        S2 = log(A+1); S = sqrt(S2);
        mu = log(abs(D)) - S2/2;
        L = exp(-1/2*((log(Y) - mu)./S).^2) ./ ...
            (sqrt(2*pi).*Y.*S);
    case 'gamma'
        mu = abs(D); A = mu.^2./sigma.^2; B = sigma.^2./mu;
        Z = Y./B;
        U = (A-1).*log(Z) - Z - gammaln(A);
        L = exp(U)./B;
    case 'tnormal'
        mu = max(D,eps); a = 0; b = inf;
        L = tnormpdf(Y,mu,sigma,a,b);
end
% Return BMA log-likelihood
ell = - sum(log(L * beta'));

end