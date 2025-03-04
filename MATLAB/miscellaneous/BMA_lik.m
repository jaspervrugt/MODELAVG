function [ell,L] = BMA_lik(x,D,Y,options)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% This function calculates log-likelihood of BMA mixture distribution     %
%                                                                         %
% SYNOPSIS: [ell,L] = BMA_lik(x,D,Y,options)                              %
%  where                                                                  %
%   x         [input] 1xd vector with values of BMA parameters            %
%   D         [input] nxK matrix with forecasts of ensemble members       %
%   Y         [input] nxK matrix with verifying data (= K copies)         %
%   options   [input] structure with BMA algorithmic variables            %
%   ell       [outpt] scalar log-likelihood BMA mixture distribution      %
%   L         [outpt] nxK matrix likelihoods of ensemble members          %
%                                                                         %
% (c) Written by Jasper A. Vrugt, Feb 2012                                %
% University of California Irvine                                         %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

if nargin < 4
    error('MODELAVG:BMA_calc:TooFewInputs',['Requires at least ' ...
        'four input arguments.']);
end

[n,K] = size(D);            % # forecasts and # ensemble members
beta = x(1:K)';             % Unpack weights of member's conditional pdf

switch options.VAR  % VARIANCE OPTION -> (n x K)-matrix "S" 
                    % with forecast standard deviations
    case {'1'} % 1: common constant variance
        S = x(K+1) * ones(n,K);
    case {'2'} % 2: individual constant variance
        S = bsxfun(@times,x(K+1:2*K),ones(n,K));
    case {'3'} % 3: common non-constant variance
        c = x(K+1); S = c * abs(D); 
    case {'4'} % 4: individual non-constant variance
        c = x(K+1:2*K); S = bsxfun(@times,c,abs(D));
    otherwise
        error('MODELAVG:BMA_calc','Unknown variance option');
end
S = max(S,eps); % each element (n x K)-matrix sigma >= 2.22e-16

switch options.PDF % CONDITIONAL DISTRIBUTION
                    
    case {'normal'}     % NORMAL: mean "mu = D" and std. "S"
        if strcmp(options.CPU,'yes')    % --> fastest computation
            L = exp(-1/2*(( Y - D)./S).^2)./(sqrt(2*pi).*S);
        else                            % --> built-in computation
            L = pdf('normal',Y,D,S);
        end
    
    case {'lognormal'}  % LOGNORMAL: mean "mu" and std. "S"
        logn_impl = 2;
        switch logn_impl
            case 1 % original implementation: works only constant variance
                mu = log(abs(D)) - S.^2/2; % --> mean is now equal to "D"
            case 2 % TEST
                sigma2_fxs = S.^2;
                % --> we define directly sigma2_fxs not S2.
                % A = sigma2_fxs./(D.^2);
                % ( exp(S2) - 1 ) = A;
                % exp(S2) = A + 1;
                % S2 = log(A+1);
                % Mu = log(abs(D)) - S2/2;
                % Mu = log(abs(D)) - log(A+1)/2;
                % Back-substitute
                % sigma2_fxs = exp(2*Mu + S.^2) * ...   % Variance each fk
                %                  (exp(S.^2) - 1); 
                % sigma2_fxs = exp(2*(log(abs(D)) - ...
                %                  log(A+1)/2) + log(A+1)) * ...
                %                  (exp(log(A+1)) - 1);          
                % sigma2_fxs = D.^2 * A;                     
                A = sigma2_fxs./(D.^2);
                S2 = log(A+1); S = sqrt(S2);
                mu = log(abs(D)) - S2/2;
        end        
        % Compute BMA likelihood 
        if strcmp(options.CPU,'yes')    % --> fastest computation
            L = exp(-1/2*((log(Y) - mu)./S).^2) ./ (sqrt(2*pi).*Y.*S);
        else                            % --> built-in computation
            L = pdf('logn',Y,mu,S);
        end   

    case {'tnormal'}    % TRUNCATED NORMAL: mean µ = |D| & std. S 
                        % given a = 0 and b = inf
        % Cannot set D equal to mean of truncated normal 
        % --> no closed form solution for mu
        % Rather we set D equal to mode (peak of tnormal); then mu = abs(D) 
        mu = abs(D); a = 0; b = inf; 
        L = tnormpdf(Y,mu,S,a,b);
        
    case {'gen_normal'} % GENERALIZED NORMAL: mean µ = D, std. S, shape tau
        d = size(x,2); % Number of members of vector x
        switch options.TAU
            case {'1'} % 1: common value of tau
                tau = x(d);
            case {'2'} % 2: individual tau
                tau = x(d-K+1:d);
        end
        % Note: if A: (N x p) and b: (1 x p); A.*b allowed in MATLAB 2016
        % Note: if A: (N x p) and b: (1 x p); A.^b allowed in MATLAB 2016
        % Note: if A: (N x p) and b: (1 x p); b./A allowed in MATLAB 2016
        calc_method = 2;
        switch calc_method
            case 1
                L = tau./(2*S.*repmat(gamma(1./tau),n,1)) .* ...
                    exp(-(abs(Y - D)./S).^tau);
            case 2
                Tau = repmat(tau,n,1);
                L = Tau./(2*S.*gamma(1./Tau)) .* ...
                    exp(-(abs(Y - D)./S).^Tau);
            case 3
                L = nan(n,K); gTau = gamma(1./tau);
                for Kk = 1:K
                    L(1:n,Kk) = tau(Kk)./(2*S(1:n,Kk).*gTau(Kk)) .* ...
                        exp(-(abs(Y(1:n,1) - D(1:n,Kk))./ ...
                        S(1:n,Kk)).^tau(Kk));
                end
        end

    case {'gamma'}      % GAMMA: shape "A" and scale "B"
        mu = abs(D); A = mu.^2./S.^2; B = S.^2./mu;
        if strcmp(options.CPU,'yes')    % --> fastest computation
            z = Y./B;
            u = (A-1).*log(z) - z - gammaln(A);
            L = exp(u)./B;
        else                            % --> built-in computation
            L = pdf('gamma',Y,A,B);
        end

    case {'weibull'}    % WEIBULL: wblpdf(Lambda,Kk) 
                        % Scale "Lambda", shape "Kk"
        wbl_impl = 2;
        switch wbl_impl
            case 1 % Shape parameter treated sigma & Lambda (scale paramtr)
                % derived so that mean of Weibull PDFs equals D
                Kk = S; Lambda = abs(D)./gamma(1+1./Kk);
                % Lambda cannot be zero
                Lambda = max(Lambda,realmin);
            case 2  % Scale parameter (Lambda) treated as sigma (= corrct!)
                    % and shape parameter (Kk) set so that PDF has µ = D
                Lambda = S; X = abs(D)./Lambda;
                c = sqrt(2*pi)/exp(1) - gamma(1.461632);
                Lx = @(x) log((x+c)/sqrt(2*pi)); A = Lx(X);
                B = A ./ real(lambertw(A/exp(1))) + 1/2;
                Kk = 1./(B - 1);
                % --> this implementation is not stable yet; requires the
                % inverse of gamma function: y = gam(x); x = gam^{-1}(y)
                % can try linear interpolation:
                % xx = linspace(eps,150,1e5); 
                % fxx = gamma(xx); Y = interp1(fxx,xx,X);
                % Kk = 1./(Y - 1);
        end
        if strcmp(options.CPU,'yes')    % --> fastest computation
            L = Kk./Lambda .* (Y./Lambda).^(Kk-1) .* exp(-(Y./Lambda).^Kk);
        else
            L = pdf('wbl',Y,Lambda,Kk);
        end

    case {'gev'}    % GEV distribution: gevdf(xi,S,mu), 
                    % shape xi, scale S, location mu
        d = size(x,2); % Number of members of vector x
        switch options.TAU
            case {'1'} % 1: common value of tau
                xi = x(d) * ones(n,K);
            case {'2'} % 2: individual tau
                xi = repmat(x(d-K+1:d),n,1);
        end
        g1 = gamma(1 - xi);
        mu = D - (g1-1) .* S./xi;
        % If xi == 0 --> Mu = abs(D) - S.* Eul_constant
        kk = find(abs(xi(1,1:K)) < eps);
        Eul_constant = 0.577215664901532860606512090082;
        mu(1:n,kk) = D(1:n,kk) - Eul_constant * S(1:n,kk);
        L = pdf('gev',Y,xi,S,mu);

    case {'gpareto'}    % Generalized Pareto: gppdf(kk,S,Theta)
                        % shape kk, scale S, location mu
        d = size(x,2);  % # members vector x
        switch options.TAU
            case {'1'} % 1: common value of tau
                kk = x(d) * ones(n,K);
            case {'2'} % 2: individual tau
                kk = repmat(x(d-K+1:d),n,1);
        end
        Theta = abs(D) - S./(1-kk); % Scale
        if strcmp(options.CPU,'yes')    % --> fastest computation
            L = zeros(n,K);
            z = (Y - Theta)./S; z(z<0) = Inf; % max drops NaNs
            j = (abs(kk) < eps); 
            L(j) = exp(-z(j));
            t = z.*kk;
            jj = (t<=-1 & kk < -eps);
            % t(jj) = 0; % temporarily silence warnings from log1p
            j = ~j;
            % L(j) = exp((-1 - 1./kk(j)).*log1p(t(j))); 
            % (1 + z.*k).^(-1 - 1./k)
            L(j) = 1./S(j) .* ( 1 + t(j) ).^ ( - 1./kk(j) - 1 );               
            L(jj) = 0;
        else
            L = pdf('gp',Y,kk,S,Theta);
        end

end

% Replace likelihood with zero if D < 0
%idx = D < 0; L(idx) = 0;    % Always true for gamma (not needed), 
% %                             % but with this statement also imposed for 
% %                             % the normal & gen_normal conditional PDFs
lik = L*beta + realmin;     % (nx1)-vector of likelihoods BMA model at y
ell = sum(log(lik));        % Return log-likelihood of BMA model

end