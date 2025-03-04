function [rnd,mu_mix,sigma2_mix] = BMA_draw(x,D,N,options)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% This function calculates log-likelihood of BMA mixture distribution     %
%                                                                         % 
% SYNOPSIS: [rnd,mu_mix,sigma2_mix] = BMA_draw(x,D,N,options)             %
%  where                                                                  %
%   x          [input] 1xd vector with max. likelihood BMA parameters     %
%   D          [input] nxK matrix with forecasts of ensemble members      %
%   N          [input] number samples drawn BMA mixture distribution      %
%   options    [input] structure with BMA algorithmic variables           %
%   rnd        [outpt] nxN matrix samples drawn BMA mixture dstribution   %
%   mu_mix     [outpt] nx1 vector mean of BMA mixture distribution        %
%   sigma2_mix [outpt] nx1 vector variance BMA mixture distribution       %
%                                                                         %
% (c) Written by Jasper A. Vrugt, Feb 2022                                %
% University of California Irvine                                         %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

if nargin < 4
    error('MODELAVG:BMA_sample:TooFewInputs',['Requires at least four ' ...
        'input arguments.']);
end

[n,K] = size(D);            % Number of forecasts and number of ensemble members
beta = x(1:K)';             % Unpack weights of member's conditional pdf
mu_mix = D * beta;          % Now determine mean of the mixture (= exact)
%mu_mix = abs(D) * beta;     % Should we take absolute value?
rnd = nan(n,N);             % Initialize the random draw from predictive PDF

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
        error('MODELAVG:BMA_sample','Unknown variance option');
end
S = max(S,eps);     % each element (n x K)-matrix S >= 2.22e-16
S2 = S.^2;          % variance (n x K)-matrix

% Mixture variance: See Equation 3.129 of
% https://www.value-at-risk.net/mixtures-of-distributions/
% sigma2_mix = (sigma2 + Mu.^2) * beta - mu_mix.^2;

% Compute parameters and exact variance of BMA forecast density
switch options.PDF
    case {'normal'}     % NORMAL: mean "Mu = D" and variance "S2"
        Mu = D;
        sigma2_mix = (S2 + Mu.^2) * beta - mu_mix.^2;

    case {'lognormal'}  % LOGNORMAL: mean "Mu" and variance "S2"
        mu_fxs = D;                                         % Mean each fk
        logn_impl = 2;
        switch logn_impl
            case 1 % lognrnd(Mu,S): original implementation: works only for a constant variance
                Mu = log(abs(D)) - S2/2; % --> mean is now equal to "D"
                % M = exp(Mu + S^2/2)
                % V = exp(2*Mu + S^2) * (exp(S^2) - 1)
                % sigma2_fxs = exp(2*Mu + S.^2)*(exp(S.^2) - 1);                         % Variance each fk
                % sigma2_fxs = exp(2*(log(abs(D)) - S.^2/2) + S.^2)*(exp(S.^2) - 1);     % Variance each fk
                % sigma2_fxs = exp(2*log(abs(D)) )*(exp(S.^2) - 1);                      % Variance each fk
                % sigma2_fxs = D.^2 .* (exp(S.^2) - 1);                                  % Variance each fk
                sigma2_fxs = D.^2 .* ( exp(S2) - 1 );                                   % Variance each fk
                sigma2_mix = (sigma2_fxs + mu_fxs.^2) * beta - mu_mix.^2; % CORRECT
            case 2 % TEST lognrnd(Mu,Sigma)
                sigma2_fxs = S2; % Also called "V" for variance 
                % Given by MATLAB:
                % S = sqrt(log(V/M^2 + 1)); V = variance of distribution
                % Mu = log(M^2 / sqrt(V+M^2)); M = mean of distribution
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
                Mu = log(abs(D)) - S2/2;
                sigma2_mix = (sigma2_fxs + mu_fxs.^2) * beta - mu_mix.^2; % CORRECT
        end

    case {'gamma'}      % GAMMA: shape "A" and scale "B"
        Mu = abs(D); A = Mu.^2./S2; B = S2./Mu;
        sigma2_mix = (S2 + Mu.^2) * beta - mu_mix.^2;       % CORRECT

    case {'gen_normal'} % GENERALIZED NORMAL: mean "Mu = D", std. "S", shape "tau"
        d = size(x,2);
        switch options.TAU
            case {'1'} % 1: common value of tau
                tau = x(d);
            case {'2'} % 2: individual tau
                tau = x(d-K+1:d);
        end
        Mu = D;
        sigma2 = S2 .* gamma(3./repmat(tau,n,1)) ./ ...
            gamma(1./repmat(tau,n,1));                           % Variance of fk's
        sigma2_mix = (sigma2 + Mu.^2) * beta - mu_mix.^2;        % CORRECT

    case {'weibull'}    % WEIBULL: scale "lambda > 0" and shape "k > 0"
        wbl_impl = 2;
        switch wbl_impl
            case 1 % The shape parameter is treated as sigma and Lambda (scale parameter)
                % derived so that mean of Weibull PDFs equals D
                Kk = S; Lambda = abs(D)./gamma(1+1./Kk);
                % Lambda cannot be zero
                Lambda = max(Lambda,realmin);
            case 2 % The scale parameter (Lambda) is treated as sigma (= correct!)
                % and the shape parameter (Kk) is set so that PDF has mean equal to D
                Lambda = S; X = abs(D)./Lambda;
                c = sqrt(2*pi)/exp(1) - gamma(1.461632);
                Lx = @(x) log((x+c)/sqrt(2*pi)); A = Lx(X);
                B = A ./ real(lambertw(A/exp(1))) + 1/2;
                Kk = 1./(B - 1);
                % --> this implementation is not stable yet
        end
        Mu = abs(D);
        sigma2 = Lambda.^2 .* (gamma(1 + 2./Kk) - gamma(1+1./Kk).^2);   % Variance of fk's
        sigma2_mix = (sigma2 + Mu.^2) * beta - mu_mix.^2;               % CORRECT?
    
    case {'gev'}        % GEV(xi,S,Mu)
        % scale parameter "S > 0", and location parameter "Mu"
        d = size(x,2); % Number of members of vector x
        switch options.TAU
            case {'1'} % 1: common value of tau
                xi = x(d) * ones(n,K);
            case {'2'} % 2: individual tau
                xi = repmat(x(d-K+1:d),n,1);
        end
        g = @(a,xi) gamma(1 - a*xi);
        Mu = D - (g(1,xi) - 1) .* S./xi;
        % If xi == 0 --> Mu = abs(D) - S.* Eul_constant
        kk = find(abs(xi(1,1:K)) < eps);
        Eul_constant = 0.577215664901532860606512090082;
        Mu(1:n,kk) = D(1:n,kk) - S(1:n,kk) * Eul_constant;
        % built-in function to get mean and variance of each GEV component
        [mu_fxs,sigma2_fxs] = gevstat(xi,S,Mu);
        % mu_fxs is equal to abs(D), as requested by setting mean
        sigma2_mix = (sigma2_fxs + mu_fxs.^2) * beta - mu_mix.^2;       % CORRECT

    case {'gpareto'}     % GP(xi,S,Theta)
        d = size(x,2); % Number of members of vector x
        switch options.TAU
            case {'1'} % 1: common value of tau
                kk = x(d) * ones(n,K);
            case {'2'} % 2: individual tau
                kk = repmat(x(d-K+1:d),n,1);
        end
        Theta = abs(D) - S./(1-kk);
        % built-in function to get mean and variance of each GP component
        [mu_fxs,sigma2_fxs] = gpstat(kk,S,Theta);
        % mu_fxs is equal to abs(D), as requested by setting mean
        sigma2_mix = (sigma2_fxs + mu_fxs.^2) * beta - mu_mix.^2;       % CORRECT

    case {'tnormal'}
        % Cannot set D equal to mean, no closed form solution for mu
        % Rather we set D equal to mode 
        Mu = abs(D); % a = 0; b = inf; 
        % Compute mean of distributions
        npdf = @(x) 1/sqrt(2*pi) * exp(-1/2*x.^2);
        ncdf = @(x) 1/2 + 1/2*erf(x/sqrt(2)); 
        a = 0; b = inf; alfa = (a - Mu)./S; 
        Z = 1 - ncdf(alfa);
        mu_fxs = Mu + npdf(alfa) .* S ./ Z;
        sigma2_fxs = S.^2 .* (1 + alfa.*ncdf(alfa)./ Z - ( ncdf(alfa)./Z ).^2 );
        % mu_fxs is not equal to mu_fxs as abs(D) is set equal to mode!
        sigma2_mix = (sigma2_fxs + mu_fxs.^2) * beta - mu_mix.^2;       % CORRECT

end

% draw random numbers
r = rand(N,1); w = cumsum(beta');
% now determine which component of the mixture we need to draw from
id_mix = K - sum(r<w,2) + 1;
% Find elements for each component - and number of points each component
id = cell(1,K); n_id = nan(1,K);
for k = 1:K
    id{k} = find(id_mix==k); n_id(k) = numel(id{k});
end
% Now sample
for t = 1:n
    for k = 1:K
        switch options.PDF % CONDITIONAL DISTRIBUTION
            case {'normal'}     % NORMAL: mean "Mu" and std. "S"
                rnd(t,id{k}) = normrnd(Mu(t,k),S(t,k),[1 n_id(k)]);
            case {'lognormal'}  % LOGNORMAL: mean "Mu" and std. "S"
                rnd(t,id{k}) = lognrnd(Mu(t,k),S(t,k),[1 n_id(k)]);
            case {'gamma'}      % GAMMA: shape "A" and scale "B"
                rnd(t,id{k}) = gamrnd(A(t,k),B(t,k),[1 n_id(k)]);
            case {'weibull'}    % WEIBULL: shape "Kk" and scale "Lambda"
                rnd(t,id{k}) = wblrnd(Lambda(t,k),Kk(t,k),[1 n_id(k)]);
            case {'gen_normal'} % GENERALIZED NORMAL: mean "Mu", std. "S", shape "tau"
                rnd(t,id{k}) = gnormrnd(Mu(t,k),S(t,k),tau(k),[1 n_id(k)]);
            case {'gev'}        % GENERALIZED EXTREME VALUE: mean "Mu", std. "S", shape "xi"
                rnd(t,id{k}) = gevrnd(xi(t,k),S(t,k),Mu(t,k),[1 n_id(k)]);
            case {'gpareto'}    % GENERALIZED PARETO: mean "Mu", std. "S", shape "xi"
                rnd(t,id{k}) = gprnd(kk(t,k),S(t,k),Theta(t,k),[1 n_id(k)]);
            case {'tnormal'}    % TRUNCATED NORMAL: mean "Mu", std. "S", interval "a" and "b"
                rnd(t,id{k}) = tnormrnd(Mu(t,k),S(t,k),a,b,[1 n_id(k)]);
                
        end
    end
end

% %     case {'power_law'}  % LOGNORMAL POWER-LAW DISTRIBUTION: mean "MU" and variance "sigma2"
% %         % ADJUST mean and sigma according to mean and variance of distribution
% %         mu = abs(D); var = sigma2; d = size(x,2); % Number of members of vector x
% %         switch options.TAU
% %             case {'1'} % 1: common value of tau
% %                 tau = x(d);
% %             case {'2'} % 2: individual tau
% %                 tau = x(d-K+1:d);
% %         end
% %         var_fxs = tau ./ (tau - 2) .* exp(2*(var + mu));
% %         var_mix = (var_fxs + mu.^2) * beta - mu_mix.^2; % Mixture variance, CHECK

% %                % CHECK mean and variance, matching?
% %                 mom_rnd = [mean(rnd(t,id)) , var(rnd(t,id))]
% %                 mom_anal = [abs(D(t,k)) sigma2(t,k)]
% %             case {'power_law'}
% %                 rnd(t,id) = powrnd(D(t,k),S(t,k),alfa(k),[1 n_id]);
