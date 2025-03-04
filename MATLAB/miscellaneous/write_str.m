function varargout = write_str(beta,D,n,K,S,mu_mix,options)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% This function returns mean, scale and location parameters BMA mixture   %
% components - and a string for the cdf of the BMA mixture                %
%                                                                         %
% SYNOPSIS: varargout = write_str(beta,D,n,K,S,mu_mix,options)            % 
%  where                                                                  %
%   beta      [input] 1xd vector with BMA weights                         %
%   D         [input] nxK matrix with forecasts of ensemble members       %
%   n         [input] number of verifying measurements                    %
%   K         [input] number of ensemble members                          %
%   S         [input] nxK matrix with standard deviation forecasts        %
%   mu_mix    [input] nx1 vector with µ BMA mean forecast                 %
%   options   [input] structure with BMA algorithmic variables            %
%   varargout [outpt] variance BMA mixture, cdf string, and shape         %
%                     mixture components BMA distribution                 %
%                                                                         %
% (c) Written by Jasper A. Vrugt, Feb 2012                                %
% University of California Irvine                                         %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

S2 = S.^2;
% Conditional distribution
switch options.PDF      
    case 'normal'       % normpdf(Mu,S)
        Mu = D;         % Mean of normal distribution ( matrix )
        CDF_str = strcat('Fx = @(x,i,gam)',['beta(1) * ' ...
            'normcdf(x,Mu(i,1),S(i,1))']);
        for k = 2:K
            CDF_str = strcat(CDF_str,'+ beta(',num2str(k),[') * ' ...
                'normcdf(x,Mu(i,'],num2str(k),'),S(i,',num2str(k),'))');
        end
        sigma2_mix = (S2 + Mu.^2) * beta - mu_mix.^2;
        varargout(1:4) = {sigma2_mix,CDF_str,Mu,S};

    case 'lognormal'    % lognpdf(Mu,S) (impl = 1); lognpdf(Mu,sigma_fxs)
        mu_fxs = D;     % Mean each fk
        logn_impl = 2;
        switch logn_impl
            case 1 % originl implmntation: works only for constant variance
                Mu = log(abs(D)) - S.^2/2; % --> mean is now equal to "D"
                sigma2_fxs = D.^2 .* ( exp(S2) - 1 );   % Variance each fk
                sigma2_mix = (sigma2_fxs + mu_fxs.^2) * beta - mu_mix.^2; % CORRECT
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
                Mu = log(abs(D)) - S2/2;
                sigma2_mix = (sigma2_fxs + mu_fxs.^2) * beta - mu_mix.^2;   % CORRECT
        end
        CDF_str = strcat('Fx = @(x,i,gam)',['beta(1) * ' ...
            'logncdf(x,Mu(i,1),S(i,1))']);
        for k = 2:K
            CDF_str = strcat(CDF_str,'+ beta(',num2str(k),[') * ' ...
                'logncdf(x,Mu(i,'],num2str(k),'),S(i,',num2str(k),'))');
        end
        varargout(1:4) = {sigma2_mix,CDF_str,Mu,S};

    case {'tnormal'}  %% tnormpdf(Mu,S,a,b)
        Mu = abs(D);            % Mode truncated normal     (matrix)
        CDF_str = strcat('Fx = @(x,i,gam)',['beta(1) * ' ...
            'tnormcdf(x,Mu(i,1),S(i,1),a,b)']);
        for k = 2:K
            CDF_str = strcat(CDF_str,'+ beta(',num2str(k),[') * ' ...
                'tnormcdf(x,Mu(i,'],num2str(k),'),S(i,', ...
                num2str(k),'),a,b)');
        end
        % Compute mean of distributions
        npdf = @(x) 1/sqrt(2*pi) * exp(-1/2*x.^2);
        ncdf = @(x) 1/2 + 1/2*erf(x/sqrt(2));
        a = 0; b = inf; alfa = (a - Mu)./S;
        Z = 1 - ncdf(alfa);
        mu_fxs = Mu + npdf(alfa) .* S ./ Z;
        sigma2_fxs = S.^2 .* (1 + alfa.*ncdf(alfa)./ ...
            Z - ( ncdf(alfa)./Z ).^2 );
        % mu_fxs is not equal to mu_fxs as abs(D) is set equal to mode!
        sigma2_mix = (sigma2_fxs + mu_fxs.^2) * beta - mu_mix.^2;           % CORRECT
        varargout(1:6) = {sigma2_mix,CDF_str,Mu,S,a,b};

    case {'gen_normal'}  % gnormpdf(Mu,S,tau)
        Mu = D; d = size(x,2); % Number of members of vector x
        switch options.TAU
            case {'1'} % 1: constant value of delta for all models
                tau = x(d) * ones(1,K);
            case {'2'} % 2: individual delta for each model
                tau = x(d-K+1:d);
        end
        CDF_str = strcat('Fx = @(x,i,gam)',['beta(1) * ' ...
            'gnormcdf(x,Mu(i,1),S(i,1),tau(1))']);
        for k = 2:K
            CDF_str = strcat(CDF_str,'+ beta(',num2str(k),[') * ' ...
                'gnormcdf(x,Mu(i,'],num2str(k),'),S(i,',num2str(k), ...
                '),tau(',num2str(k),'))');
        end
        sigma2 = S2 .* gamma(3./repmat(tau,n,1)) ./ ...
            gamma(1./repmat(tau,n,1));      % Variance of fk's
        sigma2_mix = (sigma2 + Mu.^2) * beta - mu_mix.^2;                   % CORRECT
        varargout(1:5) = {sigma2_mix,CDF_str,Mu,S,tau};

    case {'gamma'}      % gampdf(A,B)
        Mu = abs(D);        % Mean of gamma distribution    (matrix)
        A = (Mu.^2)./S.^2;  % Shape parameter               (matrix)
        B = S.^2./Mu;       % Scale parameter               (matrix)
        CDF_str = strcat('Fx = @(x,i,gam)',['beta(1) * ' ...
            'gamcdf(x,A(i,1),B(i,1))']);
        for k = 2:K
            CDF_str = strcat(CDF_str,'+ beta(',num2str(k),[') * ' ...
                'gamcdf(x,A(i,'],num2str(k),'),B(i,',num2str(k),'))');
        end
        sigma2_mix = (S2 + Mu.^2) * beta - mu_mix.^2;                       % CORRECT
        varargout(1:4) = {sigma2_mix,CDF_str,A,B};

    case {'weibull'}  % wblpdf(Lambda,Kk)
        Mu = abs(D); wbl_impl = 2;
        switch wbl_impl
            case 1  % Shape parameter treated sigma & Lambda (scale parmtr)
                % derived so that mean of Weibull PDFs equals D
                Kk = S; Lambda = Mu./gamma(1+1./Kk);
                % Lambda cannot be zero
                Lambda = max(Lambda,realmin);
            case 2  % Scale parameter (Lambda) treatd as sigma (= correct!)
                % and shape parameter (Kk) set so that PDF mean = D
                Lambda = S; X = Mu./Lambda;
                c = sqrt(2*pi)/exp(1) - gamma(1.461632);
                Lx = @(x) log((x+c)/sqrt(2*pi)); A = Lx(X);
                B = A ./ real(lambertw(A/exp(1))) + 1/2;
                Kk = 1./(B - 1);
                % --> this implementation is not stable yet; requires the
                % inverse of the gamma function
        end
        % -> does not work perfectly
        CDF_str = strcat('Fx = @(x,i,gam)',['beta(1) * wblcdf(x,' ...
            'Lambda(i,1),Kk(i,1))']);
        for k = 2:K
            CDF_str = strcat(CDF_str,'+ beta(',num2str(k),[') * ' ...
                'wblcdf(x,Lambda(i,'],num2str(k),'),Kk(i,', ...
                num2str(k),'))');
        end
        % Variance of fk's
        sigma2 = Lambda.^2 .* (gamma(1 + 2./Kk) - gamma(1+1./Kk).^2);   
        sigma2_mix = (sigma2 + Mu.^2) * beta - mu_mix.^2;                   % CORRECT?
        varargout(1:4) = {sigma2_mix,CDF_str,Lambda,Kk};

    case {'gev'}  %% gevpdf(xi,S,Mu)
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
        % Create string for BMA mixture of GEV distributions
        CDF_str = strcat('Fx = @(x,i,gam)',['beta(1) * ' ...
            'gevcdf(x,xi(i,1),S(i,1),Mu(i,1))']);
        for k = 2:K
            CDF_str = strcat(CDF_str,'+ beta(',num2str(k),[') * ' ...
                'gevcdf(x,xi(i,',num2str(k),'),S(i,',num2str(k), ...
                '),Mu(i,'],num2str(k),'))');
        end
        % built-in function to get mean and variance of each GEV component
        [mu_fxs,sigma2_fxs] = gevstat(xi,S,Mu);
        % mu_fxs is equal to abs(D), as requested by setting mean
        sigma2_mix = (sigma2_fxs + mu_fxs.^2) * beta - mu_mix.^2;           % CORRECT
        varargout(1:5) = {sigma2_mix,CDF_str,xi,S,Mu};

    case {'gpareto'}  %% gppdf(kk,S,Theta)
        d = size(x,2); % Number of members of vector x
        switch options.TAU
            case {'1'} % 1: common value of tau
                kk = x(d) * ones(n,K);
            case {'2'} % 2: individual tau
                kk = repmat(x(d-K+1:d),n,1);
        end
        Theta = abs(D) - S./(1-kk);
        % Create string for BMA mixture of GEV distributions
        CDF_str = strcat('Fx = @(x,i,gam)',['beta(1) * ' ...
            'gpcdf(x,kk(i,1),S(i,1),Theta(i,1))']);
        for k = 2:K
            CDF_str = strcat(CDF_str,'+ beta(',num2str(k),[') * ' ...
                'gpcdf(x,kk(i,',num2str(k),'),S(i,',num2str(k), ...
                '),Theta(i,'],num2str(k),'))');
        end
        % built-in function to get mean and variance of each GP component
        [mu_fxs,sigma2_fxs] = gpstat(kk,S,Theta);
        % mu_fxs is equal to abs(D), as requested by setting mean
        sigma2_mix = (sigma2_fxs + mu_fxs.^2) * beta - mu_mix.^2;           % CORRECT
        varargout(1:5) = {sigma2_mix,CDF_str,kk,S,Theta};

end

end

