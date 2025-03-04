function [X,pdf_BMA,cdf_BMA,plimit] = BMA_mixture(x,D,T,options)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% This function returns CDF & PDF of the BMA mixture distribution at T    %
%                                                                         %
% SYNOPSIS: [X,pdf,cdf,pred] = BMA_mixture(x,D,T,options)                 %
%  where                                                                  %
%   x         [input] 1xd vector with maximum likelihood BMA parameters   %
%   D         [input] nxK matrix with forecasts of ensemble members       %
%   T         [input] mx1 vector with time indices of D [= row numbers]   %
%   options   [input] structure with BMA algorithmic variables            %
%   X         [outpt] 1xm vector with entity of interest for all T        %
%   pdf_BMA   [outpt] 1xm vector with pdf of BMA distribution at X        %
%   cdf_BMA   [outpt] 1xm vector with cdf of BMA distribution at X        %
%   plimit    [outpt] lower and upper end 1-alpha prediction interval     %
%                                                                         %
% (c) Written by Jasper A. Vrugt, Feb 2006                                %
% University of California Irvine                                         %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

[n,K] = size(D);                    % Unpack number of observations
beta = x(1:K);                      %#ok Unpack weights
count = 0;                          % Set counter to zero
m = numel(T);                       % # time entries of T
N = 1e5;                            % # discretized x-values
X = linspace(0,2*max(max(D)),N);    % Create x-values
[pdf_BMA,cdf_BMA] = deal(nan(m,N)); % Initialize PDF and CDF observations

% % % Now calculate the bias-corrected forecasts
% % for k = 1 : K
% %     D(:,k) = a(k) + b(k) * D(:,k);
% % end

% Unpack alpha values - compute confidence levels and sort quantiles
g = 1-options.alpha; gam = sort([(1-g)/2 , (1-g)/2 + g]); nG = numel(gam);

switch options.VAR      % Variance option
    case {'1'}          % 1: common constant variance
        S = x(K+1) * ones(n,K);
    case {'2'}          % 2: individual constant variance
        S = repmat(x(K+1:2*K),n,1);
    case {'3'}          % 3: common non-constant variance
        c = x(K+1); S = c * abs(D); 
    case {'4'}          % 4: individual non-constant variance
        c = x(K+1:2*K); S = bsxfun(@times,c,abs(D));
    otherwise
        error('MODELAVG:BMA_quantile','Unknown variance option');
end

% sigma is positive and initialize density at measured values
S = abs(S);

switch options.PDF      % Conditional distribution
    
    case {'normal'}     % normpdf(Mu,S)
        Mu = D;         %#ok Mean of normal distribution ( matrix )
        CDF_str = strcat('Fx = @(x,i,gam)',['beta(1) * ' ...
            'normcdf(x,Mu(i,1),S(i,1))']);
        for k = 2:K
            CDF_str = strcat(CDF_str,'+ beta(',num2str(k),[') * ' ...
                'normcdf(x,Mu(i,'],num2str(k),'),S(i,',num2str(k),'))');
        end

    case {'lognormal'}  % lognpdf(Mu,S) (impl = 1); lognpdf(Mu,sigma_fxs)
        logn_impl = 2;
        switch logn_impl
            case 1 % originl implmntation: works only for constant variance
                Mu = log(abs(D)) - S.^2/2; %#ok --> mean is equal to "D"
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
                S2 = log(A+1); S = sqrt(S2);   %#ok
                Mu = log(abs(D)) - S2/2;       %#ok
        end 
        CDF_str = strcat('Fx = @(x,i,gam)',['beta(1) * ' ...
            'logncdf(x,Mu(i,1),S(i,1))']);
        for k = 2:K
            CDF_str = strcat(CDF_str,'+ beta(',num2str(k),[') * ' ...
                'logncdf(x,Mu(i,'],num2str(k),'),S(i,',num2str(k),'))');
        end

    case {'tnormal'}  %% tnormpdf(Mu,S,a,b)
        Mu = abs(D);            %#ok Mode truncated normal     (matrix)
        a = 0; b = inf;         %#ok Lower and upper end point of support    
        CDF_str = strcat('Fx = @(x,i,gam)',['beta(1) * ' ...
            'tnormcdf(x,Mu(i,1),S(i,1),a,b)']);
        for k = 2:K
            CDF_str = strcat(CDF_str,'+ beta(',num2str(k),[') * ' ...
                'tnormcdf(x,Mu(i,'],num2str(k),'),S(i,', ...
                num2str(k),'),a,b)');
        end
        
    case {'gen_normal'}  % gnormpdf(Mu,S,tau)
        Mu = D; d = size(x,2); %#ok Number of members of vector x
        switch options.TAU
            case {'1'} % 1: constant value of delta for all models
                tau = x(d) * ones(1,K);    %#ok
            case {'2'} % 2: individual delta for each model
                tau = x(d-K+1:d);          %#ok
        end
        CDF_str = strcat('Fx = @(x,i,gam)',['beta(1) * ' ...
            'gnormcdf(x,Mu(i,1),S(i,1),tau(1))']);
        for k = 2:K
            CDF_str = strcat(CDF_str,'+ beta(',num2str(k),[') * ' ...
                'gnormcdf(x,Mu(i,'],num2str(k),'),S(i,',num2str(k), ...
                '),tau(',num2str(k),'))');
        end
        
    case {'gamma'}  % gampdf(A,B)
        Mu = abs(D);        % Mean of gamma distribution    (matrix)
        A = (Mu.^2)./S.^2;  %#ok Shape parameter            (matrix)
        B = S.^2./Mu;       %#ok Scale parameter            (matrix)
        CDF_str = strcat('Fx = @(x,i,gam)',['beta(1) * ' ...
            'gamcdf(x,A(i,1),B(i,1))']);
        for k = 2:K
            CDF_str = strcat(CDF_str,'+ beta(',num2str(k),[') * ' ...
                'gamcdf(x,A(i,'],num2str(k),'),B(i,',num2str(k),'))');
        end
    
    case {'weibull'}  % wblpdf(Lambda,Kk)
        wbl_impl = 2;
        switch wbl_impl
            case 1  % Shape parameter treated sigma & Lambda (scale parmtr)
                    % derived so that mean of Weibull PDFs equals D
                Kk = S; Lambda = abs(D)./gamma(1+1./Kk);
                % Lambda cannot be zero
                Lambda = max(Lambda,realmin);     %#ok
            case 2  % Scale parameter (Lambda) treatd as sigma (= correct!)
                    % and shape parameter (Kk) set so that PDF mean = D
                Lambda = S; X = abs(D)./Lambda;
                c = sqrt(2*pi)/exp(1) - gamma(1.461632);
                Lx = @(x) log((x+c)/sqrt(2*pi)); A = Lx(X);
                B = A ./ real(lambertw(A/exp(1))) + 1/2;
                Kk = 1./(B - 1);                  %#ok
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
        Mu(1:n,kk) = D(1:n,kk) - S(1:n,kk) * Eul_constant;  %#ok
        % Create string for BMA mixture of GEV distributions
        CDF_str = strcat('Fx = @(x,i,gam)',['beta(1) * ' ...
            'gevcdf(x,xi(i,1),S(i,1),Mu(i,1))']);
        for k = 2:K
            CDF_str = strcat(CDF_str,'+ beta(',num2str(k),[') * ' ...
                'gevcdf(x,xi(i,',num2str(k),'),S(i,',num2str(k), ...
                '),Mu(i,'],num2str(k),'))');
        end

    case {'gpareto'}  %% gppdf(kk,S,Theta)
        d = size(x,2); % Number of members of vector x
        switch options.TAU
            case {'1'} % 1: common value of tau
                kk = x(d) * ones(n,K);
            case {'2'} % 2: individual tau
                kk = repmat(x(d-K+1:d),n,1);
        end
        Theta = abs(D) - S./(1-kk);                            %#ok 
        % Create string for BMA mixture of GEV distributions
        CDF_str = strcat('Fx = @(x,i,gam)',['beta(1) * ' ...
            'gpcdf(x,kk(i,1),S(i,1),Theta(i,1))']);
        for k = 2:K
            CDF_str = strcat(CDF_str,'+ beta(',num2str(k),[') * ' ...
                'gpcdf(x,kk(i,',num2str(k),'),S(i,',num2str(k), ...
                '),Theta(i,'],num2str(k),'))');
        end
        
end
% Copy mixture cdf to create mixture pdf called fx
PDF_str = CDF_str;
% Now activate the mixture cdf with alfa for root finding
Fx = strcat(CDF_str,'- gam;'); eval(char(Fx));
% Now activate the mixture pdf
PDF_str = strcat(PDF_str,';'); PDF_str = strrep(PDF_str,'cdf','pdf');
PDF_str = strrep(PDF_str,',gam',''); fx = strrep(PDF_str,'Fx','fx'); 
eval(fx);
% Compute PDF and CDF at discretized x-values
for i = 1:m
    pdf_BMA(i,1:N) = fx(X,T(i)); cdf_BMA(i,1:N) = Fx(X,T(i),0);
end

% Now loop over each observation and determine BMA quantiles
plimit = nan(m,nG);
for i = 1 : m
    % Print progress
    if mod(i,floor(m/25)) == 0
        if i > 1
            fprintf(1, repmat('\b',1,count)); % delete line before
            count = fprintf(['BMA mixture calculation, %% done: ' ...
                '%3.2f'],100*(i/m));
        end
    end
    % Get exact answer using root finding of CDF, but is slower
    for z = 1 : nG
        plimit(i,z) = fzero(@(x) Fx(x,T(i),gam(z)),abs(mean(D(T(i),1:K))));
    end
end
fprintf(1,'\n');

end
