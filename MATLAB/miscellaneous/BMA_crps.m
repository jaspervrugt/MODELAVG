function crps = BMA_crps(x,D,y,options)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% This function computes the Continuous Ranked Probability Score for      %
% BMA mixture distribution using numerical (trapezoidal) integration of   %
% the quantile formulation:                                               %
%     CRPS(P,w) = w(1 - 2F_P(w)) + 2\int_{0}^{1}tau F_P^{-1}(tau)d tau    %
%                 - 2\int_{F_P(w)}^{1} F_P^{-1}(tau)d tau                 %
%                                                                         %
% SYNOPSIS: [mean_crps,crps,num_nan] = BMA_crps(x,D,Y,options)            %
%  where                                                                  %
%   x         [input] 1xd vector with maximum likelihood BMA parameters   %
%   D         [input] nxK matrix with forecasts of ensemble members       %
%   y         [input] nx1 vector with verifying observations              %
%   options   [input] structure with BMA algorithmic variables            %
%   crps      [outpt] nx1 vector with CRPS values of BMA mixture CDF      %
%                                                                         %
% Reference:                                                              %
%                                                                         %
% Note: the CRPS is positively-oriented!! (larger is better)              %
%                                                                         %
% (c) Written by Jasper A. Vrugt, Feb 2022                                %
% University of California Irvine                                         %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

[n,K] = size(D);                    % Unpack number of observations
beta = x(1:K);                      %#ok Unpack weights
count = 0;                          % Set counter to zero
crps = nan(n,1);                    % Initialize CRPS of BMA mixture CDF

% Note: we use "p" for tau
P1 = [ 0.001 0.002 0.005 0.01 0.02 0.05:0.05:0.45 ];
P = [ P1 , 0.5 , sort(1 - P1) ]'; nP = numel(P);
%P = [0.001 0.01 :0.01:0.99 0.999]'; nP = numel(P);

% Determine the std. of each members' conditional PDF
switch options.VAR      % Variance option
    case {'1'}          % 1: common constant variance
        S = x(K+1) * ones(n,K);
    case {'2'}          % 2: individual constant variance
        S = repmat(x(K+1:2*K),n,1);
    case {'3'} % 3: common non-constant variance
        c = x(K+1); S = c * abs(D); 
    case {'4'} % 4: individual non-constant variance
        c = x(K+1:2*K); S = bsxfun(@times,c,abs(D));
    otherwise
        error('MODELAVG:BMA_quantile','Unknown variance option');
end
S = max(S,eps);     % each element (n x K)-matrix S >= 2.22e-16
S2 = S.^2;          % variance (n x K)-matrix

%% Step 1: Define CDF of mixture distribution
switch options.PDF  % Check conditional distribution

    case 'normal'   % normpdf(Mu,S)
        Mu = D;
        % Create string for BMA mixture of normal distributions
        CDF_str = strcat('F_P = @(x,i,p)',['beta(1) * ' ...
            'normcdf(x,Mu(i,1),S(i,1))']);
        for k = 2:K
            CDF_str = strcat(CDF_str,'+ beta(',num2str(k),[') * ' ...
                'normcdf(x,Mu(i,'],num2str(k),'),S(i,',num2str(k),'))');
        end

    case 'lognormal' % lognpdf(Mu,S)
        logn_impl = 2;
        switch logn_impl
            case 1 % original implementation: works only constant variance
                Mu = log(abs(D)) - S.^2/2; % --> mean is now equal to "D"
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
        end 
        % Create string for BMA mixture of lognormal distributions
        CDF_str = strcat('F_P = @(x,i,p)',['beta(1) * ' ...
            'logncdf(x,Mu(i,1),S(i,1))']);
        for k = 2:K
            CDF_str = strcat(CDF_str,'+ beta(',num2str(k),[') * ' ...
                'logncdf(x,Mu(i,'],num2str(k),'),S(i,',num2str(k),'))');
        end

    case {'tnormal'}    % tnormpdf(
        Mu = abs(D);            % Mode truncated normal         ( matrix )
        a = 0; b = inf;         % Lower and upper end point of support    
        CDF_str = strcat('F_P = @(x,i,p)',['beta(1) * ' ...
            'tnormcdf(x,Mu(i,1),S(i,1),a,b)']);
        for k = 2:K
            CDF_str = strcat(CDF_str,'+ beta(',num2str(k),[') * ' ...
                'tnormcdf(x,Mu(i,'],num2str(k),'),S(i,', ...
                num2str(k),'),a,b)');
        end
        
    case 'gen_normal' % gnormpdf(Mu,S,tau)
        d = size(x,2); 
        switch options.TAU
            case {'1'} % 1: common value of tau
                tau = x(d) * ones(1,K);
            case {'2'} % 2: individual tau
                tau = x(d-K+1:d);
        end
        Mu = D; 
        % Create string for BMA mixture of generalized normal distributions
        CDF_str = strcat('F_P = @(x,i,p)',['beta(1) * ' ...
            'gnormcdf(x,Mu(i,1),S(i,1),tau(1))']);
        for k = 2:K
            CDF_str = strcat(CDF_str,'+ beta(',num2str(k),[') * ' ...
                'gnormcdf(x,Mu(i,'],num2str(k),'),S(i,',num2str(k), ...
                '),tau(',num2str(k),'))');
        end

    case 'gamma' % gampdf(A,B)
        Mu = abs(D); A = Mu.^2./S2; B = S2./Mu; 
        % Create string for BMA mixture of gamma distributions
        CDF_str = strcat('F_P = @(x,i,p)',['beta(1) * ' ...
            'gamcdf(x,A(i,1),B(i,1))']);
        for k = 2:K
            CDF_str = strcat(CDF_str,'+ beta(',num2str(k),[') * ' ...
                'gamcdf(x,A(i,'],num2str(k),'),B(i,',num2str(k),'))');
        end

    case 'weibull' % wblpdf(Lambda,Kk)
        wbl_impl = 2;
        switch wbl_impl
            case 1  % Shape parameter treated sigma & Lambda (scale prmtr)
                    % derived so that mean of Weibull PDFs equals D
                Kk = S; Lambda = abs(D)./gamma(1+1./Kk);
                % Lambda cannot be zero
                Lambda = max(Lambda,realmin);
            case 2  % Scale parameter (Lambda) treated as sigma (= corrct!)
                    % and the shape parameter (Kk) is set so that PDF µ = D
                Lambda = S; X = abs(D)./Lambda;
                c = sqrt(2*pi)/exp(1) - gamma(1.461632);
                Lx = @(x) log((x+c)/sqrt(2*pi)); A = Lx(X);
                B = A ./ real(lambertw(A/exp(1))) + 1/2;
                Kk = 1./(B - 1);
                % --> this implementation is not stable yet; requires the
                % inverse of the gamma function
        end
        % Create string for BMA mixture of Weibull distributions
        CDF_str = strcat('F_P = @(x,i,p)',['beta(1) * wblcdf(x,' ...
            'Lambda(i,1),Kk(i,1))']);
        for k = 2:K
            CDF_str = strcat(CDF_str,'+ beta(',num2str(k),[') * ' ...
                'wblcdf(x,Lambda(i,'],num2str(k), ...
                '),Kk(i,',num2str(k),'))');
        end

    case 'gev'  % gevpdf(xi,S,Mu)
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
        CDF_str = strcat('F_P = @(x,i,p)',['beta(1) * ' ...
            'gevcdf(x,xi(i,1),S(i,1),Mu(i,1))']);
        for k = 2:K
            CDF_str = strcat(CDF_str,'+ beta(',num2str(k),[') * ' ...
                'gevcdf(x,xi(i,',num2str(k),'),S(i,',num2str(k), ...
                '),Mu(i,'],num2str(k),'))');
        end

    case {'gpareto'}  % gppdf(kk,S,Theta) --> kk in [-2.5,1/2]
        d = size(x,2); % Number of members of vector x
        switch options.TAU
            case {'1'} % 1: common value of tau
                kk = x(d) * ones(n,K);
            case {'2'} % 2: individual tau
                kk = repmat(x(d-K+1:d),n,1);
        end
        Theta = abs(D) - S./(1-kk);    %#ok
        % Create string for BMA mixture of GEV distributions
        CDF_str = strcat('F_P = @(x,i,p)',['beta(1) * ' ...
            'gpcdf(x,kk(i,1),S(i,1),Theta(i,1))']);
        for k = 2:K
            CDF_str = strcat(CDF_str,'+ beta(',num2str(k),[') * ' ...
                'gpcdf(x,kk(i,',num2str(k),'),S(i,',num2str(k), ...
                '),Theta(i,'],num2str(k),'))');
        end

    otherwise
        error('MODELAVG:BMA_crps','Unknown conditional PDF');

end
% Now activate the mixture cdf with p (tau) for root finding
F_P = strcat(CDF_str,'- p;'); eval(char(F_P))

%% Step 2: Solve for values of y so that F_P^{-1}(y) = P
%% Step 3: Compute CRPS 
FinvP = nan(nP,1);
for i = 1 : n
    % Print progress
    if mod(i,floor(n/25)) == 0
        if i > 1
            fprintf(1, repmat('\b',1,count)); % delete line before
            count = fprintf(['BMA CRPS calculation, %% done: ' ...
                '%3.2f'],100*(i/n));
        end
    end
    % get exact answer using root finding of CDF, but slower
    for z = 1 : nP
        switch z
            case 1
                FinvP(z,1) = fzero(@(y) F_P(y,i,P(z)),abs(mean(D(i,1:K))));
            otherwise
                FinvP(z,1) = fzero(@(y) F_P(y,i,P(z)),FinvP(z-1,1));
        end
    end
    % Evaluate mixture CDF at omega
    Fy = F_P(y(i,1),i,0); ii = P > Fy; 
    % Check whether we have at least 1 element of ii
    if sum(ii) == 0
        % Area right of observation is almost zero (omega close to P = 1)
% %     P_temp = 1/2*(1 + Fw); FinvP_temp = fzero(@(x) F_P(x,i,P_temp),Fw);
% %     intgral2 = - 2*trapz([Fw ; P_temp],[Y(i,1) ; FinvP_temp])
% %     % This integral is very close to zero --> discard
        % Compute CRPS using trapezoidal integration rule (left from omega)
        crps(i,1) = y(i,1)*(1 - 2*Fy) + 2*trapz(P,P.*FinvP);
    else
        % Compute CRPS using trapezoidal integration rule
        crps(i,1) = y(i,1)*(1 - 2*Fy) + 2*trapz(P,P.*FinvP) ...
            - 2*trapz([Fy ; P(ii)],[y(i,1) ; FinvP(ii)]);
    end
end
fprintf(1,'\n');

end
