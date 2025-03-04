function [lowup,mix_norm] = BMA_norm(x,D,y,options,method)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% This function calculates the α-norm of the BMA mixture density          %
%                                                                         %
% SYNOPSIS: [lowup ,mix_norm] = BMA_norm(x,D,options,method)              %
%  where                                                                  %
%   x         [input] 1xd vector with maximum likelihood BMA parameters   %
%   D         [input] nxK matrix with forecasts of ensemble members       %
%   y         [input] nx1 matrix with verifying data                      %
%   options   [input] structure with BMA algorithmic variables            %
%   method    [input] method used, 'standard' or 'robust' (default)       %
%   lowup     [outpt] lower & upper (a/2,1-a/2) pred. intervl, a = 1e-4   %
%   mix_norm  [outpt] α-norm of BMA predictive density                    %
%                                                                         %
% (c) Written by Jasper A. Vrugt, Feb 2022                                %
% University of California Irvine                                         %
%                                                                         %
% EXPLANATION: Standard and Robust                                        %
% 1. Find value y that corresponds quantles "a/2" and "(1-a/2)", a = 1e-4 %
% 2. Use CDF and root finding to find: low = F^-1(a) and up = F^-1(1-a)   %
% Standard method then continues with                                     %
%   3. Discretize [low,up] support in ny values and evaluate PDF          %
%   4. We yield the δ-norm (δ=1 and δ=2) using numerical integration      %
% Robust method continues after step 2 with                               %
%   3. Use built-in integral function "fx" on [low,up] and δ=1 and δ=2    %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

if nargin < 4
    method = 'robust';  % Get α-norm using built-in integral function
end
[n,K] = size(D);        % Unpack number of observations
beta = x(1:K);          %#ok Unpack weights
count = 0;              % Set counter to zero
a = 1e-4;               % Consider integral between cdf = a and cdf = 1-a

% Define lower and upper quantile
alpha = sort([ a , (1-a) ]); nA = numel(alpha);

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

% sigma is positive and initialize density at measured values
S = abs(S);

switch options.PDF      % Conditional distribution
    case {'normal'}
        Mu = D;         %#ok Mean of normal distribution ( matrix )
        CDF_str = strcat('Fx = @(x,i,alfa)',['beta(1) * ' ...
            'normcdf(x,Mu(i,1),S(i,1))']);
        for k = 2:K
            CDF_str = strcat(CDF_str,'+ beta(',num2str(k),[') * ' ...
                'normcdf(x,Mu(i,'],num2str(k),'),S(i,',num2str(k),'))');
        end

    case {'lognormal'}
        logn_impl = 2;
        switch logn_impl
            case 1 % original implementation: works only constant variance
                Mu = log(abs(D)) - S.^2/2; %#ok --> mean is now "D"
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
                S2 = log(A+1); S = sqrt(S2);     %#ok
                Mu = log(abs(D)) - S2/2;         %#ok
        end 
        CDF_str = strcat('Fx = @(x,i,alfa)',['beta(1) * ' ...
            'logncdf(x,Mu(i,1),S(i,1))']);
        for k = 2:K
            CDF_str = strcat(CDF_str,'+ beta(',num2str(k),[') * ' ...
                'logncdf(x,Mu(i,'],num2str(k),'),S(i,',num2str(k),'))');
        end

    case {'tnormal'}
        Mu = abs(D);     %#ok Mode truncated normal             ( matrix )
        a = 0; b = inf;  %#ok Lower and upper end point of support    
        CDF_str = strcat('Fx = @(x,i,alfa)',['beta(1) * ' ...
            'tnormcdf(x,Mu(i,1),S(i,1),a,b)']);
        for k = 2:K
            CDF_str = strcat(CDF_str,'+ beta(',num2str(k),[') * ' ...
                'tnormcdf(x,Mu(i,'],num2str(k),'),S(i,', ...
                num2str(k),')',',a,','b',')');
        end

    case {'gen_normal'}
        Mu = D; d = size(x,2); %#ok Number of members of vector x
        switch options.TAU
            case {'1'} % 1: constant value of delta for all models
                tau = x(d) * ones(1,K);   %#ok
            case {'2'} % 2: individual delta for each model
                tau = x(d-K+1:d);         %#ok
        end
        CDF_str = strcat('Fx = @(x,i,alfa)',['beta(1) * ' ...
            'gnormcdf(x,Mu(i,1),S(i,1),tau(1))']);
        for k = 2:K
            CDF_str = strcat(CDF_str,'+ beta(',num2str(k),[') * ' ...
                'gnormcdf(x,Mu(i,'],num2str(k),'),S(i,',num2str(k), ...
                '),tau(',num2str(k),'))');
        end
        
    case {'gamma'}
        Mu = abs(D);        % Mean of gamma distribution        ( matrix )
        A = (Mu.^2)./S.^2;  %#ok Shape parameter                ( matrix )
        B = S.^2./Mu;       %#ok Scale parameter                ( matrix )
        CDF_str = strcat('Fx = @(x,i,alfa)',['beta(1) * ' ...
            'gamcdf(x,A(i,1),B(i,1))']);
        for k = 2:K
            CDF_str = strcat(CDF_str,'+ beta(',num2str(k),[') * ' ...
                'gamcdf(x,A(i,'],num2str(k),'),B(i,',num2str(k),'))');
        end

    case {'weibull'}
        wbl_impl = 2;
        switch wbl_impl
            case 1
                Kk = S; Lambda = abs(D)./gamma(1+1./Kk);
                % Lambda cannot be zero
                Lambda = max(Lambda,realmin);    %#ok
            case 2
                Lambda = S; X = abs(D)./Lambda;
                c = sqrt(2*pi)/exp(1) - gamma(1.461632);
                Lx = @(x) log((x+c)/sqrt(2*pi)); A = Lx(X);
                B = A ./ real(lambertw(A/exp(1))) + 1/2;
                Kk = 1./(B - 1);                %#ok
        end
        % -> does not work for predict int.
        CDF_str = strcat('Fx = @(x,i,alfa)',['beta(1) * wblcdf(x,' ...
            'Lambda(i,1),Kk(i,1))']);
        for k = 2:K
            CDF_str = strcat(CDF_str,'+ beta(',num2str(k),[') * ' ...
                'wblcdf(x,Lambda(i,'],num2str(k),'),Kk(i,', ...
                num2str(k),'))');
        end

    case {'gev'}
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
        Mu(1:n,kk) = D(1:n,kk) - S(1:n,kk) * Eul_constant;          %#ok
        % Create string for BMA mixture of GEV distributions
        CDF_str = strcat('Fx = @(x,i,alfa)',['beta(1) * ' ...
            'gevcdf(x,xi(i,1),S(i,1),Mu(i,1))']);
        for k = 2:K
            CDF_str = strcat(CDF_str,'+ beta(',num2str(k),[') * ' ...
                'gevcdf(x,xi(i,',num2str(k),'),S(i,',num2str(k), ...
                '),Mu(i,'],num2str(k),'))');
        end

    case {'gpareto'}
        d = size(x,2); % Number of members of vector x
        switch options.TAU
            case {'1'} % 1: common value of tau
                kk = x(d) * ones(n,K);
            case {'2'} % 2: individual tau
                kk = repmat(x(d-K+1:d),n,1);
        end
        Theta = abs(D) - S./(1-kk);                                 %#ok
        % Create string for BMA mixture of GEV distributions
        CDF_str = strcat('Fx = @(x,i,alfa)',['beta(1) * ' ...
            'gpcdf(x,kk(i,1),S(i,1),Theta(i,1))']);
        for k = 2:K
            CDF_str = strcat(CDF_str,'+ beta(',num2str(k),[') * ' ...
                'gpcdf(x,kk(i,',num2str(k),'),S(i,',num2str(k), ...
                '),Theta(i,'],num2str(k),'))');
        end        

end
% Copy mixture cdf to create mixture pdf called fx
PDF_str = CDF_str;
% Now add "-alfa" to the end (root finding) and activate BMA mixture cdf
Fx = strcat(CDF_str,'- alfa;'); eval(char(Fx)); 
% Now activate the mixture pdf
PDF_str = strcat(PDF_str,';'); PDF_str = strrep(PDF_str,'cdf','pdf');
% Replace alpha (significance level) with α-norm (called delta, δ-norm)
PDF_str = strrep(PDF_str,',alfa',',delta'); PDF_str = strrep(PDF_str, ...
    'Fx','fx');
PDF_str = strrep(PDF_str,'@(x,i,delta)','@(x,i,delta) (');
switch options.PDF      
    case {'tnormal'}
        fx = strrep(PDF_str,');',') ).^delta;');
    otherwise
        fx = strrep(PDF_str,'));',')) ).^delta;');
end
% Activate the BMA mixture pdf
eval(fx);

% How many equidistant samples of BMA mixture?
ny = 1e3; lowup = nan(n,nA); mix_norm = nan(n,2);
% Determine lower and upper end of BMA mixture
for i = 1 : n
    % Print progress
    if mod(i,floor(n/25)) == 0
        if i > 1
            fprintf(1, repmat('\b',1,count)); % delete line before
            count = fprintf('BMA α-norm calculation, %% done: %3.2f', ...
                100*(i/n));
        end
    end
    % Get exact answer using root finding of CDF, but slower
    for z = 1 : nA
%        lowup(i,z) = fzero(@(y) Fx(y,i,alpha(z)),abs(mean(D(i,1:K))));
        lowup(i,z) = fzero(@(y) Fx(y,i,alpha(z)),y(i));
    end
    % Two implementations
    switch method
        case 'standard'
            % ny equally spaced values between low and up
            yi = linspace(lowup(i,1),lowup(i,2),ny)';
            % Evaluate BMA mixture density at each yi
            pdf_yi = fx(yi,i);
            % Compute δ-norm of BMA mixture: δ=1
            mix_norm(i,1) = trapz(yi,pdf_yi); % = sum(pdf_yi)*(yi(2)-yi(1))
            % Compute δ-norm of BMA mixture: δ=2 for QS and SS
            mix_norm(i,2) = sqrt(trapz(yi,pdf_yi.^2));
        case 'robust'
            % Compute δ-norm of BMA mixture
            for delta = 1:2
                mix_norm(i,delta) = integral(@(y) fx(y,i,delta), ...
                    lowup(i,1),lowup(i,2)).^(1/delta);
            end
    end
end
fprintf(1,'\n');

end
