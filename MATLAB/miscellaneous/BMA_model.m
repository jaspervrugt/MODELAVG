function varargout = BMA_model(x,D,y,options,return_arg,varargin)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%BMA_MODEL This multifunctional subroutine embeds the following functions %
% 1. [ell,L] = BMA_lik(x,D,Y,options)                                     %
% 2. [plimits,pdf_y,cdf_y] = BMA_quantile(x,D,y,options)                  %
% 3. [lowup,mix_norm] = BMA_norm(x,D,options)                             %
% 4. [X,pdf_BMA,cdf_BMA,plimits] = BMA_mixture(x,D,T,options)             %
% 5. [rnd,mu_mix,sigma2_mix] = BMA_draw(x,D,N,options)                    %
% 6. [crps] = BMA_crps(x,D,y,options)                                     %
% and computes the following items of the BMA mixture distribution        %
% 1. likelihood                                                           %
% 2. prediction limits, pdf and cdf @ y                                   %
% 3. α-norm                                                               %
% 4. BMA mixture cdf and pdf                                              %
% 5. random samples, mean and variance, and                               %
% 6. continuous ranked probability score for y                            %
%                                                                         %
% SYNOPSIS: varargout = BMA_model(x,D,y,options,func_name,varargin)       %
%  where                                                                  %
%   x         [input] 1xd vector of maximum likelihood BMA parameters     %
%   D         [input] nxK matrix of forecasts of ensemble members         %
%   y         [input] nx1 vector of verifying data                        %
%   options   [input] structure with BMA algorithmic variables            %
%   func_name [input] name (string) of BMA function                       %
%    = 'norm'   1. [lowup,mix_norm] = BMA_norm(x,D,options)               %
%    = 'quant'  2. [pred,pdf_y,cdf_y] = BMA_quantile(x,D,y,options)       %
%    = 'crps'   3. [crps] = BMA_crps(x,D,y,options)                       %
%    = 'lik'    4. [ell,L] = BMA_lik(x,D,Y,options)                       %
%    = 'mix'    5. [X,pdf,cdf,pred] = BMA_mixture(x,D,T,options)          %
%    = 'draw'   6. [rnd,mu_mix,sigma2_mix] = BMA_draw(x,D,N,options)      %
%   T         [input] indices of 1:n at which we return BMA pdf & cdf     %
%   varargout [outpt] return arguments                                    %
%    = [ell,L] if func_name = 'lik'                                       %
%    = [plimits,pdf_y,cdf_y] if func_name = 'quant'                       %
%    = [lowup,mix_norm] if func_name = 'norm';                            %
%    = [X,pdf_BMA,cdf_BMA,plimits] if func_name = 'mix'                   %
%    = [rnd,mu_mix,sigma2_mix] if func_name = 'draw'                      %
%    = [crps] if func_name = 'crps'                                       %
%                                                                         %
% EXAMPLE USAGE                                                           %
%  [ell,L] = BMA_model(x,D,y,options,'lik');                              %
%  [plimits,pdf_y,cdf_y] = BMA_model(x,D,y,options,'quant');              %
%  [lowup,mix_norm] = BMA_model(x,D,y,options,'norm');                    %
%  [X,pdf_BMA,cdf_BMA,plimits] = BMA_model(x,D,y,options,'mix',T);        %
%  [rnd,mu_mix,sigma2_mix] = BMA_model(x,D,y,options,'draw',N);           %
%  [crps] = BMA_model(x,D,y,options,'crps');                              %
%                                                                         %
% (c) Written by Jasper A. Vrugt, Feb 2006                                %
% University of California Irvine                                         %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

if nargin < 5
    error('BMA_model: method not defined. Use help function');
end
if nargin == 5
    if strcmpi(return_arg,'mix')
        error(['BMA_model:Undefined input vector T with indices of y ' ...
            'at which BMA pdf and cdf are computed'])
    end
    if strcmpi(return_arg,'draw')
        warning(['BMA_model:Undefined number of draws N of BMA ' ...
            'mixture distribution. N = 1e4 (DEFault)']);
    end
end    
if nargin == 6
    switch return_arg
        case 'mix'
            T = varargin{1};
        case 'draw'
            N = varargin{1};
        otherwise
    end
end

[n,K] = size(D);                    % Unpack number of observations
beta = x(1:K)';                     % Unpack weights
count = 0;                          % Set counter to zero
mu_mix = D * beta;                  % Mean of mixture (= exact)

% % % Now calculate the bias-corrected forecasts
% % for k = 1 : K
% %     D(:,k) = a(k) + b(k) * D(:,k);
% % end

g = 1-options.alpha;                % Unpack alpha values
gam = sort([(1-g)/2,(1-g)/2 + g]);  % Confidence levels and sort quantiles 
nG = numel(gam);                    % # confidence levels

switch options.VAR                  % Variance option
    case {'1'}                      % 1: common constant variance
        S = x(K+1) * ones(n,K);
    case {'2'}                      % 2: individual constant variance
        S = bsxfun(@times,x(K+1:2*K),ones(n,K));
    case {'3'}                      % 3: common non-constant variance
        c = x(K+1); S = c * abs(D);
    case {'4'}                      % 4: individual non-constant variance
        c = x(K+1:2*K); S = bsxfun(@times,c,abs(D));
    otherwise
        error('MODELAVG:BMA_quantile','Unknown variance option');
end

S = max(S,eps);                     % each entry nxK-matrix S >= 2.22e-16
Y = repmat(y,1,K);                  % replicate verifying data nxK matrix

switch options.PDF                  % Conditional distribution

    case {'normal'}                 % normpdf(Mu,S)
        [sigma2_mix,CDF_str,Mu,S] = ...         % Info BMA mixture PDF
            write_str(beta,D,n,K,S, ...
            mu_mix,options);
        L = pdf('normal',Y,Mu,S);               % BMA likelihood 
        
    case {'lognormal'}              % lognpdf(Mu,S)
        [sigma2_mix,CDF_str,Mu,S] = ...         % Info BMA mixture PDF
            write_str(beta,D,n,K,S, ...
            mu_mix,options);
        L = pdf('logn',Y,Mu,S);                 % BMA likelihood 
        
    case {'tnormal'}                % tnormpdf(Mu,S,a,b)
        [sigma2_mix,CDF_str,Mu,S,a,b] = ...     % Info BMA mixture PDF
            write_str(beta,D,n,K,S, ...
            mu_mix,options);
        L = tnormpdf(Y,Mu,S,a,b);               % BMA likelihood 
        
    case {'gen_normal'}             % gnormpdf(Mu,S,tau)
        [sigma2_mix,CDF_str,Mu,S,tau] = ...     % Info BMA mixture PDF
            write_str(beta,D,n,K,S, ...
            mu_mix,options);
        L = gnormpdf(Y,Mu,S,tau);               % BMA likelihood 
        
    case {'gamma'}                  % gampdf(A,B)
        [sigma2_mix,CDF_str,A,B] = ...          % Info BMA mixture PDF
            write_str(beta,D,n,K,S, ...
            mu_mix,options);
        L = pdf('gamma',Y,A,B);                 % BMA likelihood 
        
    case {'weibull'}                % wblpdf(Lambda,Kk)
        [sigma2_mix,CDF_str,Lambda,Kk] = ...    % Info BMA mixture PDF        
            write_str(beta,D,n,K,S, ...
            mu_mix,options);
        L = pdf('wbl',Y,Lambda,Kk);             % BMA likelihood 
        
    case {'gev'}                    % gevpdf(xi,S,Mu)
        [sigma2_mix,CDF_str,xi,S,Mu] = ...      % Info BMA mixture PDF
            write_str(beta,D,n,K,S, ...
            mu_mix,options);
        L = pdf('gev',Y,xi,S,Mu);               % BMA likelihood 
        
    case {'gpareto'}                % gppdf(kk,S,Theta)
        [sigma2_mix,CDF_str,kk,S,Theta] = ...   % Info BMA mixture PDF
            write_str(beta,D,n,K,S, ...
            mu_mix,options);
        L = pdf('gp',Y,kk,S,Theta);             % BMA likelihood 

end

% BMA_norm
switch return_arg

    case 'lik'  % [ell,L] = BMA_lik(x,D,Y,options)
        lik = L*beta + realmin;                 % (nx1)-vector lik. BMA @ y
        ell = sum(log(lik));                    % Log-likelihood BMA model
        varargout = {ell,L};

    case 'mix'  % [X,pdf_BMA,cdf_BMA,plimit] = BMA_mixture(x,D,T,options)
        m = numel(T);                           % # time entries of T
        N = 1e5;                                % # discretized x-values
        X = linspace(0,2*max(max(D)),N);        % Create x-values
        [pdf_BMA,cdf_BMA] = deal(nan(m,N));     % Init pdf/cdf obsrvations
        PDF_str = CDF_str;                      % Copy mixture pdf calld fx
        Fx = strcat(CDF_str,'- gam;'); 
        eval(char(Fx));                         % Activate mixture cdf
        PDF_str = strcat(PDF_str,';'); 
        PDF_str = strrep(PDF_str, ...           % Activate mixture pdf
            'cdf','pdf');
        PDF_str = strrep(PDF_str,',gam',''); 
        fx = strrep(PDF_str,'Fx','fx');
        eval(fx);
        for i = 1:m
            pdf_BMA(i,1:N) = fx(X,T(i));        % Compute pdf/cdf x-values
            cdf_BMA(i,1:N) = Fx(X,T(i),0);
        end
        plimit = nan(m,nG);                     % Prediction limits
        for i = 1 : m
            if mod(i,floor(m/25)) == 0          % Print progress
                if i > 1
                    fprintf(1, ...              % delete line before
                        repmat('\b',1,count)); 
                    count = fprintf(['BMA mixture calculation, ' ...
                        '%% done: %3.2f'],100*(i/m));
                end
            end
            for z = 1 : nG
                plimit(i,z) = fzero(@(x)...     % Root finding cdf exact
                    Fx(x,T(i),gam(z)), ...
                    abs(mean(D(T(i),1:K))));
            end
        end
        fprintf(1,'\n'); varargout = {X,pdf_BMA,cdf_BMA,plimit};

    case 'norm' % [lowup,mix_norm] = BMA_norm(x,D,options)
        method = 'robust';                      % α-norm built-in integrl
        PDF_str = CDF_str;                      % Copy mixture cdf 
        Fx = strcat(CDF_str,'- gam;');          % Add "-gam" to cdf
        eval(char(Fx));
        PDF_str = strcat(PDF_str,';'); 
        PDF_str = strrep(PDF_str, ...           % Activate mixture pdf
            'cdf','pdf');
        PDF_str = strrep(PDF_str, ...           % delta, δ-norm
            ',gam',',delta'); 
        PDF_str = strrep(PDF_str,'Fx','fx');
        PDF_str = strrep(PDF_str, ...
            '@(x,i,delta)','@(x,i,delta) (');
        switch options.PDF
            case {'tnormal'}
                fx = strrep(PDF_str,');',') ).^delta;');
            otherwise
                fx = strrep(PDF_str,'));',')) ).^delta;');
        end
        eval(fx);                               % Activate BMA mixture pdf
        a = 1e-4;                               % Between cdf = a & 1-a
        gam = sort([a,1-a]); nA = numel(gam);   % Lower and upper quantiles
        
        ny = 1e3; lowup = nan(n,nA);            % # equidistant samples
        mix_norm = nan(n,2);                    
        
        for i = 1 : n                           % Lower & upper end BMA mix
            if mod(i,floor(n/25)) == 0          % Print progress
                if i > 1
                    fprintf(1, repmat('\b', ... % delete line before
                        1,count)); 
                    count = fprintf(['BMA α-norm calculation, ' ...
                        '%% done: %3.2f'],100*(i/n));
                end
            end
            for z = 1 : nA
                lowup(i,z) = fzero(@(y) ...     % Exact using root finding
                    Fx(y,i,gam(z)), ...
                    abs(mean(D(i,1:K))));
            end
            switch method
                case 'standard'
                    yi = linspace( ...          % ny equally spacd low & up
                        lowup(i,1),lowup(i,2),ny)';
                    pdf_yi = fx(yi,i);          % Evaluate BMA mixtre at yi
                    mix_norm(i,1) = trapz(...   % δ-norm BMA mixture: δ=1
                        yi,pdf_yi);             %=sum(pdf_yi)*(yi(2)-yi(1))
                    mix_norm(i,2) = sqrt( ...   % δ-norm BMA mixture: δ=2 
                        trapz(yi,pdf_yi.^2));   % quadratic & sphrical rule                    
                case 'robust'
                    for delta = 1:2             % δ-norm BMA mixture
                        mix_norm(i,delta) = ...
                            integral(@(y) ...
                            fx(y,i,delta), ...
                            lowup(i,1), ...
                            lowup(i,2)).^(1/delta);
                    end
            end
        end
        fprintf(1,'\n'); varargout = {lowup,mix_norm};

    case 'draw' % [rnd,mu_mix,sigma2_mix] = BMA_draw(x,D,N,options)
        
        rnd = nan(n,N);                         % Initialze random draw pdf
        r = rand(N,1); w = cumsum(beta');       % draw random numbers
        id_mix = K - sum(r<w,2) + 1;            % component to draw from
        id = cell(1,K); n_id = nan(1,K);        % # points each component
        for k = 1:K
            id{k} = find(id_mix==k); 
            n_id(k) = numel(id{k});
        end
        for t = 1:n                             % Now sample
            for k = 1:K
                switch options.PDF              % Conditional distiribution
                    case {'normal'}             % µ Mu & std. S
                        rnd(t,id{k}) = ...
                            normrnd(Mu(t,k), ...
                            S(t,k),[1 n_id(k)]);
                    case {'lognormal'}          % µ Mu & std. S
                        rnd(t,id{k}) = ...
                            lognrnd(Mu(t,k), ...
                            S(t,k),[1 n_id(k)]);
                    case {'gamma'}              % shape A & scale B
                        rnd(t,id{k}) = ...
                            gamrnd(A(t,k), ...
                            B(t,k),[1 n_id(k)]);
                    case {'weibull'}            % shape Kk & scale Lambda
                        rnd(t,id{k}) = ...
                            wblrnd(Lambda(t,k), ...
                            Kk(t,k),[1 n_id(k)]);
                    case {'gen_normal'}         % µ Mu, std. S, shape tau
                        rnd(t,id{k}) = ...
                            gnormrnd(Mu(t,k), ...
                            S(t,k),tau(k), ...
                            [1 n_id(k)]);
                    case {'gev'}                % µ Mu, std. S, shape xi
                        rnd(t,id{k}) = ...
                            gevrnd(xi(t,k), ...
                            S(t,k),Mu(t,k), ...
                            [1 n_id(k)]);
                    case {'gpareto'}            % µ Mu, std. S, shape xi
                        rnd(t,id{k}) = ...
                            gprnd(kk(t,k), ...
                            S(t,k),Theta(t,k), ...
                            [1 n_id(k)]);
                    case {'tnormal'}            % µ Mu, std. S, int. a & b
                        rnd(t,id{k}) = ...
                            tnormrnd(Mu(t,k), ...
                            S(t,k),a,b,[1 n_id(k)]);
                end
            end
        end
        varargout = {rnd,mu_mix,sigma2_mix};

    case 'crps' % crps = BMA_crps(x,D,y,options)
        crps = nan(n,1);                        % Initilze CRPS BMA mixture
        P1 = [0.001 0.002 0.005 0.01 ...        % Note: we use "p" for tau
            0.02 0.05:0.05:0.45];
        P = [P1,0.5,sort(1 - P1)]';         
        nP = numel(P);
        F_P = strcat(CDF_str,'- gam;'); 
        eval(char(F_P))                         % Mixture cdf with p (tau)
        % Step 2: Solve for values of y so that F_P^{-1}(y) = P
        % Step 3: Compute CRPS
        FinvP = nan(nP,1);
        for i = 1 : n
            if mod(i,floor(n/25)) == 0          % Print progress
                if i > 1
                    fprintf(1, repmat('\b', ... % delete line before
                        1,count)); 
                    count = fprintf(['BMA CRPS calculation, %% done: ' ...
                        '%3.2f'],100*(i/n));
                end
            end
            for z = 1 : nP                      % Root finding [exact cdf]
                switch z
                    case 1      
                        FinvP(z,1) = ...        % F_P(a,b,c) --> Fx(a,b,c)
                            fzero(@(y) ...
                            Fx(y,i,P(z)), ...
                            abs(mean(D(i,1:K))));                   %#ok
                    otherwise  
                        FinvP(z,1) = ...        % F_P(a,b,c) --> Fx(a,b,c)    
                            fzero(@(y) ...
                            Fx(y,i,P(z)), ...
                            FinvP(z-1,1));                          %#ok
                end
            end
            Fy = Fx(y(i,1),i,0); ii = P > Fy;   % Evaluate mixture at omega: 
                                                % F_P(a,b,c) --> Fx(a,b,c)
            if sum(ii) == 0                     % At least 1 element of ii
                % Area right of observation 
                % almost 0 (omega close to P = 1)
                %  P_temp = .5*(1+Fw); 
                % FinvP_temp = fzero(@(x) ...
                %     F_P(x,i,P_temp),Fw);
                %  intgral2 = - 2*trapz( ...
                %     [Fw ; P_temp],[Y(i,1) ; ...
                %     FinvP_temp])
                % This integral is very close 
                % to zero --> discard
                % Compute CRPS using trapezoidal 
                % rule (left from omega)
                crps(i,1) = y(i,1)*(1 - 2*Fy) + 2*trapz(P,P.*FinvP);
            else
                % Compute CRPS using trapezoidal integration rule
                crps(i,1) = y(i,1)*(1 - 2*Fy) + 2*trapz(P,P.*FinvP) ...
                    - 2*trapz([Fy ; P(ii)],[y(i,1) ; FinvP(ii)]);
            end
        end
        fprintf(1,'\n'); varargout = {crps};

    case 'quant' % [plimit,pdf_y,cdf_y] = BMA_quantile(x,D,y,options)
        g = 1-options.alpha;                    % Unpack alpha values
        gam = sort([(1-g)/2,(1-g)/2 + g]);      % Sort confidence levels
        nA = numel(gam);                        % # gamma values?
        [pdf_y,cdf_y] = deal(nan(n,1));         % Initialize pdf/cdf at y
        PDF_str = CDF_str;                      % Copy to create mixtre pdf
        Fx = strcat(CDF_str,'- gam;'); 
        eval(char(Fx));                         % Activate mixture cdf 
                                                % with alfa for root findng
        PDF_str = strcat(PDF_str,';'); 
        PDF_str = strrep(PDF_str, ...
            'cdf','pdf');
        PDF_str = strrep(PDF_str,',gam',''); 
        fx = strrep(PDF_str, ...
            'Fx','fx'); eval(fx);               % Activate mixture pdf
        for i = 1:n
            pdf_y(i) = fx(y(i),i);              % pdf @ measured values, y
            cdf_y(i) = Fx(y(i),i,0);            % cdf @ measured values, y
        end
        plimit = nan(n,nA);                     % Determine BMA quantiles
        for i = 1 : n
            if mod(i,floor(n/25)) == 0          % Print progress
                if i > 1
                    fprintf(1, repmat('\b', ... % delete line before
                        1,count)); 
                    count = fprintf(['BMA quantile calculation, ' ...
                        '%% done: %3.2f'],100*(i/n));
                end
            end
            for z = 1 : nA                      % Root finding [exact cdf]
                plimit(i,z) = fzero(@(x) Fx(x,i,gam(z)), ...
                    abs(mean(D(i,1:K))));
            end
        end
        fprintf(1,'\n'); varargout = {plimit,pdf_y,cdf_y};
end

end
