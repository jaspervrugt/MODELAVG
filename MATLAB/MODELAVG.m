function varargout = MODELAVG(method,D,y,options)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%                                                                       %%
%% MM      MM  OOOOOO  DDDDDDD  EEEEEEE LL        AAAA   VV   VV GGGGGGG %%
%% MMM     MM OOOOOOOO DDDDDDDD EEEEEEE LL       AA  AA  VV   VV GG   GG %%
%% MMMM  MMMM OO    OO DD    DD EE      LL       AA  AA  VV   VV GG   GG %%
%% MM MMMM MM OO    OO DD    DD EEEE    LL       AA  AA  VV   VV GGGGGGG %%
%% MM  MM  MM OO    OO DD    DD EEEE    LL       AAAAAA   VV VV  GGGGGGG %%
%% MM      MM OO    OO DD    DD EE      LL      AA    AA   VVV        GG %%
%% MM      MM OOOOOOOO DDDDDDDD EEEEEEE LLLLLLL AA    AA   VVV       GGG %%
%% MM      MM  OOOOOO  DDDDDDD  EEEEEEE LLLLLLL AA    AA    V    GGGGGGG %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%                                                                       %%
%% Multi-model averaging is widely used in various scientific and        %%
%% engineering disciplines to post-process forecast ensembles and/or     %%
%% quantify conceptual model uncertainty. This MATLAB toolbox, called    %%
%% MODELAVG, implements a suite of different model averaging techniques, %%
%% including (amongst others) equal weights averaging (EWA),             %%
%% Bates-Granger model averaging (BGA), Bayesian model averaging (BMA),  %%
%% Mallows model averaging (MMA), and Granger-Ramanathan averaging (GRA) %%
%% For BMA the user can select among different options for the           %%
%% conditional forecast distribution of each individual member of the    %%
%% ensemble. Options include a normal distribution with homoscedastic    %%
%% and heteroscedastic varaiance, a gamma distribution with constant or  %%
%% heteroscedastic variance, and a generalized normal distribution with  %%
%% homoscedastic or non-constant variance. The toolbox returns the       %%
%% optimal values posterior distribution of the weights and related      %%
%% parameters of each model averaging method, along with graphical       %%
%% output of the results. For MMA and BMA the user has access to the     %%
%% entire posterior distribution of the weights and/or variances derived %%
%% from MCMC simulation with the DREAM algorithm. Three case studies     %%
%% involving forecast ensembles of hydrologic and meteorologic models    %%
%% are used to illustrate the capabilities and functionalities of the    %%
%% MODELAVG toolbox                                                      %%
%%                                                                       %%
%% SYNOPSIS: [beta] = MODELAVG(method,D,y)                               %%
%%           [beta,output] = MODELAVG(method,D,y,options)                %%
%%  where                                                                %%
%%   method    [input] String with model averaging used                  %%
%%   D         [input] nxK matrix forecasts ensemble members             %%
%%   y         [input] nx1 vector with verifying observations            %%
%%   options   [input] structure BMA algorithmic variables               %%
%%         .PDF        string: conditional PDF for BMA method (MANUAL)   %%
%%                      = 'normal'     = normal distribution             %%
%%                      = 'lognormal'  = lognormal distribution          %%
%%                      = 'tnormal'    = truncated normal ( > 0)         %%
%%                      = 'gen_normal' = generalized normal              %%
%%                      = 'gamma'      = gamma distribution              %%
%%                      = 'weibull'    = Weibull distribution [!]        %%
%%                      = 'gev'        = generalized extreme value       %%
%%                      = 'gpareto'    = generalized Pareto              %%
%%         .VAR        string: variance treatment BMA method (MANUAL)    %%
%%                      =  '1'          = constant group variance        %%
%%                      =  '2'          = constant individual variance   %%
%%                      =  '3'          = nonconstant group variance     %%
%%                      =  '4'          = nonconstant ind. variance      %%
%%         .TAU        string: trtmnt shape par gnrlzd nrml BMA method   %%
%%                      =  '1'          = common (= group) value of tau  %%
%%                      =  '2'          = individual value of tau        %%
%%         .alpha      1xna vector (1-sign. level) conf & pred limits    %%
%%                      = [0.5 0.7 0.9 0.95 0.99]                        %%
%%         .print      string: print results to screen or not            %%
%%                      = 'yes'                                          %%
%%                      = 'no'                                           %%
%%         .CPU        string: CPU-efficient solution or not             %%
%%                      = 'yes'                                          %%
%%                      = 'no'                                           %%
%%         .p          1 x K vector: with # parameters of models         %%
%%                      = [10 5 2 9 18 ...]                              %%
%%         .postproc   string: postprocess results BMA method or not     %%
%%                      = 'yes'                                          %%
%%                      = 'no' (= faster for ML parameters BMA method)   %%
%%   beta      [outpt] 1xK vector weights/array chain samples BMA/MMA    %%
%%   output    [outpt] structure of results DREAM BMA/MMA                %%
%%         .AR         qx2 matrix sample # & DREAM_ZS acceptance rate    %%
%%         .R_stat     qx(d+1) matrix sample # & R-hat conv. diag.       %%
%%         .MR_stat    qx2 matrix multivariate R-hat conv. diag.         %%
%%         .RunTime    Run time in seconds of DREAM_ZS algorithm         %%
%%         .CORR       Correlation matrix  posterior parameter samples   %%
%%         .loglik     BMA log-likelihood maximum likelihood parameters  %%
%%         .std        Posterior standard deviations of BMA parameters   %%
%%         .par_unc    nx2na matrix lower/upper confidence limits alpha  %%
%%         .pred       nx2na matrix lower/upper prediction limits alpha  %%
%%         .pdf_Y      nx1 vector with BMA mixture pdf at y              %%
%%         .cdf_Y      nx1 vector with BMA mixture cdf at y              %%
%%         .coverage   Coverage (%) 100alpha pred. intervals BMA model   %%
%%         .spread     Spread of 100alpha pred. intervals BMA model      %%
%%         .Ye         nx1 vector with µ BMA mixture forecast [exact]    %%
%%         .R2         coefficient detrm. weighted-average BMA forecast  %%
%%         .RMSE       Root mean square err weighted-average BMA forcst  %%
%%         .R          Pearson coef. meas & weighted-average BMA forcst  %%
%%         .ML         maximum likelihood BMA weights & shape parametrs  %%
%%         .CRPS       nx1 vector continuous ranked probability score    %%
%%         .QS         nx1 vector quadratic score BMA forecst dstrbtion  %%
%%         .LS         nx1 vector logrithmc score BMA forecst dstrbtion  %%
%%         .SS         nx1 vector spherical score BMA forecst dstrbtion  %%
%%         .mRLBL_anal reliability of BMA mixture distribution           %%
%%         .mCv_anal   time-averaged coef. var. BMA mixture [= exact]    %%
%%         .mQS        time-averaged quadratic score BMA forecast dis.   %%
%%         .mLS        time-averaged logrithmic score BMA forecast dis.  %%
%%         .mSS        time-averaged spherical score BMA forecast dis.   %%
%%         .mCRPS      time-averaged continuous rankd probability score  %%
%%         .mix_norm   nx2 matrix 1-norm & 2-norm BMA mixture PDF        %%
%%         .mu_mix     nx1 vector mean BMA mixture forecast [= exact]    %%
%%         .var_mix    nx1 vector variance BMA mixture forcst [= exact]  %%
%%         .KGE        Kling-Gupta eff. weighted-average BMA forecast    %%
%%         .str_table  1 x d cell names of BMA parameters (Latex)        %%
%%         .RMSE_mod   Root mean square err. models ensemble frecasts D  %%
%%                                                                       %%
%% This general purpose MATLAB code is designed to postprocess forecast  %%
%% ensembles of dynamic simulation models. The user can select among 7   %%
%% diferent model averaging methods including                            %%
%%      1. Equal weights averaging                                       %%
%%      2. Bates-Granger averaging                                       %%
%%      3. Information criterion averaging                               %%
%%         Akaike (AIC) & Bayes (BIC)                                    %%
%%      4. Granger-Ramanathan averaging                                  %%
%%      5. Bayesian model averaging                                      %%
%%      6. Mallows model averaging (with/without weights unit Simplex)   %%
%%                                                                       %%
%% For BMA the user can select among different options for the           %%
%% conditional forecast distribution of the ensemble members. This       %%
%% includes the                                                          %%
%%          a. normal distribution                                       %%
%%          b. lognormal distribution                                    %%
%%          c. generalized normal distribution                           %%
%%          d. truncated normal distribution                             %%
%%          e. gamma distribution                                        %%
%%          f. Weibull distribution                                      %%
%%          g. generalized extreme value distribution                    %%
%%          h. generalized Pareto distribution                           %%
%% For each conditional PDF the user can determine whether to implement  %%
%% homoscedastic (= constant) or heteroscedastic (= nonconstant)         %%
%% variance with a single or group treatment                             %%
%%                                                                       %%
%% THIS TOOLBOX HAS BEEN DESCRIBED IN                                    %%
%%   Vrugt, J.A. (2023), MODELAVG: A MATLAB toolbox for postprocessing   %%
%%       of model ensembles, Manual, Version 2.0, pp. 1 - XX, 2023       %%
%% FOR MORE INFORMATION, PLEASE READ                                     %%
%%   Vrugt, J.A. (2024), Distribution-Based Model Evaluation and         %%
%%       Diagnostics: Elicitability, Propriety, and Scoring Rules for    %%
%%       Hydrograph Functionals, Water Resources Research, 60,           %%
%%       e2023WR036710, https://doi.org/10.1029/2023WR036710             %%
%%   Vrugt, J.A. (2015), Markov chain Monte Carlo simulation using the   %%
%%       DREAM software package: Theory, concepts, and MATLAB            %%
%%       implementation, Environmental Modeling and Software, 75,        %%
%%       pp. 273-316                                                     %%
%%   Diks, C.G.H., and J.A. Vrugt (2010), Comparison of point forecast   %%
%%       accuracy of model averaging methods in hydrologic applications, %%
%%       Stochastic Environmental Research and Risk Assessment, 24(6),   %%
%%       809-820, https://doi.org/10.1007/s00477-010-0378-z              %%
%%   Vrugt, J.A., C.G.H. Diks, and M.P. Clark (2008), Ensemble Bayesian  %%
%%       model averaging using Markov chain Monte Carlo sampling,        %%
%%       Environmental Fluid Mechanics, 8(5-6), 579-595,                 %%
%%           https://doi.org/10.1007/s10652-008-9106-3                   %%
%%   Vrugt, J.A., C.J.F. ter Braak, M.P. Clark, J.M. Hyman and B.A.      %%
%%       Robinson (2008), Treatment of input uncertainty in hydrologic   %%
%%       modeling: Doing hydrology backward with Markov chain Monte      %%
%%       Carlo simulation, 44 (12), https://doi.org/10.1029/2007WR006720 %%
%%   Vrugt, J.A., and B.A. Robinson (2007), Treatment of uncertainty     %%
%%       using ensemble methods: Comparison of sequential data           %%
%%       assimilation and Bayesian model averaging, Water Resources      %%
%%       Research, 43, W01411, https://doi.org/10.1029/2005WR004838      %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%                                                                       %%
%% COPYRIGHT (c) 2015  the author                                        %%
%%                                                                       %%
%%   This program is free software: you can modify it under the terms    %%
%%   of the GNU General Public License as published by the Free Software %%
%%   Foundation, either version 3 of the License, or (at your option)    %%
%%   any later version                                                   %%
%%                                                                       %%
%%   This program is distributed in the hope that it will be useful, but %%
%%   WITHOUT ANY WARRANTY; without even the implied warranty of          %%
%%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU    %%
%%   General Public License for more details                             %%
%%                                                                       %%
%%   (c) Written by Jasper A. Vrugt, Feb 2012                            %%
%%   University of California Irvine                                     %%
%%                                                                       %%
%%   Version 1: March 2015                                               %%
%%   Version 2: March 2022                                               %%
%%   Version 3: June 2024                                                %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%                                                                       %%
%% BUILT-IN CASE STUDIES                                                 %%
%%   Example 1   24-hour forecasts of river discharge                    %%
%%   Example 2   48-forecasts of sea surface temperature                 %%
%%   Example 3   48-forecasts of sea surface pressure                    %%
%%   Example 11  Forecasts of water levels                               %%
%%   Example 12  Hydrologic modeling                                     %%
%%   Example 13  Flood modeling                                          %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

%% Check:  http://faculty.sites.uci.edu/jasper
%% Papers: http://faculty.sites.uci.edu/jasper/publications/
%% Google Scholar: https://scholar.google.com/citations?user=zkNXecUAAAAJ&hl=nl

% If less than three input arguments --> produce error
if nargin < 3, error('MODELAVG:TooFewInputs', ['Requires at least ' ...
        'three input arguments.']); end
% If less than 3 input arguments --> define DE as structure
if nargin < 4, options = struct; end
% Now call MODELAVG_check to evaluate setup user
[method,options] = MODELAVG_check(method,D,y,options);
% Now call MODELAVG_setup to define computational framework
[beta,chain,n,K,sigma_2,log_L,par,count,options] = ...
    MODELAVG_setup(method,D,y,options);
% Now start the timer ( to calculate run time )
tic;

% Now switch between the different model averaging methods
switch method
    
    case 'ewa'              % equal weights averaging
        beta = 1/K * ones(1,K);
    case 'bga'              % Bates-Granger averaging
        beta = (1./sigma_2) ./ sum( 1./sigma_2 );
    case 'aica'             % information criterion averaging - AIC
        % Calculate information criterion
        I = -2*log_L + 2 * options.p;
        % Now adjust for underflow
        I = I - min(I);
        % Calculate weights
        beta = exp( -I/2 ) / sum( exp( -I/2 ) );
    case 'bica'             % information criterion averaging - BIC
        % Calculate information criterion
        I = -2 * log_L + options.p * log( n );
        % Now adjust for underflow
        I = I - min(I);
        % Calculate weights
        beta = exp( -I/2 ) ./ sum( exp( -I/2 ) );
    case 'gra'              % Granger - Ramanathan averaging
        % Calculate the weights (least squares type approach)
        %beta = inv( D' * D ) * D' * y
        beta = D'*D\D'*y;
        % Now tranpose to horizontal vector
        beta = beta';
    case {'bma'}            % Bayesian model averaging
        % Now check implementation of conditional PDF
        switch options.VAR
            case '1'    % 1: common constant variance
                d = K + 1;                 
                if strcmp(options.PDF,{'lognormal'})
                    par.max = [ ones(1,K) std(y) ];
                else
                    par.max = [ ones(1,K) 10*std(y) ];
                end
            case '2'    % 2: individual constant variance
                d = 2 * K; 
                if strcmp(options.PDF,{'lognormal'})
                    par.max = [ ones(1,K) std(y)*ones(1,K) ];
                else
                    par.max = [ ones(1,K) 10*std(y)*ones(1,K) ];
                end
            case {'3'}  % 3: common non-constant variance
                d = K + 1; 
                if strcmp(options.PDF,{'weibull'})
                    par.max = [ ones(1,K) 25 ];
                elseif sum(strcmp(options.PDF,{'lognormal','gev', ...
                        'gpareto'}))
                    par.max = [ones(1,K) 1];
                else
                    par.max = [ones(1,K) 2];
                end
            case {'4'}  % 4: individual non-constant variance
                d = 2 * K; 
                if strcmp(options.PDF,{'weibull'})
                    par.max = [ ones(1,K) 25*ones(1,K) ];
                elseif sum(strcmp(options.PDF,{'lognormal','gev', ...
                        'gpareto'}))
                    par.max = [ones(1,K) ones(1,K)];
                else
                    par.max = [ones(1,K) 2*ones(1,K)];
                end
        end
        if strcmp(options.PDF,{'gen_normal'})
            switch options.TAU
                case {'1'}
                    par.max = [par.max 20]; d = d + 1;
                case {'2'}
                    par.max = [par.max 20*ones(1,K)]; d = d + K;
            end
        end
        if sum(strcmp(options.PDF,{'gev','gpareto'}))
            switch options.TAU
                case {'1'} % xi cannot exceed 1/2 --> variance not defined
                    par.max = [par.max 0.5]; d = d + 1; l = 0;
                case {'2'}
                    par.max = [par.max ones(1,K)/2]; d = d + K; l = K - 1;
            end
        end
        % Weights must be on simplex, otherwise pred. density meaningless
        par.min = zeros(1,d);
        % Overwrite values of zero of lower bound of xi
        if strcmp(options.PDF,'gpareto')
            par.min(d-l:d) = -1.0;
        end
        % For GEV: zero implies that skew cannot become negative
        % This avoids numerical instabilities and convergence problems for
        % some of the shape parameters [= otherwise very difficult to
        % solve]
        Y = repmat(y,1,K); % copy the verifying data K times
        % Problem settings defined by user
        log_pdf = @(x) BMA_lik(x,D,Y,options); % Function handle for DREAM
        % Run basic implementation of the DREAM_ZS algorithm
        code_to_use = 'dream_zs';
        switch code_to_use
            case 'dream_zs'
                T = d * 2500; %e4;              % # generations
                N = 3;                          % # chains
                [chain,output] = MODELAVG_dream_zs(par,log_pdf,N,T,d,K,1);
            case 'dream'
                T = d * 1e3;                    % # generations
                N = 20;                         % # chains
                [chain,output] = MODELAVG_dream(par,log_pdf,N,T,d,K,1);
        end
    case {'mma','mma-s'}    % Mallows model averaging
        d = K;                                  % # parameters
        T = d * 1e4;                            % # generations
        N = 3;                                  % # chains
        % Calculate error variance of each model
        sigma_2 = sum(bsxfun(@minus,D,y).^2)/n;
        % Most complex model of ensemble?
        [~,idx] = max(options.p); var_err = sigma_2( idx(1) );
        % Function handle for DREAM
        log_pdf = @(x) MMA_lik(x,D,y,var_err,options.p);     
        if strcmp(method,'mma')
            % Define parameter ranges ( only for initial sampling )
            par.min = -ones(1,K); par.max = ones(1,K);
            % Run basic implementation of the DREAM_ZS algorithm
            [chain,output] = MODELAVG_dream_zs(par,log_pdf,N,T,d,K,0);
        elseif strcmp(method,'mma-s')
            % Define parameter ranges ( weights on Simplex )
            par.min = zeros(1,K); par.max = ones(1,K);
            % Run basic implementation of the DREAM_ZS algorithm
            [chain,output] = MODELAVG_dream_zs(par,log_pdf,N,T,d,K,1);
        end
end

% Calculate wall time of MODELAVG code
output.RunTime = toc;
% Store the RMSE of the forecast error of each model
output.RMSE_mod = sqrt(sigma_2);

% Now create the output that is printed in postprocessor
[beta,output,str_plot] = MODELAVG_end(method,D,y,beta,chain, ...
    options,output);
% Print progress
if strcmp(options.postproc,'yes')
    for t = 2:3
        if t == 2, print_name = '........'; end
        if t == 3, print_name = '........ done'; end
        if ( t > 1 )
            fprintf(1, repmat('\b',1,count)); % Delete line before
            count = fprintf('MODELAVG postprocessor, please wait %s',...
                print_name );
            if t == 2
                MODELAVG_postproc(method,D,y,beta,sigma_2,par,chain,...
                    options,output,str_plot);
            end
        end
    end
    fprintf(1,'\n');
end

% Return output argument
varargout = { beta output };