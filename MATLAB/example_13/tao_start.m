%% function varargout = MODELAVG(method,D,Y,options)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%                                                                                    %%
%% MM       MM  OOOOOOO  DDDDDDDD  EEEEEEEE LL           AAA     VV       VV GGGGGGGG %%
%% MMM      MM OOOOOOOOO DDDDDDDDD EEEEEEEE LL          AA AA    VV       VV GG    GG %%
%% MMMM   MMMM OO     OO DD     DD EE       LL          AA AA     VV     VV  GG    GG %%
%% MM MM MM MM OO     OO DD     DD EEEEE    LL         AA   AA    VV     VV  GGGGGGGG %%
%% MM  MMM  MM OO     OO DD     DD EEEEE    LL        AAAAAAAAA    VV   VV   GGGGGGGG %%
%% MM       MM OO     OO DD     DD EE       LL        AA     AA     VV VV          GG %%
%% MM       MM OOOOOOOOO DDDDDDDDD EEEEEEEE LLLLLLLL AA       AA    VV VV         GGG %%
%% MM       MM  OOOOOOO  DDDDDDDD  EEEEEEEE LLLLLLLL AA       AA     VVV     GGGGGGGG %%
%%                                                                                    %%
%% SYNOPSIS: [beta] = MODELAVG(method,D,Y)                                            %%
%%           [beta,output] = MODELAVG(method,D,Y,options)                             %%
%%  where                                                                             %%
%%   method    [input]  REQUIRED: 1 x d vector with maximum likelihood BMA parameters %%
%%   D         [input]  REQUIRED: n x K matrix with forecasts of ensemble members     %%
%%   Y         [input]  REQUIRED: n x 1 vector with verifying observations            %%
%%   options   [input]  REQUIRED: structure with BMA algorithmic variables            %%
%%   beta      [output] 1 x K vector with weights / array chain trajectories BMA/MMA  %%
%%   output    [output] structure with results DREAM for BMA and MMA (fewer fields)   %%
%%         .AR          q x 2 matrix with sample # and DREAM_ZS acceptance rate       %%     
%%         .R_stat      q x d+1 matrix with sample # and R-hat scale reduction factor %%
%%         .MR_stat     q x 2 matrix with multivariate R-hat scale reduction factor   %%
%%         .RunTime     Run time in seconds of the DREAM_ZS algorithm                 %%
%%         .CORR        Correlation matrix of posterior parameter samples             %%
%%         .loglik      BMA log-likelihood of maximum likelihood parameters           %%
%%         .std         Posterior standard deviations of BMA parameters               %%
%%         .pred        n x 2na matrix with lower/upper prediction quantiles of alpha %%
%%         .pdf_Y       n x 1 vector probability density BMA mixture at measurements  %%
%%         .cdf_Y       n x 1 vector cumulative distribution BMA mixture at ...       %%
%%         .coverage    Coverage (%) 100alpha prediction intervals BMA distribution   %%
%%         .spread      Spread of 100alpha prediction intervals BMA distribution      %%
%%         .Ye          n x 1 vector with the mean of BMA mixture forecast [= exact]  %% 
%%         .R2          coefficient of determination weighted-average BMA forecast    %%
%%         .RMSE        Root mean squared error of weighted-average BMA forecast      %%
%%         .R           Pearson coeff. measured and weighted-average BMA forecast     %%
%%         .ML          maximum likelihood values of BMA weights and shape parameters %%
%%         .CRPS        n x 1 vector with continuous ranked probability score ...     %%
%%         .QS          n x 1 vector with quadratic score BMA forecast distribution   %%
%%         .LS          n x 1 vector with logarithmic score BMA forecast distribution %%
%%         .SS          n x 1 vector with spherical score BMA forecast distribution   %%
%%         .mRLBL_anal  reliability of BMA mixture distribution                       %%
%%         .mCv_anal    time-averaged coefficient of variation BMA mixture [= exact]  %%
%%         .mQS         time-averaged quadratic score BMA forecast distribution       %%
%%         .mLS         time-averaged logarithmic score BMA forecast distribution     %%
%%         .mSS         time-averaged spherical score BMA forecast distribution       %%
%%         .mCRPS       time-averaged continuous ranked probability score ...         %%
%%         .mix_norm    n x 2 matrix with 1-norm and 2-norm of BMA mixture PDF        %%
%%         .mu_mix      n x 1 vector with the mean of BMA mixture forecast [= exact]  %% 
%%         .var_mix     n x 1 vector with variance of BMA mixture forecast [= exact]  %%
%%         .KGE         Kling-Gupta efficiency of weighted-average BMA forecast       %%
%%         .str_table   1 x d cell with names of the BMA parameters (Latex)           %%
%%         .RMSE_mod    Root mean squared errors of models of ensemble forecasts D    %%
%%                                                                                    %%
%% This general purpose MATLAB code is designed to postprocess forecast ensembles of  %%
%% dynamic simulation models. The user can select among 7 diferent model averaging    %%
%% methods including                                                                  %%
%%      1. Equal weights averaging                                                    %%
%%      2. Bates-Granger averaging                                                    %%
%%      3. Information criterion averaging with akaike (AIC) and bayes (BIC) IC       %%
%%      4. Granger-Ramanathan averaging                                               %%
%%      5. Bayesian model averaging                                                   %%
%%      6. Mallows model averaging (with and without weights restricted to Simplex)   %%
%%                                                                                    %%
%% For BMA the user can select among different options for the conditional forecast   %% 
%% distribution of the ensemble members. This includes                                %% 
%%          a. normal distribution                                                    %% 
%%          b. lognormal distribution                                                 %%     
%%          c. generalized normal distribution                                        %% 
%%          d. truncated normal distribution                                          %% 
%%          e. gamma distribution                                                     %%  
%%          f. Weibull distribution                                                   %% 
%%          g. generalized extreme value distribution                                 %% 
%%          h. generalized Pareto distribution                                        %%
%% For each conditional PDF the user can determine whether to implement a             %%
%% homoscedastic (= constant) or heteroscedastic (= nonconstant) variance with a      %%
%% single or group treatment                                                          %%
%%                                                                                    %%
%% This algorithm has been described in:                                              %%
%%   Vrugt, J.A. (2023), MODELAVG: A MATLAB toolbox for postprocessing of model       %%
%%       ensembles, Manual, Version 2.0, pp. 1 - XX, 2023.                            %%
%% For more information please read:                                                  %%
%%   Vrugt, J.A. (2023), Distribution-based Model Evaluation and Diagnostics:         %%
%%       Elicitability, Propriety and Scoring Rules for Hydrograph Functionals,       %%
%%       Water Resources Research, In Review.                                         %%
%%   Vrugt, J.A. (2015), Markov chain Monte Carlo simulation using the DREAM software %%
%%       package: Theory, concepts, and MATLAB implementation, Environmental Modeling %%
%%       and Software, XX, XX, doi:XX/XX.                                             %%
%%   Diks, C.G.H (2010), and J.A. Vrugt (2010), Comparison of point forecast accuracy %%
%%       of model averaging methods in hydrologic applications, Stochastic            %%
%%       Environmental Research and Risk Assessment, 24(6), 809-820,                  %%
%%       doi:10.1007/s00477-010-0378-z.                                               %%
%%   Vrugt, J.A., C.G.H. Diks, and M.P. Clark (2008), Ensemble Bayesian model         %%
%%       averaging using Markov chain Monte Carlo sampling, Environmental Fluid       %%
%%       Mechanics, 8(5-6), 579-595, doi:10.1007/s10652-008-9106-3.                   %%
%%   Vrugt, J.A., and B.A. Robinson (2007), Treatment of uncertainty using ensemble   %%
%%       methods: Comparison of sequential data assimilation and Bayesian model       %%
%%       averaging, Water Resources Research, 43, W01411, doi:10.1029/2005WR004838.   %%
%%                                                                                    %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%                                                                                    %%
%% Copyright (C) 2015  the author                                                     %%
%%                                                                                    %%
%%   This program is free software: you can modify it under the terms of the GNU      %%
%%   General Public License as published by the Free Software Foundation, either      %%
%%   version 3 of the License, or (at your option) any later version.                 %%
%%                                                                                    %%
%%   This program is distributed in the hope that it will be useful, but WITHOUT ANY  %%
%%   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS        %%
%%   FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.   %%
%%                                                                                    %%
%%   (c) Written by Jasper A. Vrugt, Feb 2012                                         %%
%%   University of California Irvine 			        	                          %%
%%                                                                                    %%
%%   Version 1: March 2015                                                            %%
%%   Version 2: March 2022                                                            %%
%%                                                                                    %%
%% ---------------------------------------------------------------------------------- %%

% Tao, I copy here only the call to BMA
options.VAR % any of '1' , '2' , '3' , '4'
options.PDF % any of 'normal','lognormal','gen_normal','tnormal','gamma','weibull','gev','gpareto'
% --> see manual
% Y is a n x 1 vector with training data
% D is a n x K matrix with predictions of K models (bias-corrected)

% Now check implementation of conditional PDF
switch options.VAR
    case '1'    % 1: common constant variance
        if strcmp(options.PDF,{'lognormal'})
            d = K + 1; par.max = [ ones(1,K) std(Y) ];
        else
            d = K + 1; par.max = [ ones(1,K) 10*std(Y) ];
        end
    case '2'    % 2: individual constant variance
        if strcmp(options.PDF,{'lognormal'})
            d = 2 * K; par.max = [ ones(1,K) std(Y)*ones(1,K) ];
        else
            d = 2 * K; par.max = [ ones(1,K) 10*std(Y)*ones(1,K) ];
        end
    case {'3'}  % 3: common non-constant variance
        if strcmp(options.PDF,{'weibull'})
            d = K + 1; par.max = [ ones(1,K) 25 ];
        elseif sum(strcmp(options.PDF,{'lognormal','gev','gpareto'}))
            d = K + 1; par.max = [ ones(1,K) 1 ];
        else
            d = K + 1; par.max = [ ones(1,K) 2 ];
        end
    case {'4'}  % 4: individual non-constant variance
        if strcmp(options.PDF,{'weibull'})
            d = 2 * K; par.max = [ ones(1,K) 25*ones(1,K) ];
        elseif sum(strcmp(options.PDF,{'lognormal','gev','gpareto'}))
            d = 2 * K; par.max = [ ones(1,K) ones(1,K) ];
        else
            d = 2 * K; par.max = [ ones(1,K) 2*ones(1,K) ];
        end
end
if strcmp(options.PDF,{'gen_normal'})
    switch options.TAU
        case {'1'}
            par.max = [ par.max 20 ]; d = d + 1;
        case {'2'}
            par.max = [ par.max 20*ones(1,K) ]; d = d + K;
    end
end
if sum(strcmp(options.PDF,{'gev','gpareto'}))
    switch options.TAU
        case {'1'} % xi cannot exceed 1/2 --> variance not defined
            par.max = [ par.max 0.5 ]; d = d + 1; l = 0;
        case {'2'}
            par.max = [ par.max ones(1,K)/2 ]; d = d + K; l = K - 1;
    end
end
% Weights must be on simplex - otherwise predictive density is meaningless
par.min = zeros(1,d);
% Overwrite values of zero of lower bound of xi
if strcmp(options.PDF,'gpareto')
    par.min(d-l:d) = -1.0;
end
% For GEV: zero implies that skew cannot become negative
% This avoids numerical instabilities and convergence problems for
% some of the shape parameters [= otherwise very difficult to
% solve]

% Problem settings defined by user
log_pdf = @(x) BMA_calc(x,D,Y,options); % Function handle for DREAM
% Run basic implementation of the DREAM_ZS algorithm
code_to_use = 'dream_zs';
switch code_to_use
    case 'dream_zs'
        T = d * 2500; %e4;                      % Number of generations
        N = 3;                                  % Number of chains
        [chain,output] = MODELAVG_dream_zs(par,log_pdf,N,T,d,K,1);
    case 'dream'
        T = d * 1e3;                            % Number of generations
        N = 20;                                 % Number of chains
        [chain,output] = MODELAVG_dream(par,log_pdf,N,T,d,K,1);
end