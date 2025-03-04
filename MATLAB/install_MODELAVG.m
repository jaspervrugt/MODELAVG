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
%%       pp. 273-316, 2016                                               %%
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
%%                                                                       %%
%% http://faculty.sites.uci.edu/jasper                                   %%
%% http://faculty.sites.uci.edu/jasper/publications/                     %%
%% https://scholar.google.com/citations?user=zkNXecUAAAAJ&hl=nl          %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Go to main DREAM directory 
addpath(pwd,[pwd '/postprocessing'],[pwd '/miscellaneous']);
% Now go to example directory; say example_1
cd example_1
% Now execute this example by typing in command prompt: "example_1" 