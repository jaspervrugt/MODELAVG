# ####################################################################### #
#                                                                         #
#  MM      MM  OOOOOO  DDDDDDD  EEEEEEE LL        AAAA   VV   VV GGGGGGG  #
#  MMM     MM OOOOOOOO DDDDDDDD EEEEEEE LL       AA  AA  VV   VV GG   GG  #
#  MMMM  MMMM OO    OO DD    DD EE      LL       AA  AA  VV   VV GG   GG  #
#  MM MMMM MM OO    OO DD    DD EEEE    LL       AA  AA  VV   VV GGGGGGG  #
#  MM  MM  MM OO    OO DD    DD EEEE    LL       AAAAAA   VV VV  GGGGGGG  #
#  MM      MM OO    OO DD    DD EE      LL      AA    AA   VVV        GG  #
#  MM      MM OOOOOOOO DDDDDDDD EEEEEEE LLLLLLL AA    AA   VVV       GGG  #
#  MM      MM  OOOOOO  DDDDDDD  EEEEEEE LLLLLLL AA    AA    V    GGGGGGG  #
#                                                                         #
# ####################################################################### #
#                                                                         #
# Multi-model averaging is widely used in various scientific and          #
# engineering disciplines to post-process forecast ensembles and/or       #
# quantify conceptual model uncertainty. This PYTHON toolbox, called      #
# MODELAVG, implements a suite of different model averaging techniques,   #
# including (amongst others) equal weights averaging (EWA),               #
# Bates-Granger model averaging (BGA), Bayesian model averaging (BMA),    #
# Mallows model averaging (MMA), and Granger-Ramanathan averaging (GRA)   #
# For BMA the user can select among different options for the             #
# conditional forecast distribution of each individual member of the      #
# ensemble. Options include a normal distribution with homoscedastic      #
# and heteroscedastic varaiance, a gamma distribution with constant or    #
# heteroscedastic variance, and a generalized normal distribution with    #
# homoscedastic or non-constant variance. The toolbox returns the         #
# optimal values posterior distribution of the weights and related        #
# parameters of each model averaging method, along with graphical         #
# output of the results. For MMA and BMA the user has access to the       #
# entire posterior distribution of the weights and/or variances derived   #
# from MCMC simulation with the DREAM algorithm. Three case studies       #
# involving forecast ensembles of hydrologic and meteorologic models      #
# are used to illustrate the capabilities and functionalities of the      #
# MODELAVG toolbox                                                        #
#                                                                         #
# ####################################################################### #
# Evaluates performance for evaluation period                             #
#                                                                         #
# SYNOPSIS: [val,D_eval] = MODELAVG_eval(method,D_eval,y_eval, ...        #
#               options,a,b,output)                                       #
#  where                                                                  #
#   method    [input] String with model averaging used                    #
#   D_eval    [input] mxK matrix forecasts ensemble members               #
#   y_eval    [input] mx1 vector with verifying observations              #
#   options   [input] structure BMA algorithmic variables                 #
#         .PDF        string: conditional PDF for BMA method (MANUAL)     #
#                      = 'normal'     = normal distribution               #
#                      = 'lognormal'  = lognormal distribution            #
#                      = 'tnormal'    = truncated normal ( > 0)           #
#                      = 'gen_normal' = generalized normal                #
#                      = 'gamma'      = gamma distribution                #
#                      = 'weibull'    = Weibull distribution [!]          #
#                      = 'gev'        = generalized extreme value         #
#                      = 'gpareto'    = generalized Pareto                #
#         .VAR        string: variance treatment BMA method (MANUAL)      #
#                      =  '1'          = constant group variance          #
#                      =  '2'          = constant individual variance     #
#                      =  '3'          = nonconstant group variance       #
#                      =  '4'          = nonconstant ind. variance        #
#         .TAU        string: trtmnt shape par gnrlzd nrml BMA method     #
#                      =  '1'          = common (= group) value of tau    #
#                      =  '2'          = individual value of tau          #
#         .alpha      1xna vector (1-sign. level) conf & pred limits      #
#                      = [0.5 0.7 0.9 0.95 0.99]                          #
#         .print      string: print results to screen or not              #
#                      = 'yes'                                            #
#                      = 'no'                                             #
#         .CPU        string: CPU-efficient solution or not               #
#                      = 'yes'                                            #
#                      = 'no'                                             #
#         .p          1 x K vector: with # parameters of models           #
#                      = [10 5 2 9 18 ...]                                #
#         .postproc   string: postprocess results BMA method or not       #
#                      = 'yes'                                            #
#                      = 'no' (= faster for ML parameters BMA method)     #
#   a         [input] 1xK vector of slopes linear bias correction         #
#   b         [input] 1xK vector of intercepts linear bias correction     #
#   output    [input] structure of results MODELAVG training period       #
#         .AR         qx2 matrix sample # & DREAM_ZS acceptance rate      #
#         .R_stat     qx(d+1) matrix sample # & R-hat conv. diag.         #
#         .MR_stat    qx2 matrix multivariate R-hat conv. diag.           #
#         .RunTime    Run time in seconds of DREAM_ZS algorithm           #
#         .CORR       Correlation matrix  posterior parameter samples     #
#         .loglik     BMA log-likelihood maximum likelihood parameters    #
#         .std        Posterior standard deviations of BMA parameters     #
#         .par_unc    nx2na matrix lower/upper confidence limits alpha    #
#         .pred       nx2na matrix lower/upper prediction limits alpha    #
#         .pdf_Y      nx1 vector with BMA mixture pdf at y                #
#         .cdf_Y      nx1 vector with BMA mixture cdf at y                #
#         .coverage   Coverage (#) 100alpha pred. intervals BMA model     #
#         .spread     Spread of 100alpha pred. intervals BMA model        #
#         .Ye         nx1 vector with µ BMA mixture forecast [exact]      # 
#         .R2         coefficient detrm. weighted-average BMA forecast    #
#         .RMSE       Root mean square err weighted-average BMA forcst    #
#         .R          Pearson coef. meas & weighted-average BMA forcst    #
#         .ML         maximum likelihood BMA weights & shape parametrs    #
#         .CRPS       nx1 vector continuous ranked probability score      #
#         .QS         nx1 vector quadratic score BMA forecst dstrbtion    #
#         .LS         nx1 vector logrithmc score BMA forecst dstrbtion    #
#         .SS         nx1 vector spherical score BMA forecst dstrbtion    #
#         .mRLBL_anal reliability of BMA mixture distribution             #
#         .mCv_anal   time-averaged coef. var. BMA mixture [= exact]      #
#         .mQS        time-averaged quadratic score BMA forecast dis.     #
#         .mLS        time-averaged logrithmic score BMA forecast dis.    #
#         .mSS        time-averaged spherical score BMA forecast dis.     #
#         .mCRPS      time-averaged continuous rankd probability score    #
#         .mix_norm   nx2 matrix 1-norm & 2-norm BMA mixture PDF          #
#         .mu_mix     nx1 vector mean BMA mixture forecast [= exact]      #
#         .var_mix    nx1 vector variance BMA mixture forcst [= exact]    #
#         .KGE        Kling-Gupta eff. weighted-average BMA forecast      #
#         .str_table  1 x d cell names of BMA parameters (Latex)          #
#         .RMSE_mod   Root mean square err. models ensemble frecasts D    #
#                                                                         #
# THIS TOOLBOX HAS BEEN DESCRIBED IN                                      #
#   Vrugt, J.A. (2023), MODELAVG: A MATLAB toolbox for postprocessing     #
#       of model ensembles, Manual, Version 2.0, pp. 1 - XX, 2023         #
# FOR MORE INFORMATION, PLEASE READ                                       #
#   Vrugt, J.A. (2024), Distribution-Based Model Evaluation and           #
#       Diagnostics: Elicitability, Propriety, and Scoring Rules for      #
#       Hydrograph Functionals, Water Resources Research, 60,             #
#       e2023WR036710, https://doi.org/10.1029/2023WR036710               #
#   Vrugt, J.A. (2015), Markov chain Monte Carlo simulation using the     # 
#       DREAM software package: Theory, concepts, and MATLAB              #
#       implementation, Environmental Modeling and Software, 75,          #
#       pp. 273-316                                                       #
#   Diks, C.G.H., and J.A. Vrugt (2010), Comparison of point forecast     #
#       accuracy of model averaging methods in hydrologic applications,   #
#       Stochastic Environmental Research and Risk Assessment, 24(6),     # 
#       809-820, https://doi.org/10.1007/s00477-010-0378-z                #
#   Vrugt, J.A., C.G.H. Diks, and M.P. Clark (2008), Ensemble Bayesian    #
#       model averaging using Markov chain Monte Carlo sampling,          #
#       Environmental Fluid Mechanics, 8(5-6), 579-595,                   #
#           https://doi.org/10.1007/s10652-008-9106-3                     #
#   Vrugt, J.A., C.J.F. ter Braak, M.P. Clark, J.M. Hyman and B.A.        #
#       Robinson (2008), Treatment of input uncertainty in hydrologic     #
#       modeling: Doing hydrology backward with Markov chain Monte        #
#       Carlo simulation, 44 (12), https://doi.org/10.1029/2007WR006720   #
#   Vrugt, J.A., and B.A. Robinson (2007), Treatment of uncertainty       #
#       using ensemble methods: Comparison of sequential data             #
#       assimilation and Bayesian model averaging, Water Resources        #
#       Research, 43, W01411, https://doi.org/10.1029/2005WR004838        #
#                                                                         #
# ####################################################################### #
#                                                                         #
# BUILT-IN CASE STUDIES                                                   #
#   Example 1   24-hour forecasts of river discharge                      #
#   Example 2   48-forecasts of sea surface temperature                   #
#   Example 3   48-forecasts of sea surface pressure                      #
#   Example 11  Forecasts of water levels                                 #
#   Example 12  Hydrologic modeling                                       #
#   Example 13  Flood modeling                                            #
#                                                                         #
# ####################################################################### #
#                                                                         #
# COPYRIGHT (c) 2015  the author                                          #
#                                                                         #
#   This program is free software: you can modify it under the terms      #
#   of the GNU General Public License as published by the Free Software   #
#   Foundation, either version 3 of the License, or (at your option)      #
#   any later version                                                     #
#                                                                         #
#   This program is distributed in the hope that it will be useful, but   #
#   WITHOUT ANY WARRANTY; without even the implied warranty of            #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU      #
#   General Public License for more details                               #
#                                                                         #
# ####################################################################### #
#                                                                         #
#  PYTHON CODE:                                                           #
#  © Written by Jasper A. Vrugt using GPT-4 OpenAI's language model       # 
#    University of California Irvine                                      #
#  Version 2.0    Dec 2024                                                #
#                                                                         #
# ####################################################################### #

# Check:  http://faculty.sites.uci.edu/jasper
# Papers: http://faculty.sites.uci.edu/jasper/publications/
# Google Scholar: https://scholar.google.com/citations?user=zkNXecUAAAAJ&hl=nl

import numpy as np
import time, os, sys
import scipy.special    # for the gamma function!

# Add the present directory to the Python path
module_path = os.getcwd()
if module_path not in sys.path:
    sys.path.append(module_path)

# Add related functions to Python path
parent_directory = os.path.join(module_path, 'miscellaneous')
sys.path.append(parent_directory)

from MODELAVG_functions import *

def MODELAVG(method, D, y, options=None):
    # If fewer than 3 arguments --> raise error
    if len(locals()) < 3:
        raise ValueError("MODELAVG: TooFewInputs", "Requires at least three input arguments.")
    
    # If options is not provided, initialize it as an empty dictionary
    if options is None:
        options = {}

    # Check if method is valid and make sure model is integer, JAV Python
    method = method.lower();

    # Now call MODELAVG_check to evaluate the user setup
    method, options = MODELAVG_check(method, D, y, options)
    
    # Now call MODELAVG_setup to define the computational framework
    beta, chain, n, K, sigma_2, log_L, par, count, options = MODELAVG_setup(method, D, y, options)
    
    # Start the timer
    start_time = time.time()

    # Use setdefault to initialize 'output' if it does not exist
    # output = {}  # options.setdefault('output', {})
    # output = options.setdefault('output', {})

    # Switch between different model averaging methods
    if method == 'ewa':  # equal weights averaging
        beta = np.array([np.ones(K) / K]).reshape(-1,1).flatten()
    elif method == 'bga':  # Bates-Granger averaging
        beta = (1 / sigma_2) / np.sum(1 / sigma_2)
    elif method == 'aica':  # information criterion averaging - AIC
        I = -2 * log_L + 2 * options['p']
        I = I - np.min(I)  # Adjust for underflow
        beta = np.exp(-I / 2) / np.sum(np.exp(-I / 2))
        print(beta,beta.shape) # OK
    elif method == 'bica':  # information criterion averaging - BIC
        I = -2 * log_L + options['p'] * np.log(n)
        I = I - np.min(I)  # Adjust for underflow
        beta = np.exp(-I / 2) / np.sum(np.exp(-I / 2))
    elif method == 'gra':  # Granger-Ramanathan averaging
        beta = np.linalg.pinv(D.T @ D) @ D.T @ y
        beta = beta.T
    elif method == 'bma':  # Bayesian model averaging
        par = {}; std_y = np.std(y)
        # Check implementation of conditional PDF
        if options['VAR'] == '1':  # Common constant variance
            d = K + 1
            if options['PDF'] == 'lognormal':
                par['max'] = np.concatenate([np.ones(K), std_y * np.ones(1)])
            else:
                par['max'] = np.concatenate([np.ones(K), 10 * std_y * np.ones(1)])
        elif options['VAR'] == '2':  # Individual constant variance
            d = 2 * K
            if options['PDF'] == 'lognormal':
                par['max'] = np.concatenate([np.ones(K), std_y * np.ones(K)])
            else:
                par['max'] = np.concatenate([np.ones(K), 10 * std_y * np.ones(K)])
        elif options['VAR'] == '3':  # Common non-constant variance
            d = K + 1
            if options['PDF'] == 'weibull':
                par['max'] = np.concatenate([np.ones(K), 25 * np.ones(1)])
            elif options['PDF'] in ['lognormal', 'gev', 'gpareto']:
                par['max'] = np.concatenate([np.ones(K), np.ones(1)])
            else:
                par['max'] = np.concatenate([np.ones(K), 2 * np.ones(1)])
        elif options['VAR'] == '4':  # Individual non-constant variance
            d = 2 * K
            if options['PDF'] == 'weibull':
                par['max'] = np.concatenate([np.ones(K), 25 * np.ones(K)])
            elif options['PDF'] in ['lognormal', 'gev', 'gpareto']:
                par['max'] = np.ones(2*K)
            else:
                par['max'] = np.concatenate([np.ones(K), 2 * np.ones(K)])
        if options['PDF'] == 'gen_normal':
            if options['TAU'] == '1':
                d += 1
                par['max'] = np.concatenate([par['max'], 20 * np.ones(1)])
            elif options['TAU'] == '2':
                d += K
                par['max'] = np.concatenate([par['max'], 20 * np.ones(K)])
        if options['PDF'] in ['gev', 'gpareto']:
            if options['TAU'] == '1':
                d += 1; l = 0
                par['max'] = np.concatenate([par['max'], 0.5 * np.ones(1)])
            elif options['TAU'] == '2':
                d += K; l = K - 1
                par['max'] = np.concatenate([par['max'], 0.5 * np.ones(K)])
        
        # Ensure weights are on simplex
        par['min'] = np.zeros(d)
        if options['PDF'] == 'gpareto':
            par['min'][-l:] = -1.0
        Y = np.tile(y, (K, 1)).T  # Copy the verifying data K times
        # Define the log-likelihood function
        log_pdf = lambda x: BMA_lik(x, D, Y, options)  # Define the likelihood function for DREAM
        # Call the DREAM_ZS function (you need to implement this or use an existing library)
        T = int(d * 2500)
        N = int(3)
        chain, output = MODELAVG_dream_zs(par, log_pdf, N, T, d, K, 1)
    elif method in ['mma', 'mma-s']:  # Mallows model averaging
        par = {}
        d = K
        T = int(d * 10000)
        N = int(3)
#        sigma_2 = np.sum((D - y) ** 2, axis=0) / n
        sigma_2 = np.sum((D - y[:,np.newaxis]) ** 2, axis=0) / n
        idx = np.argmax(options['p'])
        var_err = sigma_2[idx]
        
        log_pdf = lambda x: MMA_lik(x, D, y, var_err, options['p'])
        
        if method == 'mma':
            par['min'] = -np.ones(K)
            par['max'] = np.ones(K)
            chain, output = MODELAVG_dream_zs(par, log_pdf, N, T, d, K, 0)

        elif method == 'mma-s':
            par['min'] = np.zeros(K)
            par['max'] = np.ones(K)
            chain, output = MODELAVG_dream_zs(par, log_pdf, N, T, d, K, 1)

    if 'RunTime' not in options:
        # Calculate wall time
        output['RunTime'] = time.time() - start_time
        # Store RMSE of the forecast error for each model
        output['RMSE_mod'] = np.sqrt(sigma_2)

    # Now create the output that is printed in postprocessor
    beta, output, str_plot = MODELAVG_end(method, D, y, beta, chain, options, output)
    # Print progress
    if options['postproc'] == 'yes':
        for t in range(2, 4):
            if t == 2:
                print_name = '........'
            if t == 3:
                print_name = '........ done'
            if t > 1:
                print(f"\r\nMODELAVG postprocessor, please wait {print_name}", end="")
                if t == 2:
                    MODELAVG_postproc(method, D, y, beta, sigma_2, par, chain, options, output, str_plot)

    # Return output
    return beta, output
