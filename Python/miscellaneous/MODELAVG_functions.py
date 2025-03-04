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
import scipy.stats as stats
from scipy.stats import norm, lognorm, truncnorm, pearsonr, gamma, weibull_min, genextreme, genpareto, t 
from scipy.special import erf, lambertw, gammaln, psi, gammainc
from scipy.optimize import fsolve
from scipy.integrate import quad
import scipy.special    # for the gamma function!
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from statsmodels.graphics.tsaplots import plot_acf
from math import pi
import array, os
from datetime import datetime
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages # type: ignore
from screeninfo import get_monitors # type: ignore
from scipy import integrate
import warnings

def MODELAVG_eval(method, D_eval, y_eval, options, a, b, output):
    """
    Evaluates performance for evaluation period for the BMA model or other methods.
    
    Parameters:
        method (str): The model averaging method, e.g., 'bma'.
        D_eval (ndarray): The ensemble forecasts for the evaluation period.
        y_eval (ndarray): The verifying observations for the evaluation period.
        options (dict): A dictionary with options for model averaging.
        a (ndarray): Intercept values for bias correction.
        b (ndarray): Slope values for bias correction.
        output (dict): The output from the model averaging procedure containing model parameters.
    
    Returns:
        val (dict): A dictionary containing the evaluation results.
        D_eval (ndarray): The bias-corrected ensemble forecasts.
    """
    # Bias correction evaluation ensemble using intercept "a" and slope "b"
    m, K = D_eval.shape
    for k in range(K):
        D_eval[:, k] = a[k] + b[k] * D_eval[:, k]
    
    print(f"----- Evaluate {method.upper()} for the evaluation period -----\n")
    
    # Now check setup -> make sure that structure options is correct
    method, options = MODELAVG_check(method, D_eval, y_eval, options)
    # Now run setup function of toolbox (This is assumed to be a separate function)
    _, _, _, _, _, _, _, _, options = MODELAVG_setup(method, D_eval, y_eval, options, 1)
    # Extract maximum likelihood values weights, etc.
    beta_opt = output['ML']; val = {}
    
    if method == 'bma':  # Derive quantiles and scoring rules
        # Now derive the prediction uncertainty ranges of BMA model
        val['pred'], val['pdf_Y'], val['cdf_Y'] = BMA_quantile(beta_opt, D_eval, y_eval, options)
        # Set loglikelihood equal to the training period (assuming it's provided in output)
        val['loglik'] = np.sum(np.log(val['pdf_Y']))
        
        # number of alpha values
        na = len(options['alpha'])
	
        # And now the % observations enclosed - for each alpha value
        val['coverage'] = np.zeros(na)
        val['spread'] = np.zeros(na)
        for i in range(na):
            # Compute coverage percentage
            val['coverage'][i] = 100 * np.sum((val['pred'][:, i] < y_eval) & (y_eval < val['pred'][:, 2 * na - (i + 1)])) / m
            # Compute the spread of the prediction interval
            val['spread'][i] = np.mean(val['pred'][:, 2 * na - (i + 1)] - val['pred'][:, i])
        
        # Compute scoring rules of BMA distribution forecasts
        val = BMA_score(beta_opt, D_eval, y_eval, options, val)
    
    # Calculate point forecast, RMSE, and Pearson correlation coefficient (R)
    val['Ye'] = np.dot(D_eval, beta_opt[:K])
    # Compute the mean of the data
    Ym_eval = np.mean(y_eval)
    # Compute the total sum of squares (SStot_eval)
    SStot_eval = np.sum((y_eval - Ym_eval) ** 2)
    # Compute the residuals
    e_eval = y_eval - val['Ye']
    # Compute the residual sum of squares (SSres_eval)
    SSres_eval = np.sum(e_eval ** 2)
    # Compute coefficient of determination (R2)
    val['R2'] = 1 - SSres_eval / SStot_eval
    # Compute RMSE (Root Mean Square Error)
    val['RMSE'] = np.sqrt(np.sum(e_eval ** 2) / m)  # sqrt(sigma2)?
    # Compute Pearson correlation coefficient (R)
    val['R'] = np.corrcoef(val['Ye'], y_eval)[0, 1]
    # Compute KGE (Kling-Gupta Efficiency) of weighted average BMA forecast
    val['KGE'] = compKGE(y_eval, val['Ye'])
    
    return val, D_eval


def MODELAVG_check(method, D, y, options):
    # suppress warnings
    warnings.filterwarnings("ignore")

    # Derive current time and set deadline
    deadline = datetime(2025, 2, 28)
    # Check if trial version has ended
    if (deadline - datetime.now()).days < 0:
        raise ValueError('MODELAVG ERROR: Trial version of MODELAVG V2.0 has ended')
    
    # Open an output file with warnings
    with open('warning_file.txt', 'w+') as fid:
        fid.write('-------------- MODELAVG warning file --------------\n')

    with open('warning_file.txt', 'a') as file:
        # Check method
        if not method:
            raise ValueError('MODELAVG ERROR: Input argument method should not be empty but a string enclosed between quotes with model averaging method to use')
        if not isinstance(method, str):
            raise ValueError('MODELAVG ERROR: Input argument method should be a string enclosed between quotes')
        
        # Check D
        if D is None or len(D) == 0:
            raise ValueError('MODELAVG ERROR: Input argument D should not be empty but a n x K matrix with the n forecasts of K models')
        if not isinstance(D, np.ndarray):
            raise ValueError('MODELAVG ERROR: Input argument D should be a n x K matrix ( = numerical values )')
        
        # Check y
        if y is None or len(y) == 0:
            raise ValueError('MODELAVG ERROR: Input argument y should not be empty but a n x 1 vector with the verifying observations of the training data set')
        if not isinstance(y, np.ndarray):
            raise ValueError('MODELAVG ERROR: Input argument y should be a n x 1-vector ( = numerical values )')

        # Check options
        if not isinstance(options, dict):
            raise ValueError('MODELAVG ERROR: Input argument options should be a structure with fields')

        n, K = D.shape

        if n == 1:
            raise ValueError('MODELAVG ERROR: No point to use model averaging if matrix D stores only a single prediction of each of the K models, that is n = 1')

        if K > n:
            raise ValueError('MODELAVG ERROR: Ill-determined problem as you have more models in the ensemble than their corresponding predictions/training observations; K > n')

        if K == 1:
            raise ValueError('MODELAVG ERROR: No point to use model averaging if only a single model is used ( matrix D only has one column )')
        
        if y.shape[0] != n:
            raise ValueError('MODELAVG ERROR: Dimension mismatch. The number of rows, n of y should match number of rows of D')
        
        # Ensure y is a column vector
        if y.shape[0] > 1:
            print(f'MODELAVG WARNING: Vector y with verifying observations of training data set should be column vector; transposed in code')
            file.write(f'MODELAVG WARNING: Vector y with verifying observations of training data set should be column vector; transposed in code\n')
            y = y.flatten()  # Ensure y is a column vector
        
        # if m != n:
        #   raise ValueError('MODELAVG ERROR: Number of rows, n, of n x K matrix D and number of elements ( = training data points ) of vector y do not match')
        
        # Convert method to lowercase
        method = method.lower()

        valid_methods = ['ewa', 'bga', 'aica', 'bica', 'gra', 'bma', 'mma', 'mma-s']
        if method not in valid_methods:
            raise ValueError(f'MODELAVG ERROR: Unknown model averaging method: method should be one of {", ".join(valid_methods)}')

        # List of specific field names to capitalize
        fields_to_capitalize = ['pdf', 'var', 'tau', 'cpu']
        for field in fields_to_capitalize:
            for key in list(options.keys()):  # Use list() to avoid modifying dict during iteration
                if key.lower() == field.lower():
                    options[key.upper()] = options.pop(key)

        # List of specific field names to use lower case letters
        fields_to_not_capitalize = ['print','p','postproc','alpha']
        for field in fields_to_not_capitalize:
            for key in list(options.keys()):  # Use list() to avoid modifying dict during iteration
                if key.lower() == field.lower():
                    options[key.lower()] = options.pop(key)
        
        # List of specific field values that should be lower case
        fields_to_modify = ['PDF','CPU','postproc','print']
        for field in fields_to_modify:
            for key in list(options.keys()):  # Use list() to avoid modifying dict during iteration
                if key.lower() == field.lower():
                    # Check if the value is a string, and if so, convert it to lowercase
                    if isinstance(options[key], str):
                        options[key] = options[key].lower()
                        # If the value is a list, apply lowercase to string elements (if any)
                    elif isinstance(options[key], list):
                        options[key] = [str(item).lower() if isinstance(item, str) else item for item in options[key]]

        # Check 'p' field in options
        if method in ['aica', 'bica', 'mma', 'mma-s']:
            if 'p' not in options:
                raise ValueError(f'MODELAVG ERROR: Unknown complexity of each model for {method} method --> Define field "p" of structure options')
            elif isinstance(options['p'], array.array):
                options['p'] = np.array(options['p']).reshape(1, -1)
            else:
                p = np.array(options['p'])
                # p = options.get('p').flatten()
            if p is None or len(p) == 0:
                raise ValueError('MODELAVG ERROR: Field "p" of structure options should not be empty but a (K x 1)-vector with number of "free" parameters of each model')
            if not isinstance(p, np.ndarray):
                raise ValueError('MODELAVG ERROR: Field "p" of structure options should be a (K x 1)-vector with number of "free" parameters of each model')

            if p.shape[0] != K:
                raise ValueError('MODELAVG ERROR: Number of elements of field "p" of structure options does not match K, the number of models of the ensemble in matrix D')
            
            if method in ['aica', 'bica']:
                options['p'] = p.flatten()
            elif method in ['mma', 'mma-s']:
                options['p'] = p

        elif method in ['ewa', 'bga', 'gra', 'bma']:
            if 'p' in options:
                file.write(f'MODELAVG WARNING: Field "p" of structure options defined but not used by {method} method\n')
                del options['p']
                # options['p'] = np.array(options['p']).reshape(1, -1)

        # Check 'print' field in options
        if 'print' not in options:
            file.write(f'MODELAVG WARNING: Field "print" of structure options not defined as "yes" or "no" (default of options.print = "yes" used)\n')
        
        if 'print' in options:
            print_value = options['print']
            if print_value not in ['yes', 'no']:
                raise ValueError('MODELAVG ERROR: Field "print" of structure options should equal "yes" or "no"')

        if 'alpha' not in options:
            evalstr = "MODELAVG WARNING: Field 'alpha' of structure options not defined -> default of alpha = 0.05 is used"
            print(evalstr)
            file.write(evalstr + '\n')
            options['alpha'] = np.array([0.05])
        else:

            if options['alpha'] is None:
                raise ValueError("MODELAVG ERROR: Field 'alpha' of structure "
                                 "options should not be empty but should be "
                                 "a numeric value (options: (0-1))")
            elif isinstance(options['alpha'], array.array):
                options['alpha'] = np.array(options['alpha']).reshape(1, -1)

            if isinstance(options['alpha'], list):
                print("Field alpha of options is provided as a list, please make this a numpy array.")

            options['alpha'] = np.array(options['alpha'])
            if np.any(np.array(options['alpha']) < 0):
                raise ValueError("MODELAVG ERROR: Field 'alpha' of structure "
                             "options cannot have negative values (e.g. use 0.10 or 0.05 or 0.01)")
            if np.any(options['alpha'] > 1):
                raise ValueError("MODELAVG ERROR: Field 'alpha' of structure "
                                 "options cannot be larger than one (e.g. use 0.10 or 0.05 or 0.01)")
            if np.any(options['alpha'] > 3/4):
                evalstr = "MODELAVG WARNING: Field 'alpha' of structure options set rather unusual (e.g. use 0.10 or 0.05 or 0.01)"
                print(evalstr)
                file.write(evalstr + '\n')
    
        # Now check whether BMA method is used
        if method == 'bma':
            # Which percentage of ensemble forecasts is negative
            pct = 100 * np.sum(D < 0) / D.size
            # Which percentage of training data (y) is negative
            pct_2 = 100 * np.sum(y < 0) / y.size
    
            # Warning if BMA conditional forecast PDF is not defined
            if 'PDF' not in options:
                raise ValueError("MODELAVG ERROR: BMA method is used but conditional "
                                "distribution is not defined in field 'PDF' of structure "
                                 "options")

            if 'PDF' in options:
                if options['PDF'] is None:
                    raise ValueError("MODELAVG ERROR: Field 'PDF' of structure options "
                                     "should not be empty but should be a string (options: "
                                     "'normal' or 'gamma' or 'gen_normal')")
                if not isinstance(options['PDF'], str):
                   raise ValueError("MODELAVG ERROR: Field 'PDF' of structure options "
                                     "should be a string (options: 'normal' or 'gamma' "
                                     "or 'gen_normal')")
                if options['PDF'] not in ['normal', 'gamma', 'gen_normal', 'lognormal', 
                                          'weibull', 'power_law', 'tnormal', 'gev', 'gpareto']:
                    raise ValueError("MODELAVG ERROR: Unknown conditional PDF of BMA method "
                                     "(use 'normal' or 'lognormal' or 'gen_normal' "
                                     "or 'gamma' or 'weibull' or 'gev' or 'gpareto')")
        
            if options['PDF'] == 'gamma':
                # Provide warning
                if pct > 0:
                    evalstr = f"MODELAVG WARNING: Gamma distribution is used but {round(100 * pct) / 100}% of the forecasts of the ensemble are negative"
                    print(evalstr)
                    file.write(evalstr + '\n')
                # Provide warning
                if pct_2 > 0:
                    evalstr = f"MODELAVG WARNING: Gamma distribution is used but {round(100 * pct_2) / 100}% of the verifying observations are negative"
                    print(evalstr)
                    file.write(evalstr + '\n')

            if options['PDF'] in ['gen_normal', 'gev']:
                if 'TAU' not in options:
                    raise ValueError("MODELAVG ERROR: Field 'TAU' of structure "
                                     "options should be specified for generalized "
                                     "normal conditional PDF (options: '1' or '2')")
                if options['TAU'] is None:
                    raise ValueError("MODELAVG ERROR: Field 'TAU' of structure "
                                     "options should not be empty but contain a "
                                     "string (options: '1' or '2')")
                if not isinstance(options['TAU'], str):
                    raise ValueError("MODELAVG ERROR: Field 'TAU' of structure "
                                     "options should be a string (options: '1' or '2')")
                if options['TAU'] not in ['1', '2']:
                    raise ValueError("MODELAVG ERROR: Unknown value of field "
                                     "'TAU' of structure options (use '1' or '2')")
    
            if 'VAR' not in options:
                evalstr = "MODELAVG WARNING: BMA method is used but field 'VAR' of structure options not defined (default of options.VAR = '1' is used)"
                print(evalstr)
                file.write(evalstr + '\n')

            if 'VAR' in options:
                if options['VAR'] is None:
                    raise ValueError("MODELAVG ERROR: Field 'VAR' of structure options "
                                     "should not be empty but contain a string (options: "
                                     "'1' or '2' or '3' or '4')")
            if not isinstance(options['VAR'], str):
                    raise ValueError("MODELAVG ERROR: Field 'VAR' of structure options "
                                     "should be a string (options: '1' or '2' or '3' "
                                     "or '4')")
            if options['VAR'] not in ['1', '2', '3', '4']:
                    raise ValueError("MODELAVG ERROR: Unknown value of field 'VAR' "
                                     "of structure options (use '1' or '2' or '3' "
                                     "or '4')")

            # Now check sigma ahead of time
            if options['VAR'] in ['3', '4']:
                if pct > 0:
                    evalstr = f"MODELAVG WARNING: Variance option {options['VAR']} used but {round(100 * pct) / 100}% of forecasts of the ensemble are negative"
                    print(evalstr)
                    file.write(evalstr + '\n')

                if options['VAR'] == '3':
                    evalstr = "MODELAVG WARNING: As 'sigma(:,k) = c*D(:,k)' (c > 0) then sigma(:,k) can be < 0 --> code uses 'sigma(:,k) = c*abs(D(:,k))'; k in [1...K]"
                    print(evalstr)
                    file.write(evalstr + '\n')

                if options['VAR'] == '4':
                    evalstr = "MODELAVG WARNING: As 'sigma(:,k) = c(k)*D(:,k)' (c(k) > 0) then sigma(:,k) can be < 0 --> code uses 'sigma(:,k) = c(k)*abs(D(:,k))'; k in [1...K]"
                    print(evalstr)
                    file.write(evalstr + '\n')

            # If method is not 'bma'
        else:
            # Remove PDF field
            if 'PDF' in options:
                del options['PDF']
            # Remove VAR field
            if 'VAR' in options:
                del options['VAR']
            # Remove TAU field
            if 'TAU' in options:
                del options['TAU']
            # Remove alpha field (commented out in original MATLAB code)
            # if 'alpha' in options:
            #     del options['alpha']


    return method, options


def MODELAVG_setup(method, D, Y, options, flag=0):
    """
    Defines computational framework for the MODELAVG method.
    
    Parameters:
    - method: str, method type, e.g., 'bma'
    - D: np.array, ensemble forecast matrix (n x K)
    - Y: np.array, observation values (n,)
    - options: dict, options for the setup
    - flag: int, flag to control print behavior, default is 0

    Returns:
    - beta: empty array, placeholder
    - chain: empty array, placeholder
    - n: number of observations
    - K: number of ensemble members
    - sigma_2: error variance of each model
    - log_L: log-likelihood of each model
    - par: empty array, placeholder
    - count: 0, placeholder
    - options: updated options dictionary
    """
    # Assign default settings if not specified
    if method == 'bma':
        if 'VAR' not in options:
            options['VAR'] = '1'  # Set variance treatment to '1'
        if 'alpha' not in options:
            options['alpha'] = np.array([0.05])

    if 'print' not in options:
        options['print'] = 'yes'
    if 'postproc' not in options:
        options['postproc'] = 'yes'
    if 'CPU' not in options:
        options['CPU'] = 'no'

    # Initialize beta and chain to be empty
    beta = []
    chain = []

    # Determine the size of the ensemble forecast
    n, K = D.shape

    # Calculate error variance of each model
    sigma_2 = np.sum((D - Y[:, np.newaxis]) ** 2, axis=0) / n

    # Calculate log-likelihood of each model
    log_L = -0.5 * (n * sigma_2 + n)

    # Scale to avoid numerical underflow
    log_L -= np.max(log_L)

    # Initialize count and structure par
    par = []
    count = 0

    # Print header information
    print('  -------------------------------------------------------------------------------------------            ')
    print('  MMM        MMM    OOOOOO    DDDDDDDD   EEEEEEEEE LLL           AAA     VVV    VVV GGGGGGGGG            ')
    print('  MMMM      MMMM   OOOOOOOO   DDDDDDDDD  EEEEEEEEE LLL          AAAAA    VVV    VVV GGGGGGGGG            ')
    print('  MMMMM    MMMMM  OOO    OOO  DDD    DDD EEE       LLL         AAA AAA   VVV    VVV GGG   GGG            ')
    print('  MMMMMM  MMMMMM OOO      OOO DDD    DDD EEE       LLL        AAA   AAA  VVV    VVV GGG   GGG            ')
    print('  MMM MMMMMM MMM OOO      OOO DDD    DDD EEEEE     LLL       AAA     AAA VVV    VVV GGG   GGG     /^ ^\  ')
    print('  MMM  MMMM  MMM OOO      OOO DDD    DDD EEEEE     LLL       AAAAAAAAAAA VVV    VVV  GGGGGGGG    / 0 0 \ ')
    print('  MMM   MM   MMM OOO      OOO DDD    DDD EEE       LLL       AAA     AAA  VVV  VVV        GGG    V\ Y /V ')
    print('  MMM        MMM  OOO    OOO  DDD    DDD EEE       LLL       AAA     AAA   VVVVVV         GGG     / - \  ')
    print('  MMM        MMM   OOOOOOOO   DDDDDDDDD  EEEEEEEEE LLLLLLLLL AAA     AAA    VVVV     GGGGGGGG    /     | ')
    print('  MMM        MMM    OOOOOO    DDDDDDDD   EEEEEEEEE LLLLLLLLL AAA     AAA     VV     GGGGGGGGG    V__) || ')
    print('  -------------------------------------------------------------------------------------------            ')
    print('  © Jasper A. Vrugt, University of California Irvine & GPT-4 OpenAI''s language model')
    print('    ________________________________________________________________________')
    print('    Version 2.0, Dec. 2024, Beta-release: MATLAB implementation is benchmark')
    print('\n')

    # Print summary of the settings if flag is 0
    if flag == 0:
        print(f"\nModel averaging method used: {method.upper()}")
        print("\n---------- Summary of the main settings used: options -----------------")
        for key, value in options.items():
            if isinstance(value, np.ndarray):                                           # Check if value is a NumPy array
                flattened_values = value.flatten()                                      # Flatten to 1D array
                formatted_values = ', '.join([f"{v:.4f}" for v in flattened_values])    # Format each element
                print(f"{key:<8}: {formatted_values}")
            else:
                print(f"{key:<8}: {value:<6}")
        print("\n-----------------------------------------------------------------------")

    return beta, chain, n, K, sigma_2, log_L, par, count, options


def MODELAVG_dream_zs(par, log_pdf, N, T, d, K, bound_check):
    # Default parameters shared across DREAM algorithm
    Delta, c, c_star, n_CR, p_ug = 3, 5e-2, 1e-12, 3, 0.2
    k, p_s, m0 = 10, 0.1, max(N, 20 * d)
    m = m0
    chain = np.full((T, d + 1, N), np.nan)              # Preallocate memory for Markov chains
    gamma_RWM = lambda a, b: 2.38 / np.sqrt(2 * a * b)  # Function handle with RWM-based jump rate

    count, ct, accept, t_old = 0, 1, 0, 1
    steps = T // 50
    MR_stat, AR = np.full((T // steps, 2), np.nan), np.full((T // steps, 2), np.nan)
    R_stat = np.full((T // steps, d + 1), np.nan)
    MR_stat[0,0], R_stat[0,0], AR[0,0] = 1, 1, 1 

    # fixed by JAV
    parmin = np.tile(par['min'], (N, 1))
    parmax = np.tile(par['max'], (N, 1))

    # Initialize DREAM algorithm
    p_CR = np.ones(n_CR) / n_CR     # Crossover: initial n_CR selection probabilities
    J = np.ones(n_CR)               # Crossover: jump distance and frequency of use
    n_id = np.ones(n_CR)            # Crossover: frequency of use
    # Initialize Z matrix
    Z = LH_sampling(parmin[0, :], parmax[0, :], m0 + N)  # Initial population creation
    if bound_check:
        Z[:m0 + N, :K] = Z[:m0 + N, :K] / Z[:m0 + N, :K].sum(axis=1, keepdims=True)  # Weights sum to one
    
    # Compute log-likelihood of initial chain states
    ell = np.full((m0 + N,1),np.nan)
    for i in range(m0, m0 + N):
        ell[i,0] = log_pdf(Z[i, :d])

    # combine Z and ell
    Z = np.hstack([Z, ell])
    X = Z[m0:m0 + N, :d + 1]            # Extract initial states from archive
    std_Z = np.std(Z[:m0, :d], axis=0)  # Standard deviation of each dimension in archive
    chain[0, :d + 1, :] = X.T           # Store in chain the N initial states
 
    flag = 0
    Xp = np.full((N, d+1), np.nan)      # Initialize flag and candidate population

    # DYNAMIC PART (EVOLVE N CHAINS T DIFFERENT STEPS)
    for t in range(1, T):
        method = np.random.choice([1, 2], p=[1 - p_s, p_s])         # Jump method in current generation
        eps = c_star * np.random.randn(N, d)                        # Ergodicity perturbation
        lambda_ = 1 + c * (2 * np.random.rand(N, d) - 1)            # Draw lambda ~ U[-c,c] + add 1
        delta = np.random.randint(1, Delta + 1)                     # Draw delta from [1,...,Delta]
        R = np.random.choice(m, size=(3 * delta, N), replace=True)  # Sample 3*delta*N integers
        dX = np.zeros((N, d))                                       # Set N jump vectors to zero
        log_alfa_J = np.zeros(N)                                    # snooker correction
        if method == 1:  # Parallel Direction Move
            # id = np.random.choice(n_CR, N, p=p_CR)                # Select N crossover indexes
            id = 1 + np.random.choice(n_CR, N, p=p_CR)              # Select N crossover indexes, python fix to match matlab
            gamma = np.random.choice([0, 1], N, p=[1 - p_ug, p_ug]) # Gamma equal to 0 or 1
            U = np.random.rand(N, d)                                # Random labels between 0 and 1

            for i in range(N):                                      # Create proposal in each chain
                a = R[:delta, i]
                b = R[delta:2 * delta, i]
                if gamma[i] == 0:                                   # Subspace sampling
                    A = np.where(U[i, :] < (id[i] / n_CR))[0]
                    d_star = len(A)
                    if d_star == 0:
                        A = np.argmin(U[i, :])
                        d_star = 1
                    gamma_d = 1 / 2 * gamma_RWM(delta, d_star)      # Compute default jump rate
                    A = np.array(A).flatten()                       # Ensure 'A' is 1D
                else:                                               # All dimensions update
                    A = np.arange(d)
                    gamma_d = 1
                    a = a[0]
                    b = b[0]
                    id[i] = n_CR               

                    a = np.array(a).flatten()  # Ensure 'a' is 1D
                    b = np.array(b).flatten()  # Ensure 'b' is 1D
                    A = np.array(A).flatten()  # Ensure 'A' is 1D

                dX[i, A] = lambda_[i, A] * gamma_d * np.sum(Z[np.ix_(a, A)] - Z[np.ix_(b, A)],axis=0) + eps[i, A]  # Compute jump vector

            Xp[:, :d] = X[:, :d] + dX                               # Compute candidate points
            
        elif method == 2:  # Snooker Move
            a = R[0, :]
            b = R[1, :]
            c_sn = R[2, :]
            id = n_CR * np.ones(N)        
            
            for i in range(N):  # Create proposal in each chain
                F = X[i, :d] - Z[c_sn[i], :d]
                FF = np.maximum(F @ F.T, np.finfo(float).eps)
                # Extra lines
                # zP = F * (np.sum((Z[a[i], :d] - Z[b[i], :d]) * F, axis=1) / FF)
                # Correction in python
                zP = F * (np.sum((Z[a[i], :d] - Z[b[i], :d]) * F) / FF)
                gamma_s = 1.2 + np.random.rand()
                dX[i, :] = lambda_[i, :] * gamma_s * zP + eps[i, :]
                # Vectorized implementation according to eDREAMPackage  
                # Xp[i, :] = X[i, :d] + dX[i, :]
                # XpZ = np.linalg.norm(Xp[i, :d] - Z[c_sn[i], :d], 2) + np.finfo(float).eps
                # XZ = np.linalg.norm(X[i, :d] - Z[c_sn[i], :d], 2) + np.finfo(float).eps
                # log_alfa_J[i] = (d - 1) * np.log(XpZ / XZ)

            Xp[:, :d] = X[:, :d] + dX  # Compute candidate point
            log_alfa_J = (d - 1) / 2 * (np.log(np.sum((Xp[0:N, 0:d] - Z[c_sn, 0:d]) ** 2, axis=1)) - np.log(np.sum((X[0:N, 0:d] - Z[c_sn, 0:d]) ** 2, axis=1) + np.finfo(float).eps) )

        if bound_check:
            Xp[0:N,0:d] = bnd_check(Xp[0:N,0:d], X[0:N,0:d], K, parmin, parmax)  # Boundary constraints

        # Compute log-likelihood for each proposal
        for i in range(N):
            Xp[i,d] = log_pdf(Xp[i, :d])

        # alfa_L = np.exp(Xp[:, d] - X[:, d])                # Likelihood ratio of proposals
        # p_acc = np.exp(log_alfa_J) * alfa_L                # Acceptance probability of proposals
        # p_acc = np.exp(log_alfa_J + Xp[:, d] - X[:, d])    # Acceptance probability of proposals
        # diff = log_alfa_J + Xp[:, d] - X[:, d]
        min_exp_value = -700; max_exp_value = 700   # Reasonable lower and upper bounds for float64
        diff = np.clip(log_alfa_J + (Xp[:, d] - X[:, d]), min_exp_value, max_exp_value)
        p_acc = np.exp(diff)                        # Acceptance probability of proposals
        u = np.random.rand(N)                       # Draw N labels from U[0, 1]

        for i in range(N):
            if p_acc[i] > u[i]:
                X[i, :d + 1] = Xp[i, :d + 1]
                accept += 1

        chain[t, :d + 1, :] = X.T  # Add current states to chain

        # Update crossover probabilities after burn-in
        if t < T / 10:
            for i in range(N):
                # J[int(id[i])] += np.sum((dX[i, 0:d] / std_Z) ** 2)
                J[int(id[i])-1] += np.sum((dX[i, 0:d] / std_Z) ** 2)
                # n_id[int(id[i])] += 1
                n_id[int(id[i])-1] += 1
        elif t >= T // 10 and flag == 0:
            p_CR = J / n_id
            p_CR /= np.sum(p_CR)
            flag = 1

        # Append X to archive Z
        if t % k == 0:
            # Z[m:m + N, :d + 1] = X
            Z = np.vstack((Z, X))
            m += N
            std_Z = np.std(Z[m0 + 1:m, :d], axis=0)

        if t % steps == 0:
            # Print progress
            if t > 1:
                print(f"\rDREAM_(ZS) progress, % done: {100 * (t / T):.2f}", end='')

            # Convergence diagnostics (Placeholder function)
            GR_uniR, GR_multiR = MODELAVG_gelman(chain[t // 2 + 1: t, :d, :], t)
            R_stat[ct, 0:d+1] = np.hstack((t, GR_uniR))
            MR_stat[ct, 0:2] = [t, GR_multiR]
            AR[ct, 0:2] = [t, 100 * (accept / (N * (t - t_old)))]
            accept = 0
            t_old = t
            ct += 1

    print('\n')
    output = {'AR': AR, 'R_stat': R_stat, 'MR_stat': MR_stat}
    return chain, output


# Placeholder function for bnd_check
def bnd_check(Xp, X, K, parmin, parmax):
    # Apply boundary checking as needed (adapt from the original code)
    method = 'reflection'
    ii_low = np.where(Xp < parmin)
    ii_up = np.where(Xp > parmax)

    if method == 'set to bound':        # not recommended
        Xp[ii_low] = parmin[ii_low]
        Xp[ii_up] = parmax[ii_up]
    elif method == 'reflection':        # optimization methods, ok but not reversible
        Xp[ii_low] = 2 * parmin[ii_low] - Xp[ii_low]
        Xp[ii_up] = 2 * parmax[ii_up] - Xp[ii_up]
    elif method == 'folding':           # satisfied detailed balance: see Vrugt and Ter Braak, 2011: HESS
        Xp[ii_low] = parmax[ii_low] - (parmin[ii_low] - Xp[ii_low])
        Xp[ii_up] = parmin[ii_up] + (Xp[ii_up] - parmax[ii_up])
    elif method == 'alternative':
        Xp[ii_low] = np.random.uniform(parmin[ii_low], X[ii_low])
        Xp[ii_up] = np.random.uniform(X[ii_up], parmax[ii_up])

    # Check all steps as reflection & folding can still go outside the prior ranges
    ii_low = np.where(Xp < parmin)
    ii_up = np.where(Xp > parmax)
    Xp[ii_low] = parmin[ii_low]
    Xp[ii_up] = parmax[ii_up]

    # Now make sure that first K parameter values [= weights] add up to one
    Xp[:, :K] = Xp[:, :K] / Xp[:, :K].sum(axis=1, keepdims=True)

    return Xp


def MODELAVG_gelman(chain, t):
    """
    Computes the Gelman-Rubin convergence diagnostics (\hat{R}) for MCMC chains.
    
    Parameters:
    - chain: A 3D numpy array with shape (n, d, N) representing the Markov chains.
      n = number of iterations, d = number of parameters, N = number of chains.
    - t: The current iteration number (used for logging warnings).
    
    Returns:
    - hatR: Univariate Gelman-Rubin diagnostic for each parameter (1D array of size d).
    - hatRd: Multivariate Gelman-Rubin diagnostic (scalar).
    """
    n, d, N = chain.shape
    
    # Early exit if there are fewer than 10 iterations
    if n < 10:
        return np.nan * np.ones(d), np.nan

    # STEP 0: Compute the chain means and store in a N x d matrix
    mu_chains = np.mean(chain, axis=0).T  # N x d
    # STEP 1: Compute the N within-chain variances
#    s2_chains = np.array([np.var(chain[:, :, i], axis=0) for i in range(N)])  # N x d
    s2_chains = np.array([np.var(chain[:, :, i], axis=0, ddof = 1) for i in range(N)])  # N x d
    # STEP 2: Compute the N within-chain covariances
    cov_chains = np.array([np.cov(chain[:, :, i], rowvar=False) for i in range(N)])  # N x d x d
    
    # Univariate hatR diagnostics
    # STEP 1: Compute variance B of N chain means
    B = n * np.var(mu_chains, axis=0)  # d
    # STEP 2: Compute 1xd vector W with mean of within-chain variances
    W = np.mean(s2_chains, axis=0)  # d
    # STEP 3: Estimate target variance = sum of within- and between-chain s2's
    sigma2 = ((n - 1) / n) * W + (1 / n) * B
    # STEP 4: Compute univariate hatR diagnostic for each parameter
    hatR = np.sqrt((N + 1) / N * (sigma2 / W) - (n - 1) / (N * n))

    # Multivariate hatRd diagnostic
    # STEP 1: Compute dxd matrix with mean W of within-chain covariances
    W_cov = np.mean(cov_chains, axis=0) + np.finfo(float).eps * np.eye(d)
    # STEP 2: Compute covariance B of N chain means
    B_cov = np.cov(mu_chains.T) + np.finfo(float).eps * np.eye(d)
    # STEP 3: Compute multivariate scale reduction factor, hatRd
    hatRd = np.sqrt((N + 1) / N * np.max(np.abs(np.linalg.eigvals(np.linalg.inv(W_cov) @ B_cov))) + (n - 1) / n)
    
    # Step #1: Generate a matrix of small values (eps * identity matrix)
    # identity_matrix = np.eye(d, d, dtype=np.float32)
    # small_matrix = np.finfo(np.float32).eps * np.random.randn(d, d).astype(np.float32) * identity_matrix
    # Step #2: Generate an 8x1 matrix of random values from the standard normal distribution
    # W_inv_B = np.linalg.solve(W_cov + small_matrix, B_cov)
    # eig_values = np.linalg.eigvals(W_inv_B)
    # hatRd = np.sqrt(((N + 1) / N) * np.max(np.abs(eig_values)) + (N - 1) / N)

    # Check if any warnings were issued during computation
    if np.any(np.isnan(hatR)) or np.any(np.isnan(hatRd)):
        with open("warning_file.txt", "a+") as fid:
            warning_message = (f"MODELAVG WARNING: NaN detected in multivariate "
                               f"R-statistic of Brooks and Gelman at iteration {t}.\n")
            fid.write(warning_message)

    return hatR, hatRd


# Latin Hypercube Sampling function
def LH_sampling(mn, mx, N):
    """
    Latin Hypercube Sampling.
    
    Args:
        mn: Lower bound vector
        mx: Upper bound vector
        N: Number of samples to generate

    Returns:
        N x d matrix of Latin Hypercube samples
    """
    d = len(mn)  # Number of parameters
    rng = np.array(mx) - np.array(mn)  # 1 x d vector with parameter ranges
    y =  np.random.rand(N, d)  # N x d matrix with uniform random labels
    # really important change below so that X stays in bound! as list is from 0 - N-1 rather than 1 to N
    id_matrix = 1 + np.argsort(np.random.rand(N, d), axis=0)  # Random sort (1:N without replacement)
    M = (id_matrix - y) / N  # Multiplier matrix (y introduces randomness)
    R = np.add(np.multiply(M, rng), mn)  # N x d matrix of stratified LH samples
    
    return R


def MODELAVG_end(method, D, y, beta, chain, options, output):
    """
    Postprocess results of MODELAVG toolbox (performance of BMA model, scoring rules, etc).

    Arguments:
    method : str
        The method used for the model (e.g., 'bma', 'mma', etc.)
    D : ndarray
        Forecast matrix of ensemble members (n x K).
    y : ndarray
        Verifying data (ground truth).
    beta : ndarray
        Parameter estimates.
    chain : ndarray
        Markov chain of posterior samples (for MCMC methods).
    options : dict
        Configuration options including alpha, postproc, etc.
    output : dict
        Output dictionary to store results.

    Returns:
    beta : ndarray
        The optimal parameter estimates after processing.
    output : dict
        Updated output dictionary with additional results.
    str_plot : str
        Plot string for further analysis or visualization.
    """
    # Compute size of ensemble and number of significance levels
    n, K = D.shape
    
    if len(options['alpha'].shape) == 1:
        options['alpha'] = options['alpha'].reshape(1, -1)  # Convert 1D array to 2D (1 row, multiple columns)
    _, na = options['alpha'].shape
    
    # Open a file to store warnings
    with open('warning_file.txt', 'a+') as fid:

        # Derive estimate of beta from posterior values
        if method in ['bma', 'mma', 'mma-s']:
            # First assemble all chain into one matrix
            beta_matrix = genparset(chain)  # Function 'genparset' needs to be defined elsewhere
            N, d_all = beta_matrix.shape
            d = d_all - 1

            # Check if DREAM converged
            if not np.all(output['R_stat'][-1, 1:d + 1] < 1.2):
                evalstr = 'MODELAVG WARNING: DREAM DID NOT CONVERGE FORMALLY ACCORDING TO MULTIVARIATE R-STATISTIC\n'
                print(evalstr)
                fid.write(evalstr)

            # Take the last 25% of the posterior samples
            beta = beta_matrix[int(3/4 * N):, :d_all]

            # Find the maximum a posteriori (MAP) parameter values
            idx = np.argmax(beta[:, -1])
            beta_opt = beta[idx, :d]

            if 'CORR' not in output:
                # Initialize 'par_unc' with the correct shape, for example, based on y_opt and D
                output['CORR'] =  np.full((d, d), np.nan)  # Or adjust the shape as needed
            # Calculate posterior parameter correlation matrix (R-values)
            output['CORR'] = np.corrcoef(beta[:, :d], rowvar=False)

            if 'loglik' not in output:
                output['loglik'] =  np.full((1, 1), np.nan)  # Or adjust the shape as needed
            # Maximum log likelihood value
            output['loglik'] = beta[idx, -1]

            if 'loglik' not in output:
                output['loglik'] =  np.full((d, 1), np.nan)  # Or adjust the shape as needed
            # Standard deviation of DREAM
            output['std'] = np.std(beta[:, :d], axis=0)

        else:
            beta_opt = beta

        # Get the confidence level(s) from significance level(s)
        g = 1 - options['alpha']
        gam = np.sort(np.concatenate([ (1 - g) / 2, (1 - g) / 2 + g ]))
        # added by JAV
        gam = gam.reshape(1,-1)
        ng = g.size
        output['par_unc'] =  np.full((n, 2*ng), np.nan)  
        output['pred'] =  np.full((n, 2*ng), np.nan)
        output['coverage'] = np.full((1,na), np.nan)
        output['spread'] = np.full((1,na), np.nan)
        # end added by JAV

        if method == 'bma':
            # Extract posterior draws of weights
            beta_weights = beta[:, :K]
            # Calculate the gamma% posterior simulation uncertainty due to parameter uncertainty
 #           output['par_unc'] = np.percentile(beta_weights @ D.T, 100 * gam, axis=1).T
            output['par_unc'] = np.percentile(beta_weights @ D.T, 100 * gam, axis=0).T
            # must reshape nxngx1 into nxng
            output['par_unc'] = output['par_unc'].squeeze()
        else:
            # Compute variance of averaged model
            y_opt = D @ beta_opt[:K]
            df = n - K
            sigma2_avg = np.sum((y - y_opt) ** 2) / df
            s_avg = np.sqrt(sigma2_avg)

            # Covariance matrix of weights
            C = sigma2_avg * np.linalg.inv(D.T @ D)
            t_crit = t.ppf(gam, df)
            
            for z in range(n):   # changed t to z in loop as to avoid conflict with built-in t distribution 
                output['par_unc'][z, :] = y_opt[z] + t_crit * np.sqrt(D[z, :] @ C @ D[z, :].T)
                output['pred'][z, :] = y_opt[z] + t_crit * s_avg

            output['std'] = np.sqrt(np.diag(C))

        # Calculate point forecast, RMSE, R (Pearson correlation coefficient)
        output['Ye'] = D @ beta_opt[:K]
        Ym = np.mean(y)
        SStot = np.sum((y - Ym) ** 2)
        e = y - output['Ye']
        SSres = np.sum(e ** 2)
        output['R2'] = 1 - SSres / SStot
        output['RMSE'] = np.sqrt(np.sum(e ** 2) / n)
        output['R'] = pearsonr(output['Ye'], y)[0]
        output['KGE'] = compKGE(y, output['Ye'])  # Function 'compKGE' needs to be defined elsewhere
        output['ML'] = beta_opt

        # Derive printing string for postprocessing
        str_table, str_plot = create_str_plot(method, options, K)  # Function 'create_str_plot' needs to be defined elsewhere
        output['str_table'] = str_table

        # Check whether we compute quantiles
        if options['postproc'] == 'no':
            return beta, output, str_plot

        if method == 'bma':
            # Derive quantiles for BMA
            output['pred'], output['pdf_Y'], output['cdf_Y'] = BMA_quantile(beta_opt, D, y, options)  # Function 'BMA_quantile' needs to be defined elsewhere
            # Compute scoring rules of BMA model
            output = BMA_score(beta_opt, D, y, options, output)  # Function 'BMA_score' needs to be defined elsewhere
        
        for i in range(na):
            # Compute coverage and spread
#            output['coverage'][i] = 100 * np.sum((output['pred'][:, i] < y) & (y < output['pred'][:, 2 * na - (i - 1)])) / n
#            output['spread'][i] = np.mean(output['pred'][:, 2 * na - (i - 1)] - output['pred'][:, i])
            output['coverage'][0,i] = 100 * np.sum((output['pred'][:, i] < y) & (y < output['pred'][:, 2 * na - i - 1])) / n
            output['spread'][0,i] = np.mean(output['pred'][:, 2 * na - i - 1] - output['pred'][:, i])

    return beta, output, str_plot


def genparset(chain):
    """
    Generates a 2D matrix `parset` from a 3D array `chain`.

    Parameters:
    chain (numpy.ndarray): A 3D array with shape (T, d, N)

    Returns:
    numpy.ndarray: A 2D matrix `parset` derived from the `chain`
    """
    T, d, N = chain.shape  # Get dimensions of chain
    
    parset = []  # Initialize parset as an empty list

    # If T is 0, do nothing (parset remains empty)
    if T != 0:
        # Generate parset from all chain elements
        for n in range(N):
            # Append the values to the parset list, concatenate the third dimension
            # along with the time indices (1:T)
            parset.append(np.hstack((chain[:, :, n], np.arange(1, T + 1).reshape(-1, 1))))
        
        # Convert the list to a 2D numpy array
        parset = np.vstack(parset)
        
        # Sort by the last column (the time indices)
        parset = parset[np.argsort(parset[:, d], axis=0)]
        
        # Return only the first d columns (drop the time indices)
        parset = parset[:, :d]

    return parset


def create_str_plot(method, options, K):
    """
    Generates strings for output tables and postprocessor figures.

    Parameters:
    method (str): The method to generate the plot and table strings for.
    options (dict): A dictionary containing 'VAR', 'PDF', and 'TAU' options.
    K (int): The number of elements for the beta parameters.

    Returns:
    str_table (list): The list for the table output (non-LaTeX).
    str_plot (list): The list for the plot output (LaTeX formatted).
    """
    # First print only the beta values
    str_ = [f'beta_{{{i}}}' for i in range(1, K+1)]
    
    # Now for BMA method expand the vector
    if method == 'bma':
        if options['VAR'] == '1':
            str_.append('sigma')
        elif options['VAR'] == '2':
            str_.extend([f'sigma_{{{i}}}' for i in range(1, K+1)])
        elif options['VAR'] == '3':
            str_.append('c')
        elif options['VAR'] == '4':
            str_.extend([f'c_{{{i}}}' for i in range(1, K+1)])

        # Additional parameters based on PDF and TAU options
        if options['PDF'] in ['gen_normal', 'gev', 'gpareto']:
            if options['PDF'] == 'gen_normal':
                if options['TAU'] == '1':
                    str_.append('tau')
                elif options['TAU'] == '2':
                    str_.extend([f'tau_{{{i}}}' for i in range(1, K+1)])

            elif options['PDF'] == 'gev':
                if options['TAU'] == '1':
                    str_.append('xi')
                elif options['TAU'] == '2':
                    str_.extend([f'xi_{{{i}}}' for i in range(1, K+1)])

            elif options['PDF'] == 'gpareto':
                if options['TAU'] == '1':
                    str_.append('k')
                elif options['TAU'] == '2':
                    str_.extend([f'k_{{{i}}}' for i in range(1, K+1)])

    # Table does not need LaTeX print, so return as it is
    str_table = str_

    # LaTeX formatting for plotting
    str_plot = [s.replace('beta', '\\beta')
                .replace('sigma', '\\sigma')
                .replace('tau', '\\tau')
                .replace('xi', '\\xi') for s in str_]
    
    # Enclose each string in LaTeX math mode '$...$'
    str_plot = [f'${item}$' for item in str_plot]

    return str_table, str_plot


def BMA_lik(x, D, Y, options):
    """
    Calculate log-likelihood of BMA mixture distribution
    
    Arguments:
    x : array-like
        1D vector with values of BMA parameters (weights of members and variance/shape parameters)
    D : array-like
        2D matrix (n x K) with forecasts of ensemble members
    Y : array-like
        2D matrix (n x K) with verifying data (ground truth)
    options : dict
        Dictionary with BMA algorithmic variables:
        - 'VAR': Variance model (1: common, 2: individual, 3: non-constant, 4: individual non-constant)
        - 'PDF': Conditional distribution type (e.g., 'normal', 'lognormal', 'gamma', etc.)
        - 'CPU': Whether to use faster computation (e.g., 'yes')
        - 'TAU': For 'gen_normal', 'gev', and 'gpareto', the shape parameter method ('1' for common tau, '2' for individual tau)
    
    Returns:
    ell : float
        Scalar log-likelihood of BMA mixture distribution
    L : ndarray
        2D matrix of likelihoods of ensemble members (n x K)
    """
    n, K = D.shape  # number of forecasts and ensemble members
    beta = x[:K]    # Weights for ensemble members
    
    # Variance calculation
    if options['VAR'] == '1':  # Common constant variance
        S = x[K] * np.ones((n, K))
    elif options['VAR'] == '2':  # Individual constant variance
        S = np.multiply(x[K:2*K], np.ones((n, K)))
    elif options['VAR'] == '3':  # Common non-constant variance
        c = x[K]
        S = c * np.abs(D)
    elif options['VAR'] == '4':  # Individual non-constant variance
        c = x[K:2*K]
        S = np.multiply(c, np.abs(D))
    else:
        raise ValueError("Unknown variance option")
    
    S = np.maximum(S, np.finfo(float).eps)  # Ensure no zero variance

    # Conditional distribution
    if options['PDF'] == 'normal':  # Normal distribution
        if options['CPU'] == 'yes':  # Fast computation
            L = np.exp(-0.5 * ((Y - D) / S) ** 2) / (np.sqrt(2 * np.pi) * S)
        else:  # Built-in computation
            L = norm.pdf(Y, D, S)
    
    elif options['PDF'] == 'lognormal':  # Lognormal distribution
        logn_impl = 2
        if logn_impl == 2:
            sigma2_fxs = S ** 2
            A = sigma2_fxs / (D ** 2)
            S2 = np.log(A + 1)
            S = np.sqrt(S2)
            mu = np.log(np.abs(D)) - S2 / 2
        
        if options['CPU'] == 'yes':  # Fast computation
            L = np.exp(-0.5 * ((np.log(Y) - mu) / S) ** 2) / (np.sqrt(2 * np.pi) * Y * S)
        else:  # Built-in computation
            L = lognorm.pdf(Y, S, scale=np.exp(mu))
    
    elif options['PDF'] == 'tnormal':  # Truncated Normal distribution
        mu = np.abs(D)
        a, b = 0, np.inf  # truncation limits
        Alfa = (a - mu) / S
        Beta = (b - mu) / S
        # L = truncnorm.pdf(Y, a, b, loc=mu, scale=S)
        L = truncnorm.pdf(Y, Alfa, Beta, loc=mu, scale=S)

    elif options['PDF'] == 'gen_normal':  # Generalized Normal distribution
        #tau = x[K] if options['TAU'] == '1' else x[K-K+1:K]
        d = len(x)
        if options['TAU'] == '1':
            tau = x[d] * np.ones((n, K))
        elif options['TAU'] == '2':
            tau = np.tile(x[d-K:d], (n, 1))
        L = tau / (2 * S * np.repeat(scipy.special.gamma(1 / tau), n, axis=0)) * np.exp(- (np.abs(Y - D) / S) ** tau)
    
    elif options['PDF'] == 'gamma':  # Gamma distribution
        mu = np.abs(D)
        A = mu ** 2 / S ** 2
        B = S ** 2 / mu
        if options['CPU'] == 'yes':  # Fast computation
            z = Y / B
            u = (A - 1) * np.log(z) - z - np.gammaln(A)
            L = np.exp(u) / B
        else:  # Built-in computation
            L = gamma.pdf(Y, A, scale=B)
    
    elif options['PDF'] == 'weibull':  # Weibull distribution
        wbl_impl = 2
        if wbl_impl == 1:
            Kk = S
            Lambda = np.abs(D) / scipy.special.gamma(1 + 1 / Kk)
            Lambda = np.maximum(Lambda, np.finfo(float).eps)  # Avoid zero Lambda
        elif wbl_impl == 2:
            Lambda = S
            X = np.abs(D) / Lambda
            c = np.sqrt(2 * np.pi) / np.exp(1) - scipy.special.gamma(1.461632)
            Lx = lambda x: np.log((x + c) / np.sqrt(2 * np.pi))
            A = Lx(X)
            B = A / np.real(lambertw(A / np.exp(1))) + 1 / 2
            Kk = 1 / (B - 1)
        if options['CPU'] == 'yes':  # Fast computation
            L = Kk/Lambda * (Y/Lambda)^(Kk-1) * np.exp(-(Y/Lambda)^Kk);
        else:  # Built-in computation
            L = weibull_min.pdf(Y, Kk, scale=Lambda)
    
    elif options['PDF'] == 'gev':  # Generalized Extreme Value (GEV) distribution
        # xi = np.repeat(x[K], n, axis=0) if options['TAU'] == '1' else np.repeat(x[K-K+1:K], n, axis=0)
        d = len(x)
        if options['TAU'] == '1':
            xi = x[d] * np.ones((n, K))
        elif options['TAU'] == '2':
            xi = np.tile(x[d-K:d], (n, 1))
        g1 = scipy.special.gamma(1 - xi)
        mu = D - (g1 - 1) * S / xi
        # Handle case where xi is very close to zero
        Eul_constant = 0.577215664901532860606512090082
        mu[np.abs(xi) < np.finfo(float).eps] = D[np.abs(xi) < np.finfo(float).eps] - S[np.abs(xi) < np.finfo(float).eps] * Eul_constant
        L = genextreme.pdf(Y, xi, loc=mu, scale=S)
    
    elif options['PDF'] == 'gpareto':  # Generalized Pareto distribution
        # kk = x[K] if options['TAU'] == '1' else np.repeat(x[K-K+1:K], n, axis=0)
        d = len(x)
        if options['TAU'] == '1':
            kk = x[d] * np.ones((n, K))
        elif options['TAU'] == '2':
            kk = np.tile(x[d-K:d], (n, 1))
        Theta = np.abs(D) - S / (1 - kk)
        L = genpareto.pdf(Y, kk, loc=Theta, scale=S)
    
    else:
        raise ValueError("Unknown PDF option")
    
    # Replace likelihood with zero if D < 0 (only needed for certain distributions)
    # L[D < 0] = 0

    # Calculate the likelihood weighted by beta
    lik = np.dot(L, beta) + np.finfo(float).tiny
    ell = np.sum(np.log(lik))  # Log-likelihood of the BMA model
    
    return ell #, L


def BMA_quantile(x, D, y, options):
    n, K = D.shape  # Unpack number of observations
    beta = x[:K]  # Unpack weights
    count = 0  # Set counter to zero
    pdf_y = np.nan * np.ones(n)  # Initialize PDF observations
    cdf_y = np.nan * np.ones(n)  # Initialize CDF observations

    # Unpack alpha values
    g = 1 - np.array(options['alpha'])
    gam = np.sort([(1 - g) / 2, (1 - g) / 2 + g])  # Quantile levels
    # added by JAV to make sure that gam is a row vector with nA elements
    gam = gam.reshape(-1)
    nA = len(gam)

    # Set variance based on the options.VAR setting
    if options['VAR'] == '1':  # Common constant variance
        S = x[K] * np.ones((n, K))
    elif options['VAR'] == '2':  # Individual constant variance
        S = np.tile(x[K:2*K], (n, 1))
    elif options['VAR'] == '3':  # Common non-constant variance
        c = x[K]
        S = c * np.abs(D)
    elif options['VAR'] == '4':  # Individual non-constant variance
        c = x[K:2*K]
        S = np.multiply(c[:, None], np.abs(D))
    else:
        raise ValueError("Unknown variance option")

    # Ensure S is positive
    S = np.abs(S)

    # Conditional distribution case
    if options['PDF'] == 'normal':
        Mu = D  # Mean of normal distribution

        def PDF(x, i):
            result = beta[0] * norm.pdf(x, Mu[i, 0], S[i, 0])
            for k in range(1, K):
                result += beta[k] * norm.pdf(x, Mu[i, k], S[i, k])
            return result

        def CDF(x, i):
            result = beta[0] * norm.cdf(x, Mu[i, 0], S[i, 0])
            for k in range(1, K):
                result += beta[k] * norm.cdf(x, Mu[i, k], S[i, k])
            return result

    elif options['PDF'] == 'lognormal':
        logn_impl = 2
        if logn_impl == 1:
            Mu = np.log(np.abs(D)) - S**2 / 2  # Mean for lognormal distribution
        elif logn_impl == 2:
            sigma2_fxs = S**2
            A = sigma2_fxs / (D**2)
            S2 = np.log(A + 1)
            S = np.sqrt(S2)
            Mu = np.log(np.abs(D)) - S2 / 2
        
        def PDF(x, i):
            result = beta[0] * lognorm.pdf(x, S[i, 0], scale=np.exp(Mu[i, 0]))
            for k in range(1, K):
                result += beta[k] * lognorm.pdf(x, S[i, k], scale=np.exp(Mu[i, k]))
            return result

        def CDF(x, i):
            result = beta[0] * lognorm.cdf(x, S[i, 0], scale=np.exp(Mu[i, 0]))
            for k in range(1, K):
                result += beta[k] * lognorm.cdf(x, S[i, k], scale=np.exp(Mu[i, k]))
            return result

    elif options['PDF'] == 'tnormal':
        Mu = np.abs(D)  # Mode truncated normal
        a = 0  # Lower end point
        b = np.inf  # Upper end point
        Alfa = (a - Mu) / S
        Beta = (b - Mu) / S

        def PDF(x, i):
            result = beta[0] * truncnorm.pdf(x, Alfa[i,0], Beta[i,0], loc=Mu[i, 0], scale=S[i, 0])
            for k in range(1, K):
                result += beta[k] * truncnorm.pdf(x, Alfa[i,k], Beta[i,k], loc=Mu[i, k], scale=S[i, k])
            return result

        def CDF(x, i):
            result = beta[0] * truncnorm.cdf(x, Alfa[i,0], Beta[i,0], loc=Mu[i, 0], scale=S[i, 0])
            for k in range(1, K):
                result += beta[k] * truncnorm.cdf(x, Alfa[i,k], Beta[i,k], loc=Mu[i, k], scale=S[i, k])
            return result

    elif options['PDF'] == 'gen_normal':
        Mu = D
        d = len(x)
        if options['TAU'] == '1':
            tau = x[d] * np.ones((n, K))
        elif options['TAU'] == '2':
            tau = np.tile(x[d-K:d], (n, 1))

        def PDF(x, i):
            result = beta[0] * gnormpdf(x, Mu[i, 0], S[i, 0], tau[i,0])
            for k in range(1, K):
                result += beta[k] * gnormpdf(x, Mu[i, k], S[i, k], tau[i,k])
            return result

        def CDF(x, i):
            result = beta[0] * gnormcdf(x, Mu[i, 0], S[i, 0], tau[i,0])
            for k in range(1, K):
                result += beta[k] * gnormcdf(x, Mu[i, k], S[i, k], tau[i,k])
            return result

    elif options['PDF'] == 'gamma':
        Mu = np.abs(D)  # Mean of gamma distribution
        A = (Mu**2) / S**2  # Shape parameter
        B = S**2 / Mu  # Scale parameter

        def PDF(x, i):
            result = beta[0] * gamma.pdf(x, A[i, 0], scale=B[i, 0])
            for k in range(1, K):
                result += beta[k] * gamma.pdf(x, A[i, k], scale=B[i, k])
            return result

        def CDF(x, i):
            result = beta[0] * gamma.cdf(x, A[i, 0], scale=B[i, 0])
            for k in range(1, K):
                result += beta[k] * gamma.cdf(x, A[i, k], scale=B[i, k])
            return result

    elif options['PDF'] == 'weibull':
        wbl_impl = 2
        if wbl_impl == 1:
            Kk = S
            Lambda = np.abs(D) / scipy.special.gamma(1 + 1 / Kk)
            Lambda = np.maximum(Lambda, np.finfo(float).eps)  # Avoid zero Lambda
        elif wbl_impl == 2:
            Lambda = S
            X = np.abs(D) / Lambda
            c = np.sqrt(2 * np.pi) / np.exp(1) - scipy.special.gamma(1.461632)
            Lx = lambda x: np.log((x + c) / np.sqrt(2 * np.pi))
            A = Lx(X)
            B = A / np.real(lambertw(A / np.exp(1))) + 1 / 2
            Kk = 1 / (B - 1)

        def PDF(x, i):
            result = beta[0] * weibull_min.pdf(x, Kk[i, 0], scale=Lambda[i, 0])
            for k in range(1, K):
                result += beta[k] * weibull_min.pdf(x, Kk[i, k], scale=Lambda[i, k])
            return result

        def CDF(x, i):
            result = beta[0] * weibull_min.cdf(x, Kk[i, 0], scale=Lambda[i, 0])
            for k in range(1, K):
                result += beta[k] * weibull_min.cdf(x, Kk[i, k], scale=Lambda[i, k])
            return result

    elif options['PDF'] == 'gev':
        d = len(x)
        if options['TAU'] == '1':
            xi = x[d] * np.ones((n, K))
        elif options['TAU'] == '2':
            xi = np.tile(x[d-K:d], (n, 1))

        g = lambda a, xi: scipy.special.gamma(1 - a * xi)
        Mu = D - (g(1, xi) - 1) * S / xi

        # Handle case where xi is very close to zero
        Eul_constant = 0.577215664901532860606512090082
        Mu[np.abs(xi) < np.finfo(float).eps] = D[np.abs(xi) < np.finfo(float).eps] - S[np.abs(xi) < np.finfo(float).eps] * Eul_constant

        def PDF(x, i):
            result = beta[0] * genextreme.pdf(x, xi[i, 0], loc=Mu[i, 0], scale=S[i, 0])
            for k in range(1, K):
                result += beta[k] * genextreme.pdf(x, xi[i, k], loc=Mu[i, k], scale=S[i, k])
            return result

        def CDF(x, i):
            result = beta[0] * genextreme.cdf(x, xi[i, 0], loc=Mu[i, 0], scale=S[i, 0])
            for k in range(1, K):
                result += beta[k] * genextreme.cdf(x, xi[i, k], loc=Mu[i, k], scale=S[i, k])
            return result

    elif options['PDF'] == 'gpareto':
        d = len(x)
        if options['TAU'] == '1':
            kk = x[d] * np.ones((n, K))
        elif options['TAU'] == '2':
            kk = np.tile(x[d-K:d], (n, 1))

        Theta = np.abs(D) - S / (1 - kk)

        def PDF(x, i):
            result = beta[0] * genpareto.pdf(x, kk[i,0], loc=Theta[i, 0], scale=S[i, 0])
            for k in range(1, K):
                result += beta[k] * genpareto.pdf(x, kk[i, k], loc=Theta[i, k], scale=S[i, k])
            return result
        
        def CDF(x, i):
            result = beta[0] * genpareto.cdf(x, kk[i,0], loc=Theta[i, 0], scale=S[i, 0])
            for k in range(1, K):
                result += beta[k] * genpareto.cdf(x, kk[i, k], loc=Theta[i, k], scale=S[i, k])
            return result

    # Compute PDF and CDF at measured values
    for i in range(n):
        pdf_y[i] = PDF(y[i], i)
        cdf_y[i] = CDF(y[i], i)

    # Now activate the mixture CDF with alpha for root finding
    def Fx(x, i, gam):
      #  if x < 0:
      #      return 1e10
      #  else:
        return CDF(x, i) - gam

    # Now loop over each observation and determine BMA quantiles
    plimit = np.nan * np.ones((n, nA))
    x = np.linspace(min(0, 2 * np.min(np.min(D), axis=0)), 2 * np.max(y), int(1e5))  # Define x values
    calc_method = 1  # Set calculation method (1 or 2)
    for i in range(n):
        if i % (n // 25) == 0:
            if i > 0:
                print(f'BMA quantile calculation, {100*(i/n):.2f}% done', end='\r')

        if calc_method == 1:  # Exact answer using root finding of CDF, but slower
            for z in range(nA):
                # Find root using the fzero equivalent in Python (optimize.root_scalar)
                # plimit[i, z] = fsolve(lambda x: Fx(x, i, gam[z]), abs(np.mean(D[i, :])))
                plimit[i, z], info, ier, msg = fsolve(lambda x: Fx(x, i, gam[z]), y[i], full_output=True)
                if ier != 1:
                    plimit[i, z] = fsolve(lambda x: Fx(x, i, gam[z]), abs(np.mean(D[i, :])))
                # be careful => initial value must make sense, poor convergence shows in the quantile plot
                # MATLAB appears much more robust here, meaning the initial value of np.mean(abs(D[i, :])) works well across forecast PDFs and different variables/time series
                # Must check whether fsolve converged properly [= reported in output argument] and then use linear interpolation as alternative
                # plimit[i, z], info, ier, msg = fsolve(lambda x: Fx(x, i, gam[z]), abs(np.mean(D[i, :])), full_output=True)
                # check ier and do linear interpolation!
        elif calc_method == 2:  # Approximate answer using linear interpolation
            cdf_x = Fx(x, i, 0)  # Get the CDF values
            ii = np.diff(cdf_x) > 0  # Find where the CDF is increasing
            ii = np.concatenate(([True], ii))
            plimit[i, :nA] = np.interp(gam, cdf_x[ii], x[ii])  # Linear interpolation
    
    print("\n BMA quantile calculation complete.")

    return plimit, pdf_y, cdf_y


def BMA_crps(x, D, y, options):
    n, K = D.shape      # Unpack the number of observations
    beta = x[:K]        # Unpack weights
    count = 0           # Set counter to zero
    crps = np.full(n, np.nan)  # Initialize CRPS of BMA mixture CDF

    # Define tau points
    P1 = np.concatenate([np.array([0.001, 0.002, 0.005, 0.01, 0.02, 0.05]) , np.arange(0.05, 0.45, 0.05)])
    P = np.concatenate([P1, [0.5], 1 - P1[::-1]])  # tau points for integration
    nP = len(P)

    # Step 1: Determine the variance based on options.VAR
    if options['VAR'] == '1':  # common constant variance
        S = x[K] * np.ones((n, K))
    elif options['VAR'] == '2':  # individual constant variance
        S = np.tile(x[K:2*K], (n, 1))
    elif options['VAR'] == '3':  # common non-constant variance
        c = x[K]
        S = c * np.abs(D)
    elif options['VAR'] == '4':  # individual non-constant variance
        c = x[K:2*K]
        S = c[:, np.newaxis] * np.abs(D)
    else:
        raise ValueError("Unknown variance option")

    S = np.maximum(S, np.finfo(float).eps)  # Ensure S >= eps
    S2 = S ** 2  # variance matrix

    # Step 2: Define the CDF of the mixture distribution
    if options['PDF'] == 'normal':          # normal PDF
        Mu = D
        # Create CDF function for normal mixture
        def F_P(x, i, p):
            F = beta[0] * norm.cdf(x, Mu[i, 0], S[i, 0])
            for k in range(1, K):
                F += beta[k] * norm.cdf(x, Mu[i, k], S[i, k])
            return F - p

    elif options['PDF'] == 'lognormal':     # lognormal PDF
        Mu = np.log(np.abs(D)) - S2 / 2
        # Create CDF function for lognormal mixture
        def F_P(x, i, p):
            F = beta[0] * lognorm.cdf(x, S[i, 0], scale=np.exp(Mu[i, 0]))
            for k in range(1, K):
                F += beta[k] * lognorm.cdf(x, S[i, k], scale=np.exp(Mu[i, k]))
            return F - p

    elif options['PDF'] == 'tnormal':       # truncated normal PDF
        Mu = np.abs(D)
        a, b = 0, np.inf  # Lower and upper bounds for truncation
        Alfa = (a - Mu) / S
        Beta = (b - Mu) / S
        # Create CDF function for truncated normal mixture
        def F_P(x, i, p):
#            F = beta[0] * tnormcdf(x, Mu[i, 0], S[i, 0], a, b)
            F = beta[0] * truncnorm.cdf(x, Alfa[i,0], Beta[i,0], loc=Mu[i, 0], scale=S[i, 0])
            for k in range(1, K):
#                F += beta[k] * tnormcdf(x, Mu[i, k], S[i, k], a, b)
                F += beta[k] * truncnorm.cdf(x, Alfa[i,k], Beta[i,k], loc=Mu[i, k], scale=S[i, k])
            return F - p

    elif options['PDF'] == 'gen_normal':       # generalized normal PDF
        d = len(x)
        if options['TAU'] == '1':
            tau = x[d] * np.ones((n, K))
        elif options['TAU'] == '2':
            tau = np.tile(x[d-K:d], (n, 1))
        Mu = D
        # Create CDF function for generalized normal mixture
        def F_P(x, i, p):
            F = beta[0] * gnormcdf(x, Mu[i, 0], S[i, 0], tau[i,0])
            for k in range(1, K):
                F += beta[k] * gnormcdf(x, Mu[i, k], S[i, k], tau[i,k])
            return F - p

    elif options['PDF'] == 'gamma':         # gamma PDF
        Mu = np.abs(D)
        A = Mu**2 / S2
        B = S2 / Mu
        # Create CDF function for gamma mixture
        def F_P(x, i, p):
            F = beta[0] * gamma.cdf(x, A[i, 0], scale=B[i, 0])
            for k in range(1, K):
                F += beta[k] * gamma.cdf(x, A[i, k], scale=B[i, k])
            return F - p

    elif options['PDF'] == 'weibull':       # Weibull distribution parameters
        wbl_impl = 2
        if wbl_impl == 1:
            Kk = S
            Lambda = np.abs(D) / scipy.special.gamma(1 + 1 / Kk)
            Lambda = np.maximum(Lambda, np.finfo(float).eps)
        elif wbl_impl == 2:
            Lambda = S
            X = np.abs(D) / Lambda
            c = np.sqrt(2 * np.pi) / np.exp(1) - scipy.special.gamma(1.461632)
            A = np.log((X + c) / np.sqrt(2 * np.pi))
            B = A / np.real(np.log(A / np.exp(1))) + 1 / 2
            Kk = 1 / (B - 1)

        # Create CDF function for Weibull mixture
        def F_P(x, i, p):
            F = beta[0] * weibull_min.cdf(x, Kk[i, 0], scale=Lambda[i, 0])
            for k in range(1, K):
                F += beta[k] * weibull_min.cdf(x, Kk[i, k], scale=Lambda[i, k])
            return F - p

    elif options['PDF'] == 'gev':           # Generalized Extreme Value (GEV) distribution
        d = len(x)
        if options['TAU'] == '1':  # Common tau
            xi = x[d] * np.ones((n, K))
        elif options['TAU'] == '2':  # Individual tau
            xi = np.tile(x[d-K:d], (n, 1))
        
        g = lambda a, xi: scipy.special.gamma(1 - a * xi)
        Mu = D - (g(1, xi) - 1) * S / xi
        
        # Handle case where xi is very close to zero
        Eul_constant = 0.577215664901532860606512090082
        Mu[np.abs(xi) < np.finfo(float).eps] = D[np.abs(xi) < np.finfo(float).eps] - S[np.abs(xi) < np.finfo(float).eps] * Eul_constant

        # Create CDF function for GEV mixture
        def F_P(x, i, p):
            F = beta[0] * genextreme.cdf(x, xi[i, 0], loc=Mu[i, 0], scale=S[i, 0])
            for k in range(1, K):
                F += beta[k] * genextreme.cdf(x, xi[i, k], loc=Mu[i, k], scale=S[i, k])
            return F - p        

    elif options['PDF'] == 'gpareto':     # generalized Pareto PDF
        d = len(x)
        if options['TAU'] == '1':    # Common tau
            kk = x[d] * np.ones((n, K))
        elif options['TAU'] == '2':  # Individual tau
            kk = np.tile(x[d-K:d], (n, 1))

        Theta = np.abs(D) - S / (1 - kk)

        # Create CDF function for generalized Pareto mixture
        def F_P(x, i, p):
            F = beta[0] * genpareto.cdf(x, kk[i,0], loc=Theta[i, 0], scale=S[i, 0])
            for k in range(1, K):
                F += beta[k] * genpareto.cdf(x, kk[i, k], loc=Theta[i, k], scale=S[i, k])
            return F - p

    # Step 3: Compute CRPS for each observation
    FinvP = np.nan * np.ones(nP)
    for i in range(n):
        if i % (n // 25) == 0:
            if i > 0:
                print(f'BMA CRPS calculation, {100*(i/n):.2f}% done', end='\r')

        # Solve for values of y so that F_P^-1(y) = P using fsolve
        for z in range(nP):
            if z == 0:
#                FinvP[z] = fsolve(lambda y: F_P(y, i, P[z]), abs(np.mean(D[i, :]))) # initial value OK in MATLAB
                FinvP[z], info, ier, msg = fsolve(lambda y: F_P(y, i, P[z]), abs(np.mean(D[i, :])), full_output=True) # initial value OK in MATLAB 
            else:
#                FinvP[z] = fsolve(lambda y: F_P(y, i, P[z]), FinvP[z-1])
                FinvP[z], info, ier, msg = fsolve(lambda y: F_P(y, i, P[z]), FinvP[z-1], full_output=True) 
            if ier != 1:  # back-up if not converged properly
                FinvP[z] = fsolve(lambda y: F_P(y, i, P[z]), y[i])

        # Evaluate mixture CDF at omega
        Fy = F_P(y[i], i, 0)  # F_P evaluated at observation y(i)
        ii = P > Fy
        if np.sum(ii) == 0:
            crps[i] = y[i] * (1 - 2 * Fy) + 2 * np.trapezoid(P * FinvP,P)
        else:
#            crps[i] = y[i] * (1 - 2 * Fy) + 2 * np.trapezoid(P * FinvP,P) - 2 * np.trapezoid(np.concatenate(y[i],FinvP[ii]),np.concatenate(Fy,P[ii]))
            crps[i] = y[i] * (1 - 2 * Fy) + 2 * np.trapezoid(P * FinvP,P) - 2 * np.trapezoid(np.concatenate([y[i].flatten(), FinvP[ii].flatten()]),np.concatenate([Fy.flatten(), P[ii].flatten()]))

    print("\n BMA CRPS calculation complete.")
    return crps


def BMA_score(x, D, y, options, output):
    # Draw m samples BMA mixture PDF (= Phi) + return exact µ/var. BMA mixture
    m = 1000
    Phi, mu_mix, var_mix = BMA_draw(x, D, m, options)

    # Compute CRPS - quantile integration using BMA mixture CDF
    output['CRPS'] = BMA_crps(x, D, y, options)

    # Compute α-norm of BMA mixture distribution (robust option)
    lowup, mix_norm = BMA_norm(x, D, y, options, 'robust')

    # Compute some scoring rules - numerically
    nY = D.shape[0]
    output['QS'] = 2 * output['pdf_Y'] - mix_norm[:nY, 1]**2
    output['LS'] = np.log2(output['pdf_Y'])
    output['SS'] = output['pdf_Y'] / mix_norm[:nY, 1]

    # Compute Energy Score
#    output['ES'] = energy_score(Phi, y, 2)

#    output['ES'] = np.full(nY, np.nan)
#    count = 0
#    for t in range(nY):
#        # Print progress
#        if t % (nY // 25) == 0:
#            if t > 0:
#                print(f'\rEnergy score calculation, % done: {100 * (t / nY):3.2f}', end='')
#            else:
#                count = print(f'Energy score calculation, % done: {100 * (t / nY):3.2f}', end='')

        # Compute Energy score
#        output['ES'][t] = energy_score(Phi[t, :m], y[t], 2)
    
#    print('\n')

    # Compute exact reliability for BMA model
    output['mRLBL_anal'] = BMA_rlbl(output['cdf_Y'])

    # Compute coefficient of variation exactly
    Cv_anal = np.sqrt(var_mix) / mu_mix
    output['mCv_anal'] = np.mean(Cv_anal)

    # Compute time-averaged scores
    id1 = np.arange(nY)
    output['mQS'] = np.mean(output['QS'][id1])
    output['mLS'] = np.mean(output['LS'][id1])
    output['mSS'] = np.mean(output['SS'][id1])
    output['mCRPS'] = np.mean(output['CRPS'][id1])
  #  output['mES'] = np.mean(output['ES'][id1])

    # Store mixture normalizations and distributions
    output['mix_norm'] = mix_norm
    output['mu_mix'] = mu_mix
    output['var_mix'] = var_mix

    return output


def energy_score(fcst, obs, beta=1):
    """
    The ENERGY SCORE is a generalization of the continuous rank probability score.

    Args:
    - fcst (numpy.ndarray): n x m matrix of ensemble forecasts
    - obs (numpy.ndarray): n x 1 vector of measured data (observations)
    - beta (float): scalar with value of beta (default: 1)

    Returns:
    - mean_ES (float): Mean of non-missing energy score values
    - ES_value (numpy.ndarray): n x 1 vector with energy score values
    - num_nan (int): Number of missing (NaN) values of energy score
    """
    # Check dimensions of input
    n, m = fcst.shape  # n measurement times, m ensemble members
    
    if obs.shape[0] != n:
        raise ValueError("The length of the observation vector does not match the number of rows in the forecast matrix.")
    #if obs.shape[1] != 1:
    #    raise ValueError("The observation vector should have one column only.")
    
    ES_value = np.full(n, np.nan)  # Initialize energy score values
    ys = np.sort(fcst, axis=1)     # Sort forecast matrix rows in increasing order
    
    # Approach B1: For empirical CDF calculation
    for t in range(n):
        second_term = (1 / m) * np.sum(np.abs(ys[t, :m] - obs[t])**beta)
        first_term = 0
        for i in range(m):
            first_term += np.sum(np.abs(ys[t, :m] - ys[t, i])**beta)
        
        ES_value[t] = (1 / 2) * (1 / m**2) * first_term - second_term
    
    mean_ES = np.nanmean(ES_value)  # Compute mean of non-missing energy scores
    num_nan = np.sum(np.isnan(ES_value))  # Number of NaN values
    
    return mean_ES, ES_value, num_nan

# Example usage:
#fcst = np.random.normal(0, 1, (1000, 1000))  # Example ensemble forecasts (1000 observations, 1000 members)
#obs = np.random.rand(1000, 1)  # Example observations (1000 measurements)
#mean_ES, ES_value, num_nan = energy_score(fcst, obs)


def BMA_draw(x, D, N, options):
    """
    This function calculates random draws from a BMA (Bayesian Model Averaging) mixture distribution.

    Parameters:
        x (numpy array): 1D vector with max. likelihood BMA parameters
        D (numpy array): nxK matrix with forecasts of ensemble members
        N (int): Number of samples to draw from the BMA mixture distribution
        options (dict): Dictionary with BMA algorithmic variables, including variance and PDF options

    Returns:
        rnd (numpy array): nxN matrix of samples drawn from the BMA mixture distribution
        mu_mix (numpy array): nx1 vector mean of the BMA mixture distribution
        sigma2_mix (numpy array): nx1 vector variance of the BMA mixture distribution
    """

    n, K = D.shape  # Number of forecasts (n) and number of ensemble members (K)
    beta = x[:K]  # Unpack weights of member's conditional pdf
    mu_mix = D.dot(beta)  # Compute mean of the mixture

    rnd = np.full((n, N), np.nan)  # Initialize the random draw from the predictive PDF

    # Handle variance options based on input `options`
    if options['VAR'] == '1':  # Common constant variance
        S = x[K] * np.ones((n, K))
    elif options['VAR'] == '2':  # Individual constant variance
        S = np.multiply(x[K:2*K], np.ones((n, K)))
    elif options['VAR'] == '3':  # Common non-constant variance
        c = x[K]
        S = c * np.abs(D)
    elif options['VAR'] == '4':  # Individual non-constant variance
        c = x[K:2*K]
        S = np.multiply(c, np.abs(D))
    else:
        raise ValueError('Unknown variance option')

    S = np.maximum(S, np.finfo(float).eps)  # Ensure S is >= 2.22e-16
    S2 = S ** 2  # Compute variance matrix

    # Calculate mixture variance based on PDF options
    if options['PDF'] == 'normal':  # Normal distribution
        Mu = D
        sigma2_mix = (S2 + Mu ** 2) @ beta - mu_mix ** 2

    elif options['PDF'] == 'lognormal':  # Lognormal distribution
        mu_fxs = D
        logn_impl = 2
        if logn_impl == 1:
            Mu = np.log(np.abs(D)) - S2 / 2
            sigma2_fxs = D ** 2 * (np.exp(S2) - 1)
            sigma2_mix = (sigma2_fxs + mu_fxs ** 2) @ beta - mu_mix ** 2
        elif logn_impl == 2:
            sigma2_fxs = S2
            A = sigma2_fxs / (D ** 2)
            S2 = np.log(A + 1)
            S = np.sqrt(S2)
            Mu = np.log(np.abs(D)) - S2 / 2
            sigma2_mix = (sigma2_fxs + mu_fxs ** 2) @ beta - mu_mix ** 2

    elif options['PDF'] == 'tnormal':  # Truncated normal distribution
        Mu = np.abs(D)
        a = 0
        b = np.inf  # 'b' is unused, but kept for clarity
        Alfa = (a - Mu) / S    # Python implementation
        Beta = (b - Mu) / S    # Python implementation
        Z = 1 - norm.cdf(Alfa)
        mu_fxs = Mu + norm.pdf(Alfa) * S / Z
        sigma2_fxs = S**2 * (1 + Alfa * norm.cdf(Alfa) / Z - (norm.cdf(Alfa) / Z)**2)
        sigma2_mix = (sigma2_fxs + mu_fxs**2) @ beta - mu_mix ** 2

    elif options['PDF'] == 'gamma':  # Gamma distribution
        Mu = np.abs(D)
        A = Mu ** 2 / S2
        B = S2 / Mu
        sigma2_mix = (S2 + Mu ** 2) @ beta - mu_mix ** 2

    elif options['PDF'] == 'weibull':  # Weibull distribution
        Lambda = S
        Kk = np.abs(D) / Lambda
        Mu = np.abs(D)
        sigma2 = Lambda ** 2 * (scipy.special.gamma(1 + 2 / Kk) - scipy.special.gamma(1 + 1 / Kk) ** 2)
        sigma2_mix = (sigma2 + Mu ** 2) @ beta - mu_mix ** 2

    elif options['PDF'] == 'gen_normal':  # Generalized normal
        #tau = x[-1] if options['TAU'] == '1' else x[-K:]
        d = len(x)
        if options['TAU'] == '1':
            tau = x[d] * np.ones((n, K))
        elif options['TAU'] == '2':
            tau = np.tile(x[d-K:d], (n, 1))
        Mu = D
        sigma2 = S2 * scipy.special.gamma(3 / np.repeat(tau, n)) / scipy.special.gamma(1 / np.repeat(tau, n))
        sigma2_mix = (sigma2 + Mu ** 2) @ beta - mu_mix ** 2

    elif options['PDF'] == 'gev':  # Generalized extreme value (GEV) distribution
        d = len(x)
        if options['TAU'] == '1':  # Common tau
            xi = x[d] * np.ones((n, K))
        elif options['TAU'] == '2':  # Individual tau
            xi = np.tile(x[d-K:d], (n, 1))
#        Mu = D - (scipy.special.gamma(1, xi) - 1) * S / xi
        g = lambda a, xi: scipy.special.gamma(1 - a * xi)
        Mu = D - (g(1, xi) - 1) * S / xi
        mu_fxs, sigma2_fxs = gevstat(xi, S, Mu)
        sigma2_mix = (sigma2_fxs + mu_fxs ** 2) @ beta - mu_mix ** 2

    elif options['PDF'] == 'gpareto':     # Generalized Pareto distribution
        d = len(x)
        if options['TAU'] == '1':    # Common tau
            kk = x[d] * np.ones((n, K))
        elif options['TAU'] == '2':  # Individual tau
            kk = np.tile(x[d-K:d], (n, 1))

        Theta = np.abs(D) - S / (1 - kk)
        # Python built-in function
        mu_fxs, sigma2_fxs = gpstat(kk, S, Theta)
        sigma2_mix = (sigma2_fxs + mu_fxs ** 2) @ beta - mu_mix ** 2

    else:
        raise ValueError('Unknown PDF option')

    # Draw random numbers based on mixture
    r = np.random.rand(N)
    w = np.cumsum(beta)
    id_mix = K - np.sum(r[:, None] < w, axis=1) + 1

    # Assign elements to each mixture component and count the number of points per component
    id = [np.where(id_mix == k)[0] for k in range(1, K + 1)]
    n_id = [len(idx) for idx in id]

    # Sample from each mixture component
    for t in range(n):
        for k in range(K):
            if options['PDF'] == 'normal':
                rnd[t, id[k]] = norm.rvs(loc=Mu[t, k], scale=S[t, k], size=n_id[k])
            elif options['PDF'] == 'lognormal':
                rnd[t, id[k]] = lognorm.rvs(sigma=S[t, k], scale=np.exp(Mu[t, k]), size=n_id[k])
            elif options['PDF'] == 'gamma':
                rnd[t, id[k]] = gamma.rvs(A[t, k], scale=B[t, k], size=n_id[k])
            elif options['PDF'] == 'weibull':
                rnd[t, id[k]] = weibull_min.rvs(Kk[t, k], scale=Lambda[t, k], size=n_id[k])
            elif options['PDF'] == 'gen_normal':
                rnd[t, id[k]] = gnormrnd(Mu[t, k], S[t, k], tau[t,k], size=n_id[k])
            elif options['PDF'] == 'gev': # loc, scale
                rnd[t, id[k]] = genextreme.rvs(xi[t, k], loc=Mu[t, k], scale=S[t, k], size=n_id[k])
            elif options['PDF'] == 'gpareto':  # loc, scale
                rnd[t, id[k]] = genpareto.rvs(kk[t, k], loc=Theta[t, k], scale=S[t, k], size=n_id[k])
            elif options['PDF'] == 'tnormal': # loc, scale
                rnd[t, id[k]] = truncnorm.rvs(Alfa[t, k], Beta[t, k], loc=Mu[t, k], scale=S[t, k], size=n_id[k])
                #rnd[t, id[k]] = tnormal.rvs(xi[t, k], loc=Mu[t, k], scale=S[t, k], size=n_id[k])

    return rnd, mu_mix, sigma2_mix


def BMA_norm(x, D, y, options, method='robust'):
    # Unpack number of observations and ensemble members
    n, K = D.shape
    beta = x[:K]  # Mixture weights (first K elements)
    a = 1e-4  # Significance level for quantiles

    # Variance calculation based on options
    if options['VAR'] == '1':  # Common constant variance
        S = np.full_like(D, x[K], dtype=np.float64)
    elif options['VAR'] == '2':  # Individual constant variance
        S = np.tile(x[K:2*K], (n, 1))
    elif options['VAR'] == '3':  # Common non-constant variance
        c = x[K]
        S = c * np.abs(D)
    elif options['VAR'] == '4':  # Individual non-constant variance
        c = x[K:2*K]
        S = c * np.abs(D)
    else:
        raise ValueError('Unknown variance option')

    S = np.abs(S)  # Ensure positive variances

    # Initialize lowup and mix_norm arrays
    lowup = np.nan * np.ones((n, 2))
    mix_norm = np.nan * np.ones((n, 2))

    # Conditional distribution case
    if options['PDF'] == 'normal':
        Mu = D  # Mean of normal distribution

        def PDF(x, i):
            result = beta[0] * norm.pdf(x, Mu[i, 0], S[i, 0])
            for k in range(1, K):
                result += beta[k] * norm.pdf(x, Mu[i, k], S[i, k])
            return result

        def CDF(x, i):
            result = beta[0] * norm.cdf(x, Mu[i, 0], S[i, 0])
            for k in range(1, K):
                result += beta[k] * norm.cdf(x, Mu[i, k], S[i, k])
            return result

    elif options['PDF'] == 'lognormal':
        logn_impl = 2
        if logn_impl == 1:
            Mu = np.log(np.abs(D)) - S**2 / 2  # Mean for lognormal distribution
        elif logn_impl == 2:
            sigma2_fxs = S**2
            A = sigma2_fxs / (D**2)
            S2 = np.log(A + 1)
            S = np.sqrt(S2)
            Mu = np.log(np.abs(D)) - S2 / 2
        
        def PDF(x, i):
            result = beta[0] * lognorm.pdf(x, S[i, 0], scale=np.exp(Mu[i, 0]))
            for k in range(1, K):
                result += beta[k] * lognorm.pdf(x, S[i, k], scale=np.exp(Mu[i, k]))
            return result

        def CDF(x, i):
            result = beta[0] * lognorm.cdf(x, S[i, 0], scale=np.exp(Mu[i, 0]))
            for k in range(1, K):
                result += beta[k] * lognorm.cdf(x, S[i, k], scale=np.exp(Mu[i, k]))
            return result

    elif options['PDF'] == 'tnormal':
        Mu = np.abs(D)  # Mode truncated normal
        a = 0  # Lower end point
        b = np.inf  # Upper end point
        Alfa = (a - Mu) / S
        Beta = (b - Mu) / S

        def PDF(x, i):
            result = beta[0] * truncnorm.pdf(x, Alfa[i,0], Beta[i,0], loc=Mu[i, 0], scale=S[i, 0])
            for k in range(1, K):
                result += beta[k] * truncnorm.pdf(x, Alfa[i,k], Beta[i,k], loc=Mu[i, k], scale=S[i, k])
            return result

        def CDF(x, i):
            result = beta[0] * truncnorm.cdf(x, Alfa[i,0], Beta[i,0], loc=Mu[i, 0], scale=S[i, 0])
            for k in range(1, K):
                result += beta[k] * truncnorm.cdf(x, Alfa[i,k], Beta[i,k], loc=Mu[i, k], scale=S[i, k])
            return result

    elif options['PDF'] == 'gen_normal':
        Mu = D
        d = len(x)
        if options['TAU'] == '1':
            tau = x[d] * np.ones((n, K))
        elif options['TAU'] == '2':
            tau = np.tile(x[d-K:d], (n, 1))

        def PDF(x, i):
            result = beta[0] * gnormpdf(x, Mu[i, 0], S[i, 0], tau[i,0])
            for k in range(1, K):
                result += beta[k] * gnormpdf(x, Mu[i, k], S[i, k], tau[i,k])
            return result

        def CDF(x, i):
            result = beta[0] * gnormcdf(x, Mu[i, 0], S[i, 0], tau[i,0])
            for k in range(1, K):
                result += beta[k] * gnormcdf(x, Mu[i, k], S[i, k], tau[i,k])
            return result

    elif options['PDF'] == 'gamma':
        Mu = np.abs(D)  # Mean of gamma distribution
        A = (Mu**2) / S**2  # Shape parameter
        B = S**2 / Mu  # Scale parameter

        def PDF(x, i):
            result = beta[0] * gamma.pdf(x, A[i, 0], scale=B[i, 0])
            for k in range(1, K):
                result += beta[k] * gamma.pdf(x, A[i, k], scale=B[i, k])
            return result

        def CDF(x, i):
            result = beta[0] * gamma.cdf(x, A[i, 0], scale=B[i, 0])
            for k in range(1, K):
                result += beta[k] * gamma.cdf(x, A[i, k], scale=B[i, k])
            return result

    elif options['PDF'] == 'weibull':
        wbl_impl = 2
        if wbl_impl == 1:
            Kk = S
            Lambda = np.abs(D) / scipy.special.gamma(1 + 1 / Kk)
            Lambda = np.maximum(Lambda, np.finfo(float).eps)  # Avoid zero Lambda
        elif wbl_impl == 2:
            Lambda = S
            X = np.abs(D) / Lambda
            c = np.sqrt(2 * np.pi) / np.exp(1) - scipy.special.gamma(1.461632)
            Lx = lambda x: np.log((x + c) / np.sqrt(2 * np.pi))
            A = Lx(X)
            B = A / np.real(lambertw(A / np.exp(1))) + 1 / 2
            Kk = 1 / (B - 1)

        def PDF(x, i):
            result = beta[0] * weibull_min.pdf(x, Kk[i, 0], scale=Lambda[i, 0])
            for k in range(1, K):
                result += beta[k] * weibull_min.pdf(x, Kk[i, k], scale=Lambda[i, k])
            return result

        def CDF(x, i):
            result = beta[0] * weibull_min.cdf(x, Kk[i, 0], scale=Lambda[i, 0])
            for k in range(1, K):
                result += beta[k] * weibull_min.cdf(x, Kk[i, k], scale=Lambda[i, k])
            return result

    elif options['PDF'] == 'gev':
        d = len(x)
        if options['TAU'] == '1':
            xi = x[d] * np.ones((n, K))
        elif options['TAU'] == '2':
            xi = np.tile(x[d-K:d], (n, 1))

        g = lambda a, xi: scipy.special.gamma(1 - a * xi)
        Mu = D - (g(1, xi) - 1) * S / xi

        # Handle case where xi is very close to zero
        Eul_constant = 0.577215664901532860606512090082
        Mu[np.abs(xi) < np.finfo(float).eps] = D[np.abs(xi) < np.finfo(float).eps] - S[np.abs(xi) < np.finfo(float).eps] * Eul_constant

        def PDF(x, i):
            result = beta[0] * genextreme.pdf(x, xi[i, 0], loc=Mu[i, 0], scale=S[i, 0])
            for k in range(1, K):
                result += beta[k] * genextreme.pdf(x, xi[i, k], loc=Mu[i, k], scale=S[i, k])
            return result

        def CDF(x, i):
            result = beta[0] * genextreme.cdf(x, xi[i, 0], loc=Mu[i, 0], scale=S[i, 0])
            for k in range(1, K):
                result += beta[k] * genextreme.cdf(x, xi[i, k], loc=Mu[i, k], scale=S[i, k])
            return result

    elif options['PDF'] == 'gpareto':
        d = len(x)
        if options['TAU'] == '1':
            kk = x[d] * np.ones((n, K))
        elif options['TAU'] == '2':
            kk = np.tile(x[d-K:d], (n, 1))

        Theta = np.abs(D) - S / (1 - kk)

        def PDF(x, i):
            result = beta[0] * genpareto.pdf(x, kk[i,0], loc=Theta[i, 0], scale=S[i, 0])
            for k in range(1, K):
                result += beta[k] * genpareto.pdf(x, kk[i, k], loc=Theta[i, k], scale=S[i, k])
            return result
        
        def CDF(x, i):
            result = beta[0] * genpareto.cdf(x, kk[i,0], loc=Theta[i, 0], scale=S[i, 0])
            for k in range(1, K):
                result += beta[k] * genpareto.cdf(x, kk[i, k], loc=Theta[i, k], scale=S[i, k])
            return result

    # Now activate the mixture CDF with alpha for root finding
    def Fx(x, i, gam):
      #  if x < 0:
      #      return 1e10
      #  else:
        return CDF(x, i) - gam

    # Root finding to determine low and up quantiles (alpha/2 and 1-alpha/2)
    alpha = [a, (1 - a)]
    for i in range(n):
        for z in range(2):
            # Use fsolve to find the quantile corresponding to a or (1-a)
#            lowup[i, z] = fsolve(lambda y: CDF(y, i) - alpha[z], abs(np.mean(D[i, :])))
#            lowup[i, z] = fsolve(lambda y: CDF(y, i) - alpha[z], y[i])
            lowup[i, z], info, ier, msg = fsolve(lambda y: CDF(y, i) - alpha[z], abs(np.mean(D[i, :])), full_output=True) # initial value OK in MATLAB 
            if ier != 1:  # back-up if not converged properly
                lowup[i, z] = fsolve(lambda y: CDF(y, i) - alpha[z], y[i])

    #method = 'standard' #    method = 'robust'
    ny = 10000  # Number of discretization points for integration
    # Calculation of the α-norm using either 'standard' or 'robust' method
    for i in range(n):
        if i % (n // 25) == 0:
            if i > 0:
                print(f'BMA α-norm calculation, {100*(i/n):.2f}% done', end='\r')

        if method == 'standard':
            # Discretize between low and up quantiles
            yi = np.linspace(lowup[i, 0], lowup[i, 1], ny)
            # Evaluate the BMA mixture density at each yi
            pdf_yi = PDF(yi, i)
            # Compute δ-norm of BMA mixture: δ=1
            mix_norm[i, 0] = np.trapezoid(pdf_yi, yi); 
            # Compute δ-norm of BMA mixture: δ=2 for QS and SS
            mix_norm[i, 1] = np.sqrt(np.trapezoid(pdf_yi **2, yi))
        elif method == 'robust':
            for delta in [1, 2]:  # delta=1 and delta=2 for δ-norm
                # Robust method uses numerical integration with quad (more accurate)
                # out = quad(lambda y: PDF(y, i), lowup[i, 0], lowup[i, 1])[0]
                # Perform integration using scipy.integrate.quad (equivalent to MATLAB's integral)
                integral_result, _ = integrate.quad(lambda y: PDF(y, i) ** delta, lowup[i, 0], lowup[i, 1])
                mix_norm[i, delta - 1] = integral_result ** (1/delta)
            
    print("\n BMA α-norm calculation complete.")
                
    return lowup, mix_norm


def BMA_rlbl(cdf_y):
    """
    This function computes reliability of BMA forecast distribution using the formulation of Renard (2011).

    Parameters:
    cdf_y : array-like
        nx1 vector with CDF BMA mixture at verifying data.

    Returns:
    rlbl : float
        Scalar with reliability according to Renard (2011).
    """
    n = len(cdf_y)                              # Number of observations
    eCDF = np.sort(cdf_y)                       # Empirical CDF (CDF in ascending order)
    uCDF = (np.arange(1, n + 1) / n)            # CDF of uniform distribution
    rlbl = 2 / n * np.sum(np.abs(uCDF - eCDF))  # Reliability: negatively oriented, smaller = better

    return rlbl


def MODELAVG_postproc(method, D, y, beta, sigma_2, par, chain, options, output, str_plot):
    """
    Visualize results of MODELAVG toolbox: marginal distribution of parameters,
    chain convergence, autocorrelation functions, confidence/prediction limits
    of BMA model confidence intervals of other model averaging methods, etc.
    
    Arguments:
        method -- model averaging method (e.g., 'bma', 'mma', 'mma-s', etc.)
        D -- matrix of ensemble forecasts (N_meas x K)
        y -- observed data
        beta -- parameter estimates
        sigma_2 -- variance estimates
        par -- parameter configurations
        chain -- MCMC chain
        options -- options dictionary
        output -- output dictionary from the model
        str_plot -- list of strings for plotting legends
    """
    
    # Set the default color scheme
    N_meas, K = D.shape

    if K <= 8:
        colororder = np.array([
            [0, 0, 1],    # Blue
            [0, 1, 1],    # Cyan
            [0, 0.5, 0],  # Dark Green
            [1, 0.5, 0.25],  # Orange
            [1, 0, 1],    # Pale Magenta
            [0.6, 0.5, 0.4],  # Dark Brown
            [0, 1, 0],    # Pale Green
            [0, 0, 0]     # Black
        ])
    else:
        clrs = plt.cm.viridis(np.linspace(0, 1, 256))
        colororder = clrs[::int(len(clrs) / K), :3]

    maxbins = 25
    g = 1 - options['alpha']
    gam = np.sort(np.concatenate([ (1 - g) / 2, (1 - g) / 2 + g ]))

    if 'p' not in options:
        options['p'] = np.zeros(K)
    
    # Derive overall min and max for data
    max_Y = max(np.max(y), np.max(D))
    min_Y = min(np.min(y), np.min(D))
    
#    _,n_pars = beta.shape
    if len(beta.shape) == 1:
        # In the case of a 1D array, just set n_pars to the length of beta
        n_pars = beta.shape[0]
    else:
        # In the case of a 2D array, unpack the shape
        _, n_pars = beta.shape    
    
    if method in ['mma', 'mma-s', 'bma']:
        n_pars -= 1  # Remove log-likelihood values

    # Open the output file
    with open('MODELAVG_output.txt', 'w') as fid:
        if method == 'bma':
            fid.write(f"TABLE 1: MODEL AVERAGING METHOD: {method.upper()} WITH {options['PDF'].upper()} CONDITIONAL PDF (MANUAL FOR DETAILS) \n")
        else:
            fid.write(f"TABLE 1: MODEL AVERAGING METHOD: {method.upper()} (MANUAL FOR DETAILS)\n")
        
        fid.write("=" * 38 + "   " + "=" * 20 + "\n")
        
        # Printing results per model
        if method in ['mma', 'mma-s', 'bma', 'gra']:
            fid.write('MODEL  PARAMETER     OPT       STD*       COMPLEXITY   RMSE\n')
        else:
            fid.write('MODEL  PARAMETER     OPT       STD        COMPLEXITY   RMSE\n')
        
        fid.write("-" * 38 + "   " + "-" * 20 + "\n")
        
        # Loop over models to print
        fmt_1 = '%3d \t %-7s %8.3f  %8.3f                 %8.3f\n'
        fmt_2 = '%3d \t %-7s %8.3f  %8.3f      %7u    %8.3f\n'
        fmt_3 = '%3d \t %-7s %8.3f                           %8.3f\n'
        fmt_4 = '%3d \t %-7s %8.3f                %7u    %8.3f\n'

        # Added by JAV to avoid error with printing as str_table and ML are a set and not a list
        if isinstance(output['str_table'], set):
            output['str_table'] = list(output['str_table'])  # Convert to list
        if isinstance(output['ML'], set):
            output['ML'] = list(output['ML'])  # Convert to list
        # CHeck whether end result is still OK; order of weights

        if isinstance(sigma_2, set):
            sigma_2 = list(sigma_2)  # Convert to list

        for i in range(K):
            if method in ['mma', 'mma-s', 'bma', 'gra']:
                if options['p'][i] == 0:
                    fid.write(fmt_1 % (i+1, str(output['str_table'][i]), safe_float(output['ML'][i]), safe_float(output['std'][i]), np.sqrt(sigma_2[i])))
                else:
                    fid.write(fmt_2 % (i+1, str(output['str_table'][i]), safe_float(output['ML'][i]), safe_float(output['std'][i]), safe_float(options['p'][i]), np.sqrt(sigma_2[i])))
            else:
                if options['p'][i] == 0:
                    fid.write(fmt_3 % (i+1, str(output['str_table'][i]), safe_float(output['ML'][i]), np.sqrt(sigma_2[i])))
                else:
                    fid.write(fmt_4 % (i+1, str(output['str_table'][i]), safe_float(output['ML'][i]), safe_float(options['p'][i]), np.sqrt(sigma_2[i])))

        # Table for performance (RMSE and R)
        fid.write("=" * 38 + "\n")
        fid.write("TABLE 2: AVERAGED MODEL PERFORMANCE: TRAINING DATA SET\n")
        fid.write("=" * 38 + "\n")
        fid.write(" METHOD       RMSE      R\n")
        fid.write("-" * 38 + "\n")
        fid.write(f" {method.upper()} \t  {output['RMSE']:.3f}  {output['R']:.3f}\n")
        fid.write("=" * 38 + "\n")
        
        # BMA specific statistics
        if method == 'bma':
            alf_prct = -100 * np.sort(-np.array(options['alpha']))
            _,na = alf_prct.shape
            fid.write(f"TABLE 3: BMA MODEL STATISTICS: TRAINING DATA SET\n")
            fid.write("=" * 61 + "\n")
            fid.write(f" BMA CONDITIONAL PDF      Coverage   Spread   Log-likelihood  \n")
            fid.write("-" * 61 + "\n")
            fmt_1 = '%20s \t                       \t %8.3f\n'
            fmt_2 = '                  %5f%% \t   %5.3f    %5.3f      \n'
            if options['PDF'] == 'gamma':
                fid.write(fmt_1 % ('Gamma Distribution', safe_float(output['loglik'])))
            elif options['PDF'] == 'normal':
                fid.write(fmt_1 % ('Normal Distribution', safe_float(output['loglik'])))
            elif options['PDF'] == 'gen_normal':
                fid.write(fmt_1 % ('Generalized Normal', safe_float(output['loglik'])))
            for i in range(na):
                fid.write(fmt_2 % (safe_float(alf_prct[0,i]), safe_float(output['coverage'][0,i]), safe_float(output['spread'][0,i])))
            fid.write("=" * 61 + "\n")

        # If method is MMA or MMA-S, output specific information
        if method in ['mma', 'mma-s']:
            fid.write(f"TABLE 3: {method.upper()} MODEL STATISTICS: TRAINING DATA SET\n")
            fid.write("=" * 29 + "\n")
            fid.write(f" {method.upper()}   Log-likelihood  \n")
            fid.write("-" * 29 + "\n")
            fmt_1 = '          %8.3f\n'
            fid.write(fmt_1 % safe_float(output['loglik']))
            fid.write("=" * 29 + "\n")
        
        # Correlation coefficients for posterior samples
        if method in ['bma', 'mma', 'mma-s']:
            fid.write("TABLE 4: CORRELATION COEFFICIENTS OF POSTERIOR SAMPLES\n")
            prt = "=========="
            prt = prt * (n_pars - 1) + "========\n"
            fid.write(prt)
            fid.write(f'        {output["str_table"][0]:<10}')
            for i in range(1, n_pars-1):
                fid.write(f'{output["str_table"][i]:<10}')
            fid.write(f'{output["str_table"][n_pars-1]}\n')
            fmt_1 = '%-8s'
            for i in range(n_pars):
                fmt_1 += ' %9.3f'
            fmt_1 += '\n'
            for i in range(n_pars):
                fid.write(fmt_1 % (output['str_table'][i], *output['CORR'][i, :]))
            fid.write(prt)
    
        fid.write("------------------------ end MODELAVG output file --------------------\n")

    # Now plot figures if user wants screen output
    if options['print'] == 'no':
        return
    else:
        pass    # Continue with figure output # plt.ion()

    # ----------------------------------------------------------------------- %
    # Now plot empty figure for PDF file
    # ----------------------------------------------------------------------- %

    # Print wait statement to the screen
    print("\n MODELAVG plotting: Please wait ...")

    # Define name of program
    n_program = 'MODELAVG'
    
    # Define name of figures file
    file_name = f'{n_program}_figures.pdf'
    
    # Determine screen size (using matplotlib to get screen dimensions)
    monitor = get_monitors()[0]
    screen_width = monitor.width
    screen_height = monitor.height
    x_mult = screen_width / 1920
    y_mult = screen_height / 1080
    t_mult = min(x_mult, y_mult)

    # Define fontsize for figures
    fontsize_xylabel = 16 * t_mult
    fontsize_axis = 16 * t_mult
    fontsize_legend = 14 * t_mult
    fontsize_text = 14 * t_mult
    fontsize_title = 18 * t_mult
    fontsize_titlepage = 20 * t_mult
    markersize_symbols = 5

    with PdfPages(file_name) as pdf:

        ## --------------------------------------
        ## Plot Empty Figure for PDF
        ## --------------------------------------
        plt.figure(figsize=(12, 6))
        plt.plot([], [], 'ro')  # Empty plot
        plt.axis([0, 1, 0, 1])
        plt.gca().set_facecolor('w')
        plt.gca().set_xticklabels([])
        plt.gca().set_yticklabels([])        
        plt.gca().set_xticks([])
        plt.gca().set_yticks([])
        plt.text(0.26 * x_mult, 0.6 * y_mult, r'${\rm Visual \; results \; of \; MODELAVG \; toolbox}$', fontsize = fontsize_titlepage) #, ha = 'center', va = 'center')
        plt.text(0.27 * x_mult, 0.5 * y_mult, r'$\;\;\;\;\;{\rm Tables \; are \; not \; printed \; to \; PDF \; file}$', fontsize = fontsize_titlepage) #, ha = 'center', va = 'center') #, fontweight='bold')
        ax = plt.gca()  # Get current axis
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        pdf.savefig()    
        plt.show()

        # PLOT univariate scale reduction factor
        if 'R_stat' in output:
            fig, ax = plt.subplots(figsize=(10, 6))
            for i in range(n_pars):
                ax.plot(output['R_stat'][:, 0], output['R_stat'][:, i + 1], label=str_plot[i], linewidth=1.5)
        
            ax.set_xlabel('Number of generations', fontsize=fontsize_xylabel, labelpad = 10)
            ax.set_ylabel('Gelman-Rubin diagnostic', fontsize=fontsize_xylabel, labelpad = 10)
            ax.set_title(r"$\text{DREAM}_{\rm (ZS)}: \text{ Evolution of } \hat{R} \text{ convergence diagnostic for }$" + f'{method.upper()}' + r"$\text{ method}$", fontsize=fontsize_title)
            ax.axhline(y=1.2, color='r', linestyle='--', linewidth=2)
 
            x_max = output['R_stat'][-1, 0]
            # Check if the value is valid
            if np.isnan(x_max) or np.isinf(x_max):
                print("Invalid value for axis limit:", x_max)
                x_max = 10000  # Set to a default valid value (e.g., 1)

            ax.axis([1, x_max, 0.8, 5])
            ax.legend(fontsize=fontsize_legend)
            pdf.savefig()
            plt.show # (block=False)
#            plt.close()

        # PLOT multivariate scale reduction factor
        if 'MR_stat' in output:
            # Create a new figure for the MR_stat plot
            fig, ax = plt.subplots(figsize=(10, 6))
            # Plot the Gelman and Rubin R-statistic for each parameter
            ax.semilogy(output['MR_stat'][:, 0], output['MR_stat'][:, 1], 'r', linewidth=1.5)
            # Add labels
            ax.set_xlabel('Number of generations', fontsize=fontsize_xylabel, labelpad = 10)
            ax.set_ylabel(r'$\hat{R}^{d} - \text{convergence diagnostic}$', fontsize=fontsize_xylabel, labelpad = 10)
            # Add title
            ax.set_title(r"$\text{DREAM}_{\rm (ZS)}: \text{ Evolution of } \hat{R}^{d} \text{ convergence diagnostic for }$" + f'{method.upper()}' + r"$\text{ method}$", fontsize=fontsize_title)            
            # Plot the theoretical convergence value of 1.2 as a horizontal line
            x_max = output['MR_stat'][-1, 0]
            # Check if the value is valid
            if np.isnan(x_max) or np.isinf(x_max):
                print("Invalid value for axis limit:", x_max)
                x_max = 10000  # Set to a default valid value (e.g., 1)
            ax.plot([1, x_max], [1.2, 1.2], '--', linewidth=2, color=[0.5, 0.5, 0.5])
            # Set axis limits
            ax.axis([0, x_max, 0.8, 20])
            pdf.savefig()
            plt.show # (block=False)
#            plt.close()

        # PLOT Acceptance Rate 
        if 'AR' in output:
            # Create a new figure for the AR plot
            fig, ax = plt.subplots(figsize=(10, 6))
            # Plot the acceptance rate
            ax.plot(output['AR'][:, 0], output['AR'][:, 1], 'r', linewidth=1.5)
            # Add labels
            ax.set_xlabel('Number of generations', fontsize=fontsize_xylabel, labelpad = 10)
            ax.set_ylabel('Acceptance rate', fontsize=fontsize_xylabel, labelpad = 10)
            # Add title
            ax.set_title(r"$\text{DREAM}_{\rm (ZS)}: \text{ Evolution of acceptance rate for }$" + f'{method.upper()}' + r"$\text{ method}$", fontsize=fontsize_title)
            # Calculate maximum y-axis value for the plot
            y_max = min(100, 1.2 * np.max(output['AR'][min(5, output['AR'].shape[0]):, 1]))
            # Set axis limits
            x_max = output['AR'][-1, 0]
            # Check if the value is valid
            if np.isnan(x_max) or np.isinf(x_max):
                print("Invalid value for axis limit:", x_max)
                x_max = 10000  # Set to a default valid value (e.g., 1)
            ax.axis([1, x_max, 0, y_max])
            pdf.savefig()
            plt.show # (block=False)
#            plt.close()
        
        # PLOT histograms of marginal densities of the parameters
        if method in ['bma', 'mma', 'mma-s']:
            # Create a new figure for the histograms
            row, col = 2, 4
            idx_y_label = np.arange(0, n_pars + 1 , col)
            # Calculate the number of bins for each parameter
            Nbins = np.array([calcnbins(beta[:, i]) for i in range(n_pars)])
            nbins = min(np.min(Nbins), maxbins)
            nbins = max(5, round(nbins / 2))  # Adjust number of bins
            counter = 0
            while counter <= n_pars-1:
                # Create a row x col grid of subplots
                fig, axs = plt.subplots(row, col, figsize=(12, 8), zorder = 1)
                # Adjust the spacing between subplots
                plt.subplots_adjust(wspace=0.4, hspace=0.4)  # Increase the space horizontally and vertically
                # Loop over each subplot
                for ax in axs.flat:
                    # Compute histogram for parameter j
                    M, edges = np.histogram(beta[:, counter], bins=nbins, density=True)
                    X = 0.5 * (edges[:-1] + edges[1:])
                    ax.bar(X, M / np.max(M), width=(edges[1] - edges[0]), color='gray', alpha=0.7,zorder=2)
                    # Add labels
                    ax.set_xlabel(str_plot[counter], fontsize=fontsize_xylabel, labelpad = 10)
                    # Adjust axis limits (tight)
                    yticks = np.arange(0, 1.02, 0.2)
                    ax.set_yticks(yticks)  # Set x-tick positions
                    # ax.set_yticklabels([str(tick) for tick in yticks])
                    ax.set_yticklabels([str(round(tick,1)) for tick in yticks])
                    ax.set_ylim(0, 1.02)                # Adjust y limits with padding
                    # Add y-label
                    if counter in idx_y_label:
                        ax.set_ylabel('Empirical density', fontsize=fontsize_xylabel, labelpad=10)
                    # Plot the MAP value
                    ax.plot(output['ML'][counter], 1.02, 'bx', markersize=12, markeredgewidth=3, linewidth=3, zorder=3, clip_on=False)
                    counter += 1
                    if counter == n_pars:
                        break

                # Set the title of the figure
                fig.suptitle(r"$\text{DREAM}_{\rm (ZS)}: \text{ Sampled marginal distribution of the }$" + f'{method.upper()}' + r"$\text{ method}$", fontsize=fontsize_title)
                # Optionally adjust spacing to avoid overlap with subplots
                fig.tight_layout(rect=[0, 0, 1, 0.95])
                pdf.savefig()
                plt.show # (block=False)
                #plt.close()

        # Check for the existence of the CORR field in output
        if 'CORR' in output:
            # Plot the correlation matrix
            fig, ax = plt.subplots(figsize=(12, 8))
            cax = ax.imshow(output['CORR'], cmap='viridis', aspect='auto')
            cbar = fig.colorbar(cax, ax=ax, orientation='vertical', pad = 0.01)
            cbar.set_label('r(x,y)', fontsize=fontsize_legend)
            ax.set_xticks(np.arange(n_pars), str_plot,fontsize=fontsize_xylabel)
            ax.set_yticks(np.arange(n_pars), str_plot,fontsize=fontsize_xylabel)
            # Title of the plot
            ax.set_title(f"Map of posterior correlation coefficients of {method.upper()} parameters", fontsize=fontsize_title)
            # Show the plot
            pdf.savefig()
            plt.show # (block=False)
#            plt.close()

        # Correlation plots of the posterior parameter samples
        if 'CORR' in output and n_pars in range(2, 15):
            # Create a matrix plot for the marginal distributions and bivariate scatter plots
            fig, axs = plt.subplots(n_pars, n_pars, figsize=(15, 15))
            # Calculate the number of bins for each parameter
            Nbins = np.array([calcnbins(beta[:, i]) for i in range(n_pars)])
            nbins = min(np.min(Nbins), maxbins)
            nbins = max(5, round(nbins / 2))  # Adjust number of bins
            for i in range(n_pars):
                for j in range(n_pars):
                    if i != j:
                        axs[i, j].scatter(beta[:, j], beta[:, i], color='gray', s=12)
                        # Add a least-squares line for off-diagonal plots
                        # You can use numpy's polyfit for fitting a line if necessary
                        fit = np.polyfit(beta[:, j], beta[:, i], 1)
                        axs[i, j].plot(beta[:, j], np.polyval(fit, beta[:, j]), 'b--', linewidth=1)
                        # adjust the axes
                        axs[i, j].set_xlim([min(beta[:, j]), max(beta[:, j])])
                        axs[i, j].set_ylim([min(beta[:, i]), max(beta[:, i])])
                    if i == j:
                        # make a histogram
                        axs[i, j].hist(beta[:, i], bins=nbins, density=True, alpha=0.6, color='gray', edgecolor='black')
                        axs[i, j].set_xlim([min(beta[:, i]), max(beta[:, i])])

                    # Set custom x-ticks and x-tick labels
                    x_min, x_max = axs[i, j].get_xlim()
                    dx = x_max - x_min
                    xticks = np.array([x_min + 1/12*dx, x_min + 6/12*dx, x_min + 11/12*dx])
                    axs[i, j].set_xticks(xticks)  # Set x-tick positions first, then labels, otherwise warning
                    axs[i, j].set_xticklabels([str(round(tick,2)) for tick in xticks])
                    y_min, y_max = axs[i, j].get_ylim()
                    dy = y_max - y_min
                    yticks = np.array([y_min + 1/12*dy, y_min + 6/12*dy, y_min + 11/12*dy])
                    axs[i, j].set_yticks(yticks)  # Set y-tick positions first, then labels, otherwise warning
                    axs[i, j].set_yticklabels([str(round(tick,2)) for tick in yticks])
                    # Add values and labels to the axes
                    if i == n_pars - 1:
                        axs[i, j].set_xlabel(str_plot[j], fontsize=fontsize_xylabel)
                    else:
                        axs[i, j].set_xticklabels([])

                    # Add values and labels to the axes
                    if j == 0:
                        axs[i, j].set_ylabel(str_plot[i], fontsize=fontsize_xylabel)
                    else:
                        axs[i, j].set_yticklabels([])

            # Title of the figure
            fig.suptitle(f"Marginal distribution and bivariate scatter plots of posterior samples of {method.upper()} method", fontsize=fontsize_title)
            pdf.savefig()
            plt.show # (block=False)
#            plt.close()

        # Plot Traceplots for Parameters
        if method in ['bma', 'mma', 'mma-s']:
            T, _, N = chain.shape
            symbols = ['cs', 'rx', 'g+', 'ko', 'm<']
            leg_string = ['chain 1']
            for i in range(1, min(5, N)):
                leg_string.append(f'chain {i + 1}')
            leg_string.append('ML')

            counter = 0
            while counter <= n_pars - 1:
                # Create a row x col grid of subplots
                fig, axs = plt.subplots(2, 1, figsize=(12, 8), zorder=1)
                # Adjust the spacing between subplots
                plt.subplots_adjust(wspace=0.4, hspace=0.4)  # Increase the space horizontally and vertically
                # Loop over subplot
                for ax in axs.flat:
                    for i in range(min(N, 5)):
                        ax.plot(np.arange(1, T + 1), chain[:, counter, i], symbols[i], markersize=markersize_symbols, linewidth=3, zorder=2)
                    # Adjust axis limits (tight)
                    ax.set_xlim(1, T)                                                          # Set x-axis limits
                    y_min, y_max = ax.get_ylim()
                    ax.set_ylim(min(0.9 * y_min, 1.1 * y_min), max(0.9 * y_max, 1.1 * y_max))  # Adjust y limits with padding
                    # Set custom x-ticks and x-tick labels: 1, T/5, 2*T/5, ... , T
                    xticks = np.concatenate([[1], np.arange(T//5, T+1, T//5)])
                    ax.set_xticks(xticks)  # Set x-tick positions
                    ax.set_xticklabels([str(int(tick)) for tick in xticks])  # Set x-tick labels
                    # Plot the MAP value
                    ax.plot(T,output['ML'][counter], 'bx', markersize=12, markeredgewidth=3, linewidth=5, zorder=3, clip_on=False)
                    ax.set_xlabel('Sample number of chain', fontsize=fontsize_xylabel, labelpad = 10)
                    ax.set_ylabel(f'Parameter {str_plot[counter]}', fontsize=fontsize_xylabel, labelpad = 10)
                    ax.set_title(f"DREAM"r"$_\text{(ZS)}\;$" f"sampled traceplot of {method.upper()} parameter {str_plot[counter]} ", fontsize=fontsize_title)
                    ax.legend(leg_string, loc='best', fontsize=fontsize_legend)
                    ax.tick_params(axis='both', labelsize=fontsize_axis)
                    counter += 1                
                    if counter == n_pars:
                        break
                #plt.tight_layout()
                pdf.savefig()                
                plt.show # (block=False)

            # Plot the log-likelihood function
            fig, ax = plt.subplots(1, 1, figsize=(12, 8))
            for i in range(min(N, 5)):
                ax.plot(np.arange(1, T + 1), chain[:, n_pars, i], symbols[i], markersize=markersize_symbols, linewidth=3)

            # Adjust axis limits (tight)
            ax.set_xlim(1, T)                                                          # Set x-axis limits
            y_min, y_max = ax.get_ylim()
            ax.set_ylim(min(0.9 * y_min, 1.1 * y_min), max(0.9 * y_max, 1.1 * y_max))  # Adjust y limits with padding
            # Set custom x-ticks and x-tick labels: 1, T/5, 2*T/5, ... , T
            xticks = np.concatenate([[1], np.arange(T//5, T+1, T//5)])
            ax.set_xticks(xticks)  # Set x-tick positions
            ax.set_xticklabels([str(int(tick)) for tick in xticks])  # Set x-tick labels
            ax.set_xlabel('Sample number of chain', fontsize=fontsize_xylabel, labelpad = 10)
            ax.set_ylabel(r'$\mathcal{L}(\mathbf{x} \vert \widetilde{\mathbf{y}})$', fontsize=18, loc='center', labelpad = 10)
            ax.legend(leg_string[0:-1], loc='best', fontsize=fontsize_legend)
            ax.set_title(f"DREAM"r"$_\text{(ZS)}\;$" f"sampled traceplot of log-likelihood using the {method.upper()} method", fontsize=fontsize_title)
            pdf.savefig()
            plt.show # (block=False)
            # plt.tight_layout()
            # plt.close()

        # PLOT Data, Ensemble and Mean Forecast
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))  # Equivalent to 'units', 'normalized', 'outerposition'
        # Plot the ensemble
        HA = ax.plot(D, linewidth=1.35)  # Plot the ensemble
        if K < len(colororder):  # Check if K is within the bounds of colororder
            for i in range(K):
                HA[i].set_color(colororder[i][:3])  # Set the color for each model
        # Plot verifying observations
        ax.plot(np.arange(1, N_meas+1), y, 'ro', markersize=6, linewidth=1, markerfacecolor='w')
        # Plot the averaged prediction (mean forecast)
        ax.plot(np.arange(1, N_meas + 1),output['Ye'], 'k', linewidth=1.5)
        # Define the axis limits
        delta_Y = 0.07 * (max_Y - min_Y)
        min_X = 1
        max_X = min(500, N_meas)
        # Adjust axis limits
        ax.axis([min_X, max_X, min_Y, max_Y + delta_Y])
        # Add labels
        ax.set_xlabel(f'Row (index) of training data set', fontsize=fontsize_xylabel, verticalalignment='top', horizontalalignment='center', labelpad=10)
        ax.set_ylabel(f'Forecast', fontsize=fontsize_xylabel, verticalalignment='top', horizontalalignment='center', labelpad=20)
        # Add custom legend
        x_loc = min_X + 0.05 * (max_X - min_X)
        y_loc = min_Y + 1 * (max_Y - min_Y)
        dx = 0.03 * (max_X - min_X)
        dy = 0.04 * (max_Y - min_Y)
        x1 = x_loc
        x2 = x1 + dx
        # Plot model entries
        for k in range(K):
            ax.plot([x1, x2], [y_loc, y_loc], color=colororder[k][:3], linewidth=4)
            ax.text(x1 + 1.4 * dx, y_loc, f'Model {k+1}', fontsize=fontsize_legend, color=colororder[k][:3], verticalalignment='center')
            y_loc -= dy
        # Plot the averaged prediction line
        ax.plot([x1, x2], [y_loc, y_loc], color='k', linewidth=4)
        ax.text(x1 + 1.4 * dx, y_loc, r'$g_{j}^{\bullet}$: Weighted-average prediction', fontsize=fontsize_legend, verticalalignment='center')
        # Plot the verifying observations marker
        y_loc -= dy
        ax.plot(x1 + dx / 2, y_loc, 'ro', markersize=8, linewidth=2, markerfacecolor='w')
        ax.text(x1 + 1.4 * dx, y_loc, 'Verifying observations', fontsize=fontsize_legend, color='r', verticalalignment='center')
        # Create title
        fig.suptitle(f'Snapshot of model ensemble (colored lines), {method.upper()} mean forecast (black line) and verifying data (red circles)', fontsize=fontsize_title, horizontalalignment='center')
        pdf.savefig()
        # Show the plot
        plt.show # (block=False)

        # PLOT confidence and prediction intervals
        _,nG = gam.shape      
        clr = np.full((nG, 3), np.nan)  # Initialize the color array with NaNs
        # Open new figure with normalized size and outer position
        fig, ax = plt.subplots(1,1, figsize=(16, 9), dpi=80)
        # Start with the gamma% total uncertainty: first entry of fill_between is larger of the two
        for i in range(nG):
            clr[i, :] = np.ones(3) - (i+1) * 1 / (2 * nG)  # Create the color gradient
            # Fill the ranges (replace 'Fill_Ranges' with fill_between)
            y1 = output['pred'][:, 2 * nG - i - 1]; y2 = output['pred'][:, i]
            ax.fill_between(np.arange(1, N_meas+1), y1, y2, where=(y1 > y2), color=clr[i, :], alpha = 0.7, interpolate=True)
        
        # Add gamma% uncertainty due to parameter uncertainty (dark gray)
        y1 = output['par_unc'][:, -1]; y2 = output['par_unc'][:, 0]
        ax.fill_between(np.arange(1, N_meas + 1), y1, y2, where=(y1 > y2), color=[0.2, 0.2, 0.2], alpha = 0.7, interpolate=True)
        # Plot the verifying observations (red circles)
        ax.plot(np.arange(1, N_meas + 1), y, 'ro', markersize=6, linewidth=1, markerfacecolor='w')
        # Plot the averaged prediction (black line)
        ax.plot(np.arange(1, N_meas + 1),output['Ye'], 'k', linewidth=1.5)
        # Define axis limits
        delta_Y = 0.07 * (max_Y - min_Y)
        min_X = 1
        max_X = min(500, N_meas)
        ax.axis([min_X, max_X, min_Y, max_Y + delta_Y])
        # Add labels
        ax.set_xlabel('Row (index) of training data set', fontsize=fontsize_xylabel, labelpad=10)
        ax.set_ylabel('Forecast', fontsize=fontsize_xylabel, labelpad=10)
        # Create title
        fig.suptitle(f"Snapshot of confidence (dark gray) and prediction (light gray) intervals: {method.upper()} method", fontsize=fontsize_title)
        # Define legend location
        x_loc = min_X + 0.05 * (max_X - min_X)
        y_loc = min_Y + 1 * (max_Y - min_Y)
        dx = 0.03 * (max_X - min_X)
        dy = 0.07 * (max_Y - min_Y)
        # Add axes to the figure manually (with [left, bottom, width, height] in normalized figure coordinates)
        values = [int(value) for value in 100*g[0,:]]
        pred_valstr = '/'.join(map(str, values)) + '%'
        pred_str = 'Prediction interval'
        conf_valstr = str(values[0]) + '%'
        conf_str = 'Confindence interval'
        # First entry (Prediction interval)
        for i in range(nG):
            x1 = x_loc + i * dx / nG
            x2 = x1 + dx / nG
            # Add a rectangle patch to highlight a specific region
            rect = patches.Rectangle((x1, y_loc), dx / nG, dy / 2, linewidth=2, edgecolor='k', facecolor=clr[i, :3])
            ax.add_patch(rect)

        # Add text for the prediction interval
        ax.text(x_loc + 1.4 * dx, y_loc + dy / 5, f"{pred_valstr} {pred_str}", fontsize=fontsize_text, color='k', verticalalignment='center', horizontalalignment='left')
        # Update the location for the next legend entry
        y_loc -= dy / 1.3
        # Second entry (Confidence interval)
        rect = patches.Rectangle((x_loc, y_loc), dx, dy/2 , linewidth=2, edgecolor='k', facecolor=[0.2, 0.2, 0.2])
        ax.add_patch(rect)
        ax.text(x_loc + 1.4 * dx, y_loc + dy / 5, f"{conf_valstr} {conf_str}", fontsize=fontsize_text, color=[0.2, 0.2, 0.2], verticalalignment='center', horizontalalignment='left')
        # Third entry (Weighted average prediction)
        y_loc -= dy / 1.7
        ax.plot([x_loc, x_loc + dx], [y_loc, y_loc], color='k', linewidth=4)
        ax.text(x_loc + 1.4 * dx, y_loc, r'$g_{j}^{\bullet}$: Weighted-average prediction', fontsize=fontsize_text, verticalalignment='center', horizontalalignment='left')
        # Last entry (Verifying observations)
        y_loc -= dy / 1.3
        ax.plot(x_loc + dx / 2, y_loc, 'ro', markersize=8, linewidth=3, markerfacecolor='w')
        ax.text(x_loc + 1.4 * dx, y_loc, 'Verifying observations', fontsize=fontsize_text, color='r', verticalalignment='center', horizontalalignment='left')
        # Add label with % contained (Coverage of the prediction intervals)
        x_loc = min_X + 0.6 * (max_X - min_X)
        y_loc = min_Y + 1.05 * (max_Y - min_Y)
        dx = 0.02 * (max_X - min_X)
        dy = 0.04 * (max_Y - min_Y)
        # Make a white patch for the coverage information
        x1 = x_loc - dx
        x2 = x_loc + 16 * dx
        rect = patches.Rectangle((x_loc, y_loc - (nG + 1) * dy), 15*dx, (nG + 1) * dy , linewidth=1, edgecolor='k', facecolor='w')
        ax.add_patch(rect)
        # Update location for the coverage information
        y_loc = min_Y + 1.01 * (max_Y - min_Y)
        # Display coverage information for each prediction interval
        for i in range(nG):
            coverage_value = int(100 * g[0, nG - 1 - i])  # Get the g value and convert it to a percentage
            coverage_percent = round(100 * output['coverage'][0, nG - 1 - i]) / 100  # Get the coverage percentage
            ax.text(x_loc + dx/2, y_loc - (i) * dy, f'Coverage of {coverage_value}% pred interval is {coverage_percent}%', fontsize=fontsize_text, color=clr[i, :3], verticalalignment='center', horizontalalignment='left')
        ax.tick_params(axis='both', labelsize=fontsize_axis)
        plt.tight_layout()
        pdf.savefig()
        plt.show # (block=False)

        # Plot: Residual Analysis: ACF
        res = output['Ye'] - y
        fig, ax = plt.subplots(1, 1, figsize=(12, 8))
        plot_acf(res, lags=40, alpha=0.05, ax=ax)      
        ax.set_title(f'Autocorrelation function of residuals of deterministic (mean) forecast of {method.upper()} method', fontsize=fontsize_title)
        ax.set_xlabel('Lag', fontsize=fontsize_xylabel,labelpad = 10)
        ax.set_ylabel('Autocorrelation', fontsize=fontsize_xylabel,labelpad = 10)
        plt.tight_layout()
        pdf.savefig()
        plt.show # (block=False)
        #plt.close()
    
        # QQ Plot
        fig, ax = plt.subplots(1, 1, figsize=(8, 8))
        stats.probplot(res, dist="norm", plot=ax)
        ax.set_title(f'Quantile-quantile plot of residuals of deterministic (mean) forecast of {method.upper()} method', fontsize=fontsize_title)
        ax.set_xlabel('Standard normal quantiles', fontsize=fontsize_xylabel,labelpad = 10)
        ax.set_ylabel('Quantiles of point forecast', fontsize=fontsize_xylabel,labelpad = 10)
        plt.tight_layout()
        pdf.savefig()
        plt.show # (block=False)
        #plt.close()

    # Open the final PDFs
    os.startfile(file_name)


def Bias_correction(D, y, intcpt=1):
    """
    This function provides a linear correction of the ensemble members.

    SYNOPSIS: [D_bc, a, b] = Bias_correction(D, y)
              [D_bc, a, b] = Bias_correction(D, y, intcpt)
    where:
    D         [input] nxK matrix with forecasts of ensemble members
    y         [input] nx1 vector with verifying observations
    intcpt    [input] optional: intercept (default: 1) or without (0)
    D_bc      [output] nxK matrix with bias-corrected forecasts
    a         [output] 1xK vector with intercept bias-corrected forecasts
    b         [output] 1xK vector with slope bias-corrected forecasts
    """
    
#    D = np.array(D)  # Ensure that D is a NumPy array
#    y = np.array(y)  # Ensure y is a NumPy array
    # Get the shape of D
    n, K = D.shape
    D_bc = np.full((n, K), np.nan)  # Initialize with NaNs
    
    # Initialize intercepts and slopes of linear bias correction functions
    a = np.zeros(K)
    b = np.zeros(K)
    
    # Perform linear regression for each ensemble member
    for k in range(K):
        if intcpt == 1:
            # Linear regression with intercept
            X = np.column_stack((np.ones(n), D[:, k]))  # Add intercept column
            ab = np.linalg.lstsq(X, y, rcond=None)[0]  # Solve for intercept and slope
            a[k], b[k] = ab
        else:
            # Linear regression without intercept
            b[k] = np.linalg.lstsq(D[:, k].reshape(-1, 1), y, rcond=None)[0][0]
        
        # Bias-corrected ensemble forecasts
        D_bc[:, k] = a[k] + b[k] * D[:, k]
    
    return D_bc, a, b


# Note the following two functions may hickup at sigma[sigma <= 0] = np.nan if sigma is not np.float64 (or np.float32)
def gevstat(k, sigma, mu):
    """
    GEVSTAT Mean and variance of the generalized extreme value distribution.
    Returns the mean (m) and variance (v) of the generalized extreme value (GEV) distribution 
    with shape parameter k, scale parameter sigma, and location parameter mu.
    
    Args:
        k (float or np.ndarray): Shape parameter (K) of the GEV distribution.
        sigma (float or np.ndarray): Scale parameter (SIGMA) of the GEV distribution.
        mu (float or np.ndarray): Location parameter (MU) of the GEV distribution.
    
    Returns:
        m (np.ndarray): Mean of the GEV distribution.
        v (np.ndarray): Variance of the GEV distribution.
    """

    # Handle scalar inputs by broadcasting them to the output size
    if np.ndim(k) == 0:
        k = np.full_like(sigma, k)
    if np.ndim(sigma) == 0:
        sigma = np.full_like(k, sigma)
    if np.ndim(mu) == 0:
        mu = np.full_like(k, mu)

    # Prepare arrays for storing the results
    m = np.full_like(k, np.nan, dtype=np.float64)
    v = np.full_like(k, np.nan, dtype=np.float64)

    # Avoid errors by returning NaN for out of range parameters
    sigma[sigma <= 0] = np.nan

    # Compute gamma values
    lngam1 = np.full_like(k, np.inf, dtype=np.float64)
    lngam2 = np.full_like(k, np.inf, dtype=np.float64)

    # Calculate gamma values for different ranges of k
    lngam1[k < 1] = gammaln(1 - k[k < 1])
    lngam2[k < 0.5] = gammaln(1 - 2 * k[k < 0.5])

    # Compute mean and variance
    jm = np.abs(k) < 1e-8
    m[jm] = -psi(1)  # Euler's constant
    jv = np.abs(k) < 5e-6
    v[jv] = pi**2 / 6  # psi(1, 1)

    # For k != 0 and k < 1, compute mean
    jm = ~jm
    jj = jm & (k < 1)
    m[jj] = np.expm1(lngam1[jj]) / k[jj]  # (gamma(1 - k) - 1) / k
    m[k >= 1] = np.inf

    # For k != 0 and k < 1/2, compute variance
    jv = ~jv
    jj = jv & (k < 1/2)
    v[jj] = (np.expm1(lngam2[jj]) - np.expm1(2 * lngam1[jj])) / (k[jj]**2)

    v[k >= 1/2] = np.inf

    # Scale the mean and variance with sigma and add the location mu
    m = mu + sigma * m
    v = sigma**2 * v

    return m, v


def gpstat(k, sigma, theta=0):
    """
    GPSTAT Mean and variance of the generalized Pareto distribution.
    Returns the mean (m) and variance (v) of the generalized Pareto (GP) distribution 
    with tail index (shape) parameter k, scale parameter sigma, and threshold (location) parameter theta.
    
    Args:
        k (np.ndarray or float): Tail index (shape) parameter of the GP distribution.
        sigma (np.ndarray or float): Scale parameter of the GP distribution.
        theta (np.ndarray or float, optional): Threshold (location) parameter of the GP distribution. Default is 0.
    
    Returns:
        m (np.ndarray): Mean of the GP distribution.
        v (np.ndarray): Variance of the GP distribution.
    """
    
    # Handle scalar inputs by broadcasting them to the output size
    if np.ndim(k) == 0:
        k = np.full_like(sigma, k)
    if np.ndim(sigma) == 0:
        sigma = np.full_like(k, sigma)
    if np.ndim(theta) == 0:
        theta = np.full_like(k, theta)

    # Initialize m and v with NaNs
    m = np.full_like(k, np.nan, dtype=np.float64)
    v = np.full_like(k, np.nan, dtype=np.float64)

    # Return NaN for out of range parameters
    sigma[sigma <= 0] = np.nan

    # Handle the k == 0 case
    j = np.abs(k) < np.finfo(float).eps
    m[j] = 1
    v[j] = 1

    # For k != 0 and k < 1, compute the mean
    j = ~j
    jj = j & (k < 1)
    m[jj] = 1 / (1 - k[jj])  # Mean formula for k < 1
    m[k >= 1] = np.inf  # Set mean to infinity when k >= 1

    # For k != 0 and k < 1/2, compute the variance
    jj = j & (k < 0.5)
    v[jj] = 1 / ((1 - k[jj])**2 * (1 - 2 * k[jj]))  # Variance formula for k < 1/2
    v[k >= 1/2] = np.inf

    # Scale the mean and variance with sigma and add the location mu
    m = theta + sigma * m
    v = sigma**2 * v

    return m, v


def gnormcdf(x, mu, alfa, beta):
    # Generalized normal CDF
    P = 0.5 + np.sign(x - mu) * (1 / (2 * scipy.special.gamma(1 / beta))) * gammainc(1 / beta, np.abs((x - mu) / alfa) ** beta) * scipy.special.gamma(1 / beta)
    return P


def gnormpdf(x, mu, alfa, beta):
    """
    Generalized normal probability density function (pdf).
    
    Parameters:
    - x : array_like
        Values at which the pdf should be evaluated.
    - mu : float
        Mean of the distribution.
    - alfa : float
        Scale parameter (standard deviation).
    - beta : float
        Shape (kurtosis) parameter.

    Returns:
    - Y : array_like
        The values of the generalized normal pdf evaluated at x.
    """
    # Calculate the generalized normal pdf
    Y = beta / (2 * alfa * scipy.special.gamma(1 / beta)) * np.exp(- (np.abs(x - mu) / alfa) ** beta)
    
    return Y


def gnormrnd(mu, alfa, beta, *args):
    """
    Generate random variables from the Generalized Normal Distribution.

    Parameters:
    - mu : Mean of the distribution
    - alfa : Standard deviation (scale parameter)
    - beta : Shape parameter that controls the kurtosis
    - *args : Additional arguments to specify the shape of the output array (e.g., M, N, ...)
    
    Returns:
    - R : Randomly generated numbers from the generalized normal distribution
    """
    # Check if the number of arguments is valid
    if len(args) < 1:
        raise ValueError('Too few inputs, dimensions of output must be provided')
    
    # Determine the shape of the output array
    sizeOut = args

    # Generate random uniform variables p
    p = np.random.rand(*sizeOut)

    # Generate random samples from the generalized normal distribution
    R = mu + np.sign(p - 0.5) * (alfa ** beta * gamma.ppf(2 * np.abs(p - 0.5), 1 / beta)) ** (1 / beta)
    # gamma.ppf works fine - at least when separately evaluated

    return R


def compKGE(Yobs, Ysim):
    """
    This function computes the Kling-Gupta efficiency (KGE)

    Parameters:
    Yobs (array-like): Observed values
    Ysim (array-like): Simulated values

    Returns:
    KGE (float): Kling-Gupta efficiency
    """
    method = 1

    if method == 1:  # Alternative form II
        m_X = np.mean(Ysim)
        s_X = np.std(Ysim,ddof=1)
        m_Y = np.mean(Yobs)
        s_Y = np.std(Yobs,ddof=1)
        alfa = m_X / m_Y
        beta = s_X / s_Y
        m_XY = np.mean(Ysim * Yobs)
        r = (m_XY - m_X * m_Y) / (s_X * s_Y)
        term = 3 + alfa**2 + beta**2 + r**2 - 2 * (alfa + beta + r)
        KGE = 1 - np.sqrt(term)
    
    elif method == 2:  # Alternative form I
        m_X = np.mean(Ysim)
        m2_X = m_X**2
        s_X = np.std(Ysim,ddof=1)
        s2_X = np.var(Ysim,ddof=1)
        m_Y = np.mean(Yobs)
        m2_Y = m_Y**2
        s_Y = np.std(Yobs,ddof=1)
        s2_Y = np.var(Yobs,ddof=1)
        m_XY = np.mean(Ysim * Yobs)
        m2_XY = m_XY**2
        term = (3 + m2_X/m2_Y - 2*m_X/m_Y + s2_X/s2_Y - 2*s_X/s_Y +
                (m2_XY + m2_X*m2_Y - 2*m_X*m_Y*m_XY) / (s2_X*s2_Y) -
                2*(m_XY - m_X*m_Y) / (s_X*s_Y))
        KGE = 1 - np.sqrt(term)
    
    else:  # Old form
        a = Yobs
        b = Ysim
        s_meas = np.std(a,ddof=1)
        m_meas = np.mean(a)
        s_sim = np.std(b,ddof=1)
        m_sim = np.mean(b)
        alpha = s_sim / s_meas
        gamma = m_sim / m_meas
        r = (np.mean(Ysim * Yobs) - np.mean(Ysim) * np.mean(Yobs)) / (np.std(Ysim) * np.std(Yobs))
        KGE = 1 - np.sqrt((r - 1)**2 + (alpha - 1)**2 + (gamma - 1)**2)

    return KGE


def MMA_lik(beta, D, Y, var_err, p):
    """
    This function calculates the log likelihood of MMA.
    Function of MODELAVG toolbox, V1.0
    
    B.C. Hansen, "Least Squares Model Averaging", Econometrica, vol. 75, 
    no. 4, pp. 1175-1189, 2007
    """
    if len(beta) == 0 or len(D) == 0 or len(Y) == 0:
        raise ValueError("Requires at least five input arguments.")
    
    # MMA deterministic point forecast
    G = np.dot(D, beta)  # Equivalent to D * beta' in MATLAB
    
    # Mallows criterion (Equation (11))
    Cn = np.sum((Y - G)**2) + 2 * var_err * np.dot(beta, p)
    
    # Log-likelihood of x = { MMA weights }
    ell = -0.5 * Cn + np.finfo(float).eps  # realmin equivalent in Python
    
    return ell


def calcnbins(x, method='middle', minb=1, maxb=np.inf):
    """
    Compute the "ideal" number of bins using different methods for histogram bin calculation.
    
    Parameters:
    - x: vector of data points
    - method: string with choice of method for calculating bins. Default is 'middle'.
        Options: 'fd', 'scott', 'sturges', 'middle', 'all'
    - minb: smallest acceptable number of bins (default: 1)
    - maxb: largest acceptable number of bins (default: np.inf)
    
    Returns:
    - nbins: The calculated number of bins based on the selected method.
    """
    
    # Input checking
    if not isinstance(x, (np.ndarray, list, np.generic)):
        raise ValueError('The x argument must be numeric or logical.')

    x = np.asarray(x)
    
    # Ensure the array is real, discard imaginary part
    if np.iscomplexobj(x):
        x = np.real(x)
        print('Warning: Imaginary parts of x will be ignored.')
    
    # Ensure x is a vector (1D array)
    if x.ndim != 1:
        x = x.flatten()
        print('Warning: x will be coerced to a vector.')
    
    # Remove NaN values
    x = x[~np.isnan(x)]
    if len(x) == 0:
        raise ValueError("x must contain at least one valid number.")
    
    # Choose method if not specified
    valid_methods = ['fd', 'scott', 'sturges', 'all', 'middle']
    if method not in valid_methods:
        raise ValueError(f"Unknown method: {method}")
    
    # Method selection
    if method == 'fd':
        nbins = calc_fd(x)
    elif method == 'scott':
        nbins = calc_scott(x)
    elif method == 'sturges':
        nbins = calc_sturges(x)
    elif method == 'middle':
        nbins = [calc_fd(x), calc_scott(x), calc_sturges(x)]
        nbins = np.median(nbins)
    elif method == 'all':
        nbins = {
            'fd': calc_fd(x),
            'scott': calc_scott(x),
            'sturges': calc_sturges(x)
        }
    
    # Confine number of bins to the acceptable range
    nbins = confine_to_range(nbins, minb, maxb)
    
    return nbins


def calc_fd(x):
    """Freedman-Diaconis rule"""
    h = np.subtract(*np.percentile(x, [75, 25]))  # Interquartile range (IQR)
    if h == 0:
        h = 2 * np.median(np.abs(x - np.median(x)))  # Median absolute deviation (MAD)
    
    if h > 0:
        nbins = np.ceil((np.max(x) - np.min(x)) / (2 * h * len(x) ** (-1/3)))
    else:
        nbins = 1
    return nbins


def calc_scott(x):
    """Scott's method"""
    h = 3.5 * np.std(x) * len(x) ** (-1/3)
    if h > 0:
        nbins = np.ceil((np.max(x) - np.min(x)) / h)
    else:
        nbins = 1
    return nbins


def calc_sturges(x):
    """Sturges' method"""
    nbins = np.ceil(np.log2(len(x)) + 1)
    return nbins


def confine_to_range(x, lower, upper):
    """Ensure bin count is within the specified range"""
    x = np.maximum(x, lower)
    x = np.minimum(x, upper)
    return np.floor(x)


def safe_int(value):
    return int(value) if value else 0  # or return a default value if empty


def safe_float(value):
    return float(value) if value else 0.0  # or return a default value if empty