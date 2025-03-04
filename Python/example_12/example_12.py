# ---------------------------------------------------------------------- #
# MM   MM   OOOOO   DDDDDD   EEEEEEE  LL         AAA    VV   VV   GGGGG  #
# MM   MM  OOOOOOO  DDDDDDD  EEEEEEE  LL        AA AA   VV   VV  GG   GG #
# MMM MMM  OO   OO  DD   DD  EE       LL        AA AA   VV   VV  GG   GG #
# MM M MM  OO   OO  DD   DD  EEEEE    LL       AA   AA  VV   VV  GGGGGGG #
# MM   MM  OO   OO  DD   DD  EEEEE    LL       AAAAAAA  VV   VV   GGGGGG #
# MM   MM  OO   OO  DD   DD  EE       LL       AA   AA   VV VV        GG #
# MM   MM  OOOOOOO  DDDDDDD  EEEEEEE  LLLLLLL  AA   AA    VVV     GGGGGG #
# MM   MM   OOOOO   DDDDDD   EEEEEEE  LLLLLLL  AA   AA     V     GGGGGGG # 
# ---------------------------------------------------------------------- #

# About the discharge ensemble of watershed models
#   Vrugt, J.A., and B.A. Robinson (2007), Treatment of uncertainty using 
#        ensemble methods: Comparison of sequential data assimilation and 
#        Bayesian model averaging, Water Resources Research, 43, W01411, 
#        doi:10.1029/2005WR004838.
# About the scoring rules
#   Vrugt, J.A. (2024), Distribution-Based Model Evaluation and 
#        Diagnostics: Elicitability, Propriety, and Scoring Rules for 
#        Hydrograph Functionals, Water Resources Research, 60, 
#        e2023WR036710, https://doi.org/10.1029/2023WR036710  

import sys
import os
import numpy as np
from MODELAVG import MODELAVG

# Get the current working directory
current_directory = os.getcwd()
# Go up one import numpy as np
from joblib import Parallel, delayed


# DEFINE MODEL AVERAGING METHOD
method = 'bma'  # 'ewa'/'bga'/'aica'/'bica'/'gra'/'bma'/'mma'/'mma-s'

# FIXED SETTINGS
options = {
    'alpha': np.array([0.01, 0.05, 0.1, 0.5]),  # significance levels for BMA model
    'print': 'no',                              # print output (figures, tables) to screen
    'CPU': 'yes',                               # efficient CPU implementation
}

# NOW LOAD DATA
S = np.loadtxt('discharge.txt')  # daily discharge forecasts (mm/day) of models and verifying data
idx_cal = np.arange(0, 3000)     # start/end training period
idx_eval = np.arange(5000, 7000) # evaluation period

# DEFINE ENSEMBLE AND VECTOR OF VERIFYING OBSERVATIONS
D = S[idx_cal, 0:8]  # Forecast ensemble (models 1-8)
y_cal = S[idx_cal, 8]  # Verifying observations (column 9)
D_eval = S[idx_eval, 0:8]  # Evaluation ensemble
y_eval = S[idx_eval, 8]  # Verifying evaluation observations

# APPLY LINEAR BIAS CORRECTION TO ENSEMBLE (UP TO USER)
D_cal, a_LB, b_LB = Bias_correction(D, y_cal)

# BMA CONDITIONAL DISTRIBUTIONS TO BE USED
PDF = ['normal', 'lognormal', 'gen_normal', 'tnormal', 'gamma', 'gev', 'gpareto']
V = ['1', '2', '3', '4']  	# Homoscedastic variance ('1' group, '2' single), 
				# Heteroscedastic variance ('3' group, '4' single)

# BMA -> CONDITIONAL DISTRIBUTION NEEDS TO BE DEFINED
options['TAU'] = '2'    # Gen. normal: different τ for each ensemble member
                        # Gen. Pareto: different xi (shape) for each ensemble member
                        # Gen. extreme value: different xi (shape) for each ensemble member

which_period = 'evaluation'  # 'training'

# Prepare to store names for conditional distributions
names_pdf = [
    	'1: normal: const group', 
	'2: lnormal: const group', 
	'3: gnormal: const group',
    	'4: tnormal: const group', 
	'5: gamma: const group', 
	'6: gev: const group', 
	'7: gp: const group',
    	'8: normal: const ind', 
	'9: lnormal: const ind', 
	'10: gnormal: const ind', 
	'11: tnormal: const ind',
    	'12: gamma: const ind', 
	'13: gev: const ind', 
	'14: gp: const ind', 
	'15: normal: nonconst group',
    	'16: lnormal: nonconst group', 
	'17: gnormal: nonconst group', 
	'18: tnormal: nonconst group',
    	'19: gamma: nonconst group', 
	'20: gev: nonconst group', 
	'21: gp: nonconst group',
    	'22: normal: nonconst ind', 
	'23: lnormal: nonconst ind', 
	'24: gnormal: nonconst ind',
    	'25: tnormal: nonconst ind', 
	'26: gamma: nonconst ind', 
	'27: gev: nonconst ind', 
	'28: gp: nonconst ind'
]

# INITIALIZE VARIABLES/OUTPUT ARRAYS
nV = len(V)
nPDF = len(PDF)

# Set the length of the vector based on the selected period
if which_period == 'training':
    nY = len(y_cal)
elif which_period == 'evaluation':
    nY = len(y_eval)

cal = [None] * (nV * nPDF)
val = [None] * (nV * nPDF)

# Prepare combinations of variances and PDFs for parallel processing
VPDF = []
for v in V:
    for p in PDF:
        VPDF.append([v, p])

# Define the parallel script
def parallel_script(ii, VPDF, method, D_cal, y_cal, options, which_period, a_LB, b_LB, D_eval, y_eval):
    # Assign variance and PDF to options
    options['VAR'] = VPDF[0]
    options['PDF'] = VPDF[1]
    
    # Run the BMA method
    beta, output = MODELAVG(method, D_cal, y_cal, options)
    AR = output['AR']
    
    # Store information of calibration period
    table_res_cal = [
        output['mQS'], output['mLS'], output['mSS'], output['mCRPS'],
        output['mES'], output['mRLBL_anal'], output['mCv_anal'],
        output['coverage'][0] / 100, output['spread'][0], output['loglik'],
        output['RMSE'], output['R2'], output['KGE']
    ]
    
    # Print log-likelihood
    print(output['loglik'])
    
    # Compute information for the evaluation period
    if which_period == 'evaluation':
        # Evaluate maximum likelihood BMA forecast distribution for evaluation period
        val, D_eval = MODELAVG_eval(method, D_eval, y_eval, options, a_LB, b_LB, output)
        table_res_eval = [
            val['mQS'], val['mLS'], val['mSS'], val['mCRPS'], val['mES'],
            val['mRLBL_anal'], val['mCv_anal'], val['coverage'][0] / 100, val['spread'][0],
            val['loglik'], val['RMSE'], val['R2'], val['KGE']
        ]
    else:
        val = []
        table_res_eval = ii
    
    return cal, beta, AR, table_res_cal, val, table_res_eval

if __name__ == '__main__':
	# Run BMA method for all combinations of variances/conditional PDFs in parallel
	results = Parallel(n_jobs=-1)(
    		delayed(parallel_script)(ii, VPDF[ii], method, D_cal, y_cal, options, which_period, a_LB, b_LB, D_eval, y_eval)
    		for ii in range(nV * nPDF)
		)

	# Save results
	np.save(f'results_{which_period}.npz', results)
