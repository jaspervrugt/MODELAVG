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

# CASE STUDY I: RAINFALL-RUNOFF TRANSFORMATION - ENSEMBLE OF CALIBRATED WATERSHED MODELS
# CHECK: J.A. VRUGT AND B.A. ROBINSON, WRR, 43, W01411, doi:10.1029/2005WR004838, 2007

import sys
import os
import numpy as np

# Get the current working directory
current_directory = os.getcwd()
# Go up one directory
parent_directory = os.path.abspath(os.path.join(current_directory, '..'))
# add this to path
sys.path.append(parent_directory)
# Add another directory
misc_directory = os.path.abspath(os.path.join(parent_directory, 'miscellaneous'))
# add this to path
sys.path.append(misc_directory)

from MODELAVG import MODELAVG
from MODELAVG_functions import Bias_correction, MODELAVG_eval

# DEFINE MODEL AVERAGING METHOD
method = 'bma'  # 'ewa'/'bga'/'aica'/'bica'/'gra'/'bma'/'mma'/'mma-s'

# BMA -> CONDITIONAL DISTRIBUTION NEEDS TO BE DEFINED
options = {
    'PDF': 'normal',                            # normal conditional pdf
    'VAR': '1',                                 # individual constant variance
    'TAU': '2',                                 # for generalized normal (gen_normal) - own tau for each member individually
    'alpha': np.array([0.01, 0.05, 0.1, 0.5]),  # significance levels: 1-alpha conf/pred intervals of BMA model
    'print': 'yes',                             # print output (figures, tables) to screen
    'postproc': 'yes'}                          # no postprocessing

# NOW LOAD DATA
S = np.loadtxt('discharge_data.txt')  	# daily discharge forecasts (mm/day) of models and verifying data 
idx_cal = slice(0, 500)           	    # start/end training period
idx_eval = slice(5000, 5500)      	    # evaluation period

# DEFINE ENSEMBLE AND VECTOR OF VERIFYING OBSERVATIONS
D = S[idx_cal, :8]  	# Ensemble of models
y = S[idx_cal, 8]   	# Verifying observations

# APPLY LINEAR BIAS CORRECTION TO ENSEMBLE (UP TO USER)
# Bias_correction function needs to be defined (or imported if available)
D, a, b = Bias_correction(D, y)

# NUMBER OF PARAMETERS OF EACH MODEL (ABC/GR4J/HYMOD/TOPMO/AWBM/NAM/HBV/SACSMA)
options['p'] = [3, 4, 5, 8, 8, 9, 9, 13]   # (only used for AICA, BICA, MMA, MMA-S)

if __name__ == '__main__':
	# Run MODELAVG toolbox
	phi, output = MODELAVG(method, D, y, options)

	# NOW DEFINE EVALUATION DATA
	D_eval = S[idx_eval, :8]
	y_eval = S[idx_eval, 8]

	# Now evaluate performance for evaluation period
	val = MODELAVG_eval(method, D_eval, y_eval, options, a, b, output)[0]
