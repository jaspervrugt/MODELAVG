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

# CASE STUDY: MODEL AVERAGING WITH WATER LEVELS

# I received this data from someone who was using the MODELAVG toolbox. 
# The data record is not long enough, but a fast practice for the 
# different methods 

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
from MODELAVG_functions import *

# DEFINE MODEL AVERAGING METHOD
method = 'bma'  # 'ewa'/'bga'/'aica'/'bica'/'gra'/'bma'/'mma'/'mma-s'

# BMA -> CONDITIONAL DISTRIBUTION NEEDS TO BE DEFINED
options = {
    'PDF': 'gamma',                                             # Gamma conditional pdf
    'VAR': 3,                                                   # Group nonconstant variance
    'alpha': np.array([0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]),   # Significance levels for BMA model
    'print': 'yes',                                             # Print output (figures, tables) to screen
}

# NOW LOAD DATA
W = np.loadtxt('water_levels.txt')  # Daily discharge forecasts (mm/day) and verifying data
idx_cal = slice(0, 25)              # Indices of training period (1-25 in MATLAB corresponds to 0-24 in Python)
idx_eval = slice(25, W.shape[0])    # Indices of evaluation data period

# DEFINE TRAINING ENSEMBLE AND VECTOR OF VERIFYING OBSERVATIONS
D_cal = W[idx_cal, 0:3]  # Ensemble (models 1-3)
y_cal = W[idx_cal, 3]    # Verifying data (column 4)

# APPLY LINEAR BIAS CORRECTION TO ENSEMBLE (UP TO USER)
D, a, b = Bias_correction(D_cal, y_cal)

# NUMBER OF PARAMETERS OF EACH MODEL
options['p'] = [3, 3, 4]  # Only used for AICA, BICA, MMA, MMA-S

if __name__ == '__main__':
	# Run MODELAVG toolbox
	phi, output = MODELAVG(method, D, y_cal, options)

	# DEFINE EVALUATION ENSEMBLE AND VECTOR OF VERIFYING OBSERVATIONS
	D_eval = W[idx_eval, 0:3]  # Ensemble for evaluation period
	y_eval = W[idx_eval, 3]    # Verifying data for evaluation period

	# Now evaluate BMA model
	val = MODELAVG_eval(method, D_eval, y_eval, options, a, b, output)
