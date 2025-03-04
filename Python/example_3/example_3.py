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

# CASE STUDY II: 48-FORECASTS OF SEA SURFACE PRESSURE
# CHECK: A.E. RAFTERY ET AL., MWR, 133, pp. 1155-1174, 2005.

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
    'PDF': 'gamma',  # gamma distribution
    'VAR': 1,        # individual non-constant variance
    'alpha': np.array([0.01, 0.05, 0.1]),  # significance levels: 1-alpha conf/pred intervals of BMA model
    'print': 'yes',  # print output (figures, tables) to screen
    'p': [10, 10, 10, 20, 30]  # some additional parameters (not specified in detail)
}

# NOW LOAD DATA
P = np.loadtxt('pressure_data.txt')  	# 48-hour forecasts of air-pressure (mbar) and verifying data 

# DEFINE ENSEMBLE AND VECTOR OF VERIFYING OBSERVATIONS (APRIL 16 TO JUNE 9, 2000)
start_idx = np.where((P[:, 0] == 2000) & (P[:, 1] == 4) & (P[:, 2] == 16))[0][0]
end_idx = np.where((P[:, 0] == 2000) & (P[:, 1] == 6) & (P[:, 2] == 9))[0][-1]
D = P[start_idx:end_idx + 1, 4:9]  # Ensemble data (columns 5-9)
y = P[start_idx:end_idx + 1, 3]    # Verifying data (column 4)

# APPLY LINEAR BIAS CORRECTION TO ENSEMBLE (UP TO USER)
D, a, b = Bias_correction(D, y)

if __name__ == '__main__':
	# Run MODELAVG toolbox
	phi, output = MODELAVG(method, D, y, options)

	# Evaluate performance evaluation period
	val = MODELAVG_eval(method, D_eval, y_eval, options, a, b, output)

