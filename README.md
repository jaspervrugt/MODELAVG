# MODELAVG: General-purpose model averaging toolbox in MATLAB and Python

## Description

Multi-model averaging is widely used in various scientific and engineering disciplines to post-process forecast ensembles (sharpness and skill) and/or quantify conceptual model uncertainty. Here, I present a MATLAB and Python toolbox called MODELAVG which implements a suite of different model averaging techniques, including (among others) equal weights averaging (EWA), Bates-Granger model averaging (BGA), Bayesian model averaging (BMA), Mallows model averaging (MMA), and Granger-Ramanathan averaging (GRA). Some of these methods return a point forecast only, whereas others (BMA) produce a forecast distribution of the variable(s) of interest. MCMC simulation with the DREAM$_{(ZS)}$ algorithm is used for averaging methods that do not have a direct closed-form solution for their weights (and/or shape parameters). For BMA the user can select among the normal, lognormal, generalized normal, truncated normal, gamma, Weibull, generalized extreme value and generalized Pareto distributions as conditional forecast PDFs for the ensemble members. Furthermore, they can choice among a constant and non-constance variance of these PDFs. The toolbox returns optimal values of the weights and related parameters of each model averaging method. For MMA and BMA the user has access to the entire posterior distribution of the weights and/or variances derived from MCMC simulation with the DREAM algorithm. The postprocessor summarizes the results in different figures, which are written to a PDF file. This includes traceplots of the sampled parameters, histograms of the marginal distributions of the BMA model parameters, convergence diagnostics of the sampled chains, confidence and predictive uncertainty intervals and scoring rules and performance metrics of the BMA predictive distribution. Built-in case studies with forecast ensembles of hydrologic and meteorologic models illustrate the capabilities and functionalities of the MODELAVG toolbox. 

## Getting Started

### Installing: MATLAB

* Download and unzip the zip file 'MATLAB_code_MODELAVG_V2.0.zip' in a directory 'MODELAVG'
* Add the toolbox to your MATLAB search path by running the script 'install_MODELAVG.m' available in the root directory
* You are ready to run the examples

### Executing program

* After intalling, you can simply direct to each example folder and execute the local 'example_X.m' file
* Please Make sure you read carefully the instructions (i.e., green comments) in 'install_MODELAVG.m' and the manual !!!  

### Installing: Python

* Download and unzip the zip file 'Python_code_MODELAVG_V2.0.zip' to a directory called 'MODELAVG'

### Executing program

* Go to Command Prompt and directory of example_X in the root of MODELAVG
* Now you can execute this example by typing 'python example_X.py'.
* Instructions can be found in the file 'MODELAVG.py' and in the manual !!!  

## Authors

* Vrugt, Jasper A. (jasper@uci.edu) 

## Literature
1. Vrugt, J.A. (2023), MODELAVG: A MATLAB toolbox for postprocessing of model ensembles, Manual, Version 2.0, pp. 1 - XX, 2023
2. Vrugt, J.A. (2024), Distribution-Based Model Evaluation and Diagnostics: Elicitability, Propriety, and Scoring Rules for Hydrograph Functionals, _Water Resources Research_, 60, e2023WR036710, https://doi.org/10.1029/2023WR036710
3. Vrugt, J.A. (2015), Markov chain Monte Carlo simulation using the DREAM software package: Theory, concepts, and MATLAB implementation, _Environmental Modeling and Software_, 75, pp. 273-316
4. Diks, C.G.H., and J.A. Vrugt (2010), Comparison of point forecast accuracy of model averaging methods in hydrologic applications, _Stochastic Environmental Research and Risk Assessment_, 24(6), pp. 809-820, https://doi.org/10.1007/s00477-010-0378-z
5. Vrugt, J.A., C.G.H. Diks, and M.P. Clark (2008), Ensemble Bayesian model averaging using Markov chain Monte Carlo sampling, _Environmental Fluid Mechanics_, 8(5-6), 579-595, https://doi.org/10.1007/s10652-008-9106-3
6. Vrugt, J.A., C.J.F. ter Braak, M.P. Clark, J.M. Hyman and B.A. Robinson (2008), Treatment of input uncertainty in hydrologic modeling: Doing hydrology backward with Markov chain Monte Carlo simulation, 44 (12), _Water Resources Research_, https://doi.org/10.1029/2007WR006720
7. Vrugt, J.A., and B.A. Robinson (2007), Treatment of uncertainty using ensemble methods: Comparison of sequential data assimilation and Bayesian model averaging, _Water Resources Research_, 43, W01411, https://doi.org/10.1029/2005WR004838

## Version History

* 1.0
    * Initial Release
* 2.0
    * New capabilities of BMA method (conditional PDF), added scoring rules, performance metrics and computation for evaluation period
    * Python implementation
    * Source code in MATLAB and Python

## Built-in Case Studies
1. Example 1: 24-hour forecasts of river discharge
2. Example 2: 48-forecasts of sea surface temperature
3. Example 3: 48-forecasts of sea surface pressure
4. Example 11: Forecasts of water levels
5. Example 12: Hydrologic modeling
6. Example 13: Flood modeling
   
## Acknowledgments
