# MODELAVG: General-purpose model averaging toolbox in MATLAB

## Description

**Model averaging** is statistical method that is widely used to quantify the conceptual uncertainty of environmental system models and to improve the sharpness and skill of forecast ensembles of multi-model prediction systems. Here, I present a MATLAB toolbox for postprocessing of forecast ensembles. This toolbox, called MODELAVG implements many different model averaging techniques, including methods that provide point forecasts only, and methods that produce a forecast distribution of the variable(s) of interest. MCMC simulation with the DREAM_{(ZS)} algorithm is used for averaging methods without a direct closed-form solution of their point forecasts. The toolbox returns to the user (among others) a vector (or matrix with posterior samples) of weights and (if appropriate) standard deviation(s) of the members' forecast distribution, a vector of averaged forecasts (and performance metrics thereof), and (if appropriate) estimates of the width and coverage of the forecast distribution, and convergence diagnostics of the DREAM algorithm. The toolbox also creates many different figures with the results of each method. Three case studies illustrate the capabilities of the MODELAVG toolbox.

V2.0 is an updated version of the **Model averaging** toolbox. This has new options for the conditional PDF of the ensemble members. Users can now also choose the Weibull, lognormal, truncated normal, generalized extreme value and generalized Pareto distribution as predictive PDF of each ensemble member. Furthermore, V2.0 also returns values of strictly proper scoring rules (quadratic, logarithmic, spherical and continuous ranked probability scores) and performance metrics such as the reliability and coefficient of variation of the BMA distribution forecasts. These statistics are also computed for an independent evaluation data period.

## Getting Started

### Installing

* Download and unzip the zip file 'MATLAB_pcode_MODELAVG_V2.0.zip'.
* Add the toolbox to your MATLAB search path by running the script 'install_MODELAVG.m' available in 'MATLAB_pcode_MODELAVG_V2.0'
* You are ready to run the examples.

### Executing program

* After intalling, you can simply direct to each example folder and execute the local 'example_X.m' file.
* Please Make sure you read carefully the instructions (i.e., green comments) in 'install_MODELAVG.m' and the manual in PDF !!!  

## Authors

* Vrugt, Jasper A. (jasper@uci.edu) 

## Version History

* 1.0
    * Initial Release
* 2.0
    * New Release

## Acknowledgments
