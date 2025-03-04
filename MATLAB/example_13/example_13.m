% ----------------------------------------------------------------------- %
%  MM   MM   OOOOO   DDDDDD   EEEEEEE  LL         AAA    VV   VV   GGGGG  %
%  MM   MM  OOOOOOO  DDDDDDD  EEEEEEE  LL        AA AA   VV   VV  GG   GG %
%  MMM MMM  OO   OO  DD   DD  EE       LL        AA AA   VV   VV  GG   GG %
%  MM M MM  OO   OO  DD   DD  EEEEE    LL       AA   AA  VV   VV  GGGGGGG %
%  MM   MM  OO   OO  DD   DD  EEEEE    LL       AAAAAAA  VV   VV   GGGGGG %
%  MM   MM  OO   OO  DD   DD  EE       LL       AA   AA   VV VV        GG %
%  MM   MM  OOOOOOO  DDDDDDD  EEEEEEE  LLLLLLL  AA   AA    VVV     GGGGGG %
%  MM   MM   OOOOO   DDDDDD   EEEEEEE  LLLLLLL  AA   AA     V     GGGGGGG %
% ----------------------------------------------------------------------- %

% This case study is based on the following paper
%   Vrugt, J.A. (2024), Comment on "Improving Bayesian Model Averaging for
%        Ensemble Flood Modeling Using Multiple Markov Chains Monte Carlo
%        Sampling," by Huang and Merwade, published in Water Resources
%        Research, 59, e2023WR034947, https://doi.org/10.1029/2023WR034947
% Please also check
%   Vrugt, J.A. (2023), MODELAVG: A MATLAB toolbox for postprocessing of
%        model ensembles, Manual, Version 2.0, pp. 1 - XX, 2023.
%   Vrugt, J.A. (2024), Distribution-Based Model Evaluation and
%        Diagnostics: Elicitability, Propriety, and Scoring Rules for
%        Hydrograph Functionals, Water Resources Research, 60,
%        e2023WR036710, https://doi.org/10.1029/2023WR036710
%

%% DEFINE MODEL AVERAGING METHOD
method = 'bma';             % 'ewa'/'bga'/'aica'/'bica'/'gra'/'bma'/'mma'/'mma-s'

%% FIXED SETTINGS
options.alpha = [0.01 0.05 0.1 0.5];    % significance levels
options.print = 'yes';                  % print output (figures, tables)
options.CPU = 'yes';                    % CPU-efficient computation

%% NOW LOAD DATA
data_pred = table2array(...
    readtable('PredictionData_feet_Numerical_Experiments.csv',...
    Range='A2:J101'));
data_meas = table2array(...
    readtable('ObservedData_feet_Numerical_Experiments.csv',...
    Range='A1:A100'));
idx_cal = 1:1:100;                      % start/end training period [= all data]

%% DEFINE ENSEMBLE AND VECTOR OF VERYFYING OBSERVATIONS
D = data_pred(idx_cal,1:10); y_cal = data_meas(idx_cal,1);
D_eval = []; y_eval = [];
%% APPLY LINEAR BIAS CORRECTION TO ENSEMBLE (= RECOMMENDED)
[ D_cal , a_LB , b_LB ] = Bias_correction ( D , y_cal );

%% BMA CONDITIONAL DISTRIBUTIONS TO BE USED
PDF = {'gen_normal'}    % Normal distribution
V = {'2'};          % homoscedastic variance; group: '1' ; single: '2'
options.TAU = '2';  % Generalized normal: each member its own τ (xi or k)
which_period = 'training';
name_pdf = strcat(PDF,V);
switch which_period
    case 'training'
        nY = numel(y_cal);
    case 'evaluation'
        nY = numel(y_eval);
end
VPDF = [ V , PDF ]; ii = 1;

% Call function parallel_script: for parallel implementation
[cal,beta,AR,table_res_cal,val,...
    table_res_val] = parallel_script(ii,VPDF,...
    method,D_cal,y_cal,options,which_period,a_LB,b_LB,D_eval,y_eval);
% Save the figures
figHandles = findall(0,'Type','figure'); nfig = numel(figHandles);
file_figs = char(strcat(name_pdf,'.pdf'));
print(figHandles(nfig),'-dpdf',file_figs,'-fillpage');
for zz = nfig-1:-1:1
    figure(figHandles(zz)); set(gcf,'color','w');
    exportgraphics(figHandles(zz), file_figs, 'Append', true);
end