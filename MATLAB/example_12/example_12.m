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

% About the discharge ensemble of watershed models
%   Vrugt, J.A., and B.A. Robinson (2007), Treatment of uncertainty using 
%        ensemble methods: Comparison of sequential data assimilation and 
%        Bayesian model averaging, Water Resources Research, 43, W01411, 
%        doi:10.1029/2005WR004838.
% About the scoring rules
%   Vrugt, J.A. (2024), Distribution-Based Model Evaluation and 
%        Diagnostics: Elicitability, Propriety, and Scoring Rules for 
%        Hydrograph Functionals, Water Resources Research, 60, 
%        e2023WR036710, https://doi.org/10.1029/2023WR036710  

%% DEFINE MODEL AVERAGING METHOD
method = 'bma';             % 'ewa'/'bga'/'aica'/'bica'/'gra'/'bma'/'mma'/'mma-s'

%% FIXED SETTINGS
options.alpha = [0.01 0.05 0.1 0.5];    % significance levels: 1-alpha conf/pred intervals of BMA model
options.print = 'no';                   % print output (figures, tables) to screen
options.CPU = 'yes';                    % efficient CPU implementation

%% NOW LOAD DATA
S = load('discharge.txt');  % daily discharge forecasts (mm/day) of models and verifying data
idx_cal = 1:1:3000;         % start/end training period
idx_eval = 5001:7000;       % evaluation period
%% DEFINE ENSEMBLE AND VECTOR OF VERYFYING OBSERVATIONS
D = S(idx_cal,1:8); y_cal = S(idx_cal,9); D_eval = S(idx_eval,1:8); y_eval = S(idx_eval,9);
%% APPLY LINEAR BIAS CORRECTION TO ENSEMBLE ( UP TO USER )
[ D_cal , a_LB , b_LB ] = Bias_correction ( D , y_cal );
%% NUMBER OF PARAMETERS OF EACH MODEL (ABC/GR4J/HYMOD/TOPMO/AWBM/NAM/HBV/SACSMA)
%% options.p = [ 3 4 5 8 8 9 9 13 ];   % ( only used for AICA, BICA, MMA, MMA-S)

%% BMA CONDITIONAL DISTRIBUTIONS TO BE USED
PDF = {'normal', ...        % Normal distribution
    'lognormal', ...        % Lognormal distribution
    'gen_normal', ...       % Generalized normal distribution
    'tnormal', ...          % Truncated normal distribution
    'gamma', ...            % Gamma distribution
    'gev', ...              % Generalized Pareto distribution
    'gpareto'};             % Generalized extreme value distribution
V = {'1','2','3','4'};      % homoscedastic variance; group: '1' ; single: '2'
% heteroscedastic variance; group: '3' ; single: '4'

%% BMA -> CONDITIONAL DISTRIBUTION NEEDS TO BE DEFINED
options.TAU = '2';          % Gen. normal: different τ each ensemble member
                            % Gen. Pareto: different xi (shape) par. each ensemble member
                            % Gen. extreme value: different xi (shape) par. each ensemble member
which_period = 'evaluation'; % 'training';
% Table formulation
% 'normal','lognormal','gen_normal','tnormal','gamma','gev','gpareto'
%    1          2           3           4        5      6      7     Var = '1'
%    8          9          10          11       12     13     14     Var = '2'
%   15         16          17          18       19     20     21     Var = '3'
%   22         23          24          25       26     27     28     Var = '4'
names_pdf = {'1: normal: const group',...
    '2: lnormal: const group', ...
    '3: gnormal: const group', ...
    '4: tnormal: const group', ...
    '5: gamma: const group', ...
    '6: gev: const group', ...
    '7: gp: const group', ...
    '8: normal: const ind',...
    '9: lnormal: const ind', ...
    '10: gnormal: const ind', ...
    '11: tnormal: const ind', ...
    '12: gamma: const ind', ...
    '13: gev: const ind', ...
    '14: gp: const ind', ...
    '15: normal: nonconst group',...
    '16: lnormal: nonconst group', ...
    '17: gnormal: nonconst group', ...
    '18: tnormal: nonconst group', ...
    '19: gamma: nonconst group', ...
    '20: gev: nonconst group', ...
    '21: gp: nonconst group', ...
    '22: normal: nonconst ind',...
    '23: lnormal: nonconst ind', ...
    '24: gnormal: nonconst ind', ...
    '25: tnormal: nonconst ind', ...
    '26: gamma: nonconst ind', ...
    '27: gev: nonconst ind', ...
    '28: gp: nonconst ind'};

%% INITIALIZE VARIABLES/OUTPUT ARRAYS
nV = numel(V); nPDF = numel(PDF);
switch which_period
    case 'training'
        nY = numel(y_cal);
    case 'evaluation'
        nY = numel(y_eval);
end
[cal,val] = deal(cell(1,nV*nPDF));
% Prepare to run in parallel
ct = 1; VPDF = cell(nV*nPDF,2);
for i = 1:nV
    for j = 1:nPDF
        VPDF(ct,1:2) = [V(i) , PDF(j)];
        ct = ct + 1;
    end
end

%% Now do BMA method for all combinations of variances/conditional PDFs
parfor ii = 1:nV*nPDF
    ii
    % Call function below - parallel implementation
    [cal{ii},beta{ii},AR{ii},table_res_cal(:,ii),val{ii},table_res_val(:,ii)] = ...
        parallel_script(ii,VPDF(ii,:),method,D_cal,y_cal,options,...
        which_period,a_LB,b_LB,D_eval,y_eval);
end
evalstr = strcat('save results_',which_period); eval(char(evalstr));

%% print some tables
print_tables

%% print some figures
plot_figures