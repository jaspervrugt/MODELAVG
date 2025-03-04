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

%% CASE STUDY I: RAINFALL-RUNOFF TRANSFORMION - ENSEMBLE OF CALIBRATED WATERSHED MODELS
%% CHECK: J.A. VRUGT AND B.A. ROBINSON, WRR, 43, W01411, doi:10.1029/2005WR004838, 2007

%% DEFINE MODEL AVERAGING METHOD
method = 'bma';             % 'ewa'/'bga'/'aica'/'bica'/'gra'/'bma'/'mma'/'mma-s'

%% BMA -> CONDITIONAL DISTRIBUTION NEEDS TO BE DEFINED
options.PDF = 'normal';                 % normal conditional pdf
options.VAR = '1';                      % individual constant variance
options.TAU = '2';                      % for generalized normal (gen_normal) - own tau for each member individually
options.alpha = [0.01 0.05 0.1 0.5];    % significance levels: 1-alpha conf/pred intervals of BMA model
options.print = 'yes';                  % print output (figures, tables) to screen
%options.CPU = 'yes';                    % CPU efficient or not
options.postproc = 'yes';                % no postprocessing

%% NOW LOAD DATA
S = load('discharge_data.txt'); % daily discharge forecasts (mm/day) of models and verifying data 
idx_cal = 1:500;                % start/end training period
idx_eval = 5001:5500;           % evaluation period

%% DEFINE ENSEMBLE AND VECTOR OF VERYFYING OBSERVATIONS
D = S(idx_cal,1:8); y = S(idx_cal,9);

%% APPLY LINEAR BIAS CORRECTION TO ENSEMBLE ( UP TO USER )
[ D , a , b ] = Bias_correction ( D , y );

%% NUMBER OF PARAMETERS OF EACH MODEL (ABC/GR4J/HYMOD/TOPMO/AWBM/NAM/HBV/SACSMA)
options.p = [ 3 4 5 8 8 9 9 13 ];   % ( only used for AICA, BICA, MMA, MMA-S)

%% Run MODELAVG toolbox with two outputs
[ phi , output ] = MODELAVG ( method , D , y , options );

%% NOW DEFINE EVALUATION DATA
D_eval = S(idx_eval,1:8); y_eval = S(idx_eval,9);

%% NOW CALCULATE STATISTICS FOR EVALUATION PERIOD
val = MODELAVG_eval ( method , D_eval , y_eval , options , a , b , output )