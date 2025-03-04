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

%% I received this data from someone who was using the MODELAVG toolbox. 
%% The data record is not long enough, but a fast practice for the 
%% different methods 

%% DEFINE MODEL AVERAGING METHOD
method = 'bma';             % 'ewa'/'bga'/'aica'/'bica'/'gra'/'bma'/'mma'/'mma-s'

%% BMA -> CONDITIONAL DISTRIBUTION NEEDS TO BE DEFINED
options.PDF = 'gamma';                  % gamma conditional pdf
%options.TAU = '2';                     % variable tau for all models - only for gen_normal PDF
options.VAR = '3';                      % group nonconstant variance
options.alpha = [0.01 0.05 0.1 0.2 0.3 0.4 0.5];    % significance levels: 1-alpha conf/pred intervals of BMA model
options.print = 'yes';                  % print output (figures, tables) to screen

%% NOW LOAD DATA
W = load('water_levels.txt');   % daily discharge forecasts (mm/day) of models and verifying data 
idx_cal = 1:25;                 % indices of training period
idx_eval = 26:size(W,1);        % indices of evaluation data period

%% DEFINE TRAINING ENSEMBLE AND VECTOR OF VERYFYING OBSERVATIONS
D_cal = W(idx_cal,1:3); y_cal = W(idx_cal,4); 

%% APPLY LINEAR BIAS CORRECTION TO ENSEMBLE ( UP TO USER )
[ D , a , b ] = Bias_correction ( D_cal , y_cal );

%% NUMBER OF PARAMETERS OF EACH MODEL (ABC/GR4J/HYMOD/TOPMO/AWBM/NAM/HBV/SACSMA)
options.p = [ 3 3 4 ];   % ( only used for AICA, BICA, MMA, MMA-S)

%% Run MODELAVG toolbox with two outputs
[ phi , output ] = MODELAVG ( method , D , y_cal , options );

%% DEFINE EVALUATION ENSEMBLE AND VECTOR OF VERYFYING OBSERVATIONS
D_eval = W(idx_eval,1:3); y_eval = W(idx_eval,4);

%% Now evaluate BMA model
val = MODELAVG_eval ( method , D_eval , y_eval , options , a , b , output );