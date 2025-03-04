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

%% CASE STUDY II: 48-FORECASTS OF SEA SURFACE TEMPERATURE
%% CHECK: A.E. RAFTERY ET AL., MWR, 133, pp. 1155-1174, 2005.

%% DEFINE MODEL AVERAGING METHOD
method = 'bma';             % 'ewa'/'bga'/'aica'/'bica'/'gra'/'bma'/'mma'/'mma-s'

%% BMA -> CONDITIONAL DISTRIBUTION NEEDS TO BE DEFINED
options.PDF = 'gen_normal';             % pdf predictor: normal/gamma/gen_normal
options.TAU = '2';                      % each model has its own tau value
options.VAR = '2';                      % individual constant variance
options.alpha = [0.01 0.05 0.1 0.5];    % significance levels: 1-alpha conf/pred intervals of BMA model
options.print = 'yes';                  % print output (figures, tables) to screen

%% NOW LOAD DATA
T = load('temp_data.txt');              % 48-hour forecasts temperature (Kelvin) and verifying data 

%% DEFINE ENSEMBLE AND VECTOR OF VERYFYING OBSERVATIONS ( APRIL 16 TO JUNE 9, 2000 )
idx = find(T(:,1) == 2000 & T(:,2) == 4 & T(:,3) == 16); start_idx = idx(1);
idx = find(T(:,1) == 2000 & T(:,2) == 6 & T(:,3) == 9); end_idx = idx(1);
D = T(start_idx:end_idx,5:9); y = T(start_idx:end_idx,4);

%% APPLY LINEAR BIAS CORRECTION TO ENSEMBLE ( UP TO USER )
[ D , a , b ] = Bias_correction ( D , y );

% Run MODELAVG toolbox with two outputs
[ phi , output ] = MODELAVG ( method , D , y , options ); 

% Can evaluate and print performance evaluation period using
% val = MODELAVG_eval ( method , D_eval , y_eval , options , a , b , output );