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

%% CASE STUDY III: 48-FORECASTS OF SEA SURFACE PRESSURE
%% CHECK: A.E. RAFTERY ET AL., MWR, 133, pp. 1155-1174, 2005.

%% DEFINE MODEL AVERAGING METHOD
method = 'bma';             % 'ewa'/'bga'/'aica'/'bica'/'gra'/'bma'/'mma'/'mma-s'

%% BMA -> CONDITIONAL DISTRIBUTION NEEDS TO BE DEFINED
options.PDF = 'gamma';              % gamma distribution
options.VAR = '1';                  % individual non-constant variance
options.alpha = [0.01 0.05 0.1];    % significance levels: 1-alpha conf/pred intervals of BMA model
options.print = 'yes';              % print output (figures, tables) to screen
options.p = [10 10 10 20 30];
%% NOW LOAD DATA
P = load('pressure_data.txt');      % 48-hour forecasts of air-pressure (mbar) and verifying data

%% DEFINE ENSEMBLE AND VECTOR OF VERYFYING OBSERVATIONS ( APRIL 16 TO JUNE 9, 2000 )
idx = find(P(:,1) == 2000 & P(:,2) == 4 & P(:,3) == 16); start_idx = idx(1);
idx = find(P(:,1) == 2000 & P(:,2) == 6 & P(:,3) == 9); end_idx = idx(end);
D = P(start_idx:end_idx,5:9); y = P(start_idx:end_idx,4);
 
%% APPLY LINEAR BIAS CORRECTION TO ENSEMBLE ( UP TO USER )
[ D , a , b ] = Bias_correction ( D , y );

% Run MODELAVG toolbox with two outputs
[ phi , output ] = MODELAVG ( method , D , y , options ); 

% Can evaluate and print performance evaluation period using
% val = MODELAVG_eval ( method , D_eval , y_eval , options , a , b , output );