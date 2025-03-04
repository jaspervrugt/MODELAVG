function [beta,chain,n,K,sigma_2,log_L,par,count,options] = ...
    MODELAVG_setup(method,D,Y,options,flag)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%                                                                       %%
%% MM      MM  OOOOOO  DDDDDDD  EEEEEEE LL        AAAA   VV   VV GGGGGGG %%
%% MMM     MM OOOOOOOO DDDDDDDD EEEEEEE LL       AA  AA  VV   VV GG   GG %%
%% MMMM  MMMM OO    OO DD    DD EE      LL       AA  AA  VV   VV GG   GG %%
%% MM MMMM MM OO    OO DD    DD EEEE    LL       AA  AA  VV   VV GGGGGGG %%
%% MM  MM  MM OO    OO DD    DD EEEE    LL       AAAAAA   VV VV  GGGGGGG %%
%% MM      MM OO    OO DD    DD EE      LL      AA    AA   VVV        GG %%
%% MM      MM OOOOOOOO DDDDDDDD EEEEEEE LLLLLLL AA    AA   VVV       GGG %%
%% MM      MM  OOOOOO  DDDDDDD  EEEEEEE LLLLLLL AA    AA    V    GGGGGGG %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Defines computational framework of MODELAVG                           %%
%%                                                                       %%
%% SYNOPSIS: [beta,chain,n,K,sigma_2,log_L,par,count,options] = ...      %%
%%               MODELAVG_setup(method,D,Y,options,flag)                 %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%                                                                       %%
%%   (c) Written by Jasper A. Vrugt, Feb 2012                            %%
%%   University of California Irvine                                     %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

if nargin < 5, flag = 0; end

% Assign default settings in case not specified
if strcmp(method,'bma')
    if ~isfield(options,'VAR')
        % Set variance treatment of BMA to '1'
        options.VAR = '1';
    end
    if ~isfield(options,'alpha')
        % Set alpha to 95% intervals
        options.alpha = 0.05;
    end   
end

if ~isfield(options,'print')
    % Set print to yes
    options.print = 'yes';
end

if ~isfield(options,'postproc')
    options.postproc = 'yes';
end

% Initialize beta and chain to be empty
beta = []; chain = [];

% Determine size of ensemble forecast
[ n , K ] = size ( D );

% Calculate error variance of each model
sigma_2 = sum ( bsxfun(@minus,D,Y).^2 )/n;

% Calculate log-likelihood of each model
log_L = -1/2 * ( n * sigma_2 + n );

% Scale to avoid numerical underflow
log_L = log_L - max(log_L);

% Initialize count and structure par
par = []; count = 0;

% disp('---------- Summary of the main settings used: method ------------------');
evalstr = char(strcat('Model averaging method used:',{' '},...
    upper(method),'\n')); fprintf(evalstr);
% disp('-----------------------------------------------------------------------');
% disp('---------- Summary of the main settings used: options -----------------');
if flag == 0
    % print options
    disp(options)
end
% disp('-----------------------------------------------------------------------');