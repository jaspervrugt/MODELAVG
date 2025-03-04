function [D_bc,a,b] = Bias_correction(D,y,intcpt)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% This function provides a linear correction of the ensemble members      %
%                                                                         %
% SYNOPSIS: [D_bc,a,b] = Bias_correction(D,y)                             %
%           [D_bc,a,b] = Bias_correction(D,y,intcpt)                      %
%  where                                                                  %
%   D         [input] nxK matrix with forecasts of ensemble members       %
%   y         [input] nx1 vector with verifying observations              %
%   intcpt    [input] OPT: intercept (default: 1) or without (0)          %
%   D_bc      [outpt] nxK matrix with bias-corrected forecasts            %
%   a         [outpt] 1xK vector with intercept bias-corrected forecasts  %
%   b         [outpt] 1xK vector with slope bias-corrected forecasts      %
%                                                                         %
% (c) Written by Jasper A. Vrugt, Feb 2022                                %
% University of California Irvine 			        	                  %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

if nargin < 3, intcpt = 1; end

% Determine the size of D, the matrix with ensemble forecasts
[n,K] = size(D); D_bc = nan(n,K);
% Initialize intercepts and slopes of linear bias correction functions
[a,b] = deal(zeros(1,K)); 
% Linear regression for bias correct ensemble forecasts
for k = 1:K
    switch intcpt
        case 1      % Linear regression function with slope and intercept
            ab = regress(y,[ones(n,1) D(1:n,k)]);
            a(k) = ab(1); b(k) = ab(2);
        otherwise   % Linear regression function with slope only
            b(k) = regress(y,D(1:n,k));
    end           
    % Bias-corrected ensemble forecasts
    D_bc(1:n,k) = a(k) + b(k) * D(1:n,k);
end

end
