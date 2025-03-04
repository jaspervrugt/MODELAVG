function rlbl = BMA_rlbl(cdf_y)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% This function computes reliability of BMA forecast distribution using   %
% formulation of Renard (2011)                                            %
%                                                                         %
% SYNOPSIS: rlbl = BMA_rlbl(cdf_y)                                        %
%  where                                                                  %
%   cdf_y     [input] nx1 vector with cdf BMA mixture at verifying data   %
%   rlbl      [outpt] scalar with reliability according to Renard, 2011   %
%                                                                         %
% (c) Written by Jasper A. Vrugt, Feb 2019                                %
% University of California Irvine                                         %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% This function computes the reliability of the BMA forecast distribution
% using the formulation of Renard et al. (2011)
n = numel(cdf_y);                   % # observations
eCDF = sort(cdf_y);                 % empirical CDF (CDF ascending order)
uCDF = (1:n)'./n;                   % cdf of uniform distribution
rlbl = 2/n * sum(abs(uCDF - eCDF)); % reliability: negatively oriented 
                                    % smaller = better

end
