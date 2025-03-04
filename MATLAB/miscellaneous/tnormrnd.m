function [R] = tnormrnd(mu,sigma,a,b,varargin)
%TNORMRND Random arrays from the truncated normal distribution.
%   R = TNORMRND(X,MU,SIGMA,A,B,[n m]) returns random draws from the 
%   truncated normal distribution with mean MU, standard deviation SIGMA, 
%   and lower and upper range, A and B.
%
%   R = TNORMRND(MU,SIGMA,A,B,M,N,...) or 
%   R = TNORMRND(MU,SIGMA,A,B,[M,N,...])
%   returns an M-by-N-by-... array.
%
%   See also NORMRND, NORMCDF, NORMFIT, NORMINV, NORMLIKE, NORMPDF, 
%   NORMSTAT, RANDOM, RANDN.
%   References:

%   Copyright 2022, Jasper A. Vrugt

if nargin < 4
    error(message('stats:gnormrnd:TooFewInputs'));
end

[err, sizeOut] = internal.stats.statsizechk(4,mu,sigma,a,b,varargin{:});
if err > 0
    error(message('stats:gnormrnd:InputSizeMismatch'));
end

ty = internal.stats.dominantType(mu,sigma,a,b);
p = rand(sizeOut,'like',ty);

% Normalize input variables
alfa = (a-mu)./sigma; beta = (b-mu)./sigma;

% CDF = p = 1/2*(1 + erf(x/sqrt(2)));
% iCDF = sqrt(2)*erfinv(2*p - 1);

CDF = @(x) 1/2*(1 + erf(x/sqrt(2))); 
iCDF = @(p) sqrt(2)*erfinv(2*p-1);

% Draw samples
R = iCDF (CDF(alfa) + p.*(CDF(beta) - CDF(alfa))) * sigma + mu;