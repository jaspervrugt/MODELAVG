function [R] = powrnd(mu,sigma,alfa,varargin)
%GNORMRND Random arrays from the lognormal power-law distribution
%   R = GNORMRND(X,MU,SIGMA,ALFA,[n m]) returns random draws from the 
%   generalized normal distribution with mean MU, standard deviation SIGMA, 
%   and shape parameter ALFA.
%
%   R = GNORMRND(MU,SIGMA,ALFA,M,N,...) or 
%   R = GNORMRND(MU,SIGMA,ALFA,[M,N,...])
%   returns an M-by-N-by-... array.
%
%   References:

%   Copyright 2022, Jasper A. Vrugt

if nargin < 3
    error(message('stats:gnormrnd:TooFewInputs'));
end

[err, sizeOut] = internal.stats.statsizechk(3,mu,sigma,alfa,varargin{:});
if err > 0
    error(message('stats:gnormrnd:InputSizeMismatch'));
end

ty = internal.stats.dominantType(mu,sigma,alfa);
p = rand(sizeOut,'like',ty); pn = randn(sizeOut,'like',ty);

% Get lognormal random variates
R = exp(mu + sigma*pn - 1/alfa * log(p));