function [R] = gnormrnd(mu,alfa,beta,varargin)
%GNORMRND Random arrays from the generalized normal distribution.
%   R = GNORMRND(X,MU,ALFA,BETA,[n m]) returns random draws from the 
%   generalized normal distribution with mean MU, standard deviation ALFA, 
%   and kurtosis parameter BETA.
%
%   R = GNORMRND(MU,ALFA,BETA,M,N,...) or 
%   R = GNORMRND(MU,ALFA,BETA,[M,N,...])
%   returns an M-by-N-by-... array.
%
%   See also NORMRND, NORMCDF, NORMFIT, NORMINV, NORMLIKE, NORMPDF, 
%   NORMSTAT, RANDOM, RANDN.
%   References:

%   Copyright 2022, Jasper A. Vrugt

if nargin < 3
    error(message('stats:gnormrnd:TooFewInputs'));
end

[err, sizeOut] = internal.stats.statsizechk(3,mu,alfa,beta,varargin{:});
if err > 0
    error(message('stats:gnormrnd:InputSizeMismatch'));
end

ty = internal.stats.dominantType(mu,alfa,beta);
p = rand(sizeOut,'like',ty);

% get draws from quantiles
R = mu + sign(p - 1/2).*(alfa^beta*gaminv(2*abs(p-1/2),1/beta)).^(1/beta); 
