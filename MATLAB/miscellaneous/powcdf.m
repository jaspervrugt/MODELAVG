function [P] = powcdf(x,mu,sigma,alfa)
%GNORMCDF Lognormal power-law cumulative distribution function (cdf).
%   P = POWCDF(X,MU,SIGMA,ALFA) returns the cdf of the lognormal power-law
%   distribution with mean MU, standard deviation SIGMA, and shape
%   parameter ALFA, evaluated at the values in X.
%   The size of P is the common size of X, MU, SIGMA and ALFA. A scalar
%   input functions as a constant matrix of the same size as the other
%   inputs.
%
%   References:

%   Copyright 2022, Jasper A. Vrugt

% Negative data would create complex values, which erfc cannot handle.
P = zeros(size(x)); id = (x > 0);

P(id) = 1/2 * erfc( - (log(x(id)) - mu)./(sqrt(2) * sigma)) ...
    - 1/2 * exp(alfa .* mu + alfa.^2 * sigma.^2/2) ...
    .* x(id).^(-alfa) ...
    .* erfc(alfa * sigma/sqrt(2) - (log(x(id)) - mu)/(sqrt(2) * sigma) );
% NOTE: WRONG IN WIKIPEDIA