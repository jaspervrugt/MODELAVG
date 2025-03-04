function [P] = gnormcdf(x,mu,alfa,beta)
%GNORMCDF Generalized normal cumulative distribution function (cdf).
%   P = GNORMCDF(X,MU,ALFA,BETA) returns the cdf of the generalized normal 
%   distribution with mean MU, standard deviation ALFA, and kurtosis 
%   parameter BETA, evaluated at the values in X.
%   The size of P is the common size of X, MU, ALFA and BETA. A scalar 
%   input functions as a constant matrix of the same size as the other 
%   inputs.
%
%   References:

%   Copyright 2022, Jasper A. Vrugt

P = 1/2 + sign(x-mu) .* 1/(2*gamma(1/beta)).* ...
    gammainc(abs((x-mu)/alfa).^beta,1/beta) * gamma(1/beta);
% Note multiplication with gamma(1/beta) is needed because we need to use
% undernormalized gammainc (default gammainc is normalized)
