function [Y] = powpdf(x,mu,sigma,alfa)
%GNORMPDF Lognormal power-law probability density function (pdf).
%   Y = POWPDF(X,MU,SIGMA,ALFA) returns the pdf of the lognormal power-law 
%   distribution with mean MU, standard deviation SIGMA, and shape ALFA, 
%   evaluated at the values in X. The size of Y is the common size of the 
%   input arguments. A scalar input functions as a constant matrix of the 
%   same size as the other inputs.

%   References:

%   Copyright 2022, Jasper A. Vrugt

% Negative data would create complex values, potentially creating spurious
% NaNi's in other elements of y.  Map them, and zeros, to the far right
% tail, whose pdf will be zero.
x(x <= 0) = Inf;

Y = alfa/2 * exp(alfa * mu + alfa^2 * sigma^2/2) * x.^(-1-alfa) .* ...
    erfc(1/sqrt(2) * (alfa*sigma - (log(x) - mu)/sigma));