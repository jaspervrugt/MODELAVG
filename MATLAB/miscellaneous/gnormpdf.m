function [Y] = gnormpdf(x,mu,alfa,beta)
%GNORMPDF Generalized normal probability density function (pdf).
%   Y = GNORMPDF(X,MU,ALFA,BETA) returns the pdf of the generalized normal 
%   distribution with mean MU, standard deviation ALFA, and kurtosis BETA, 
%   evaluated at the values in X. The size of Y is the common size of the 
%   input arguments. A scalar input functions as a constant matrix of the 
%   same size as the other inputs.

%   References:

%   Copyright 2022, Jasper A. Vrugt

Y = beta./(2*alfa .* gamma(1./beta)) .* exp(- (abs(x - mu) ./ alfa).^beta);
