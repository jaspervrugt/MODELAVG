function [Y] = tnormpdf(x,mu,sigma,a,b)
%TNORMPDF Truncated normal probability density function (pdf).
%   Y = TNORMPDF(X,MU,SIGMA,A,B) returns the pdf of the truncated normal 
%   distribution with mean MU, standard deviation SIGMA, lower and upper 
%   bounds of the interval, A and B, evaluated at the values in X. The 
%   size of Y is the common size of the input arguments. A scalar input 
%   functions as a constant matrix of the same size as the other inputs.

%   References:

%   Copyright 2022, Jasper A. Vrugt

f = @(x) 1/sqrt(2*pi) * exp(-1/2*x.^2);
g = @(x) 1/2*(1 + erf(x/sqrt(2)));

z = (x - mu)./sigma; a_n = (a - mu)./sigma; b_n = (b - mu)./sigma;

Y = 1./sigma .* f(z) ./ (g(b_n) - g(a_n));

id = x < a | x > b; Y(id) = 0;