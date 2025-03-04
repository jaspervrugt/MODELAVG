function [P] = tnormcdf(x,mu,sigma,a,b)
%TNORMCDF Truncated normal cumulative distribution function (cdf).
%   P = TNORMCDF(X,MU,SIGMA,A,B) returns the cdf of the truncated normal
%   distribution with mean MU, standard deviation SIGMA, and lower and 
%   upper bounds A and B, evaluated at the values in X.
%   The size of P is the common size of X, MU, SIGMA, A and B. A scalar
%   input functions as a constant matrix of the same size as the other
%   inputs.
%
%   References:

%   Copyright 2022, Jasper A. Vrugt

% Normalize input variables
x_n = (x-mu)./sigma; alfa = (a-mu)./sigma; beta = (b-mu)./sigma;
% Define cumulative distribution function
phi = @(x) 1/2*(1 + erf(x/sqrt(2)));
% Denominator
Z = phi(beta) - phi(alfa); 
% Compute CDF
P = (phi(x_n) - phi(alfa))./Z;

% Set to zero or one outside this interval
id = x < a; P(id) = 0; id = x > b; P(id) = 1;