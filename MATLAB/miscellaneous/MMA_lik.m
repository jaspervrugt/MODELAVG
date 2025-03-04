function ell = MMA_lik ( beta , D , Y , var_err , p )
% This function calculates the log likelihood of MMA 
% Function of MODELAVG toolbox, V1.0

% B.C. Hansen, "Least Squares Model Averaging", Econometrica, vol. 75, 
%   no. 4, pp. 1175-1189, 2007

if nargin<5
    error('MODELAVG:MMA_calc:TooFewInputs', ...
        'Requires at least five input arguments.');
end

G = D*beta';                                % MMA deterministic point forecast
Cn = sum((Y - G).^2) + 2*var_err*beta*p;    % Mallows criterion: ref Equation (11) 
ell = -1/2 * Cn + realmin;                  % Log-likelihood of x = { MMA weights }