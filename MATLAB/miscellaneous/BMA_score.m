function output = BMA_score(x,D,y,options,output)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
% This function calculates scoring rules for BMA forecast distribution    %
%                                                                         %
% SYNOPSIS: output = BMA_score(x,D,y,options,output)                      %
%  where                                                                  %
%   x         [input] 1 x d vector maximum likelihood BMA parameters      %
%   D         [input] n x K matrix forecasts of ensemble members          %
%   y         [input] n x 1 vector verifying observations                 %
%   options   [input] structure with BMA algorithmic variables            %
%   output    [input] structure with the results of MODELAVG code         %
%   output    [outpt] Scoring rules, perf. metrics BMA mixture distr.     %
%                                                                         %
% (c) Written by Jasper A. Vrugt, Feb 2022                                %
% University of California Irvine                                         %
%                                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

% Draw m samples BMA mixture PDF (= Phi) + return exact µ/var. BMA mixture
m = 1000; [Phi,mu_mix,var_mix] = BMA_draw(x,D,m,options);
% Compute CRPS - quantile integration using BMA mixture CDF
output.CRPS = BMA_crps(x,D,y,options);
% Compute α-norm of BMA mixture distribution
warning off; [lowup,mix_norm] = BMA_norm(x,D,y,options,'robust');     %#ok

% Compute some scoring rules - numerically
nY = size(D,1);
output.QS = 2*output.pdf_Y - mix_norm(1:nY,2).^2;
output.LS = log2(output.pdf_Y);
output.SS = output.pdf_Y ./ mix_norm(1:nY,2);
% Compute Energy Score
output.ES = nan(nY,1); count = 0;
for t = 1:nY
    % Print progress
    if mod(t,floor(nY/25)) == 0
        if t > 1
            fprintf(1, repmat('\b',1,count)); %delete line before
            count = fprintf('Energy score calculation, %% done: %3.2f',...
                100*(t/nY));
        end
    end
    % Compute Energy score
    output.ES(t,1) = energy_score(Phi(t,1:m),y(t),2);
end
fprintf('\n');

% Compute exact reliability for BMA model
output.mRLBL_anal = BMA_rlbl(output.cdf_Y);
% Compute coefficient of variation exactly
Cv_anal = sqrt(var_mix)./mu_mix;
output.mCv_anal = mean(Cv_anal);
% Compute time-averaged scores
id1 = 1:nY;
output.mQS = mean(output.QS(id1));
output.mLS = mean(output.LS(id1));
output.mSS = mean(output.SS(id1));
output.mCRPS = mean(output.CRPS(id1));
output.mES = mean(output.ES(id1));
output.mix_norm = mix_norm;
output.mu_mix = mu_mix;
output.var_mix = var_mix;

end