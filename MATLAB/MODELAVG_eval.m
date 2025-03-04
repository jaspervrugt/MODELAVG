function [val,D_eval] = MODELAVG_eval(method,D_eval,y_eval, ...
    options,a,b,output)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%                                                                       %%
%% MM      MM  OOOOOO  DDDDDDD  EEEEEEE LL        AAAA   VV   VV GGGGGGG %%
%% MMM     MM OOOOOOOO DDDDDDDD EEEEEEE LL       AA  AA  VV   VV GG   GG %%
%% MMMM  MMMM OO    OO DD    DD EE      LL       AA  AA  VV   VV GG   GG %%
%% MM MMMM MM OO    OO DD    DD EEEE    LL       AA  AA  VV   VV GGGGGGG %%
%% MM  MM  MM OO    OO DD    DD EEEE    LL       AAAAAA   VV VV  GGGGGGG %%
%% MM      MM OO    OO DD    DD EE      LL      AA    AA   VVV        GG %%
%% MM      MM OOOOOOOO DDDDDDDD EEEEEEE LLLLLLL AA    AA   VVV       GGG %%
%% MM      MM  OOOOOO  DDDDDDD  EEEEEEE LLLLLLL AA    AA    V    GGGGGGG %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Evaluates performance for evaluation period                           %%
%%                                                                       %%
%% SYNOPSIS: [val,D_eval] = MODELAVG_eval(method,D_eval,y_eval, ...      %%
%%               options,a,b,output)                                     %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%                                                                       %%
%%   (c) Written by Jasper A. Vrugt, Feb 2012                            %%
%%   University of California Irvine                                     %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% Bias correction evaluation ensemble using intercept "a" and slope "b"
[m,K] = size(D_eval);
for k = 1:K
    D_eval(1:m,k) = a(k) + b(k) * D_eval(1:m,k);
end
evalstr = strcat('----- Evaluate',{' '},upper(method),{' '},...
    'for the evaluation period ----- \n'); fprintf(char(evalstr));

% Now check setup -> make sure that structure options is correct
[method,options] = MODELAVG_check(method,D_eval,y_eval,options);

% Now run setup function of toolbox
[~,~,~,~,~,~,~,~,options] = MODELAVG_setup(method,D_eval,y_eval, ...
    options,1);

% Extract maximum likelihood values weights, etc.
beta_opt = output.ML;

if strcmp(method,'bma')  % Derive quantiles and scoring rules
    % Now derive the prediction uncertainty ranges of BMA model
    [val.pred,val.pdf_Y,val.cdf_Y] = BMA_quantile(beta_opt,D_eval, ...
        y_eval,options);
    % cpmpute loglikelihood of evaluation period
    val.loglik = sum(log(val.pdf_Y));
    % number of alpha values
    na = numel(options.alpha);
    % And now the % observations enclosed - for each alpha value
    for i = 1:na
        val.coverage(1,i) = 100 * ( sum( ( val.pred(:,i) < y_eval ) & ...
            ( y_eval < val.pred(:,2*na - (i-1)) ) ) / m );
        % And now compute width of the largest prediction intervals of alpha
        val.spread(1,i) = mean(val.pred(:,2*na - (i-1)) - val.pred(:,i));
    end
    % Compute scoring rules of BMA distribution forecasts
    val = BMA_score(beta_opt,D_eval,y_eval,options,val);
end

% Calculate point forecast, RMSE and R (Pearson correlation coefficient)
val.Ye = D_eval * beta_opt(1:K)';
% Compute the mean of the data
Ym_eval = mean(y_eval);
% Compute the total sum of squares
SStot_eval = sum((y_eval - Ym_eval).^2);
% Compute the residuals
e_eval = y_eval - val.Ye;
% Compute the residual sum of squares
SSres_eval = sum(e_eval.^2);
% Compute coefficient of determination
val.R2 = 1 - SSres_eval/SStot_eval;
% Compute RMSE
val.RMSE = sqrt(sum(e_eval.^2)/m); % sqrt(sigma2)?
% Pearson correlation coefficient
val.R = corr(val.Ye,y_eval,'type','Pearson');
% Compute KGE of weighted average BMA forecast
val.KGE = compKGE(y_eval,val.Ye);

end
