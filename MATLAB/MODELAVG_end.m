function [beta,output,str_plot] = MODELAVG_end(method,D,y,beta,chain, ...
    options,output)
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
%% Postprocess results MODELAVG toolbox [= performance BMA model,        %%
%% scoring rules, etc]                                                   %%
%%                                                                       %%
%% SYNOPSIS: [beta,output,str_plot] = MODELAVG_end(method,D,y,beta, ...  %%
%%               chain,options,output)                                   %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%                                                                       %%
%%   (c) Written by Jasper A. Vrugt, Feb 2012                            %%
%%   University of California Irvine                                     %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% Compute size of ensemble and number of significance levels
[n,K] = size(D); na = numel(options.alpha);

% open an output file with warnings
fid = fopen('warning_file.txt','a+','n');

% Derive estimate of beta from posterior values
if any(strcmp(method,{'bma','mma','mma-s'}))
    % First assemble all chain in one matrix
    beta_matrix = genparset(chain);
    % Now determine size of beta_matrix and number of parameters d
    [N,d_all] = size(beta_matrix); d = d_all - 1;
    % Check that DREAM converged
    if ~(sum(output.R_stat(end,2:d+1) < 1.2 ) == d)
        % Print warning that DREAM did not converge
        evalstr = char(['MODELAVG WARNING: DREAM_{(ZS)} DID NOT CONVERGE ' ...
            'FORMALLY ACCORDING TO MULTIVARIATE R-STATISTIC\n']);
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
    end
    % Take the last 25% of the posterior samples -- assume that these samples
    % are posterior samples (double check that R_stat < 1.2 for all parameters)
    beta = beta_matrix ( floor ( 3/4 * N ) : N , 1 : d_all );
    % Find the maximum aposteriori parameter values (last column of ParSet are log-density values!)
    [~,idx] = max(beta(:,end)); beta_opt = beta(idx(1),1:d);
    % Calculate the posterior parameter correlation matrix (R-values)
    output.CORR = corr(beta(:,1:d),'type','Pearson');
    % And maximum log likelihood value
    output.loglik = beta ( idx(1) , end );
    % And std of DREAM
    output.std = std(beta(:,1:d));
else
    beta_opt = beta;
end

% Get the confidence level(s) from significance level(s)
g = 1 - options.alpha; gam = sort([ (1-g)/2 , (1-g)/2 + g ]);

switch method
    case 'bma'
        % Extract posterior draws of weights (supposedly)
        beta_weights = beta(:,1:K);
        % Calculate the gamma% posterior simulation uncertainty due to parameter uncertainty
        output.par_unc = prctile(beta_weights*D',100*gam)';
    otherwise
        % Compute variance of averaged model
        y_opt = D*beta_opt'; df = n - K;
        sigma2_avg = sum((y-y_opt).^2)/df; s_avg = sqrt(sigma2_avg);
        % covariance matrix of weights
        C = sigma2_avg * inv(D'*D); %#ok
        % critical t-values
        t_crit = tinv(gam,df);
        % now get parameter uncertainty
        for t = 1:n
            output.par_unc(t,:) = y_opt(t,1) + t_crit * ...
                sqrt(D(t,:) * C * D(t,:)');
            output.pred(t,:) = y_opt(t,1) + t_crit * s_avg;
        end
        %     % Derive standard deviation from first-order approximation of model
        %     beta_opt = beta; output.std = sqrt(diag(sum((y-D*beta_opt').^2)/(n-K)*inv(D'*D)))';
end

% Write final line of warning file
fprintf(fid,'----------- End of MODELAVG warning file ----------\n');
% Close the warning file
fclose(fid);
% Now print to screen or not (not on unix/linux)
%if ( ispc || ismac ), edit warning_file.txt, end

% Calculate point forecast, RMSE and R (Pearson correlation coefficient)
output.Ye = D * beta_opt(1:K)';
% Compute the mean of the data
Ym = mean(y);
% Compute the total sum of squares
SStot = sum((y - Ym).^2);
% Compute the residuals
e = y - output.Ye;
% Compute the residual sum of squares
SSres = sum(e.^2);
% Compute coefficient of determination
output.R2 = 1 - SSres/SStot;
% Compute RMSE
output.RMSE = sqrt( sum ( e.^2 ) / n );
% Pearson correlation coefficient
output.R = corr(output.Ye,y,'type','Pearson');
% compute KGE efficiency of weighted average (= mean) BMA forecast
output.KGE = compKGE(y,output.Ye);
% Maximum likelihood weights/parameter values
output.ML = beta_opt;

% Derive printing string for postprocessing
[str_table,str_plot] = create_str_plot(method,options,K);
% Store in output
output.str_table = str_table;

% Check whether we compute quantiles
if strcmp(options.postproc,'no')
    return
end

if strcmp(method,'bma') % Derive quantiles
    % For quantile calculation --> can assume that a's are zero and b's unity
    %    a = zeros(1,K); b = ones(1,K);
    % Now derive the prediction uncertainty ranges of BMA model
    [output.pred,output.pdf_Y,output.cdf_Y] = BMA_quantile ( beta_opt , ...
        D , y , options );
    % Compute scoring rules of BMA model
    output = BMA_score(beta_opt,D,y,options,output);
end

for i = 1:na            % Compute coverage and spread 
    output.coverage(1,i) = 100 * sum( ( output.pred(:,i) < y ) & ...
        ( y < output.pred(:,2*na - (i-1)) ) ) / n;
    output.spread(1,i) = mean(output.pred(:,2*na - (i-1)) - ...
        output.pred(:,i));
end

end
