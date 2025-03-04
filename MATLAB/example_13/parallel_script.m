function [output,beta,AR,table_res_cal,val,table_res_eval] = ...
    parallel_script(ii,VPDF,method,D_cal,y_cal,options,which_period,a_LB,b_LB,D_eval,y_eval)

%% Now do BMA method for all combinations of variances/conditional PDFs
options.VAR = char(VPDF(1));       % Group/sole variance
% Define conditional pdf
options.PDF = char(VPDF(2));
% Run MODELAVG toolbox with two outputs
[ beta , output ] = MODELAVG ( method , D_cal , y_cal , options ); AR = output.AR;
% Store information of calibration period
table_res_cal = [ output.mQS output.mLS output.mSS output.mCRPS ...
    output.mES output.mRLBL_anal output.mCv_anal ...
    output.coverage(1)/100 output.spread(1) output.loglik ...
    output.RMSE output.R2 output.KGE ]';
% Compute information of evaluation period?
switch which_period
    case 'evaluation'
        % Evaluate maximum likelihood BMA forecast distribution for evaluation period
        [val,D_eval] = MODELAVG_eval ( method , D_eval , y_eval , options , a_LB , b_LB , output );
        % Return table with results
        table_res_eval = [ val.mQS val.mLS val.mSS val.mCRPS val.mES ...
            val.mRLBL_anal val.mCv_anal ...
            val.coverage(1)/100 val.spread(1) val.loglik ...
            val.RMSE val.R2 val.KGE ]';
    otherwise
        val = []; table_res_eval = ii;
end
