load results_tao.mat

% Only look at normal and gamma PDFs
idx_pdf = [1 5 8 12]; n_idx = numel(idx_pdf);
var_pdf = [1 1 2 2]; 
cond_pdf = {'normal','gamma','normal','gamma'};
M = 10;

for i = 1:n_idx
    pars = beta{idx_pdf(i)}; d = size(pars,2) - 1;
    idx_par = find(pars(:,d+1)==max(pars(:,d+1))); 
    ML_DREAM{i} = pars(idx_par(1),1:d+1);
    % Now run EM algorithm - lets do M trials and pick best
    options.VAR = num2str(var_pdf(i)); options.PDF = char(cond_pdf(i));
    ML_EM{i}(1,1:d+1) = [nan(1,d) -inf];
    for z = 1:M
        [EM_beta , EM_sigma , EM_loglik , t] = EM_bma(D_cal,y_cal,options);
        if EM_loglik > ML_EM{i}(1,d+1)
            ML_EM{i} = [EM_beta , EM_sigma , EM_loglik];
        end
    end
end
