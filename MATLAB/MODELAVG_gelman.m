function [hatR,hatRd] = MODELAVG_gelman(chain,t)
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
%% Calculates the \hat{R}-convergence diagnostics                        %%
%%                                                                       %%
%% SYNOPSIS: [hatR,hatRd] = MODELAVG_gelman(chain,t)                     %%
%%                                                                       %%
%% FOR MORE INFORMATION PLEASE RAAD                                      %%
%%   Gelman, A. and D.R. Rubin, (1992) Inference from iterative          %%
%%       simulation using multiple chain, Statistical Science, 7 (4),    %%
%%       457-472                                                         %%
%%   Brooks, S.P. and A. Gelman, (1998) General methods for monitoring   %%
%%       convergence of iterative simulations, Journal of Computational  %%
%%       and Graphical Statistics, 7, 434-455                            %%
%%   Note that this function is the square-root definition of R          %%
%%       (Gelman et al. 2003, Bayesian Data Analsyis, p. 297)            %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%                                                                       %%
%%   (c) Written by Jasper A. Vrugt, Feb 2012                            %%
%%   University of California Irvine                                     %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% Compute the dimensions of chain
[n,d,N] = size(chain); warning('off');

if (n < 10)
    % Set the R-statistic to a large value
    [hatR,hatRd] = deal(NaN(1,d),NaN);
else
    % STEP 0: Compute the chain means and store in a N x d matrix
    mu_chains = mean(chain); mu_chains = reshape(mu_chains,d,N)';
    % STEP 0: Compute the N within-chain variances
    s2_chains = nan(N,d);
    for ii = 1:N, s2_chains(ii,1:d) = var(chain(1:n,1:d,ii)); end
    % STEP 0: Compute the N within-chain covariances
    cov_chains = nan(d,d,N);
    for ii = 1:N, cov_chains(1:d,1:d,ii) = cov(chain(1:n,1:d,ii)); end

    %% Univariate \hat{R} diagnostics
    % STEP 1: Compute variance B of N chain means
    B = n * var(mu_chains); 
    % STEP 2: Compute 1xd vector W with mean of within-chain variances
    W = mean(s2_chains);
    % STEP 3: Estimate target variance (= Σ within- and between-chain s2's)
    sigma2 = ((n-1)/n) * W + (1/n) * B;
    % STEP 4: Compute univariate \hat{R} diagnostic for each parameter
    hatR = sqrt( (N+1)/N * (sigma2./W) - (n-1)/(N*n) );
    % For large n --> chains means become more similar, B --> 0
    % Then, (sigma2./W) --> (n-1)/n ~ 1, and 
    % sqrt( (N+1)/N * ~1 - (n-1)/(N*n) ) --> 1
    
    %% Multivariate \hat{R}^d diagnostic
    % STEP 1: Compute dxd matrix with mean W of within-chain covariances
    W = mean(cov_chains,3) + eps * eye(d);
    % STEP 2: Compute covariance B of N chain means
    B = cov(mu_chains) + eps * eye(d); 
    % STEP 3: Compute multivariate scale reduction factor, \hat{R}^d
    hatRd = sqrt( (N+1)/N * max(abs(eig(W\B))) + (n-1)/n );
    % hatRd = sqrt( (N+1)/N * max(abs(eig((W + eps*eye(d,d)*randn(d,d))\B))) + (n-1)/n )
    % For large n --> chains means become more similar, B --> 0
    % Then, (W\B) --> 0 and  
    % sqrt( (N+1)/N * ~ 0 + (n-1)/n ) --> 1
    
    % Now calculate the multivariate variant of this statistic
    % MR_stat = mpsrf_brooks(chain);
    
    % --------- Now check whether to write warning to warning_file ------
    [msgstr,msgid] = lastwarn;
    % if msgstr not empty
    if ~isempty(msgstr)
        % open file
        fid = fopen('warning_file.txt','a+');
        % now write to warning_file.txt
        evalstr = char(strcat('MODELAVG WARNING:',{' '},msgid,...
            {' '},'in multivariate R_statistic of Brooks and Gelman at',...
            {' '},num2str(t),{' '},'generations\n'));
        % Now print warning to file
        fprintf(fid,evalstr);
        % Close file
        fclose(fid);
        % Now remove warning
        lastwarn('');
    end
    % --------- End check whether to write warning to warning_file ------
end

end
