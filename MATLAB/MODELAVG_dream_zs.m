function [chain,output] = MODELAVG_dream_zs(par,log_pdf,N,T,d,K, ...
    bound_check)
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
%% DiffeRential Evolution Adaptive Metropolis with sampling from a past  %%
%% archive using parallel direction and snooker updates. This so-called  %%
%% DREAM_(ZS) algorithm is implemented specifically for BMA and MMA      %%
%% model parameter estimation using maximum likelihood estimation        %%
%% [= Bayes with uniform prior distribution]                             %%
%%                                                                       %%
%% SYNOPSIS: [chain,output] = MODELAVG_dream_zs(par,log_pdf,N,T,d,K, ... %%
%%               bound_check)                                            %%
%%                                                                       %%
%% ALGORITHM HAS BEEN DESCRIBED IN                                       %%
%%   Vrugt, J.A., C.G.H. Diks, and M.P. Clark (2008), Ensemble Bayesian  %%
%%       model averaging using Markov chain Monte Carlo sampling,        %%
%%       Environmental Fluid Mechanics, 8(5-6), 579-595,                 %%
%%           https://doi.org/10.1007/s10652-008-9106-3                   %%
%% FOR MORE INFORMATION PLEASE RAAD                                      %%
%%   Vrugt, J.A. (2016), Markov chain Monte Carlo simulation using the   %%
%%       DREAM software package: Theory, concepts, and MATLAB            %%
%%       implementation, Environmental Modeling and Software, 75,        %%
%%       pp. 273-316                                                     %%
%%   Diks, C.G.H, and J.A. Vrugt (2010), Comparison of point forecast    %%
%%       accuracy of model averaging methods in hydrologic applications, %%
%%       Stochastic Environmental Research and Risk Assessment, 24(6),   %%
%%       809-820, https://doi.org/10.1007/s00477-010-0378-z              %%
%%   Vrugt, J.A., C.J.F. ter Braak, M.P. Clark, J.M. Hyman and B.A.      %%
%%       Robinson (2008), Treatment of input uncertainty in hydrologic   %%
%%       modeling: Doing hydrology backward with Markov chain Monte      %%
%%       Carlo simulation, 44 (12), https://doi.org/10.1029/2007WR006720 %%
%%   Vrugt, J.A., and B.A. Robinson (2007), Treatment of uncertainty     %%
%%       using ensemble methods: Comparison of sequential data           %%
%%       assimilation and Bayesian model averaging, Water Resources      %%
%%       Research, 43, W01411, https://doi.org/10.1029/2005WR004838      %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%                                                                       %%
%%   (c) Written by Jasper A. Vrugt, Feb 2012                            %%
%%   University of California Irvine                                     %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

[Delta,c,c_star,n_CR,p_ug] = deal(3,5e-2,1e-12,3,0.2);  % default values shared DREAM parameters
[k,p_s,m0] = deal(10,0.1,max(N,20*d)); m = m0;          % default values DREAM_ZS algorithmic variables
chain = nan(T,d+1,N);                                   % preallocate memory for Markov chains
gamma_RWM = @(a,b) 2.38/sqrt(2*a*b);                    % function handle with RWM-based jump rate

count = 0; ct = 2; accept = 0; t_old = 1; steps = floor(T/50);
[MR_stat,AR] = deal ( nan(floor(T/steps)-1,2) ); 
R_stat = nan(floor(T/steps)-1,d+1);
R_stat(1) = 1; MR_stat(1) = 1; AR(1) = 1;
par.min = repmat(par.min,N,1); par.max = repmat(par.max,N,1);   % Copy N times min and max values
%par.max(:,K+1:d) = Inf;                                        % Overwrite upper value of non-weights -> Inf

% DREAM_ZS: INITIALIZATION
p_CR = ones(1,n_CR)/n_CR;                               % crossover: initial n_CR selection probabilities
[J,n_id] = deal(ones(1,n_CR));                          % crossover: jump distance and frequency of use
Z = LH_sampling(par.min(1,:),par.max(1,:),m0+N);        % Create initial population
% Boundary checking only applies to weights (if on simplex)
% And to lower bound of BMA estimates standard deviation
% By not checking upper bound of standard deviations we can initialize in
% smaller space; higher efficiency (MAY 31: 2018 -> CHANGED TO UPPER BOUND INITIALIZATION )
% --------------------------------------------------------------------------------------------
if bound_check
    Z(1:m0+N,1:K) = bsxfun(@rdivide,Z(1:m0+N,1:K),sum(Z(1:m0+N,1:K),2)); % Weights add up to one
end
% --------------------------------------------------------------------------------------------

for i = m0+1:m0+N
    Z(i,d+1) = log_pdf(Z(i,1:d));                       % compute log-likelihood initial chain states
end
X = Z(m0+1:m0+N,1:d+1);                                 % extract initial states from archive
std_Z = std(Z(1:m0,1:d));                               % compute std of each dimension in archive

chain(1,1:d+1,1:N) = reshape(X',1,d+1,N);               % store in chain the N initial states
%output.p_CR(1,1:n_CR+1) = [ 1 p_CR ];                   % store initial p_CR;
flag = 0; Xp = nan(N,d);                                % initialize flag and candidate population

% DYNAMIC PART (EVOLVE N CHAINS T DIFFERENT STEPS)
for t = 2:T
    method = randsample([1 2],1,'true',[1 - p_s p_s]);          % jump method in current generation
    eps = c_star * randn(N,d);                                  % calculate the ergodicity perturbation
    lambda = 1 + c * (2 * rand(N,d) - 1);                       % draw lambda ~ U[-c,c] + add 1
    delta = randsample(Delta,1);                                % draw delta from [1,...,Delta]
    R = reshape(randsample(m,3*delta*N,'true'),3*delta,N);      % sample 3*delta*N integers from [1,...,m]
    dX = zeros(N,d);                                            % set N jump vectors to zero
    switch method
        case 1                                                      % PARALLEL DIRECTION MOVE
            id = randsample(1:n_CR,N,'true',p_CR);                      % select N (crossover) indexes
            gamma = randsample([0 1],1,'true',[1 - p_ug p_ug]);         % gamma equal to 0 or 1
            U = rand(N,d);                                              % random labels between 0 and 1
            for i = 1:N                                                 % CREATE PROPOSAL IN EACH CHAIN
                a = R(1:delta,i); b = R(delta+1:2*delta,i);                 % unpack R in vectors a and b
                if ( gamma == 0 )                                           % true: subspace sampling
                    A = find( U(i,1:d) < id(i)/n_CR );                      % subset A with dimensions to sample
                    d_star = numel(A);                                      % cardinality of A
                    if d_star == 0, [~,A] = min(U(i,1:d)); d_star = 1; end  % subset A must have one element
                    gamma_d = 1/2*gamma_RWM(delta,d_star);                  % compute default jump rate
                else
                    A = 1:d; gamma_d = 1; a = a(1); b = b(1); id(i) = n_CR; % all dimensions update
                end
                dX(i,A) = lambda(i,A).*gamma_d.*sum(Z(a,A)-Z(b,A),1) + ...
                    eps(i,A);                                               % compute jump vector
            end
            Xp(1:N,1:d) = X(1:N,1:d) + dX(1:N,1:d);                     % compute candidate point ith chain
            log_alfa_J = zeros(N,1);                                    % symmetric jump distribution
        case 2                                                      % SNOOKER MOVE
            a = R(1,1:N); b = R(2,1:N); c_sn = R(3,1:N);                % three members of archive suffices
            id = n_CR*ones(1,N);                                        % per definition: full crossover
            for i = 1:N                                                 % CREATE PROPOSAL IN EACH CHAIN     
                F = X(i,1:d)-Z(c_sn(i),1:d); FF = max(F*F',realmin);        % projection: X(i,1:d)-Z(c_sn,1:d)
                zP = F * ( sum((Z(a(i),1:d)-Z(b(i),1:d)) .*F ) / FF );      % project Z(a,:)-Z(b,:) onto F
                gamma_s = 1.2 + rand;                                       % sample jump rate, g ~ U[1.2,2.2]
                dX(i,1:d) = lambda(i,1:d).*gamma_s.*zP + eps(i,1:d);        % compute jump vector
            end
            Xp(1:N,1:d) = X(1:N,1:d) + dX(1:N,1:d);                         % compute candidate point ith chain
            log_alfa_J = (d-1)/2*log(sum((Xp(1:N,1:d)-Z(c_sn,1:d)).^2,2)) - ...
                (d-1)/2*log(sum((X(1:N,1:d)-Z(c_sn,1:d)).^2,2));        % non-symmetry snooker correction
    end
    % --------------------------------------------------------------------------------------------
    if bound_check, Xp = bnd_check(Xp(1:N,1:d),X(1:N,1:d),K,par); end       % boundary constraints
    % --------------------------------------------------------------------------------------------
    for i = 1:N                                                         % usually CPU-intensive --> "parfor"
        Xp(i,d+1) = log_pdf(Xp(i,1:d));                                 % evaluate log(likelihood) of each proposal
    end
    alfa_L = exp(Xp(1:N,d+1) - X(1:N,d+1));                             % likelihood ratio of proposals
    p_acc = exp(log_alfa_J) .* alfa_L;                                  % acceptance probability proposals
    u = rand(N,1);                                                      % Draw N labels from U[0,1]
    for i = 1:N
        if p_acc(i) > u(i)                                              
            X(i,1:d+1) = Xp(i,1:d+1); accept = accept + 1;              % accept proposal if p_acc > U[0,1]
        else
            dX(i,1:d) = 0;                                              % otherwise, old state; jump is zero
        end
    end
    chain(t,1:d+1,1:N) = reshape(X',1,d+1,N);                           % add current states to chain
    if ( t < T/10 )
        for i = 1:N
            J(id(i)) = J(id(i)) + sum((dX(i,1:d)./std_Z).^2);           % update normalized squared jump distance
            n_id(id(i)) = n_id(id(i)) + 1;                              % update use of crossover frequency 
        end
    elseif ( t >= T/10 ) && flag == 0                                   % at end of burn-in period, do one time
        p_CR = J./n_id; p_CR = p_CR/sum(p_CR);                          % update selection probabilities crossover values
        flag = 1;                                                       % change flag
    end
    if (mod(t,k) == 0)
        Z(m+1:m+N,1:d+1) = X; m = m + N;                                % each k iterations append X to archive Z
        std_Z = std(Z(m0+1:m,1:d));                                     % update standard deviation of each target variable
    end
    if mod(t,steps) == 0
        % Print progress
        if ( t > 1 )
            fprintf(1, repmat('\b',1,count)); % deletes line before
            count = fprintf('DREAM_{(ZS)} progress, %% done: %3.2f',100 * ( t/T ) );
        end
        % Calculate convergence diagnostics
        [GR_uniR,GR_multiR] = MODELAVG_gelman(chain(floor(t/2) + 1 : t,1:d,1:N),t); 
        R_stat(ct,1:d+1) = [ t GR_uniR ] ; MR_stat(ct,1:2) = [ t GR_multiR ] ;
        AR(ct,1:2) = [ t 100 * ( accept/( N*(t - t_old) ) ) ]; accept = 0; t_old = t;
        ct = ct + 1;
    end
end

fprintf(1,'\n');
% Store in output, GR_uniR, GR_multiR and AR
output = struct('AR',AR,'R_stat',R_stat,'MR_stat',MR_stat);

end
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
% Secondary functions
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

% 1: bnd_check
function Xp = bnd_check(Xp,X,K,par)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %
%bnd_check: Function checks the parameter boundaries
%
% Written by Jasper A. Vrugt
% Los Alamos National Laboratory
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %

method = 'reflection'; 
% Which parameters are out of bound?
ii_low = find(Xp < par.min); ii_up = find(Xp > par.max);
% Now we modify these parameter values 
switch method
    case 'set_to_bound' % set to bound - violates detailed balance
        Xp(ii_low) = par.min(ii_low); Xp(ii_up) = par.max(ii_up); 
    case 'reflection'   % reflect in min / max - violates detailed balance
        Xp(ii_low) = 2 * par.min(ii_low) - Xp(ii_low); 
        Xp(ii_up) = 2 * par.max(ii_up) - Xp(ii_up);
        % --> used in 2024 scoring rules paper
    case 'folding'      % folding - satisfies detailed balance
        Xp(ii_low) = par.max(ii_low) - (par.min(ii_low) - Xp(ii_low));
        Xp(ii_up) = par.min(ii_up) + (Xp(ii_up) - par.max(ii_up));
    case 'alternative'  % at random between old X and bound - violates detailed balance
        Xp(ii_low) = unifrnd(par.min(ii_low),X(ii_low));
        Xp(ii_up) = unifrnd(X(ii_up),par.max(ii_up));        
%     case 'huang2024comment'
%         Xp(ii_low) = par.min(ii_low) + rand(size(ii_low))/4 .* ..
%             (X(ii_low) - par.min(ii_low)); 
%         Xp(ii_up) = par.max(ii_up) - rand(size(ii_up))/4 .* ...
%             (par.max(ii_up) - X(ii_up) ); 
end
% now correct one final time as reflection and folding can go out of bound
ii_low = find(Xp < par.min); ii_up = find(Xp > par.max);
Xp(ii_low) = par.min(ii_low); Xp(ii_up) = par.max(ii_up);

% Weights must add up to one
Xp(:,1:K) = bsxfun(@rdivide,Xp(:,1:K),sum(Xp(:,1:K),2));

end
%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
