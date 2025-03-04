function MODELAVG_postproc(method,D,y,beta,sigma_2,par,chain, ...
    options,output,str_plot)
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
%% Visualize results of MODELAVG toolbox: marginal distribution of       %%
%% parameters, chain convergence, autocorrelation functions, confidence/ %%
%% prediction limits of BMA model confidence intervals of other model    %%
%% averaging methods, etc.                                               %%
%%                                                                       %%
%% SYNOPSIS: MODELAVG_postproc(method,D,y,beta,sigma_2,par,chain, ...    %%
%%               options,output,str_plot)                                %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%                                                                       %%
%%   (c) Written by Jasper A. Vrugt, Feb 2012                            %%
%%   University of California Irvine                                     %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

set(0,'defaultTextInterpreter','latex'); %trying to set the default

% Determine the size of the ensemble forecasts
[N_meas,K] = size(D);

if K <= 8
    colororder =  [ 0   0   1   % 1 BLUE
        0   1   1                   % 4 CYAN
        0   .5  0                   % 9 GREEN (dark)
        1   .5  .25                 % 11 ORANGE
        1   0   1                   % 5 MAGENTA (pale)
        .6  .5  .4                  % 15 BROWN (dark)
        0   1   0                   % 2 GREEN (pale)
        0   0   0];                 % 7 BLACK
    %    0.8 0.7 0.6                % 15 BROWN (pale)
else
    clrs = parula; nc = size(clrs,1);
    colororder = clrs(1:floor((nc-1)/K):nc,1:3);
end
n_colors = size(colororder,1);
maxbins = 25;

% Extract confidence levels
g = 1 - options.alpha; gam = sort([ (1-g)/2 , (1-g)/2 + g ]);

% Define model complexity
if ~isfield(options,'p')
    options.p = zeros(K,1);
end

% -------------------------------------------------------------------------
% ---------------------------- PRE-PROCESSING -----------------------------
% -------------------------------------------------------------------------

% Derive overall min and overall max
max_Y = max( [ max(y) max(D) ] ); min_Y = min( [ min(y) min(D) ] );

% How many parameters?
[~,n_pars] = size(beta);

% Now check per method
if any(strcmp(method,{'mma','mma-s','bma'}))
    % Remove last column (log-likelihood values)
    n_pars = n_pars - 1;
end

% -------------------------------------------------------------------------
% ----------- WRITE TABLE WITH RESULTS OF MODELAVG TO SCREEN --------------
% -------------------------------------------------------------------------

% Get the right string
str = output.str_table;

% open an output file with results of formulation
fid = fopen('MODELAVG_output.txt','w');

% Now switch among methods
if strcmp(method,'bma')
    fprintf(fid,['TABLE 1: MODEL AVERAGING METHOD: %s WITH %s ' ...
        'CONDITIONAL PDF (MANUAL FOR DETAILS) \n'],upper(method), ...
        upper(options.PDF));
else
    fprintf(fid,['TABLE 1: MODEL AVERAGING METHOD: %s (MANUAL FOR ' ...
        'DETAILS)\n'],upper(method));
end
fprintf(fid,'======================================   ====================\n');
if any(strcmp(method,{'mma','mma-s','bma','gra'}))
    fprintf(fid,'MODEL  PARAMETER     OPT       STD*       COMPLEXITY   RMSE\n');
else
    fprintf(fid,'MODEL  PARAMETER     OPT       STD        COMPLEXITY   RMSE\n');
end
fprintf(fid,'--------------------------------------   --------------------\n');

fmt_1 = ('%3d \t %-7s %8.3f  %8.3f                 %8.3f\n');
fmt_2 = ('%3d \t %-7s %8.3f  %8.3f      %7u    %8.3f\n');
fmt_3 = ('%3d \t %-7s %8.3f                           %8.3f\n');
fmt_4 = ('%3d \t %-7s %8.3f                %7u    %8.3f\n');
% Now print
for i = 1 : K
    if any(strcmp(method,{'mma','mma-s','bma','gra'}))
        % print according to prespecified format
        if options.p(i) == 0
            fprintf(fid,fmt_1,i,char(str(i)),output.ML(i), ...
                output.std(i),sqrt(sigma_2(i)));
        else
            fprintf(fid,fmt_2,i,char(str(i)),output.ML(i), ...
                output.std(i),options.p(i),sqrt(sigma_2(i)));
        end
    else
        if options.p(i) == 0
            fprintf(fid,fmt_3,i,char(str(i)),output.ML(i),sqrt(sigma_2(i)));
        else
            fprintf(fid,fmt_4,i,char(str(i)),output.ML(i), ...
                options.p(i),sqrt(sigma_2(i)));
        end
    end
end

if any(strcmp(method,{'ewa','bga','aica','bica'}))
    fprintf(fid,'======================================   ====================\n');
    fprintf(fid,'\n');
elseif strcmp(method,'gra')
    fprintf(fid,'======================================   ====================\n');
    fprintf(fid,'*Derived from first-order approximation around listed weights\n');
    fprintf(fid,'\n');
elseif any(strcmp(method,{'mma','mma-s'}))
    fprintf(fid,'======================================   ====================\n');
    fprintf(fid,'*Derived from DREAM_ZS derived samples of the posterior weights\n');
    fprintf(fid,'\n');
else
    fprintf(fid,'--------------------------------------   --------------------\n');
    k = 1; % model number
    % Now check which parameters belong to the same model (May 2018)
    if any(strcmp(options.VAR,{'2','4'}))
        index(:,1) = K+1:2*K; index_max = 2*K;
    else
        index = zeros(K,1); index_max = K+1;
    end
    if sum(strcmp(options.PDF,{'gen_normal','gev','gpareto'}))
        if strcmp(options.TAU,'1')
            index(:,2) = zeros(K,1);
        else
            index(:,2) = index_max + 1: index_max + K;
        end
    end
    % Thus first row of index pertains to model 1, and so forth
    fmt_1 = ('%3s \t %-7s %8.3f  %8.3f      %90s\n');
    for i = K + 1 : n_pars
        pr = '';
        if ismember(i,index(k,:))
            if output.ML(k) < 1e-2
                pr = char(strcat(['NOTE: BMA density relatively ' ...
                    'insensitive to this variable as weight model'], ...
                    {' '},num2str(k),{' '},'close to zero'));
            end
        end
        fprintf(fid,fmt_1,' ---',char(str(i)),output.ML(i),output.std(i),pr);
        k = k + 1; % update model
        if k > K % reset to first model
            k = 1;
        end
    end
    fprintf(fid,'======================================   ====================\n');
    fprintf(fid,['*Posterior standard deviation from DREAM_ZS derived ' ...
        'samples of BMA parameters\n']);
    fprintf(fid,'\n');
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

fprintf(fid,'TABLE 2: AVERAGED MODEL PERFORMANCE: TRAINING DATA SET\n');
fprintf(fid,'======================================\n');
fprintf(fid,' METHOD       RMSE      R\n');
fprintf(fid,'--------------------------------------\n');
fmt_1 = (' %5s\t  %8.3f %8.3f\n');
fprintf(fid,fmt_1,upper(method),output.RMSE,output.R);
fprintf(fid,'======================================\n');
fprintf(fid,'\n');

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

if strcmp(method,'bma')
    % number of alpha values
    alf_prct = -100*sort(-options.alpha); na = numel(alf_prct);
    fprintf(fid,'TABLE 3: BMA MODEL STATISTICS: TRAINING DATA SET\n');
    fprintf(fid,'=============================================================\n');
    fprintf(fid,' BMA CONDITIONAL PDF      Coverage   Spread   Log-likelihood  \n');
    fprintf(fid,'-------------------------------------------------------------\n');
    fmt_1 = ('%20s \t                       \t %8.3f\n');
    fmt_2 = ('%20d%% \t   %5.3f    %5.3f      \n');
    if strcmp(options.PDF,'gamma')
        fprintf(fid,fmt_1,'Gamma Distribution',output.loglik);
    end
    if strcmp(options.PDF,'normal')
        fprintf(fid,fmt_1,'Normal Distribution',output.loglik);
    end
    if strcmp(options.PDF,'gen_normal')
        fprintf(fid,fmt_1,'Generalized Normal',output.loglik);
    end
    if strcmp(options.PDF,'weibull')
        fprintf(fid,fmt_1,'Weibull',output.loglik);
    end
    if strcmp(options.PDF,'lognormal')
        fprintf(fid,fmt_1,'Lognormal',output.loglik);
    end
    if strcmp(options.PDF,'gev')
        fprintf(fid,fmt_1,'Generalized EV',output.loglik);
    end
    if strcmp(options.PDF,'gpareto')
        fprintf(fid,fmt_1,'Generalized Pareto',output.loglik);
    end
    for i = 1:na
        fprintf(fid,fmt_2,alf_prct(i),output.coverage(i),output.spread(i));
    end
    fprintf(fid,'=============================================================\n');
    fprintf(fid,'\n');
end

if any(strcmp(method,{'mma','mma-s'}))
    if strcmp(method,'mma')
        fprintf(fid,'TABLE 3: MMA MODEL STATISTICS: TRAINING DATA SET\n');
    else
        fprintf(fid,'TABLE 3: MMA-S MODEL STATISTICS: TRAINING DATA SET\n');
    end
    fprintf(fid,'=============================\n');
    if strcmp(method,'mma')
        fprintf(fid,' MMA     Log-likelihood  \n');
    else
        fprintf(fid,' MMA-S   Log-likelihood  \n');
    end
    fprintf(fid,'-----------------------------\n');
    fmt_1 = ('          %8.3f\n');
    fprintf(fid,fmt_1,output.loglik);
    fprintf(fid,'=============================\n');
    fprintf(fid,'\n');
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------

% Print correlation values
if any(strcmp(method,{'bma','mma','mma-s'}))
    fprintf(fid,'TABLE 4: CORRELATION COEFFICIENTS OF POSTERIOR SAMPLES\n');
    pr = '=========='; prt = pr;
    for i = 1:n_pars-1, prt = strcat(prt,pr); end, prt = strcat(prt,'========\n');
    % Now print line to screen
    fprintf(fid,prt);
    % print header
    fprintf(fid,'        %10s',char(str(1)));
    for i = 2:n_pars-1, fprintf(fid,'%10s',char(str(i))); end
    fprintf(fid,'%10s\n',char(str(n_pars)));
    % end print header
    fmt_1 = ('%-8s'); for i = 1:n_pars, fmt_1 = strcat(fmt_1,{' '},'%9.3f'); end
    fmt_1 = char(strcat(fmt_1,'\n'));
    % now print rows of table with content
    for i = 1:n_pars
        fprintf(fid,fmt_1,char(str(i)),output.CORR(i,1:n_pars));
    end
    % Now print line to screen
    fprintf(fid,prt);

end
fprintf(fid,'\n');
fprintf(fid,['------------------------ end MODELAVG output file -----' ...
    '--------------------\n']);
fclose(fid);

% Now print to screen or not (not on unix/linux)
% if ( ispc || ismac ), edit MODELAVG_output.txt, end

% Now plot figures if user wants screen output
if strcmp(options.print,'no')
    return
end

% -------------------------------------------------------------------------
% ---------- PLOT OF UNIVARIATE R-STATISTIC OF GELMAN AND RUBIN -----------
% -------------------------------------------------------------------------

% Check if R_statistic exists or not
if isfield ( output , 'R_stat' )
    % Now plot the R_statistic for each parameter
    figure('units','normalized','outerposition',[0 0 1 1])
    % First print the R-statistic of Gelman and Rubin (each parameter a different color)
    HA = plot(output.R_stat(:,1),output.R_stat(:,2:n_pars+1),...
        'linewidth',1.5); hold on;
    % Check colors
    if K < n_colors
        switch strcmp(method,'bma')
            case 1
                switch options.VAR
                    case {'1','3'} % single variance parameter
                        index = 1:K+1;
                    case {'2','4'} % multiple variance parameters
                        index = [ 1:K , 1:K ];
                end
                if strcmp(options.PDF,'gen_normal')
                    switch options.TAU
                        case '1'
                            index = [ index , max(index) + 1 ];
                        case '2'
                            index = [ index , 1:K ];
                    end
                end
            otherwise
                index = 1:K;
        end
        for i = 1 : n_pars
            set(HA(i),'color',colororder(index(i),1:3));
        end
    else
        % Keep MATLAB colors
    end
    % Add labels
    xlabel('Number of generations','fontsize',18,'interpreter','latex');
    ylabel('$$\hat{R}$$-diagostic','interpreter','latex',...
        'fontsize',18);
    % Add title
    evalstr = strcat('DREAM$$_{\rm (ZS)}$$:',{' '},...
        'Evolution of $$\hat{R}$$ convergence diagnostic for',...
        {' '},upper(method),{' '},'model');
    title(evalstr,'fontsize',18,'interpreter','latex');

    % Now add the theoretical convergence value of 1.2 as horizontal
    % line
    HA = plot([0 output.R_stat(end,1)],[1.2 1.2],'r--','linewidth',2);
    % Change color
    set(HA,'color',[0.5 0.5 0.5]);
    % Set the the axes
    axis([0 output.R_stat(end,1) 0.8 5]);
    % Add a legend
    [~,objh,~,~] = legend(str_plot,'interpreter','latex');
    set(objh,'linewidth',2); legend boxoff;
    try set(objh,'interpreter','latex','fontsize',14), catch, end
    % Fontsize of all numbers
    set(gca,'fontsize',16);
end

% -------------------------------------------------------------------------
% ----- EVOLUTION OF MULTIVARIATE R_STATISTIC OF GELMAN AND RUBIN ---------
% -------------------------------------------------------------------------

% Check if R_statistic exists or not
if isfield ( output , 'MR_stat' )
    % Now plot the R_statistic for each parameter
    figure('units','normalized','outerposition',[0 0 1 1])
    % First print the R-statistic of Gelman and Rubin (each parameter a different color)
    semilogy(output.MR_stat(:,1),output.MR_stat(:,2),'r',...
        'linewidth',1.5); hold on;
    % Add labels
    xlabel('Number of generations','fontsize',18,'interpreter','latex');
    ylabel('$$\hat{R}^{d}-\mbox{convergence diagnostic}$$',...
        'interpreter','latex','fontsize',18);
    % Add title
    evalstr = strcat('DREAM$$_{\rm (ZS)}$$:',{' '},...
        'Evolution of $$\hat{R}^{d}$$ convergence diagnostics for',...
        {' '},upper(method),{' '},'model');
    title(evalstr,'fontsize',18,'interpreter','latex');
    % Now add the theoretical convergence value of 1.2 as horizontal line
    HA = plot([0 output.MR_stat(end,1)],[1.2 1.2],'k--','linewidth',2);
    % Change color
    set(HA,'color',[0.5 0.5 0.5]);
    % Set the the axes
    axis([0 output.MR_stat(end,1) 0.8 20]);
    % Make font size 18
    set(gca,'fontsize',16);
end

% -------------------------------------------------------------------------
% --------------------- EVOLUTION OF ACCEPTANCE RATE ----------------------
% -------------------------------------------------------------------------

% Check if R_statistic exists or not
if isfield ( output , 'AR' )
    % Now plot the R_statistic for each parameter
    figure('units','normalized','outerposition',[0 0 1 1])
    % First print the R-statistic of Gelman and Rubin (each parameter a different color)
    plot(output.AR(:,1),output.AR(:,2),'r','linewidth',1.5); hold on;
    % Add labels
    xlabel('Number of generations','fontsize',18,'interpreter','latex');
    ylabel('Acceptance rate','fontsize',18,'interpreter','latex');
    % Add title
    evalstr = strcat('DREAM$$_{\rm (ZS)}$$:',{' '},...
        'Evolution of acceptance rate for',{' '},...
        upper(method),{' '},'method');
    title(evalstr,'fontsize',18,'interpreter','latex');
    % Now check the maximum value of acceptance rate that should be plotted
    y_max = min(100, 1.2 * max(output.AR(min(5,...
        size(output.AR,1)):end,2)) );
    % Now adjust ranges
    axis([0 output.AR(end,1) 0 y_max ]);
    % Make font size 18
    set(gca,'fontsize',16);
end

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
str = str_plot;

% -------------------------------------------------------------------------
% ------------ HISTOGRAMS OF MARGINAL DENSITIES OF PARAMETERS -------------
% -------------------------------------------------------------------------

% Check if R_statistic exists or not
if any(strcmp(method,{'bma','mma','mma-s'}))

    % Plot the histograms (marginal density) of each parameter;
    % What lay out of marginal distributions is desired subplot(r,t)
    r = 3; t = 3; counter = 1; j = 1;
    % Open new figure
    figure('units','normalized','outerposition',[0 0 1 1])
    % Check which colors to use?
    if K < n_colors
        switch strcmp(method,'bma')
            case 1
                switch options.VAR
                    case {'1','3'} % single variance parameter
                        index = 1:K+1;
                    case {'2','4'} % multiple variance parameters
                        index = [ 1:K , 1:K ];
                end
                if strcmp(options.PDF,'gen_normal')
                    switch options.TAU
                        case '1'
                            index = [ index , max(index) + 1 ];
                        case '2'
                            index = [ index , 1:K ];
                    end
                end
            otherwise
                index = 1:K;
        end
    else
        % Keep MATLAB colors
    end
    % Compute number of bins based on different rules
    Nbins = nan(1,n_pars);
    for i = 1:n_pars
        Nbins(i) = calcnbins(beta(:,i));
    end
    % Take the minimum of the number of bins for each parameter
    nbins = min(min(Nbins),maxbins);
    % Added 2024 as there are too many bins
    nbins = max(5,round(nbins/2));
    % End added

    % Now plot each parameter
    while counter <= n_pars
        % Check whether to open a new figure?
        if j == (r * t) + 1
            evalstr = strcat(['DREAM$$_{\rm (ZS)}$$ sampled marginal ' ...
                'distributions of'],...
                {' '},upper(method),{' '},'parameters');
            mtit(char(evalstr),'fontsize',18,'interpreter','latex');
            % Open new figure
            figure('units','normalized','outerposition',[0 0 1 1])
            % Reset j to 1
            j = 1;
        end
        % New implementation using histcounts
        [M,edges] = histcounts(beta(:,counter),nbins,'normalization',...
            'countdensity');
        X = 1/2*(edges(1:nbins)+edges(2:nbins+1)); % midpoint of each bin
        % And plot histogram in red
        subplot(r,t,j); HA = bar(X,M./max(M),'r'); hold on;

        %         % Now create histogram
        %         [N,X] = hist(beta(:,counter)); z = N/trapz(X,N); y_s = z./max(z);
        %         % And plot histogram in red
        %         subplot(r,t,j);
        %         HA = bar(X,y_s,'r'); hold on; % --> can be scaled to 1 if using "trapz(X,N)" instead of "sum(N)"!
        % Now check this
        if K < n_colors            % Enough colors
            set(HA,'facecolor',colororder(index(counter),1:3));
        else
            set(HA,'facecolor',[0.7 0.7 0.7]);
        end
        % Make font size 14
        set(gca,'fontsize',14);
        % Add label
        xlabel(str(counter),'fontsize',18,'interpreter','latex');
        % Then add y-label (only if j == 1 or j = r;
        if j == 1 || ( min(abs(j - ([1:r]*t+1))) == 0 )
            ylabel('Marginal density','fontsize',14,'interpreter','latex');
        end
        % Now determine the min and max X values of the plot
        minX = min(X); maxX = max(X);
        % Now determine appropriate scales
        deltaX = 0.1*(maxX - minX);
        % Calculate x_min and x_max
        x_min = minX - deltaX; x_max = maxX + deltaX;
        % Lets add the MAP value
        plot(output.ML(counter),1,'rx','Markersize',15,'linewidth',3);
        % Adjust the axis
        axis([x_min x_max 0 1.02]);
        % Now update the counter
        counter = counter + 1;

        % Update j
        j = j + 1;
    end
    % Create title of figure
    evalstr = strcat(['DREAM$$_{\rm (ZS)}$$ sampled marginal ' ...
        'distributions of'],...
        {' '},upper(method),{' '},'parameters');
    p = mtit(char(evalstr),'fontsize',18,'interpreter','latex');

end

% -------------------------------------------------------------------------
% ---------------- PLOT MAP WITH CORRELATION ESTIMATES --------------------
% -------------------------------------------------------------------------

if isfield(output,'CORR')

    % Now plot the R_statistic for each parameter
    figure('units','normalized','outerposition',[0 0 1 1])
    % Plot
    imagesc(output.CORR); colorbar;
    % Now change axis
    set(gca,'fontsize',16);
    %// adjust position of ticks
    set(gca,'XTick',(1:1:n_pars+1),'xticklabel',[]);
    set(gca,'YTick',(1:1:n_pars+1),'yticklabel',[]);
    for j = 1 : n_pars
        % Now add labels as well
        h = text(j,n_pars + .5 + 0.05*n_pars,str(j),...
            'fontsize',16,'interpreter','latex'); set(h, 'rotation', 90)
        text( .5 -0.03*n_pars,j,str(j),'fontsize',16,...
            'interpreter','latex');
    end
    % set labels
    set(gca,'xtick', linspace(0.5,n_pars-0.5,n_pars), 'ytick', ...
        linspace(0.5,n_pars-0.5,n_pars));
    set(gca,'xgrid', 'on', 'ygrid', 'on', 'gridlinestyle', '-', ...
        'xcolor', 'k', 'ycolor', 'k');
    % Add title
    evalstr = strcat('Map of posterior correlation coefficients of',...
        {' '},upper(method),{' '},'parameters');
    title(evalstr,'fontsize',18,'interpreter','latex');
    % Make font size 16
    set(gca,'fontsize',16);

end

% -------------------------------------------------------------------------
% -------- CORRELATION PLOTS OF THE POSTERIOR PARAMETER SAMPLES -----------
% -------------------------------------------------------------------------

if isfield(output,'CORR')

    if ismember(n_pars,2:30)
        % Open a new plot
        figure('units','normalized','outerposition',[0 0 1 1])
        % Plot a matrix (includes unscaled marginals on main diagonal!
        [H,AX,~,P,PAx] = plotmatrix(beta(:,1:n_pars),'+r'); hold on;
        % Now increase font size of each figure - except main diagonal
        set(AX,'fontsize',12);
        evalstr = strcat(['Marginal distribution and bivariate scatter ' ...
            'plots of posterior samples of'],...
            {' '},upper(method),{' '},'method');
        title(evalstr,'fontsize',18,'interpreter','latex');
        % label the plots
        for i = 1:length(AX)
            ylabel(AX(i,1),str(i),'fontsize',16);
            xlabel(AX(end,i),str(i),'fontsize',16);
        end
        % Now determine ranges
        %    mm1 = minmax(Pars(:,1:DREAMPar.d)'); dm = mm1(:,2) - mm1(:,1);
        mm1 = [ min(beta(:,1:n_pars))' max(beta(:,1:n_pars))' ];
        dm = mm1(:,2) - mm1(:,1);
        mm(:,1) = mm1(:,1) - dm/10; mm(:,2) = mm1(:,2) + dm/10;
        rr = 10.^(2 - order_no(mm)); mm = round( rr.*mm )./rr;
        minP = mm(:,1); maxP = mm(:,2);
        idx_not = find(minP == maxP); minP(idx_not) = mm1(idx_not,1);
        maxP(idx_not) = mm1(idx_not,2);

        % Now make sure that there are not too many bins in histograms
        for i = 1:n_pars
            if verLessThan('matlab','9.1')
                set(P(i),'Facecolor',[0.5 0.5 0.5],'EdgeColor','w');
            else
                set(P(i),'NumBins',20,'Facecolor',[0.5 0.5 0.5],...
                    'EdgeColor','w','BinLimits',[minP(i) maxP(i)]);
                %,'BinWidth',[maxP(i)-minP(i)]/20);
            end
            set(PAx(i),'Xlim',[minP(i) maxP(i)]);
        end
        for i = 1:n_pars * n_pars
            set(H(i),'Marker','.','markersize',12,'color',...
                [0.5 0.5 0.5]);
        end
        % Now change all axis to lie between minP and maxP
        for i = 1 : n_pars
            for j = 1 : n_pars
                hold(AX(j,i),'on');
                if ( j ~= i )
                    if verLessThan('matlab','9.1')
                        % do nothing
                    else
                        h = lsline(AX(j,i)); set(h,'color','k',...
                            'linestyle','--','linewidth',1);
                    end
                end
                set(AX(j,i),'Xlim',[minP(i) maxP(i)],'Ylim',...
                    [minP(j) maxP(j)]); %axis tight
            end
        end
        % Now add the prior ranges - if hypercube
        if isfield(par,'min')
            for i = 1 : n_pars
                for j = 1 : n_pars
                    if i ~= j
                        hold(AX(j,i),'on');
                        % Vertical lines
                        plot(AX(j,i),[ par.min(i) par.min(i) ],...
                            [ par.min(j) par.max(j)],'b','color',...
                            [0 0.4470 0.7410]);
                        plot(AX(j,i),[ par.max(i) par.max(i) ],...
                            [ par.min(j) par.max(j)],'b','color',...
                            [0 0.4470 0.7410]);
                        % Horizontal lines
                        plot(AX(j,i),[ par.min(i) par.max(i) ],...
                            [ par.min(j) par.min(j)],'b','color',...
                            [0 0.4470 0.7410]);
                        plot(AX(j,i),[ par.min(i) par.max(i) ],...
                            [ par.max(j) par.max(j)],'b','color',...
                            [0 0.4470 0.7410]);
                    elseif (i == j)
                        %                     hold(AX(i,j),'on');
                        %                     AX(i,j); line([ Par_info.min(i), Par_info.min(i)], [0 1.5*max(P(i).Values)],'linewidth',1,'color',[0 0.4470 0.7410]);
                        %                     %set(h,'linestyle',':','linewidth',2,'color',[0.5 0.5 0.5]);
                        % %                    line([ Par_info.max(i), Par_info.max(i)], ylim,'linewidth',1,'color',[0 0.4470 0.7410]);
                        %                     %set(h,'linestyle',':','linewidth',2,'color',[0.5 0.5 0.5]);
                    end
                end
            end
        end
    end

end

% -------------------------------------------------------------------------
% ------- CONVERGENCE OF INDIVIDUAL CHAINS TO TARGET DISTRIBUTION ---------
% -------------------------------------------------------------------------

% Check if R_statistic exists or not
if any(strcmp(method,{'bma','mma','mma-s'}))

    % Size of chain
    [T,~,N] = size(chain);
    % Define colors for different chains
    symbol = {'cs','rx','g+','ko','m<'};
    % create legend
    leg_string = {'chain 1'};
    % Add legend
    for i = 2 : min(5,N)
        leg = [ 'chain ' , num2str(i)];
        leg_string = [ leg_string {leg} ];
    end
    leg_string = [ leg_string {'ML'} ];
    % Now loop over each parameter
    for j = 1:n_pars
        % Open new figure if j = 3, 5, etc.
        if ( j/2 - floor(j/2) ) > 0
            figure('units','normalized','outerposition',...
                [0 0 1 1]); id_plot = 1;
        else
            id_plot = 2;
        end
        % Now plot a number of chains
        for i = 1:min(N,5)
            subplot(2,1,id_plot),plot(1:T,chain(1:end,j,i),...
                char(symbol(i)),'markersize',3,'linewidth',3);
            if i == 1; hold on; end
        end
        % Make font size 14
        set(gca,'fontsize',16);
        % Ranges have not been defined -- need to derive them from ParSet
        axis tight; a = axis; axis([1 T min(0.9*a(3),1.1*a(3)) ...
            max(0.9*a(4),1.1*a(4)) ])
        % Set ticks
        set(gca,'xtick',[1 T/5 2*T/5:T/5:T],'xticklabels',...
            [1 T/5 2*T/5:T/5:T]);
        % Lets add the MAP value
        plot( T , output.ML(j),'bx','Markersize',15,'linewidth',3);
        % Legend
        [~,objh,~,~] = legend(leg_string,'interpreter','latex');
        legend boxoff; try set(objh,'linewidth',2, ...
                'fontsize',14); catch, end
        % Add a title
        xlabel('Sample number of chain','fontsize',18);
        % Then add a y-label
        ylabel(str(j),'fontsize',18);
        % Then add title
        evalstr_t = strcat('DREAM$$_{\rm (ZS)}$$ sampled traceplot of',...
            {' '},upper(method),{' '},'parameter',{' '},str(j));
        title(evalstr_t,'fontsize',20,'interpreter','latex');
    end

end

% -------------------------------------------------------------------------
% ------ CONVERGENCE OF LIKELIHOOD AND PRIOR OF TARGET DISTRIBUTION -------
% -------------------------------------------------------------------------

% Check which method is used
if any(strcmp(method,{'bma','mma','mma-s'}))
    % Define colors for different chains
    symbol = {'cs','rx','g+','ko','m<'};
    % Size of chain
    [T,~,N] = size(chain);
    % create legend
    leg_string = {'chain 1'};
    % Add legend
    for i = 2 : min(5,N)
        leg = [ 'chain ' , num2str(i)];
        leg_string = [ leg_string {leg} ];
    end
    leg_string = [ leg_string {'ML'} ];
    % Open figure
    figure('units','normalized','outerposition',[0 0 1 1]);
    % Now loop - only the log-likelihood - uniform prior

    % Now plot a number of chains
    for i = 1:min(N,5)
        plot_chain = chain(1:end,n_pars+1,i);
        subplot(2,1,1),plot(1:T,plot_chain,char(symbol(i)),...
            'markersize',3,'linewidth',3);
        if i == 1; hold on; end
    end
    % Make font size 14
    set(gca,'fontsize',16);
    % Ranges have not been defined -- need to derive them from ParSet
    axis tight; a = axis; da = a(4) - a(3); a_m = a(3) - da/5; 
    a_p = a(4) + da/5;
    axis([1 T a_m a_p])
    set(gca,'xtick',[1 T/5 2*T/5:T/5:T],'xticklabels',[1 T/5 2*T/5:T/5:T]);
    % Lets add the MAP value
    plot( T , output.loglik ,'bx','Markersize',15,'linewidth',3);
    % Then add a y-label
    ylabel('$$\mathcal{L}(\textbf{x}|\widetilde{\textbf{Y}})$$',...
        'interpreter','latex','fontsize',18);
    % Legend
    [~,objh,~,~] = legend(leg_string,'interpreter','latex',...
        'location','southeast'); legend boxoff;
    try set(objh,'linewidth',2,'fontsize',14); catch, end
    % Then add title
    evalstr_t = strcat(['DREAM$$_{\rm (ZS)}$$ sampled traceplot of ' ...
        'log-likelihood using the'],...
        {' '},upper(method),{' '},'method',{' '});
    title(evalstr_t,'fontsize',20,'interpreter','latex');
    % Add a title
    xlabel('Sample number of chain','fontsize',18);

end

% -------------------------------------------------------------------------
% ------------ PLOT THE DATA, ENSEMBLE AND MEAN FORECAST  -----------------
% -------------------------------------------------------------------------

% Now plot time series plot of ensemble
figure('units','normalized','outerposition',[0 0 1 1]);
% Plot the ensemble
HA = plot(D,'linewidth',1.35); hold on;
if K < size(colororder,1)
    for i = 1:K
        set(HA(i),'color',colororder(i,1:3));
    end
end

% Plot the verifying observations
plot(1:N_meas,y,'ro','markersize',6,'linewidth',1,'markerfacecolor','w');

% Plot the averaged prediction
plot(output.Ye,'k','linewidth',1.5);

% Define axis
delta_Y = 0.07 * (max_Y - min_Y); min_X = 1; max_X = min(500,N_meas);
% Adjust axis
axis([ min_X max_X min_Y max_Y + delta_Y]);
set(gca,'fontsize',16);

% Add labels
xlabel('Row (index) of training data set','fontsize',18,...
    'interpreter','latex');
ylabel('Forecast','fontsize',18,'interpreter','latex');

% Add own legend
x_loc = min_X + 0.05*(max_X - min_X); y_loc = min_Y + 1*(max_Y - min_Y);
dx = 0.03*(max_X - min_X); dy = 0.04*(max_Y - min_Y);
x1 = x_loc; x2 = x1 + dx;
% Model entries
for k = 1:K
    line([x1 x2],[y_loc y_loc],'color',colororder(k,1:3),...
        'linewidth',4);
    text(x1 + 1.4*dx,y_loc,strcat('Model',{' '},num2str(k)),...
        'interpreter','latex','fontsize',16,'color',colororder(k,1:3));
    y_loc = y_loc - dy;
end
% Second entry
line([x1 x2],[y_loc y_loc],'color','k',...
    'linewidth',4);
text(x1+1.4*dx,y_loc,'$g_{j}^{\bullet}$: Weighted-average prediction',...
    'interpreter','latex','fontsize',16);
% Last entry
y_loc = y_loc - dy;
plot(x1+dx/2,y_loc,'ro','markersize',8,'linewidth',2,...
    'markerfacecolor','w');
text(x1+1.4*dx,y_loc,'Verifying observations',...
    'interpreter','latex','fontsize',16,'color','r');

% Create title
evalstr = strcat('Snapshot of model ensemble (colored lines) and,', ...
    {' '},...
    upper(method),{' '},['mean forecast (black line) and verifying ' ...
    'data (red circles)']);
title(evalstr,'fontsize',18,'interpreter','latex');
% Fontsize of all numbers
set(gca,'fontsize',16);

% -------------------------------------------------------------------------
% --------- CONFIDENCE AND PREDICTION INTERVALS OF AVERAGED MODEL ---------
% -------------------------------------------------------------------------

% Open new figure
figure('units','normalized','outerposition',[0 0 1 1])
% We start with the gamma% total uncertainty
nG = numel(gam)/2; clr = nan(nG,3);
for i = 1:nG
    clr(i,1:3) = ones(1,3) - i * 1/(2*nG);
    Fill_Ranges(1:N_meas,output.pred(:,i),output.pred(:,2*nG+1-i),clr(i,1:3));
    if i == 1; hold on, end
end
% Now we add gamma% uncertainty due to parameter uncertainty
Fill_Ranges(1:N_meas,output.par_unc(:,1),output.par_unc(:,end),...
    [0.2 0.2 0.2]);
% Plot the verifying observations
plot(1:N_meas,y,'ro','markersize',6,'linewidth',1,'markerfacecolor','w');
% Plot the averaged prediction
plot(output.Ye,'k','linewidth',1.5);
% Define axis
delta_Y = 0.07 * (max_Y - min_Y); min_X = 1; max_X = min(500,N_meas);
% Adjust axis
axis([ min_X max_X min_Y max_Y + delta_Y]);

% Add labels
xlabel('Row (index) of training data set','fontsize',18,...
    'interpreter','latex');
ylabel('Forecast','interpreter','latex','fontsize',18);
% Create title
evalstr = strcat(['Snapshot of confidence (dark gray) and prediction ' ...
    '(light gray) intervals:'],...
    {' '},upper(method),{' '},'method');
title(evalstr,'fontsize',18,'interpreter','latex');
% Add own legend
x_loc = min_X + 0.05*(max_X - min_X); y_loc = min_Y + 1*(max_Y - min_Y);
dx = 0.03*(max_X - min_X); dy = 0.07*(max_Y - min_Y);

leg_str = [];
% First entry
for i = 1:nG
    x1 = x_loc + (i-1)*dx/nG; x2 = x1 + dx/nG;
    patch([x1 x2 x2 x1],[y_loc y_loc y_loc+dy/2 y_loc+dy/2],clr(i,1:3));
    leg_str = strcat(leg_str,num2str(100*g(i)));
    if i < nG
        leg_str = strcat(leg_str,'/'); 
    else
        leg_str = strcat(leg_str,'\%'); 
    end
end
text(x_loc+1.4*dx,y_loc+dy/5,strcat(leg_str,{' '},...
    'Prediction interval'),'interpreter','latex','fontsize',16, ...
    'color','k');
y_loc = y_loc - dy/1.3;
% Second entry
leg_str = strcat(num2str(100*g(1)),'\%');
patch([x_loc x_loc + dx x_loc + dx x_loc],[y_loc y_loc y_loc+dy/2 ...
    y_loc+dy/2],[0.2 0.2 0.2]);
text(x_loc+1.4*dx,y_loc+dy/5,strcat(leg_str,{' '},...
    'Confidence interval'),'interpreter','latex','fontsize',16,...
    'color',[0.2 0.2 0.2]);
% Third entry
y_loc = y_loc - dy/1.7;
line([x_loc x_loc+dx],[y_loc y_loc],'color','k',...
    'linewidth',4);
text(x_loc+1.4*dx,y_loc,['$g_{j}^{\bullet}$: Weighted-average ' ...
    'prediction'],...
    'interpreter','latex','fontsize',16);
% Last entry
y_loc = y_loc - dy/1.3;
plot(x_loc+dx/2,y_loc,'ro','markersize',8,'linewidth',2,...
    'markerfacecolor','w');
text(x_loc+1.4*dx,y_loc,'Verifying observations',...
    'interpreter','latex','fontsize',16,'color','r');

% Add label with % contained
x_loc = min_X + 0.6*(max_X - min_X); y_loc = min_Y + 1.05*(max_Y - min_Y);
dx = 0.02*(max_X - min_X); dy = 0.04*(max_Y - min_Y);
% Make a patch and fill with white color
x1 = x_loc - dx; x2 = x_loc + 16*dx;
patch([x1 x2 x2 x1],[y_loc - (nG+1)*dy y_loc-(nG+1)*dy y_loc y_loc],...
    'w','facealpha',0.8);
y_loc = min_Y  + 1.01*(max_Y - min_Y);
for i = 1:nG
    evalstr = strcat('Coverage of',{' '},num2str(100*g(nG+1-i)), ...
        '\% pred interval is',{' '},...
        num2str(round(100*output.coverage(nG+1-i))/100),...
        '\%');
    text(x_loc,y_loc - (i-1)*dy,evalstr,'fontsize', ...
        16,'interpreter','latex','color',clr(i,1:3));
end
set(gca,'fontsize',16);

% -------------------------------------------------------------------------
% ------------------------ NOW DO RESIDUAL ANALYSIS -----------------------
% -------------------------------------------------------------------------

% Calculate the residual
res = output.Ye - y(:);

% Now plot ACF
figure('units','normalized','outerposition',[0 0 1 1]), autocorr(res);

% Add legend
str_leg = {char(strcat('\rm Sample autocorrelation of residuals', ...
    {' '},...
    upper(method),{' '},'method')),['\rm Confidence intervals for ' ...
    'white noise']};

% --> CHECK THIS LATER
% Add legend
objh = legend(str_leg,'interpreter','latex','fontsize',14);
set(objh,'linewidth',2); legend boxoff;
try set(objh,'interpreter','latex'), catch, end

% Then add title
evalstr = strcat(['Autocorrelation function of residuals of ' ...
    'deterministic (mean) forecast of'],...
    {' '},upper(method),{' '},'method');
title(evalstr,'fontsize',18,'interpreter','latex');

% Define ylabel
ylabel('Autocorrelation','fontsize',18,'interpreter','latex');
xlabel('Lag','fontsize',18,'interpreter','latex');
% Fontsize of all numbers
set(gca,'fontsize',16);

% Now plot QQ plot
figure('units','normalized','outerposition',[0 0 1 1]), qqplot(res);

% Define xlabel
xlabel('Standard normal quantiles','fontsize',18,'interpreter','latex');

% Define ylabel
ylabel('Quantiles of point forecast','fontsize',18,'interpreter','latex');

% Create title
evalstr = strcat(['Quantile-quantile plot of residuals of ' ...
    'deterministic (mean) forecast of'],...
    {' '},upper(method),{' '},'method');
title(evalstr,'fontsize',18,'interpreter','latex');
% Fontsize of all numbers
set(gca,'fontsize',16);

% -------------------------------------------------------------------------
% -------------------------------------------------------------------------