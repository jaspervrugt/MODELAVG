% Load results of HM23 [= from example_13 ]
% load results_HM23.mat

parset = genparset(chain);
output = cal{1};

[T,d1,N] = size(chain); n_pars = d1 - 1;
Tidx = [1:10:T];

% Check if R_statistic exists or not
%% if any(strcmp(method,{'bma','mma','mma-s'}))
% Figure 1

% Define colors for different chains
%symbol = {'cs','rx','g+','ko','m<'};
symbol = {'cs','rd','go'};  Nmax = 3;
% create legend
leg_string = {'chain 1'};
% Add legend
for i = 2 : min(N,Nmax)
    leg = [ 'chain ' , num2str(i)];
    leg_string = [ leg_string {leg} ];
end
leg_string = [ leg_string {'ML'} ];
% Size of chain
[T,~,N] = size(chain);
% Name of string
str = {'$\beta_{1}$','$\beta_{2}$','$\beta_{3}$','$\beta_{4}$',...
    '$\beta_{5}$','$\beta_{6}$','$\beta_{7}$','$\beta_{8}$',...
    '$\beta_{9}$','$\beta_{10}$','$\sigma_{1}$','$\sigma_{2}$',...
    '$\sigma_{3}$','$\sigma_{4}$','$\sigma_{5}$','$\sigma_{6}$',...
    '$\sigma_{7}$','$\sigma_{8}$','$\sigma_{9}$','$\sigma_{10}$'};
% 
idx_par = [1 2 3 11 12 13]; 
ymin = [-0.02 -0.02 -0.02 -0.01 -0.01 -0.2 ];
ymax = [ 1.0    1.0   1.0   0.4   0.4   10 ];
% Now loop over each parameter
for j = 1:6
    % Open new figure if j = 3, 5, etc.
    if ( j/2 - floor(j/2) ) > 0
        figure('units','normalized','outerposition',[0 0 1 1]); id_plot = 1;
    else
        id_plot = 2;
    end
    % Now plot a number of chains
    for i = 1:min(N,Nmax)
        subplot(2,1,id_plot),plot(Tidx,chain(Tidx,idx_par(j),i),...
            char(symbol(i)),'markersize',3,'linewidth',3); 
      %  axis tight; a = axis; axis([1 T ymin(j) ymax(j) ])
      %  set(gca,'xtick',[]); set(gca,'ytick',[]); 
      %  set(gca,'XAxisLocation','top'); h=gca;
      %  h.XAxis.Visible = 'off'; h.YAxis.Visible = 'off'; set(gcf,'color','w');
      %  pause
        if i == 1; hold on; end
    end
    % Make font size 14
    set(gca,'fontsize',16);
    % Ranges have not been defined -- need to derive them from ParSet
    axis tight; a = axis; axis([1 T ymin(j) ymax(j) ])
    % Set ticks
    set(gca,'xtick',[1 T/5 2*T/5:T/5:T],'xticklabels',[1 T/5 2*T/5:T/5:T]);
    % Lets add the MAP value
    plot( T , output.ML(j),'bx','Markersize',15,'linewidth',3);
    %gca.XMinorTick = 'on'; gca.YMinorTick = 'on';
    set(gca,'tickdir','out'); set(gca,'XMinorTick','on','YMinorTick','on');
    %XAxis.MinorTickValues = 1000*[1 2 3 5 6 7 9 10 11 13 14 15 17 18 19];
    %YAxis.MinorTickValues = [0.05 0.1 0.15 0.25 0.30 0.35 0.45 0.50 0.55 0.65 0.70 0.75 0.85 0.9 0.95];
    % Legend
    [~,objh,~,~] = legend(leg_string,'interpreter','latex'); legend boxoff; try set(objh,'linewidth',2,'fontsize',14); catch ME; end
    % Add a title
    xlabel('Sample number of chain','fontsize',18);
    % Then add a y-label
    ylabel(str(idx_par(j)),'fontsize',18);
%     % Then add title
%     evalstr_t = strcat('DREAM$$_{\rm (ZS)}$$ SAMPLED TRACE PLOT OF', ...
%         {' '},upper(method),{' '},'PARAMETER',{' '},str(idx_par(j)));
%     title(evalstr_t,'fontsize',20,'interpreter','latex');
end

N_meas = numel(y_cal);

% Figure 3: Prediction uncertainty
figure('units','normalized','outerposition',[0 0 1 1])
% We start with the total uncertainty
Fill_Ranges(1:N_meas,output.pred(:,2),output.pred(:,7),[0.75 0.75 0.75]); hold on;
% Now we add alpha% uncertainty due to parameter uncertainty
Fill_Ranges(1:N_meas,output.par_unc(:,2),output.par_unc(:,7),[0.5 0.5 0.5]);
% Now add the observations
plot(1:N_meas,y_cal,'r.','markersize',15);
% Plot the verifying observations
plot(output.Ye,'k','linewidth',1.5);
% Define axis
% delta_Y = 0.1 * (max_Y - min_Y);
% Adjust axis
axis([ 1 min(500,N_meas) 0 12]);
% Add labels
xlabel('Row (index) of training data set','fontsize',18, ...
        'fontweight','bold','fontname','times','interpreter','latex');
    ylabel('Forecast','interpreter','latex','fontsize',18, ...
        'fontweight','bold','fontname','times');
    % Legend
    [~,objh,~,~] = legend(str_leg,'interpreter','latex','fontsize',14);
    set(objh,'linewidth',2); legend boxoff;
    try set(objh,'interpreter','latex'), catch ME, end
    set(gca,'fontsize',18,'tickdir','out'); 
    set(gca,'xtick',[1 10:10:100],'xticklabel',[1 10:10:100]);
    set(gca,'XMinorTick','on','YMinorTick','on');
    % -------------------------------------------------------------------------
