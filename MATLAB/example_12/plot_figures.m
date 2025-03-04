% ----------------------------------------------------------------------- %
%  MM   MM   OOOOO   DDDDDD   EEEEEEE  LL         AAA    VV   VV   GGGGG  %
%  MM   MM  OOOOOOO  DDDDDDD  EEEEEEE  LL        AA AA   VV   VV  GG   GG %
%  MMM MMM  OO   OO  DD   DD  EE       LL        AA AA   VV   VV  GG   GG %
%  MM M MM  OO   OO  DD   DD  EEEEE    LL       AA   AA  VV   VV  GGGGGGG %
%  MM   MM  OO   OO  DD   DD  EEEEE    LL       AAAAAAA  VV   VV   GGGGGG %
%  MM   MM  OO   OO  DD   DD  EE       LL       AA   AA   VV VV        GG %
%  MM   MM  OOOOOOO  DDDDDDD  EEEEEEE  LLLLLLL  AA   AA    VVV     GGGGGG %
%  MM   MM   OOOOO   DDDDDD   EEEEEEE  LLLLLLL  AA   AA     V     GGGGGGG %
% ----------------------------------------------------------------------- %

% This script plots the results of example_12.m
close all hidden;

which_fig = 3;
% -2: Scatter plot of log-likelihood, quadratic and logarithmic scores
% -1: The BMA ensemble - for portion of training period
%  0: 99/95/90/50% prediction intervals BMA model: 3x2 layout
%  1: 99/95/90/50% prediction intervals BMA model: 5x1 layout
%  2: Lognormal BMA prediction intervals and traces scoring rules
%  3: BMA mixture PDF at different times of evaluation period
%  4: Obsolete - practice

% ID: [1:7       Normal,logN,genN,truncN,gamma,gev,gpareto   group sigma [= constant variance]
%      8:14      Normal,logN,genN,truncN,gamma,gev,gpareto   ind. sigma [= constant variance]
%      15:21     Normal,logN,genN,truncN,gamma,gev,gpareto   group c [= nonconstant variance]
%      22:28     Normal,logN,genN,truncN,gamma,gev,gpareto   ind. c [= nonconstant variance]]

% Select BMA models to plot
ID = [22 23 24 25 26 27]; nID = numel(ID);
% Define maximum discharge value (y-axis) and range of x-values
Ymax = 12; YmaxSR = 8; t_id = 1650:1850; nt = numel(t_id);

% Load results of evaluation period [contains calibration results]
evalstr = strcat('load results_',which_period); eval(char(evalstr));
% Define labels and colors
fig_label = {'(a)','(b)','(c)','(d)','(e)','(f)'};
blue_color = [0 102 255]/255; red_color = [255 0 0]/255;
green_color = [51 153 51]/255; yellow_color = [255 192 0]/255;
orange_color = [122 92 113]/255; cyan_color = [0,255,255]/255;
fontsize_ax = 16; measdata_size = 3.5;
% Define model names of the ensemble
model_names = {'ABC','GR4J','HYMOD','TOPMO','AWBM','NAM','HBV','SAC-SMA'};

% Extract data
switch which_period
    case 'training'
        y = y_cal; nY = numel(y);
        [QS,LS,SS,CRPS,ES,Mu_mix] = deal(nan(nY,size(cal,2)));
        predQ = nan(nY,size(cal{1}.pred,2),size(cal,2));
        for i = 1:size(cal,2)
            QS(:,i) = cal{i}.QS;
            LS(:,i) = cal{i}.LS;
            SS(:,i) = cal{i}.SS;
            CRPS(:,i) = cal{i}.CRPS;
            ES(:,i) = cal{i}.ES;
            Mu_mix(:,i) = cal{i}.mu_mix;
            predQ(:,:,i) = cal{i}.pred;
        end
    case 'evaluation'
        y = y_eval; nY = numel(y);
        [QS,LS,SS,CRPS,ES,Mu_mix] = deal(nan(nY,size(cal,2)));
        predQ = nan(nY,size(cal{1}.pred,2),size(cal,2));
        for i = 1:size(val,2)
            QS(:,i) = val{i}.QS;
            LS(:,i) = val{i}.LS;
            SS(:,i) = val{i}.SS;
            CRPS(:,i) = val{i}.CRPS;
            ES(:,i) = val{i}.ES;
            Mu_mix(:,i) = val{i}.mu_mix;
            predQ(:,:,i) = val{i}.pred;
        end
end

% Some of the figures are large - not for a laptop screen
switch which_fig
    case -2 % Scatter plot of log-likelihood, quadratic and logarithmic scores
        measdata_size = 8;
        QSp = Table_8(1,:); LSp = Table_8(2,:); logL = Table_8(10,:);
        minL = 100*floor(min(logL)/100); maxL = 100*ceil(max(logL)/100);
        minQS = floor(min(QSp)); maxQS = ceil(max(QSp));
        minLS = floor(min(LSp)); maxLS = ceil(max(LSp));
        figure('unit','inches','PaperOrientation','portrait',...
            'position',[0.1 0.3 13.5 6.5]);
        ax1 = axes('units','inches'); axpos = [ 1.2 , 1 , 5 , 5 ];
        set(ax1,'position',axpos);
        % Now add data
        plot(ax1,logL,QSp,'bs','markerfacecolor','w','linewidth',2,...
            'markersize',measdata_size);
        % Adjust figure settings
        set(ax1,'tickdir','out','XMinorTick','on','YMinorTick','on',...
            'TickLength',[0.012 0.024],'box','off');
        set(ax1,'fontsize',fontsize_ax,'box','on');
        set(ax1,'box','off');
        axis(ax1,[minL maxL minQS maxQS]);
        xlabel(ax1,'Log-likelihood','interpreter','latex','fontsize',16);
        ylabel(ax1,'Quadratic score','interpreter','latex','fontsize',16);
        %set(ax1,'xtick',-3000:1000:0,'ytick',1.0:1.0:4.0);
        ax2 = axes('units','inches'); axpos = [ 8.0 , 1 , 5 , 5 ];
        set(ax2,'position',axpos);
        % Now add data
        plot(ax2,logL,LSp,'bs','markerfacecolor','w','linewidth',2,...
            'markersize',measdata_size);
        % Adjust figure settings
        set(ax2,'tickdir','out','XMinorTick','on','YMinorTick','on',...
            'TickLength',[0.012 0.024],'box','off');
        set(ax2,'fontsize',fontsize_ax,'box','on');
        %set(ax2,'xtick',-3000:1000:0,'ytick',-1.4:0.4:0.2);
        axis(ax2,[minL maxL minLS maxLS]);
        ax2.YAxis.MinorTickValues = [-1.6 -1.5 -1.3 -1.2 -1.1 -0.9 -0.8 ...
            -0.7 -0.5 -0.4 -0.3 -0.1 0 0.1];
        set(ax2,'box','off');
        xlabel(ax2,'Log-likelihood','interpreter','latex','fontsize',16);
        ylabel(ax2,'Logarithmic score','interpreter','latex','fontsize',16);
        
        set(gcf,'color','w');
        % print correlation Table of log-likelihood
        logL = Table_8(10,:);
        for u = 1:size(Table_8,1)
            R = corrcoef(logL,Table_8(u,:)); r(1,u) = R(1,2);
        end
        % Now print Table of Latex paper
        str = [];
        for u = 1:numel(r)
            if r(1,u) > 0
                str{u} = strcat('\hphantom{-}$',num2str(r(1,u),'%6.3f'),'$');
            else
                str{u} = strcat('-$',num2str(abs(r(1,u)),'%6.3f'),'$');
            end
        end
        fprintf('logL & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\',...
            str{1},str{2},str{3},str{4},str{5}, ...
            str{6},str{7},str{8},str{9},str{11},str{12},str{13});
        fprintf('\n');
        
    case -1 % The BMA ensemble - for portion of training period
        K = 8; blue_color = [0 102 255]/255;
        green_color = [51 153 51]/255; yellow_color = [246 252 4]/255;
        orange_color = [252 134 4]/255; cyan_color = [0,240,251]/255;
        purple_color = [141,37,219]/255; magenta_color = [240 16 208]/255;
        gray_color = [145 128 111]/255; rose_color = [235 21 210]/255;
        colors = [blue_color ; green_color ; yellow_color ; ...
            orange_color ; cyan_color ; purple_color ; ...
            rose_color ; gray_color ]; colors = flipud(colors);
        figure('unit','inches','PaperOrientation','portrait',...
            'position',[0.1 0.3 19 8]);
        t_id = 1:400; nt = numel(t_id);
        y_min = -1; y_max = 20;
        measdata_size = 7;
        ax1 = axes('units','inches'); axpos = [ 1.2 , 0.6 , 17 , 7 ];
        set(ax1,'position',axpos);
        for i = 1:K
            plot(ax1,t_id,D_cal(t_id,i),'color',colors(i,1:3),...
                'linewidth',2);
            if i == 1, hold on; end
            % add legend
            xloc1 = t_id(1) + 0.02*max(t_id); xloc2 = t_id(1) + ...
                0.05*max(t_id);
            y_loc = y_min + (0.95 - 0.06*(i-1)) * (y_max-y_min);
            line(ax1,[xloc1 xloc2],[y_loc y_loc],'color',colors(i,1:3),...
                'linewidth',3);
            text(ax1,xloc2 + 0.01*max(t_id),y_loc,char(model_names(i)),...
                'interpreter','latex','fontsize',16,'color',colors(i,1:3));
        end
        % Now add data
        plot(ax1,t_id,y_cal(t_id),'ro','markerfacecolor','w','linewidth',...
            2,'markersize',measdata_size);
        % Adjust figure settings
        set(ax1,'tickdir','out','XMinorTick','on','YMinorTick','on',...
            'TickLength',[0.012 0.024],'box','off');
        set(ax1,'fontsize',fontsize_ax,'box','on');
        ax1.YAxis.MinorTickValues = [-1 1 3 5 7 9 11 13 15 17 19 ];
        % Adjust axis
        axis([t_id(1) t_id(end) y_min y_max]);
        
    case 0 % 99/95/90/50% prediction intervals BMA model: 3x2 layout
        figure('unit','inches','PaperOrientation','portrait',...
            'position',[0.1 0.3 19 10.2]);
        y_min = [ -0.2*ones(1,nID) -YmaxSR ]; Ymax = 20;
        y_max = [ Ymax*ones(1,nID) YmaxSR ];
        cyan_color = [0,255,255]/255; row = 0; col = 1;
        for i = 1:nID
            if (i > nID/2) && (col == 1)
                col = col + 1; row = 1;
            else
                row = row + 1;
            end
            zz = ID(i);
            ax1 = axes('units','inches'); axpos = [ 1.2 + (col-1)*9 , ...
                7.3 - (row-1)*3.2 , 8 , 2.8 ]; set(ax1,'position',axpos);
            % Now we add 99% BMA uncertainty
            Fill_Ranges(t_id,predQ(t_id,1,zz),predQ(t_id,8,zz),...
                [0.9 0.9 0.9]); hold on;
            % Now we add 95% BMA uncertainty
            Fill_Ranges(t_id,predQ(t_id,2,zz),predQ(t_id,7,zz),...
                [0.7 0.7 0.7]);
            % Now we add 90% BMA uncertainty
            Fill_Ranges(t_id,predQ(t_id,3,zz),predQ(t_id,6,zz),...
                [0.5 0.5 0.5]);
            % Now we add 50% BMA uncertainty
            Fill_Ranges(t_id,predQ(t_id,4,zz),predQ(t_id,5,zz),...
                [0.3 0.3 0.3]);
            % Now add data
            plot(ax1,t_id,y(t_id),'rs','markerfacecolor','w',...
                'linewidth',1,'markersize',measdata_size);
            % And mean of BMA mixture
            plot(ax1,t_id,Mu_mix(t_id,zz),'k','linewidth',1.25);
            % Adjust axis
            axis([t_id(1) t_id(end) y_min(i) y_max(i)]);
            % Adjust figure settings
            set(ax1,'tickdir','out','XMinorTick','on','YMinorTick',...
                'on','TickLength',[0.012 0.024],'box','off');
            xtic = 1660:40:1820;
            set(ax1,'xtick',xtic,'xticklabel',xtic);
            ytic = 0:4:20;
            set(ax1,'ytick',ytic,'yticklabel',ytic);
            ax1.YAxis.MinorTickValues = [2 6 10 14 18];
            ax1.XAxis.MinorTickValues = [1650 1670 1680 1690 1710 ...
                1720 1730 1750 1760 1770 1790 1800 1810 1830 1840 1850];
            if (row ~= nID/2 && row ~= nID)
                set(ax1,'xticklabel',[]);
            end
            if col == 1
                % set ylabel
                ylabel(ax1,'${\rm Discharge,\;[mm/d]}$',...
                    'Interpreter','latex','fontsize',20);
            else
                set(ax1,'yticklabel',[]);
            end
            set(ax1,'fontsize',fontsize_ax,'box','on');
            % Add figure number (label)
            x_loc = t_id(1) + 0.005*(t_id(nt)-t_id(1)); y_loc = ...
                y_min(i) + 0.85*(y_max(i)-y_min(i));
            text(ax1,x_loc,y_loc,names_pdf(zz),'interpreter',...
                'latex','fontsize',fontsize_ax);
        end
        % Adjust figure settings
        set(ax1,'tickdir','out','XMinorTick','on','YMinorTick',...
            'on','TickLength',[0.01 0.02]);
        set(ax1,'fontsize',16,'box','on');
        % Add figure number (label)
        x_loc = t_id(1) + 0.005*(t_id(nt)-t_id(1)); y_loc = ...
            y_min(nID+1) + 0.85*(y_max(nID+1)-y_min(nID+1));
        % text(ax1,x_loc,y_loc,fig_label(nID+1),'interpreter','latex','fontsize',fontsize_ax);
        % Make all of figure white
        set(gcf,'color','w');
        
    case 1 % 99/95/90/50% prediction intervals BMA model: 5x1 layout of ID(1:5);
        figure('unit','inches','PaperOrientation','portrait',...
            'position',[0.1 0.3 12 10.5]);
        %ID = [ 1 2 3 4 5 ]; nID = numel(ID);
        fig_label = {'(a)','(b)','(c)','(d)','(e)'};
        y_min = -1/2*ones(1,nID); y_max = Ymax*ones(1,nID);
        cyan_color = [0,255,255]/255;
        t_id = 1650:1850; nt = numel(t_id);
        measdata_size = 4; blue_mean = [ 143 170 220]/255;
        for i = 1:nID - 1
            zz = ID(i);
            ax1 = axes('units','inches'); axpos = [ 1.2 , 8.56 - ...
                (i-1)*2.0 , 10.5 , 1.9 ]; set(ax1,'position',axpos);
            % Now we add 99% BMA uncertainty
            Fill_Ranges(t_id,predQ(t_id,1,zz),predQ(t_id,8,zz),...
                [0.94 0.94 0.94]); hold on;
            % Now we add 95% BMA uncertainty
            Fill_Ranges(t_id,predQ(t_id,2,zz),predQ(t_id,7,zz),...
                [0.70 0.70 0.70]);
            % Now we add 90% BMA uncertainty
            Fill_Ranges(t_id,predQ(t_id,3,zz),predQ(t_id,6,zz),...
                [0.5 0.5 0.5]);
            % Now we add 50% BMA uncertainty
            Fill_Ranges(t_id,predQ(t_id,4,zz),predQ(t_id,5,zz),...
                [0.3 0.3 0.3]);
            % Now add data
            plot(ax1,t_id,y(t_id),'ro','markerfacecolor','w',...
                'linewidth',1.25,'markersize',measdata_size);
            % And mean of BMA mixture
            plot(ax1,t_id,Mu_mix(t_id,zz),'color',blue_mean,...
                'linewidth',1.25);
            % Adjust axis
            axis([t_id(1) t_id(end) y_min(i) y_max(i)]);
            % Adjust figure settings
            set(ax1,'tickdir','out','XMinorTick','on','YMinorTick','on',...
                'TickLength',[0.012 0.024],'box','off');
            if i <= nID - 1
                set(ax1,'xticklabel',[],'xtick',[]);
            end
            % set ylabel
            %ylabel(ax1,'${\rm Discharge,\;[mm/d]}$','Interpreter','latex','fontsize',20);
            set(ax1,'fontsize',fontsize_ax);
            % Add figure number (label)
            x_loc = t_id(1) + 0.005*(t_id(nt)-t_id(1)); y_loc = y_min(i) ...
                + 0.85*(y_max(i)-y_min(i));
            text(ax1,x_loc,y_loc,names_pdf(zz),'interpreter','latex',...
                'fontsize',fontsize_ax);
        end
        % Adjust axis
        % axis([t_id(1) t_id(end) y_min(nID+1) y_max(nID+1)]);
        %set(ax1,'fontsize',16,'box','on');
        % Make all of figure white
        set(gcf,'color','w');
        
    case 2 % Figure of lognormal and traces scoring rules
        figure('unit','inches','PaperOrientation','portrait','position',...
            [0.1 0.3 12 10.5]);
        y_min = [ -1/2*ones(1,nID) -YmaxSR ]; y_max = [ Ymax*ones(1,nID) YmaxSR ];
        cyan_color = [0,255,255]/255;
        t_id = 1650:1850; nt = numel(t_id);
        blue_mean = [ 143 170 220]/255;
        % Which BMA forecast PDF + variance treatment? See index above
        zz = ID(2); measdata_size = 4;
        ax1 = axes('units','inches'); axpos = [ 1.2 , 7.4 , 10.5 , 3 ];
        set(ax1,'position',axpos);
        % ax1 = axes('units','inches'); axpos = [ 1.2 , 8.4 , 10 , 2 ]; set(ax1,'position',axpos);
        % Now we add 99% BMA uncertainty
        Fill_Ranges(t_id,predQ(t_id,1,zz),predQ(t_id,8,zz),[0.94 0.94 0.94]); hold on;
        % Now we add 95% BMA uncertainty
        Fill_Ranges(t_id,predQ(t_id,2,zz),predQ(t_id,7,zz),[0.70 0.70 0.70]);
        % Now we add 90% BMA uncertainty
        Fill_Ranges(t_id,predQ(t_id,3,zz),predQ(t_id,6,zz),[0.5 0.5 0.5]);
        % Now we add 50% BMA uncertainty
        Fill_Ranges(t_id,predQ(t_id,4,zz),predQ(t_id,5,zz),[0.3 0.3 0.3]);
        % Now add data
        plot(ax1,t_id,y(t_id),'ro','markerfacecolor','w','linewidth',...
            1.25,'markersize',measdata_size);
        % And mean of BMA mixture
        plot(ax1,t_id,Mu_mix(t_id,zz),'color',blue_mean,'linewidth',1.25);
        % Adjust axis
        axis([t_id(1) t_id(end) y_min(5) y_max(5)]);
        % Adjust figure settings
        set(ax1,'tickdir','out','XMinorTick','on','YMinorTick','on',...
            'TickLength',[0.012 0.024],'box','off');
        set(ax1,'xticklabel',[]);
        % set ylabel
        %ylabel(ax1,'${\rm Discharge,\;[mm/d]}$','Interpreter','latex','fontsize',20);
        set(ax1,'fontsize',fontsize_ax,'box','on');
        % Add figure number (label)
        x_loc = t_id(1) + 0.005*(t_id(nt)-t_id(1)); y_loc = y_min(1) + 0.85*12;
        text(ax1,x_loc,y_loc,names_pdf(zz),'interpreter','latex',...
            'fontsize',fontsize_ax);
        
        % Now we add scoring rules
        y_min = [-3 -5 -0.2 -2.2 -10] ; y_max = [3 1 1.7 0.2 1];
        for uu = 1:5
            % QS: green, LS: red, SS: blue, CRPS: yellow, ES: cyan
            ax1 = axes('units','inches');
            axpos = [ 1.2 , 7.2 - uu * 1.2 , 10.5 , 1.1 ];
            %                    axpos = [ 1.2 , 8.2 - uu * 0.9 , 10 , 0.8 ];
            set(ax1,'position',axpos);
            switch uu
                case 1
                    plot(ax1,t_id,QS(t_id,zz),'color',green_color,...
                        'linewidth',2); hold on
                    set(ax1,'ytick',[-2 0 2],'yticklabel',[-2 0 2]);
                    ax1.YAxis.MinorTickValues = [-3 -1 1 3];
                case 2
                    plot(ax1,t_id,LS(t_id,zz),'color',red_color,...
                        'linewidth',2);
                    set(ax1,'ytick',[-4 -2 0 ],'yticklabel',[-4 -2 0]);
                    ax1.YAxis.MinorTickValues = [-5 -3 -1 1 ];
                case 3
                    plot(ax1,t_id,SS(t_id,zz),'color',blue_color,...
                        'linewidth',2);
                    set(ax1,'ytick',[ 0 1 2],'yticklabel',[ 0 1 2]);
                    ax1.YAxis.MinorTickValues = [ 0.5 1.5 2.5];
                case 4
                    plot(ax1,t_id,CRPS(t_id,zz),'color',yellow_color,...
                        'linewidth',2);
                    set(ax1,'ytick',[-2 -1 0],'yticklabel',[-2 -1 0]);
                    ax1.YAxis.MinorTickValues = [-1.5 -0.5 0.5];
                case 5
                    plot(ax1,t_id,ES(t_id,zz),'color',cyan_color,...
                        'linewidth',2);
                    set(ax1,'ytick',[-8 -4 0],'yticklabel',[-8 -4 0]);
                    ax1.YAxis.MinorTickValues = [-10 -6 -2];
            end
            % Adjust axis
            axis([t_id(1) t_id(end) y_min(uu) y_max(uu)]);
            % Adjust figure settings
            set(ax1,'tickdir','out','XMinorTick','on','YMinorTick',...
                'on','TickLength',[0.01 0.02]);
            set(ax1,'fontsize',16,'box','on');
            if uu < 5
                set(ax1,'xticklabel',[],'xtick',[]);
            else
                % Add figure number (label)
                x_loc = t_id(1) + 0.005*(t_id(nt)-t_id(1)); y_loc = ...
                    y_min(end) + 0.85*(y_max(end)-y_min(end));
                text(ax1,x_loc,y_loc,names_pdf(zz),'interpreter',...
                    'latex','fontsize',fontsize_ax);
            end
        end
        % Make all of figure white
        set(gcf,'color','w');

    case 3 % BMA mixture PDF at different times of evaluation period
        % Desired times
        T = [1660 1700 1720 1740 1780 1800 1820 1840];
        for u = 1:numel(idx_col)
            idx = idx_col(u);
            % Get optimal parameter values
            a = beta{idx}; [~,ii] = max(a(:,end));
            x = a(ii(1),1:end-1);
            % Now right settings for options
            options.PDF = char(VPDF(idx,2)); options.VAR = char(VPDF(idx,1));
            % Now call function using evaluation period
            [X,pdfT{idx},cdfT{idx},predT{idx}] = BMA_mixture(x,D_eval,T,options);
        end
        % Create different colors
        clrs = parula; nc = size(clrs,1);
        clrs = clrs(1:floor((nc-1)/numel(idx_col)):nc,1:3);
        for t = 1:numel(T)
            subplot(2,4,t)
            for u = 1:numel(idx_col)
                plot(X,pdfT{idx_col(u)}(t,:),'color',clrs(u,1:3));
                if u == 1, hold on; end
            end
            xlim([0 10]), title(T(t))
            plot(y_eval(T(t)),0,'ro','linewidth',2,'markerfacecolor','w','markersize',6);
        end
        
    case 4 % Obsolete figure - for practice!
        figure('unit','inches','PaperOrientation','portrait',...
            'position',[0.1 0.3 19 10.5]);
        t_id = 1650:1850; Ymax = 18; nt = numel(t_id);
        y_min = [ 0 0 0 -12 ]; y_max = [ Ymax Ymax Ymax 12 ];
        for i = 1:3
            zz = ID(i);
            ax1 = axes('units','inches'); axpos = [ 1.2 , 8.2 - ...
                (i-1)*2.5 , 17 , 2.2 ]; set(ax1,'position',axpos);
            % Now we add 99% BMA uncertainty
            Fill_Ranges(t_id,predQ(t_id,1,zz),predQ(t_id,8,zz),...
                [0.9 0.9 0.9]); hold on;
            % Now we add 95% BMA uncertainty
            Fill_Ranges(t_id,predQ(t_id,2,zz),predQ(t_id,7,zz),...
                [0.7 0.7 0.7]);
            % Now we add 90% BMA uncertainty
            Fill_Ranges(t_id,predQ(t_id,3,zz),predQ(t_id,6,zz),...
                [0.5 0.5 0.5]);
            % Now we add 50% BMA uncertainty
            Fill_Ranges(t_id,predQ(t_id,4,zz),predQ(t_id,5,zz),...
                [0.3 0.3 0.3]);
            % Now add data
            plot(ax1,t_id,y(t_id),'rs','markerfacecolor','r','markersize',4);
            % Adjust axis
            axis([t_id(1) t_id(end) y_min(i) y_max(i)]);
            % Adjust figure settings
            set(ax1,'tickdir','out','XMinorTick','on','YMinorTick',...
                'on','TickLength',[0.01 0.02],'box','off');
            if i <= 3
                set(ax1,'xticklabel',[],'xtick',[]);
            end
            % set ylabel
            ylabel(ax1,'${\rm Discharge,\;[mm/d]}$','Interpreter',...
                'latex','fontsize',20);
            set(ax1,'fontsize',18);
            % Add figure number (label)
            x_loc = t_id(1) + 0.005*(t_id(nt)-t_id(1)); y_loc = ...
                y_min(i) + 0.9*(y_max(i)-y_min(i));
            text(ax1,x_loc,y_loc,names_pdf(zz),'interpreter','latex',...
                'fontsize',20);
        end
        % Now we add scoring rules
        ax1 = axes('units','inches'); axpos = [ 1.2 , 8.2 - 3*2.5 , ...
            17 , 2.2 ]; set(ax1,'position',axpos);
        % Now select which scoring rules to plot - of bottom ID
        zz = ID(3);
        % QS: green, LS: red, SS: blue, CRPS: yellow, ES: cyan
        plot(ax1,t_id,QS(t_id,zz),'color',green_color,'linewidth',2); hold on
        plot(ax1,t_id,LS(t_id,zz),'color',red_color,'linewidth',2);
        plot(ax1,t_id,SS(t_id,zz),'color',blue_color,'linewidth',2);
        plot(ax1,t_id,CRPS(t_id,zz),'color',yellow_color,'linewidth',2);
        plot(ax1,t_id,ES(t_id,zz),'color',cyan_color,'linewidth',2);
        % Adjust axis
        axis([t_id(1) t_id(end) -8 4.2]);
        % Adjust figure settings
        set(ax1,'tickdir','out','XMinorTick','on','YMinorTick','on',...
            'TickLength',[0.01 0.02],'box','off');
        set(ax1,'fontsize',18);
        % Add figure number (label)
        x_loc = t_id(1) + 0.005*(t_id(nt)-t_id(1)); y_loc = y_min(4) + ...
            0.9*(y_max(4)-y_min(4));
        text(ax1,x_loc,y_loc,names_pdf(zz),'interpreter','latex',...
            'fontsize',20);
        % Make all of figure white
        set(gcf,'color','w');
end